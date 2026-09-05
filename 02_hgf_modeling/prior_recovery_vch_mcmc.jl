### prior_recovery_vch_mcmc.jl  ─  VCH HGF 4-Way Prior-Based Parameter & Model Recovery
###                                 (MCMC + Bridge Sampling)
###
### MCMC-based extension of prior_recovery_vch.jl.  Instead of MAP point
### estimates, this script fits each candidate model via full MCMC (NUTS,
### MAP init, 4 chains × 1000 iterations — matching the main analysis),
### then computes marginal evidence via bridge sampling (Meng & Wong, 1996;
### Gronau et al., 2017, J Math Psych) for the confusion matrix.
###
### For each (sim_index, true_model), this script:
###   1. Samples generative parameters from the prior of the true_model
###   2. Generates a VCH stimulus sequence (same 12-block structure as the task),
###      mapping condition labels to intensities using the TRUE model's convention
###   3. Simulates binary responses via the manual forward pass
###   4. Fits ALL FOUR candidate models via MCMC (NUTS, MAP init, 4 chains ×
###      1000 iterations), each using its OWN condition-intensity mapping.
###      For each fitted model:
###        · Extracts posterior medians
###        · Checks R-hat convergence (0.9–1.1)
###        · Computes bridge sampling LME using MVN proposal in unconstrained space
###   5. Records the winner (highest bridge sampling LME) and posterior medians
###      for all 4 models
###   6. Saves one row to param_recovery/prior_based_mcmc/results/sim{i}_{true_model}.csv
###
### Bridge sampling method:
###   Uses the iterative "optimal" bridge estimator (Meng & Wong, 1996) with
###   a multivariate normal proposal fitted to the MCMC posterior samples in
###   unconstrained (bijector-transformed) parameter space.  The unnormalized
###   log posterior in unconstrained space is evaluated via a manual forward
###   pass (log prior + log likelihood + Jacobian), using Bijectors.jl for
###   the constrained→unconstrained transforms.
###
### Candidate models:
###   · 2level_empiric  — binary_2level HGF, corrected empiric intensities
###   · 3level_empiric  — binary_3level HGF, corrected empiric intensities
###   · 2level_nominal  — binary_2level HGF, raw nominal proportions
###   · 3level_nominal  — binary_3level HGF, raw nominal proportions
###
### Output columns (sim{i}_{true_model}.csv):
###   sim_index, true_model,
###   gen_nu, gen_beta, gen_omega, gen_omega3  (gen_omega3=NaN for 2-level)
###   bridge_lme_2level_empiric, converged_2level_empiric,
###     mcmc_2level_empiric_nu, mcmc_2level_empiric_beta,
###     mcmc_2level_empiric_omega
###   bridge_lme_3level_empiric, converged_3level_empiric,
###     mcmc_3level_empiric_nu, mcmc_3level_empiric_beta,
###     mcmc_3level_empiric_omega, mcmc_3level_empiric_omega3
###   bridge_lme_2level_nominal, converged_2level_nominal,
###     mcmc_2level_nominal_nu, mcmc_2level_nominal_beta,
###     mcmc_2level_nominal_omega
###   bridge_lme_3level_nominal, converged_3level_nominal,
###     mcmc_3level_nominal_nu, mcmc_3level_nominal_beta,
###     mcmc_3level_nominal_omega, mcmc_3level_nominal_omega3
###   winner  (model key with the highest bridge sampling LME)
###
### Usage:
###   julia --project=. --threads 4 prior_recovery_vch_mcmc.jl <sim_index> <true_model>
###
###   sim_index  — Simulation index (1–N_SIM); also used as the random seed
###   true_model — One of: "2level_empiric", "3level_empiric",
###                        "2level_nominal",  "3level_nominal"

################################################################################
# ─── ROLE IN THE MANUSCRIPT ───────────────────────────────────────────────────
#
# PRIOR-BASED RECOVERY, part 1 of 3.  Produces Supplementary Fig. S3a, and the
# MCMC chains that the BIC confusion matrix (Fig. S3b) is computed from.
#
#     THIS SCRIPT ──▶ prior_recovery_aic_bic.jl ──▶ prior_recovery_aic_bic.py
#               └──▶ prior_recovery_vch_mcmc.py  (recovery scatter, Fig. S3a)
#
# The question it answers
#     Before asking which model best explains real data, one must ask whether the
#     models are distinguishable at all, and whether their parameters can be
#     recovered from data they themselves generated.  This is the face-validity
#     check that precedes model comparison.
#
# Design
#     For each of 500 simulation indices × 4 true models: draw parameters from the
#     TRUE model's prior, simulate a full 360-trial session using that model's own
#     stimulus mapping, then fit ALL FOUR candidate models to those responses —
#     each with its own mapping.  Fitting uses the same NUTS settings as the main
#     analysis, so recovery is assessed under the conditions actually used.
#
#     Marginal evidence is computed here by bridge sampling.  The manuscript's
#     reported confusion matrix instead uses BIC, computed from these same chains
#     by prior_recovery_aic_bic.jl; the bridge-sampling values are retained for
#     comparison.
#
# What it showed
#     nu and beta recovered well (r > 0.8); omega acceptably (r ≈ 0.4); omega3 not
#     at all (r < 0.1, posteriors within 0.01 of the prior).  The unrecoverable
#     omega3 is the substantive reason the 3-level model was not carried forward.
#
# HPC note — read before submitting
#     Bouchet has heterogeneous CPU types and Julia's precompilation cache is
#     CPU-target specific.  Export JULIA_CPU_TARGET=generic when precompiling AND
#     when running, or tasks landing on a non-matching node will spend their whole
#     wall time recompiling and produce nothing.
################################################################################

using ActionModels, HierarchicalGaussianFiltering
using Bijectors
using CSV
using Serialization
using Turing
using DataFrames
using Distributions
using LinearAlgebra
using LogExpFunctions
using MCMCChains
using Random
using Statistics: median, mean, cov

println("Packages loaded successfully.")
flush(stdout)
# flush(stdout) is called throughout to ensure progress messages appear in HPC
# job logs in real time — without it, output may be buffered until job completion.

include("helper_functions/create_agent.jl")

################################################################################
# ─── PARSE ARGUMENTS ──────────────────────────────────────────────────────────
################################################################################

const VALID_MODELS = ("2level_empiric", "3level_empiric",
                      "2level_nominal",  "3level_nominal")

if length(ARGS) < 2
    _valid_str = join(VALID_MODELS, ", ")
    error("Usage: julia prior_recovery_vch_mcmc.jl <sim_index> <true_model>\n" *
          "  true_model: one of $(_valid_str)")
end

sim_index  = parse(Int, ARGS[1])
true_model = ARGS[2]

if true_model ∉ VALID_MODELS
    _valid_str = join(VALID_MODELS, ", ")
    error("true_model must be one of $(_valid_str), got: \"$true_model\"")
end

################################################################################
# ─── MCMC SETTINGS (match model_fitting_singleagent_forarray.jl) ──────────────
################################################################################

const N_SAMPLES = 1000   # per chain
const N_CHAINS  = 4

################################################################################
# ─── VCH TASK STRUCTURE ───────────────────────────────────────────────────────
# Fixed condition counts per block  [0.0, 0.25, 0.50, 0.75]
################################################################################

const BLOCK_COUNTS = [
    [ 2,  1,  1, 26],   # Block  1
    [ 8,  4,  4, 14],   # Block  2
    [11,  5,  6,  8],   # Block  3
    [12,  7,  6,  5],   # Block  4
    [13,  7,  6,  4],   # Block  5
    [14,  7,  6,  3],   # Block  6
    [14,  7,  7,  2],   # Block  7
    [14,  7,  7,  2],   # Block  8
    [14,  7,  7,  2],   # Block  9
    [14,  7,  7,  2],   # Block 10
    [14,  7,  7,  2],   # Block 11
    [14,  7,  7,  2],   # Block 12
]

const COND_LEVELS = [0.0, 0.25, 0.50, 0.75]

# Non-hallucinator empirical intensities for 25/50/75 % (updated 2026-06-17).
# Source: out-of-set COPE cohort, max_vh/ah_freq == 0 or NaN, n=29 QC-passing VCH.
# See out_of_set_data/README_empirical_condition_redetermination.md
const EMPIRIC_INTENSITIES = [
    0.0,                   # condition 0.0  (ground-truth: no stimulus → 0.0)
    0.4180444024563061,    # condition 0.25
    0.7115104419621175,    # condition 0.50
    0.8994252873563219,    # condition 0.75
]

const NOMINAL_INTENSITIES = [0.0, 0.25, 0.50, 0.75]

const FITTED_MODELS = ["2level_empiric", "3level_empiric",
                       "2level_nominal",  "3level_nominal"]

################################################################################
# ─── PRIOR DISTRIBUTIONS (match model_fitting_singleagent_forarray.jl) ────────
################################################################################

const PRIORS_2LEVEL = (;
    action_precision       = truncated(Normal(0.29350739^(-1), 1), lower = 0.001),
    prior_posterior_weight = truncated(Normal(0.72646851, 1),      lower = 0),
    xprob_volatility       = truncated(Normal(-5.1682685, 1),      upper = -0.5),
)
const PRIORS_3LEVEL = (;
    action_precision       = truncated(Normal(0.29350739^(-1), 1), lower = 0.001),
    prior_posterior_weight = truncated(Normal(0.72646851, 1),      lower = 0),
    xprob_volatility       = truncated(Normal(-5.1682685, 1),      upper = -0.5),
    xvol_volatility        = truncated(Normal(-6, 1),              upper = -0.5),
)

################################################################################
# ─── HELPER FUNCTIONS ─────────────────────────────────────────────────────────
################################################################################

"""Return whether a model type uses the 3-level HGF architecture."""
is_3level(model_type::String) = startswith(model_type, "3level")

"""Return the HGF string for a given model type."""
get_hgf_string(model_type::String) = is_3level(model_type) ? "binary_3level" : "binary_2level"

"""Return the prior NamedTuple for a given model type."""
get_priors(model_type::String) = is_3level(model_type) ? PRIORS_3LEVEL : PRIORS_2LEVEL

"""
    conditions_to_intensities(cond_seq, model_type) → Vector{Float64}

Map a condition-label sequence {0.0, 0.25, 0.50, 0.75} to intensity values
using the stimulus convention of the given model_type.
"""
function conditions_to_intensities(cond_seq::Vector{Float64}, model_type::String)
    intensities = endswith(model_type, "_nominal") ? NOMINAL_INTENSITIES : EMPIRIC_INTENSITIES
    mapping = Dict(zip(COND_LEVELS, intensities))
    return Float64[mapping[c] for c in cond_seq]
end

"""
    generate_vch_conditions(rng) → Vector{Float64}

Generate a 360-trial VCH condition-label sequence matching the task block
structure.  Within each of the 12 blocks of 30 trials the fixed condition
multiset is preserved; trial order within each block is shuffled.
"""
function generate_vch_conditions(rng::AbstractRNG)
    seq = Float64[]
    for block_counts in BLOCK_COUNTS
        block = Float64[]
        for (n, cond) in zip(block_counts, COND_LEVELS)
            append!(block, fill(cond, n))
        end
        shuffle!(rng, block)
        append!(seq, block)
    end
    return seq
end

"""Sample one set of parameters independently from the prior of the given model type."""
function sample_from_prior(rng::AbstractRNG, model_type::String)
    priors = get_priors(model_type)
    return NamedTuple{keys(priors)}(rand(rng, d) for d in values(priors))
end

"""
Simulate binary VCH responses via the manual forward pass.
`intensities` is the stimulus strength sequence for the true_model's convention.
"""
function simulate_responses(rng::AbstractRNG, hgf_string::String,
                             gen_params::NamedTuple, intensities::Vector{Float64})
    am    = create_agent(hgf_string)
    agent = init_agent(am)
    set_parameters!(agent, gen_params)
    reset!(agent)

    responses = Float64[]
    for input in intensities
        dist   = agent.action_model(agent.model_attributes, input)
        action = Float64(rand(rng, dist))
        push!(responses, action)
        set_actions!(agent, :action, action)
    end
    return responses
end

################################################################################
# ─── BRIDGE SAMPLING ─────────────────────────────────────────────────────────
#
# Meng & Wong (1996) optimal bridge estimator with MVN proposal in
# unconstrained parameter space.
#
# References:
#   Meng & Wong (1996). Simulating ratios of normalizing constants via a
#     simple identity. Statistica Sinica, 6, 831-860.
#   Gronau et al. (2017). A tutorial on bridge sampling. J Math Psych, 81,
#     80-97.
################################################################################

"""
    make_log_joint_fn(priors, hgf_string, intensities, responses)

Return a function `θ_u::Vector{Float64} → Float64` that evaluates the log
unnormalized posterior in *unconstrained* (bijector-transformed) parameter
space:

    log p(y, θ_u) = log p(y | θ_c) + log π(θ_c) + log |∂θ_c/∂θ_u|

where θ_c = inverse(bij)(θ_u) for each parameter's bijector.

This replaces the broken DynamicPPL.LogDensityFunction approach, which
returned -Inf for all posterior samples when the Turing model was wrapped
by ActionModels (model.model was not data-conditioned in a way that
LogDensityFunction could evaluate).

The log likelihood is computed by a manual forward pass through the agent:
for each trial, the agent produces an action distribution, we evaluate it
at the observed response, then advance the agent state via set_actions!.
"""
function make_log_joint_fn(priors::NamedTuple, hgf_string::String,
                            intensities::Vector{Float64},
                            responses::Vector{Float64})
    param_keys  = keys(priors)
    prior_dists = [priors[p] for p in param_keys]
    bijs        = [bijector(d) for d in prior_dists]
    inv_bijs    = [inverse(b) for b in bijs]

    am    = create_agent(hgf_string)
    agent = init_agent(am)

    return function log_joint_unconstrained(θ_u::Vector{Float64})
        # Unconstrained → constrained
        θ_c = [inv_bijs[j](θ_u[j]) for j in eachindex(param_keys)]

        # Log prior (in constrained space)
        log_prior = 0.0
        for j in eachindex(param_keys)
            lp = logpdf(prior_dists[j], θ_c[j])
            isfinite(lp) || return -Inf
            log_prior += lp
        end

        # Jacobian correction: log |∂θ_c/∂θ_u|
        log_jac = 0.0
        for j in eachindex(param_keys)
            lj = logabsdetjac(inv_bijs[j], θ_u[j])
            isfinite(lj) || return -Inf
            log_jac += lj
        end

        # Log likelihood via forward pass
        θ_nt = NamedTuple{param_keys}(θ_c)
        set_parameters!(agent, θ_nt)
        reset!(agent)

        log_lik = 0.0
        for (inp, resp) in zip(intensities, responses)
            dist = agent.action_model(agent.model_attributes, inp)
            lp   = logpdf(dist, resp)
            isfinite(lp) || return -Inf
            log_lik += lp
            set_actions!(agent, :action, resp)
        end

        return log_prior + log_jac + log_lik
    end
end

"""
    bridge_sampling_lme(log_density_fn, chains, priors_nt; max_iter, tol)

Estimate log marginal evidence via bridge sampling.

Arguments:
- `log_density_fn`: callable f(θ_u::Vector{Float64}) → Float64, returning
  the log unnormalized posterior in unconstrained space (log joint + Jacobian).
  Build with make_log_joint_fn.
- `chains`: MCMCChains.Chains object with posterior samples (constrained space)
- `priors_nt`: NamedTuple of prior distributions (keys = parameter names)

Returns log p(y | M), or NaN if bridge sampling fails.
"""
function bridge_sampling_lme(log_density_fn, chains, priors_nt;
                              max_iter::Int=500, tol::Float64=1e-10)
    param_keys = keys(priors_nt)
    k = length(param_keys)

    # ── Get bijectors for constrained → unconstrained transformation ──────
    prior_dists = [priors_nt[p] for p in param_keys]
    bijs = [bijector(d) for d in prior_dists]

    # ── Extract constrained posterior samples from MCMCChains ─────────────
    session_names = [Symbol("$(p).session[1]") for p in param_keys]
    chains_df = DataFrame(chains)

    n_post = nrow(chains_df)
    θ_c = Matrix{Float64}(undef, n_post, k)
    for (j, sn) in enumerate(session_names)
        θ_c[:, j] = chains_df[!, sn]
    end

    # ── Transform to unconstrained space ─────────────────────────────────
    θ_u_post = Matrix{Float64}(undef, n_post, k)
    for i in 1:n_post
        for j in 1:k
            θ_u_post[i, j] = bijs[j](θ_c[i, j])
        end
    end

    # Drop non-finite samples (rare edge cases at bijector boundaries)
    valid_mask = vec(all(isfinite, θ_u_post; dims=2))
    if sum(valid_mask) < n_post
        n_bad = n_post - sum(valid_mask)
        @warn "Dropping $n_bad non-finite unconstrained posterior samples"
        θ_u_post = θ_u_post[valid_mask, :]
        n_post = size(θ_u_post, 1)
    end

    if n_post < 50
        @warn "Too few valid posterior samples ($n_post) for bridge sampling"
        return NaN
    end

    # ── Fit MVN proposal in unconstrained space ──────────────────────────
    μ = vec(mean(θ_u_post, dims=1))
    Σ = cov(θ_u_post)
    Σ = Symmetric(Σ + 1e-6 * I(k))   # regularise for numerical stability

    if !isposdef(Σ)
        @warn "Proposal covariance not positive definite — bridge sampling failed"
        return NaN
    end

    proposal = MvNormal(μ, Σ)

    # ── Draw proposal samples ────────────────────────────────────────────
    n_prop = n_post   # match posterior count
    rng_prop = MersenneTwister(42)
    θ_u_prop = Matrix{Float64}(rand(rng_prop, proposal, n_prop)')   # n_prop × k

    # ── Evaluate log unnormalized density at all samples ──────────────────
    # g(θ_u) = log p(y, θ) in unconstrained space (includes Jacobian)
    log_g_post = Vector{Float64}(undef, n_post)
    log_g_prop = Vector{Float64}(undef, n_prop)

    for i in 1:n_post
        log_g_post[i] = try
            log_density_fn(θ_u_post[i, :])
        catch
            -Inf
        end
    end

    for i in 1:n_prop
        log_g_prop[i] = try
            log_density_fn(θ_u_prop[i, :])
        catch
            -Inf
        end
    end

    # Filter non-finite posterior evaluations (should be rare)
    post_finite = isfinite.(log_g_post)
    if sum(post_finite) < n_post
        n_bad = n_post - sum(post_finite)
        @warn "$n_bad non-finite log joint evaluations at posterior samples"
        log_g_post = log_g_post[post_finite]
        θ_u_post   = θ_u_post[post_finite, :]
        n_post     = sum(post_finite)
    end

    if n_post < 50
        @warn "Too few finite posterior evaluations ($n_post) for bridge sampling"
        return NaN
    end

    # ── Evaluate log proposal density at all samples ─────────────────────
    log_q_post = [logpdf(proposal, θ_u_post[i, :]) for i in 1:n_post]
    log_q_prop = [logpdf(proposal, θ_u_prop[i, :]) for i in 1:n_prop]

    # ── Bridge sampling iteration (Meng & Wong 1996, optimal bridge) ─────
    l1 = n_post / (n_post + n_prop)
    l2 = n_prop / (n_post + n_prop)
    log_l1 = log(l1)
    log_l2 = log(l2)

    # Initialise with importance sampling estimate
    is_terms = [log_g_prop[i] - log_q_prop[i] for i in 1:n_prop if isfinite(log_g_prop[i])]
    if isempty(is_terms)
        @warn "All importance sampling weights non-finite — bridge sampling failed"
        return NaN
    end
    log_r = logsumexp(is_terms) - log(length(is_terms))

    for iter in 1:max_iter
        log_r_old = log_r

        # Numerator: average over proposal samples
        log_num_terms = Float64[]
        for i in 1:n_prop
            if !isfinite(log_g_prop[i]); continue; end
            denom = logaddexp(log_l1 + log_g_prop[i],
                              log_l2 + log_r + log_q_prop[i])
            push!(log_num_terms, log_g_prop[i] - denom)
        end

        if isempty(log_num_terms)
            @warn "No finite numerator terms at iteration $iter"
            return NaN
        end

        # Denominator: average over posterior samples
        log_den_terms = Float64[]
        for i in 1:n_post
            denom = logaddexp(log_l1 + log_g_post[i],
                              log_l2 + log_r + log_q_post[i])
            push!(log_den_terms, log_q_post[i] - denom)
        end

        if isempty(log_den_terms)
            @warn "No finite denominator terms at iteration $iter"
            return NaN
        end

        log_r = (logsumexp(log_num_terms) - log(length(log_num_terms))) -
                (logsumexp(log_den_terms) - log(length(log_den_terms)))

        if abs(log_r - log_r_old) < tol
            println("    Bridge sampling converged after $iter iterations")
            break
        end
    end

    return log_r
end

################################################################################
# ─── FIT A SINGLE MODEL VIA MCMC + BRIDGE SAMPLING LME ──────────────────────
################################################################################

"""
    fit_mcmc_and_bridge_lme(hgf_string, priors, intensities, responses;
                             n_samples, n_chains)

Fit a model via full MCMC (NUTS, MAP init), compute posterior medians,
check convergence, and estimate marginal evidence via bridge sampling.

The log joint density function is constructed directly from the Turing
model via DynamicPPL.LogDensityFunction — no MAP required.

Returns (bridge_lme, medians_dict, converged)
where:
  bridge_lme   — bridge sampling log marginal evidence (NaN if failed)
  medians_dict — (; nu, beta, omega, omega3)
  converged    — Bool (all R-hats in [0.9, 1.1])
"""
function fit_mcmc_and_bridge_lme(hgf_string::String, priors::NamedTuple,
                                  intensities::Vector{Float64},
                                  responses::Vector{Float64};
                                  n_samples::Int=N_SAMPLES,
                                  n_chains::Int=N_CHAINS,
                                  chains_file::Union{String,Nothing}=nothing)

    # ── MCMC ──────────────────────────────────────────────────────────────
    am_mcmc = create_agent(hgf_string)
    modelfit_mcmc = create_model(am_mcmc, priors, intensities, responses;
                                  check_parameter_rejections=true, verbose=false)

    println("    Running MCMC ($n_samples samples × $n_chains chains, MAP init)...")
    flush(stdout)
    my_sampler = NUTS(adtype=AutoForwardDiff())
    chains = sample_posterior!(modelfit_mcmc, MCMCSerial();
        sampler     = my_sampler,
        init_params = :MAP,
        n_samples   = n_samples,
        n_chains    = n_chains,
        progress    = false,
    )
    println("    MCMC complete.")
    flush(stdout)

    # ── Save chains (mirrors main pipeline convention) ───────────────────
    if !isnothing(chains_file)
        serialize(chains_file, chains)
        CSV.write(replace(chains_file, r"\.jls$" => ".csv"), DataFrame(chains))
        println("    Chains saved: $chains_file")
        flush(stdout)
    end

    # ── Extract posterior medians ────────────────────────────────────────
    param_names = keys(priors)
    session_params = get_session_parameters!(modelfit_mcmc)
    session_id_sym = Symbol(session_params.session_ids[1])

    medians_dict = Dict{Symbol, Float64}()
    for pn in param_names
        medians_dict[pn] = median(session_params.value[pn][session_id_sym])
    end

    # Friendly names for output CSV
    nu    = medians_dict[:prior_posterior_weight]
    beta  = medians_dict[:action_precision]
    omega = medians_dict[:xprob_volatility]
    omega3 = get(medians_dict, :xvol_volatility, NaN)

    friendly_medians = (; nu=nu, beta=beta, omega=omega, omega3=omega3)

    # ── R-hats → convergence flag ───────────────────────────────────────
    # Threshold matches the main analysis pipeline (import_hgf_results_unified.py
    # and param_recovery_vch.py: RHAT_LOWER=0.9, RHAT_UPPER=1.1).
    rhats_df = DataFrame(rhat(chains))
    converged = true
    for row in eachrow(rhats_df)
        r = row.rhat
        if !isfinite(r) || r < 0.9 || r > 1.1
            converged = false
            break
        end
    end

    # ── Log density function: manual forward pass (replaces DynamicPPL) ────
    # DynamicPPL.LogDensityFunction returned -Inf for all posterior samples
    # when the model was wrapped by ActionModels (the internal model property
    # was not data-conditioned in a compatible way).  We instead compute the
    # log joint directly: log prior + log likelihood (forward pass) + Jacobian.
    log_joint_fn = make_log_joint_fn(priors, hgf_string, intensities, responses)

    # ── Bridge sampling LME ──────────────────────────────────────────────
    println("    Running bridge sampling...")
    flush(stdout)
    bridge_lme = try
        bridge_sampling_lme(log_joint_fn, chains, priors)
    catch e
        @warn "Bridge sampling failed for $hgf_string: $e"
        NaN
    end
    println("    Bridge sampling complete: log p(y|M) ≈ $(round(bridge_lme, digits=3))")
    flush(stdout)

    return bridge_lme, friendly_medians, converged
end

################################################################################
# ─── PATHS AND SKIP CHECK ─────────────────────────────────────────────────────
################################################################################

out_dir    = joinpath("param_recovery", "prior_based_mcmc", "results")
chains_dir = joinpath("param_recovery", "prior_based_mcmc", "chains")
mkpath(out_dir)
mkpath(chains_dir)
out_file = joinpath(out_dir, "sim$(sim_index)_$(true_model).csv")

println("=" ^ 70)
println("VCH Prior-Based Recovery — MCMC + Bridge Sampling (4-way)")
println("  sim_index  = $sim_index")
println("  true_model = $true_model")
println("  MCMC: $N_SAMPLES samples × $N_CHAINS chains (MAP init, NUTS)")
println("=" ^ 70)

if isfile(out_file)
    println("Already complete — skipping ($out_file found).")
    exit(0)
end

################################################################################
# ─── MAIN ─────────────────────────────────────────────────────────────────────
################################################################################

rng = MersenneTwister(sim_index)

# ── Step 1: Sample generative parameters from the true model's prior ─────────
gen_params = sample_from_prior(rng, true_model)
gen_nu     = gen_params.prior_posterior_weight
gen_beta   = gen_params.action_precision
gen_omega  = gen_params.xprob_volatility
gen_omega3 = hasproperty(gen_params, :xvol_volatility) ? gen_params.xvol_volatility : NaN

println("Generative parameters (sampled from prior of $true_model):")
println("  nu    = $(round(gen_nu,    digits=4))")
println("  beta  = $(round(gen_beta,  digits=4))")
println("  omega = $(round(gen_omega, digits=4))")
is_3level(true_model) && println("  omega3= $(round(gen_omega3, digits=4))")
flush(stdout)

# ── Step 2: Generate condition-label sequence ─────────────────────────────────
cond_sequence = generate_vch_conditions(rng)
println("\nGenerated $(length(cond_sequence)) condition labels (12 blocks × 30 trials)")
flush(stdout)

# ── Step 3: Simulate binary responses using the TRUE model's intensities ──────
true_intensities = conditions_to_intensities(cond_sequence, true_model)
true_hgf         = get_hgf_string(true_model)
responses        = simulate_responses(rng, true_hgf, gen_params, true_intensities)
n_detect         = Int(sum(responses))
println("Simulated $(length(responses)) responses: $n_detect detect, "
        * "$(length(responses)-n_detect) miss")
flush(stdout)

# ── Step 4: Fit all four candidate models via MCMC + bridge sampling ──────────
fit_bridge_lme = Dict{String, Float64}()
fit_medians    = Dict{String, NamedTuple}()
fit_converged  = Dict{String, Bool}()

for fitted_model in FITTED_MODELS
    println("\n" * "─" ^ 60)
    println("Fitting $fitted_model via MCMC ($N_SAMPLES × $N_CHAINS chains) ...")
    flush(stdout)
    fit_intensities = conditions_to_intensities(cond_sequence, fitted_model)
    fit_hgf         = get_hgf_string(fitted_model)
    fit_prior       = get_priors(fitted_model)

    chains_file = joinpath(chains_dir, "sim$(sim_index)_$(true_model)_$(fitted_model).jls")

    bridge_lme, medians, converged = NaN, (; nu=NaN, beta=NaN, omega=NaN, omega3=NaN), false
    try
        bridge_lme, medians, converged =
            fit_mcmc_and_bridge_lme(fit_hgf, fit_prior, fit_intensities, responses;
                                    chains_file=chains_file)
    catch e
        @warn "fit_mcmc_and_bridge_lme failed for $fitted_model: $e"
    end

    fit_bridge_lme[fitted_model] = bridge_lme
    fit_medians[fitted_model]    = medians
    fit_converged[fitted_model]  = converged

    println("  bridge_lme = $(round(bridge_lme, digits=3))")
    println("  converged  = $converged")
    flush(stdout)
end

# ── Step 5: Select winner (highest finite bridge LME) ─────────────────────────
finite_lmes = [(m, fit_bridge_lme[m]) for m in FITTED_MODELS if isfinite(fit_bridge_lme[m])]
if isempty(finite_lmes)
    winner = "none"
    @warn "All model fits failed — winner = \"none\""
else
    winner = finite_lmes[argmax(last.(finite_lmes))][1]
end
println("\n" * "=" ^ 60)
println("Winner: $winner")
flush(stdout)

# ── Step 6: Build and save output row ─────────────────────────────────────────
row = Dict{String, Any}(
    "sim_index"  => sim_index,
    "true_model" => true_model,
    "gen_nu"     => gen_nu,
    "gen_beta"   => gen_beta,
    "gen_omega"  => gen_omega,
    "gen_omega3" => gen_omega3,
    "winner"     => winner,
)

for fitted_model in FITTED_MODELS
    medians = fit_medians[fitted_model]

    row["bridge_lme_$(fitted_model)"] = fit_bridge_lme[fitted_model]
    row["converged_$(fitted_model)"]  = fit_converged[fitted_model]

    row["mcmc_$(fitted_model)_nu"]    = medians.nu
    row["mcmc_$(fitted_model)_beta"]  = medians.beta
    row["mcmc_$(fitted_model)_omega"] = medians.omega
    if is_3level(fitted_model)
        row["mcmc_$(fitted_model)_omega3"] = medians.omega3
    end
end

CSV.write(out_file, DataFrame([row]))
println("\nSaved: $out_file")
println("Done — sim_index=$sim_index, true_model=$true_model, winner=$winner")
flush(stdout)
