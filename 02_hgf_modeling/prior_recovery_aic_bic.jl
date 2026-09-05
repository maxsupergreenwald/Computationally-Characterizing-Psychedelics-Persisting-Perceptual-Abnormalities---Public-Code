### prior_recovery_aic_bic.jl
###
### Patch script: adds AIC and BIC to prior-based recovery result CSVs.
###
### For each (sim_index, true_model) this script:
###   1. Replays the RNG (same seed = sim_index, same call sequence as
###      prior_recovery_vch_mcmc.jl) to reconstruct the exact generative
###      parameters, condition sequence, and simulated responses.
###   2. Cross-checks replayed gen_params against the saved result CSV to
###      guarantee that the log-likelihood is computed for the correct
###      participant.  Any mismatch beyond 1e-8 is a fatal error.
###   3. For each of the 4 fitted models, loads the saved MCMC chains
###      (.jls) and evaluates the total log-likelihood at every posterior
###      sample via the HGF forward pass.
###   4. Computes AIC and BIC using the maximum log-likelihood over samples
###      as an approximation of the MLE log-likelihood (log L̂):
###
###         k       = number of free parameters (3 for 2-level, 4 for 3-level)
###         n       = 360 trials
###         log L̂  = max_{θ in chains} log p(responses | θ, fitted_model)
###         AIC     = 2k − 2 log L̂
###         BIC     = k log(n) − 2 log L̂
###
###   5. Appends these columns to the existing result CSV (in-place update):
###        aic_{model}, bic_{model}, log_lik_max_{model},
###        log_lik_mean_{model}, log_lik_sd_{model},
###        n_samples_lik_{model}   (for each of the 4 fitted models)
###        aic_winner, bic_winner
###
### Log-likelihood calculation:
###   Uses make_log_lik_fn, a stripped version of the bridge sampling
###   make_log_joint_fn (prior term and Jacobian removed, constrained
###   parameter space only).  This is mathematically identical to what
###   DynamicPPL.pointwise_loglikelihoods computes on the ActionModels
###   model: in the NoMissingActions path ActionModels batches all 360
###   responses into a single arraydist, so pointwise_loglikelihoods
###   returns one total log-likelihood per sample — exactly what the
###   forward pass produces.
###
### Usage:
###   julia --project=. --threads 1 prior_recovery_aic_bic.jl <sim_index> <true_model>
###
###   sim_index  — Simulation index (1–500); also used as the random seed.
###   true_model — One of: "2level_empiric", "3level_empiric",
###                        "2level_nominal",  "3level_nominal"

################################################################################
# ─── ROLE IN THE MANUSCRIPT ───────────────────────────────────────────────────
#
# PRIOR-BASED RECOVERY, part 3 of 3.  Produces the evidence values behind
# Supplementary Fig. S3b (the published model-identifiability confusion matrix).
#
#     prior_recovery_vch_mcmc.jl ──▶ THIS SCRIPT ──▶ prior_recovery_aic_bic.py
#
# What it does and why it is a separate script
#     It refits nothing.  It reopens the chains already written by
#     prior_recovery_vch_mcmc.jl and scores each of the four candidate models by
#     AIC and BIC, appending those columns to the existing per-simulation CSVs.
#     Keeping it separate means the expensive MCMC never has to be repeated in
#     order to add an information criterion.
#
#     BIC is what the manuscript reports: model identifiability is defined as the
#     proportion of the 500 generative samples for which the data-generating model
#     achieved the lowest BIC after inversion.
#
# The RNG-replay guard
#     Because the generative parameters are not re-derived but reconstructed, the
#     script replays the original random stream (seed = sim_index, same call
#     order) and then CROSS-CHECKS the replayed parameters against the values
#     stored in the saved CSV.  Any disagreement beyond 1e-8 is a fatal error, not
#     a warning — a mismatch would mean scoring one simulation's likelihood
#     against another's data, and silently producing a plausible-looking but wrong
#     confusion matrix.
################################################################################

using ActionModels, HierarchicalGaussianFiltering
using CSV
using Serialization
using DataFrames
using Distributions
using MCMCChains
using Random
using Statistics: mean, std

println("Packages loaded.")
flush(stdout)

include("helper_functions/create_agent.jl")

################################################################################
# ─── PARSE ARGUMENTS ──────────────────────────────────────────────────────────
################################################################################

const VALID_MODELS = ("2level_empiric", "3level_empiric",
                      "2level_nominal",  "3level_nominal")

if length(ARGS) < 2
    _valid_str = join(VALID_MODELS, ", ")
    error("Usage: julia prior_recovery_aic_bic.jl <sim_index> <true_model>\n" *
          "  true_model: one of $(_valid_str)")
end

sim_index  = parse(Int, ARGS[1])
true_model = ARGS[2]
true_model ∉ VALID_MODELS &&
    error("true_model must be one of $(VALID_MODELS), got: \"$true_model\"")

################################################################################
# ─── VCH TASK STRUCTURE (identical to prior_recovery_vch_mcmc.jl) ─────────────
################################################################################

const BLOCK_COUNTS = [
    [ 2,  1,  1, 26],
    [ 8,  4,  4, 14],
    [11,  5,  6,  8],
    [12,  7,  6,  5],
    [13,  7,  6,  4],
    [14,  7,  6,  3],
    [14,  7,  7,  2],
    [14,  7,  7,  2],
    [14,  7,  7,  2],
    [14,  7,  7,  2],
    [14,  7,  7,  2],
    [14,  7,  7,  2],
]

const COND_LEVELS          = [0.0, 0.25, 0.50, 0.75]
const EMPIRIC_INTENSITIES  = [0.0, 0.4180444024563061,
                               0.7115104419621175, 0.8994252873563219]
const NOMINAL_INTENSITIES  = [0.0, 0.25, 0.50, 0.75]
const FITTED_MODELS        = ["2level_empiric", "3level_empiric",
                               "2level_nominal",  "3level_nominal"]
const N_TRIALS             = 360

################################################################################
# ─── PRIOR DISTRIBUTIONS (identical to prior_recovery_vch_mcmc.jl) ────────────
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
# ─── HELPER FUNCTIONS (identical to prior_recovery_vch_mcmc.jl) ───────────────
################################################################################

is_3level(m::String)      = startswith(m, "3level")
get_hgf_string(m::String) = is_3level(m) ? "binary_3level" : "binary_2level"
get_priors(m::String)     = is_3level(m) ? PRIORS_3LEVEL : PRIORS_2LEVEL
get_n_params(m::String)   = is_3level(m) ? 4 : 3

function conditions_to_intensities(cond_seq::Vector{Float64}, model_type::String)
    intensities = endswith(model_type, "_nominal") ? NOMINAL_INTENSITIES : EMPIRIC_INTENSITIES
    mapping = Dict(zip(COND_LEVELS, intensities))
    return Float64[mapping[c] for c in cond_seq]
end

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

function sample_from_prior(rng::AbstractRNG, model_type::String)
    priors = get_priors(model_type)
    return NamedTuple{keys(priors)}(rand(rng, d) for d in values(priors))
end

function simulate_responses(rng::AbstractRNG, hgf_string::String,
                             gen_params::NamedTuple, intensities::Vector{Float64})
    am    = create_agent(hgf_string)
    agent = init_agent(am)
    set_parameters!(agent, gen_params)
    reset!(agent)
    responses = Float64[]
    for inp in intensities
        dist   = agent.action_model(agent.model_attributes, inp)
        action = Float64(rand(rng, dist))
        push!(responses, action)
        set_actions!(agent, :action, action)
    end
    return responses
end

################################################################################
# ─── LOG-LIKELIHOOD FUNCTION ──────────────────────────────────────────────────
#
# Evaluates log p(responses | θ_c, model) via the HGF forward pass.
# θ_c is a NamedTuple of constrained parameter values (same space as the
# MCMCChains output — no bijector transform needed here).
#
# Mathematically identical to what DynamicPPL.pointwise_loglikelihoods
# computes on this ActionModels model: in the NoMissingActions path,
# ActionModels batches all N_TRIALS responses into a single
#   actions ~ arraydist(distributions)
# so pointwise_loglikelihoods returns one total log-likelihood per
# posterior sample, which is the sum of per-trial log p(response | dist_t).
# The forward pass below computes that same sum directly.
################################################################################

"""
    make_log_lik_fn(hgf_string, intensities, responses)

Return a closure θ_c::NamedTuple → Float64 that computes total log-likelihood
log p(responses | θ_c) via the HGF forward pass.  Returns -Inf at the first
trial with a non-finite log-probability (matches ActionModels behaviour).
"""
function make_log_lik_fn(hgf_string::String,
                          intensities::Vector{Float64},
                          responses::Vector{Float64})
    am    = create_agent(hgf_string)
    agent = init_agent(am)

    return function log_lik_constrained(θ_c::NamedTuple)
        set_parameters!(agent, θ_c)
        reset!(agent)
        ll = 0.0
        for (inp, resp) in zip(intensities, responses)
            dist = agent.action_model(agent.model_attributes, inp)
            lp   = logpdf(dist, resp)
            isfinite(lp) || return -Inf
            ll += lp
            set_actions!(agent, :action, resp)
        end
        return ll
    end
end

################################################################################
# ─── PATHS AND SKIP CHECK ─────────────────────────────────────────────────────
################################################################################

results_dir = joinpath("param_recovery", "prior_based_mcmc", "results")
chains_dir  = joinpath("param_recovery", "prior_based_mcmc", "chains")
out_file    = joinpath(results_dir, "sim$(sim_index)_$(true_model).csv")

println("=" ^ 70)
println("VCH Prior-Based Recovery — AIC/BIC patch")
println("  sim_index  = $sim_index")
println("  true_model = $true_model")
println("=" ^ 70)

isfile(out_file) ||
    error("Result CSV not found: $out_file — run prior_recovery_vch_mcmc.jl first")

existing = CSV.read(out_file, DataFrame)

# Skip if all AIC/BIC columns already present and non-missing
if "aic_winner" in names(existing) && !ismissing(existing[1, :aic_winner])
    println("Already done — skipping ($out_file already has aic_winner)")
    exit(0)
end

################################################################################
# ─── RNG REPLAY ───────────────────────────────────────────────────────────────
#
# CRITICAL: must reproduce the EXACT same RNG sequence as
# prior_recovery_vch_mcmc.jl.  The sequence is:
#   1. MersenneTwister(sim_index)
#   2. sample_from_prior(rng, true_model)   — draws gen_params
#   3. generate_vch_conditions(rng)          — draws cond_sequence
#   4. simulate_responses(rng, ...)          — draws responses
# Any deviation = wrong participant's data → wrong log-likelihood.
################################################################################

println("\n── RNG replay ──────────────────────────────────────────────────────")
rng = MersenneTwister(sim_index)

gen_params    = sample_from_prior(rng, true_model)
cond_sequence = generate_vch_conditions(rng)

true_intensities = conditions_to_intensities(cond_sequence, true_model)
true_hgf         = get_hgf_string(true_model)
responses        = simulate_responses(rng, true_hgf, gen_params, true_intensities)

println("  Replayed $(length(responses)) responses  " *
        "($(sum(responses .== 1.0)) detect, $(sum(responses .== 0.0)) miss)")
flush(stdout)

################################################################################
# ─── CROSS-CHECK: replayed gen_params vs saved CSV ───────────────────────────
#
# Verifies that the RNG replay is faithful to the original fitting run.
# Mismatches beyond TOL indicate a change in the RNG call sequence or the
# prior/task constants — both would invalidate the log-likelihood computation.
################################################################################

const REPLAY_TOL = 1e-8

function assert_replay(param_name::String, replayed::Float64, saved_val)
    saved = Float64(saved_val)   # handle CSV integer reads
    if isnan(replayed) && isnan(saved)
        return  # both NaN is correct (omega3 for 2-level true models)
    end
    if isnan(replayed) != isnan(saved)
        error("RNG replay MISMATCH: $param_name  replayed=$replayed  saved=$saved" *
              " (one is NaN — wrong model type or prior change?)")
    end
    if abs(replayed - saved) > REPLAY_TOL
        error("RNG replay MISMATCH: $param_name  replayed=$replayed  saved=$saved" *
              "  diff=$(abs(replayed-saved))  (> tol=$REPLAY_TOL)\n" *
              "Cause: changed RNG sequence, prior, or task constants.")
    end
end

assert_replay("nu",     gen_params.prior_posterior_weight, existing[1, :gen_nu])
assert_replay("beta",   gen_params.action_precision,       existing[1, :gen_beta])
assert_replay("omega",  gen_params.xprob_volatility,       existing[1, :gen_omega])
assert_replay("omega3",
    hasproperty(gen_params, :xvol_volatility) ? gen_params.xvol_volatility : NaN,
    existing[1, :gen_omega3])

println("  gen_params cross-check PASSED (all within $(REPLAY_TOL))")
flush(stdout)

################################################################################
# ─── COMPUTE AIC/BIC FOR EACH FITTED MODEL ────────────────────────────────────
################################################################################

println("\n── AIC/BIC computation ─────────────────────────────────────────────")

aic_vals          = Dict{String, Float64}()
bic_vals          = Dict{String, Float64}()
log_lik_max_vals  = Dict{String, Float64}()
log_lik_mean_vals = Dict{String, Float64}()
log_lik_sd_vals   = Dict{String, Float64}()
n_samples_vals    = Dict{String, Int}()

for fitted_model in FITTED_MODELS
    println("\n  Fitted model: $fitted_model")
    flush(stdout)

    chains_file = joinpath(chains_dir,
                           "sim$(sim_index)_$(true_model)_$(fitted_model).jls")

    # ── Guard: missing chains file ─────────────────────────────────────────
    if !isfile(chains_file)
        @warn "  Chains file not found: $chains_file — filling with NaN"
        aic_vals[fitted_model] = bic_vals[fitted_model] = NaN
        log_lik_max_vals[fitted_model] = log_lik_mean_vals[fitted_model] = NaN
        log_lik_sd_vals[fitted_model]  = NaN
        n_samples_vals[fitted_model]   = 0
        continue
    end

    # ── Load chains and extract constrained parameter samples ──────────────
    chains    = deserialize(chains_file)
    chains_df = DataFrame(chains)
    n_total   = nrow(chains_df)
    println("  Loaded $n_total samples from chains")
    flush(stdout)

    # Column names in chains_df for constrained parameter values.
    # ActionModels/Turing stores them as "param_name.session[1]".
    fit_priors    = get_priors(fitted_model)
    param_keys_fit = keys(fit_priors)
    session_cols   = [Symbol("$(p).session[1]") for p in param_keys_fit]

    # Verify expected columns are present — catches schema drift
    for col in session_cols
        if !(string(col) in names(chains_df))
            error("  Expected chain column '$col' not found in $(basename(chains_file)).\n" *
                  "  Available: $(names(chains_df))")
        end
    end

    # ── Fitted model's intensity mapping ──────────────────────────────────
    # Each fitted model uses its OWN condition→intensity convention.
    # The responses are fixed (simulated from the true model), but the
    # inputs seen by the fitted HGF differ by mapping type.
    fit_intensities = conditions_to_intensities(cond_sequence, fitted_model)
    fit_hgf         = get_hgf_string(fitted_model)

    # ── Build log-likelihood function ──────────────────────────────────────
    log_lik_fn = make_log_lik_fn(fit_hgf, fit_intensities, responses)

    # ── Evaluate log-likelihood at each posterior sample ───────────────────
    log_liks = Vector{Float64}(undef, n_total)
    for i in 1:n_total
        θ_c = NamedTuple{param_keys_fit}(
            chains_df[i, col] for col in session_cols
        )
        log_liks[i] = try
            log_lik_fn(θ_c)
        catch
            -Inf
        end
    end

    # ── Filter non-finite evaluations (should be rare) ────────────────────
    finite_mask = isfinite.(log_liks)
    n_finite    = sum(finite_mask)
    if n_finite < n_total
        @warn "  Dropped $(n_total - n_finite) non-finite log-lik evaluations " *
              "for $fitted_model ($(round(100*(n_total-n_finite)/n_total, digits=1))%)"
    end

    # ── Guard: too few usable samples ──────────────────────────────────────
    if n_finite < 50
        @warn "  Too few finite samples ($n_finite) — filling with NaN"
        aic_vals[fitted_model] = bic_vals[fitted_model] = NaN
        log_lik_max_vals[fitted_model] = log_lik_mean_vals[fitted_model] = NaN
        log_lik_sd_vals[fitted_model]  = NaN
        n_samples_vals[fitted_model]   = n_finite
        continue
    end

    finite_liks = log_liks[finite_mask]

    # ── AIC / BIC ──────────────────────────────────────────────────────────
    k       = get_n_params(fitted_model)   # free parameters
    n       = N_TRIALS                      # observations
    log_L_hat = maximum(finite_liks)        # max log-lik over posterior samples
                                            # (approximation of MLE log-lik)
    aic = 2k       - 2 * log_L_hat
    bic = k * log(n) - 2 * log_L_hat

    aic_vals[fitted_model]          = aic
    bic_vals[fitted_model]          = bic
    log_lik_max_vals[fitted_model]  = log_L_hat
    log_lik_mean_vals[fitted_model] = mean(finite_liks)
    log_lik_sd_vals[fitted_model]   = std(finite_liks)
    n_samples_vals[fitted_model]    = n_finite

    println("    k=$(k)  n=$(n)  log_L_hat=$(round(log_L_hat, digits=3))" *
            "  AIC=$(round(aic, digits=3))  BIC=$(round(bic, digits=3))")
    flush(stdout)
end

################################################################################
# ─── DETERMINE WINNERS ────────────────────────────────────────────────────────
# AIC/BIC: LOWER value = better model fit (opposite of bridge LME).
################################################################################

finite_aic = [(m, aic_vals[m]) for m in FITTED_MODELS if isfinite(aic_vals[m])]
finite_bic = [(m, bic_vals[m]) for m in FITTED_MODELS if isfinite(bic_vals[m])]

aic_winner = isempty(finite_aic) ? "none" : finite_aic[argmin(last.(finite_aic))][1]
bic_winner = isempty(finite_bic) ? "none" : finite_bic[argmin(last.(finite_bic))][1]

println("\n── Winners ─────────────────────────────────────────────────────────")
println("  AIC winner: $aic_winner")
println("  BIC winner: $bic_winner")
flush(stdout)

################################################################################
# ─── UPDATE AND SAVE CSV ──────────────────────────────────────────────────────
# Appends new columns to the existing result CSV row in-place.
# Existing columns (bridge_lme_*, converged_*, mcmc_*) are preserved.
################################################################################

for fitted_model in FITTED_MODELS
    existing[!, "aic_$(fitted_model)"]           = [aic_vals[fitted_model]]
    existing[!, "bic_$(fitted_model)"]           = [bic_vals[fitted_model]]
    existing[!, "log_lik_max_$(fitted_model)"]   = [log_lik_max_vals[fitted_model]]
    existing[!, "log_lik_mean_$(fitted_model)"]  = [log_lik_mean_vals[fitted_model]]
    existing[!, "log_lik_sd_$(fitted_model)"]    = [log_lik_sd_vals[fitted_model]]
    existing[!, "n_samples_lik_$(fitted_model)"] = [n_samples_vals[fitted_model]]
end

existing[!, "aic_winner"] = [aic_winner]
existing[!, "bic_winner"] = [bic_winner]

CSV.write(out_file, existing)
println("\nSaved: $out_file")
println("Done — sim_index=$sim_index, true_model=$true_model  " *
        "AIC=$(aic_winner)  BIC=$(bic_winner)")
flush(stdout)
