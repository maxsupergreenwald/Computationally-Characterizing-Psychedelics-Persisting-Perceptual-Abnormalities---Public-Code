### bms_vch.jl  ─  Per-subject Laplace LME for RFX-BMS (N-model)
###
### For each participant with a staged stimulus CSV, computes the MAP estimate
### and Laplace-approximate log marginal evidence (LME) for every model in
### MODEL_TYPES.  Output is used by bms_vch.py for RFX-BMS.
###
### Laplace LME approximation (Stephan et al. 2009, NeuroImage):
###   LME ≈ log p(y,θ|M) + (k/2)·log(2π) − ½·log|H|
###   where H = ForwardDiff.hessian(map_result.f, θ_MAP)  [Hessian of neg-log-joint]
###         k = number of model parameters
###
### MAP warm-start strategy (reduces Laplace failures):
###   1. If MCMC medians file exists for this subject × model: initialize MAP
###      from posterior median (bijected to unconstrained space).
###   2. If median-init fails or no file: try up to MAX_RESTARTS random-restart MAPs.
###   3. Warn if all attempts fail — subject excluded from BMS for that model.
###
### Usage (run locally — MAP is fast, no MCMC):
###   julia --project=. --threads 1 bms_vch.jl
###
### Output: model_comparison/bms/results/<id>_lme.csv
###   Columns: record_id, lme_<model>, converged_<model>  for each model in MODEL_TYPES

################################################################################
# ─── ROLE IN THE MANUSCRIPT ───────────────────────────────────────────────────
#
# MODEL SELECTION, part 1 of 2.  Produces Supplementary Fig. S3c.
#
#     THIS SCRIPT (per-subject evidence) ──▶ bms_vch.py (group-level RFX-BMS)
#
# Why this exists
#     Two modelling choices had to be adjudicated rather than assumed: how deep
#     the hierarchy should be (2-level vs 3-level), and which stimulus convention
#     to feed the HGF (empiric vs nominal detection probabilities).  Crossing
#     those gives the four candidate models compared here.
#
#     This script supplies the per-subject evidence term.  It uses MAP + a Laplace
#     approximation rather than MCMC because model evidence is needed for every
#     subject × every candidate model, and the Laplace route is cheap enough to
#     run locally.  A subject whose Hessian is not positive definite for ANY model
#     in the set is dropped from BMS — that is a limitation of the Laplace
#     approximation at that subject's MAP, and says nothing about the validity of
#     their MCMC estimates, which are retained everywhere else.
#
# Reproducibility note
#     No random seed is set for the random-restart MAP fallback (MAX_RESTARTS,
#     below) or for Julia's default RNG generally. Re-running this script can
#     therefore land on different per-subject MAP estimates than a prior run,
#     which propagates into slightly different Ef/PXP/BOR in bms_vch.py's
#     RFX-BMS output. Set a seed before re-running if you need the result to
#     match a specific prior run exactly.
################################################################################

using ActionModels, HierarchicalGaussianFiltering
using CSV, Glob
using Turing
using Bijectors
using DataFrames
using Distributions
using ForwardDiff
using LinearAlgebra

include("helper_functions/create_agent.jl")

################################################################################
# ─── ACTIVE CONFIG ────────────────────────────────────────────────────────────
# Edit MODEL_TYPES to select which models to compare.
#
# Supported model type strings:
#   "2level_empiric"             — 2-level HGF, corrected empirical mapping (0% → 0.0)
#   "2level_nominal"             — 2-level HGF, raw condition proportions {0,0.25,0.5,0.75}
#   "3level_empiric"             — 3-level HGF, corrected empirical mapping
#   "3level_nominal"             — 3-level HGF, raw condition proportions
#
################################################################################

const MODEL_TYPES = ["2level_empiric", "3level_empiric", "2level_nominal", "3level_nominal"]

# Local base directory for MCMC result files (used for MAP warm-starting).
# Relative to 02_hgf_modeling/ (where this script runs).
# Structure: RESULTS_BASE_DIR/{model_type}/{timepoint}/{modality}/{id}_medians.csv
# Set to nothing to disable median warm-starting (random-restart fallback only).
const RESULTS_BASE_DIR = joinpath("results", "vch")

# Maximum random-restart MAP attempts after median-initialized MAP fails.
const MAX_RESTARTS = 3

################################################################################
# ─── PRIORS ───────────────────────────────────────────────────────────────────
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

# For model types without dedicated MCMC runs, borrow medians from this
# model type when warm-starting the MAP.
const MEDIAN_FALLBACK = Dict(
    "2level_nominal"             => "2level_empiric",   # borrow 2-level MCMC medians (same architecture)
    "3level_nominal"             => "3level_empiric",   # borrow 3-level MCMC medians (same architecture)
)

################################################################################
# ─── PATHS ────────────────────────────────────────────────────────────────────
################################################################################

const PROJECT   = "hppd_manuscript"
const TIMEPOINT = "hppd"
const MODALITY  = "vch"

################################################################################
# ─── HELPERS ──────────────────────────────────────────────────────────────────
################################################################################

"""Return (hgf_string, priors) for a given model_type."""
function model_config(model_type::String)
    if startswith(model_type, "2level")
        hgf = "binary_2level";  priors = PRIORS_2LEVEL
    elseif startswith(model_type, "3level")
        hgf = "binary_3level";  priors = PRIORS_3LEVEL
    else
        error("Unknown model type: $model_type")
    end
    return hgf, priors
end

"""Return the stimulus vector for a given model_type and data DataFrame."""
function get_stimulus(df::DataFrame, model_type::String)
    if endswith(model_type, "_nominal")
        return collect(Float64, df[!, :condition])
    else
        return collect(Float64, df[!, :empirical_condition])
    end
end

"""
    load_medians(record_id, model_type, priors) → NamedTuple or nothing

Load MCMC posterior medians from disk for warm-starting the MAP optimizer.
Returns a NamedTuple with the same keys as `priors`, or nothing if the file
is missing or unreadable.  Uses MEDIAN_FALLBACK for the nominal variants.
"""
function load_medians(record_id::Int, model_type::String, priors)
    RESULTS_BASE_DIR === nothing && return nothing
    mcmc_type = get(MEDIAN_FALLBACK, model_type, model_type)
    path = joinpath(RESULTS_BASE_DIR, mcmc_type, TIMEPOINT, MODALITY,
                    "$(record_id)_medians.csv")
    isfile(path) || return nothing
    try
        df          = CSV.read(path, DataFrame)
        param_names = keys(priors)
        vals        = [Float64(df[1, string(n)]) for n in param_names]
        return NamedTuple{Tuple(param_names)}(vals)
    catch
        return nothing
    end
end

"""
    medians_to_unconstrained(medians_nt, priors) → Vector{Float64}

Bijected-transform constrained posterior medians to the unconstrained space
that Turing's MAP optimizer operates in.  Uses Bijectors.jl to apply the
same transformation Turing uses internally for each truncated-Normal prior.
Values are clamped 1e-8 inside each bound before transformation to avoid
log(0) at exact boundaries.
"""
function medians_to_unconstrained(medians_nt, priors)
    θ_u = Float64[]
    for n in keys(priors)
        dist = priors[n]
        val  = Float64(medians_nt[n])
        lo, hi = minimum(dist), maximum(dist)
        isfinite(lo) && (val = max(val, lo + 1e-8))
        isfinite(hi) && (val = min(val, hi - 1e-8))
        push!(θ_u, bijector(dist)(val))
    end
    return θ_u
end

"""
    _try_lme(modelfit, k; initial_params=nothing) → (lme, converged)

Single MAP attempt.  Returns (lme::Float64, true) if the Hessian is PD,
(NaN, false) otherwise.  `initial_params` is a vector in unconstrained space.
"""
function _try_lme(modelfit, k; initial_params=nothing)
    if initial_params === nothing
        map_result = maximum_a_posteriori(modelfit.model, adtype=AutoForwardDiff())
    else
        map_result = maximum_a_posteriori(modelfit.model, adtype=AutoForwardDiff(),
                                          initial_params=initial_params)
    end
    θ_u = map_result.optim_result.u
    lp  = -map_result.f(θ_u)
    H   = ForwardDiff.hessian(map_result.f, θ_u)
    if !all(isfinite, H) || !isposdef(Symmetric(H))
        return NaN, false
    end
    lme = lp + (k / 2) * log(2π) - 0.5 * log(det(Symmetric(H)))
    return lme, true
end

"""
    compute_lme(priors, hgf_string, stimuli, responses, record_id, model_type)

Returns (lme::Float64, converged::Bool).

Attempt order:
  1. MAP warm-started from MCMC posterior medians (if file exists).
  2. Up to MAX_RESTARTS random-restart MAPs.
  3. Return (NaN, false) and warn if all attempts fail.
"""
function compute_lme(priors, hgf_string, stimuli, responses,
                     record_id::Int, model_type::String)
    am       = create_agent(hgf_string)
    modelfit = create_model(am, priors, stimuli, responses;
                            check_parameter_rejections=true, verbose=false)
    k = length(priors)

    # ── Attempt 1: warm-start from MCMC posterior medians ─────────────────────
    medians = load_medians(record_id, model_type, priors)
    if medians !== nothing
        try
            θ0       = medians_to_unconstrained(medians, priors)
            lme, conv = _try_lme(modelfit, k; initial_params=θ0)
            conv && return lme, true
        catch e
            @warn "record_id=$record_id model=$model_type: median-init MAP error — $e"
        end
    end

    # ── Fallback: random restarts ──────────────────────────────────────────────
    for attempt in 1:MAX_RESTARTS
        try
            lme, conv = _try_lme(modelfit, k)
            if conv
                note = medians === nothing ? "no medians file" : "median-init non-PD"
                @warn "record_id=$record_id model=$model_type: converged on random restart #$attempt ($note)"
                return lme, true
            end
        catch e
            @warn "record_id=$record_id model=$model_type: random restart #$attempt error — $e"
        end
    end

    n_total = 1 + MAX_RESTARTS   # 1 median-init + MAX_RESTARTS random
    @warn "record_id=$record_id model=$model_type: ALL $n_total MAP attempts failed — excluded from BMS"
    return NaN, false
end

################################################################################
# ─── MAIN ─────────────────────────────────────────────────────────────────────
################################################################################

stim_dir = joinpath("data_n_cmnds", PROJECT, TIMEPOINT, "$(MODALITY)_data")
out_dir  = joinpath("model_comparison", "bms", "results")
mkpath(out_dir)

stim_files = sort(glob("*.csv", stim_dir))
isempty(stim_files) && error("No stimulus CSVs found in $stim_dir")

_model_list = join(MODEL_TYPES, ", ")
println("=" ^ 70)
println("VCH BMS — Laplace LME  ($(length(MODEL_TYPES)) models: $(_model_list))")
println("  $(length(stim_files)) participants")
println("  MCMC warm-start: $(RESULTS_BASE_DIR === nothing ? "disabled" : RESULTS_BASE_DIR)")
println("  Random restarts: $MAX_RESTARTS")
println("=" ^ 70)

n_done = 0
n_skip = 0

for stim_file in stim_files
    record_id = parse(Int, splitext(basename(stim_file))[1])
    out_file  = joinpath(out_dir, "$(record_id)_lme.csv")

    # Determine which models still need computing.
    # If the output file exists but is missing columns for some models (e.g.
    # because MODEL_TYPES was expanded after the original run), read the existing
    # values and compute only the missing models — avoids redundant MAP calls.
    existing_row  = Dict{String, Any}()
    models_to_run = MODEL_TYPES

    if isfile(out_file)
        existing      = CSV.read(out_file, DataFrame)
        models_to_run = [mt for mt in MODEL_TYPES if "lme_$(mt)" ∉ names(existing)]
        if isempty(models_to_run)
            global n_skip += 1
            continue
        end
        # Carry over all existing column values so the merged row is complete
        for col in names(existing)
            existing_row[col] = existing[1, col]
        end
        _models_str = join(models_to_run, ", ")
        println("  →  $record_id  (partial: computing $(_models_str))")
    end

    df        = CSV.read(stim_file, DataFrame)
    responses = collect(Float64, df[!, :response])

    lme_vals  = Dict{String, Float64}()
    conv_vals = Dict{String, Bool}()

    for mt in models_to_run
        hgf_str, priors = model_config(mt)
        stimuli         = get_stimulus(df, mt)
        lme_v, conv_v   = NaN, false
        try
            lme_v, conv_v = compute_lme(priors, hgf_str, stimuli, responses,
                                        record_id, mt)
        catch e
            @warn "record_id=$record_id  model=$mt  unexpected error: $e"
        end
        lme_vals[mt]  = lme_v
        conv_vals[mt] = conv_v
    end

    # Build output row: start from existing values, overwrite/add new model columns
    row = Dict{String, Any}(existing_row)
    row["record_id"] = record_id
    for mt in models_to_run
        row["lme_$(mt)"]       = lme_vals[mt]
        row["converged_$(mt)"] = conv_vals[mt]
    end
    CSV.write(out_file, DataFrame([row]))

    status   = all(conv_vals[mt] for mt in models_to_run) ? "✓" : "!"
    lme_strs = join(["$(mt)=$(isnan(lme_vals[mt]) ? "NaN" : round(lme_vals[mt], digits=1))"
                     for mt in models_to_run], "  ")
    println("  $status  $record_id  $lme_strs")
    global n_done += 1
end

println("\n" * "=" ^ 70)
println("Done: $n_done computed, $n_skip skipped")
println("Results in: $out_dir")
