### ppc_classic_vch.jl  ─  VCH HGF Classic PPC (local, deterministic medians)
###
### For each participant in the provided medians CSV, this script:
###   1. Reads the participant's real stimulus + response data (staged CSV)
###   2. Simulates ONE synthetic response sequence at the participant's
###      posterior median parameters — no MAP, no Hessian, no covariance matrix
###   3. Computes per-condition and per-block detection rates for both the
###      empirical and simulated sequences, and saves them to a 2-row CSV
###
### Usage:
###   julia --project=. --threads 1 ppc_classic_vch.jl <model_type> <medians_csv>
###
###   model_type  — "2level_empiric", "2level_nominal", or "3level_empiric"
###   medians_csv — CSV path with columns:
###                   record_id, action_precision, prior_posterior_weight,
###                   xprob_volatility [, xvol_volatility]
###
### Output per participant:
###   param_recovery/ppc_classic/results/<model_type>/<id>_ppc_classic.csv
###   Columns: record_id, source,
###            det_rate_0.0, det_rate_0.25, det_rate_0.5, det_rate_0.75,
###            block_1 … block_12,
###            block_1_cond0 … block_12_cond0     (0%-condition trials only)
###            block_1_cond25 … block_12_cond25   (25%-condition trials only)
###            block_1_cond50 … block_12_cond50   (50%-condition trials only)
###            block_1_cond75 … block_12_cond75   (75%-condition trials only)
###            NaN for any block with zero trials of that condition

################################################################################
# ─── ROLE IN THE MANUSCRIPT ───────────────────────────────────────────────────
#
# POSTERIOR PREDICTIVE CHECKS, part 1 of 2.  Produces Supplementary Fig. S3f and
# the two PPC panels of Supplementary Fig. S2b.
#
#     THIS SCRIPT (simulate at medians) ──▶ ppc_classic_vch.py (compile + plot)
#
# What a PPC is being asked to establish here
#     Recovery analyses show the model is identifiable; they do not show it
#     resembles the data.  The manuscript set two qualitative criteria in advance:
#     simulated detection probability should (1) increase with stimulus strength
#     and (2) decay across the session.  Both were met.
#
# Why "classic"
#     One response sequence is simulated per participant, at their posterior
#     median parameters — no MAP, no Hessian, no sampling from a posterior
#     covariance.  The uncertainty shown in the figures comes from bootstrapping
#     across participants, not from within-participant parameter uncertainty.
#     That is a deliberately conservative, transparent choice, and it is the
#     variant the manuscript describes.
#
# Both stimulus conventions are run
#     model_type is a required argument, and "2level_nominal" is supported
#     alongside "2level_empiric" with its own entry in VCH_STIM_MAPPINGS.  The
#     empiric/nominal contrast in Supplementary Fig. S2b is exactly these two runs
#     placed side by side — it is the visual evidence that the nominal convention
#     departs from observed detection behaviour.
################################################################################

using ActionModels, HierarchicalGaussianFiltering
using CSV, Glob
using DataFrames
using Random
using Statistics: mean

include("helper_functions/create_agent.jl")

################################################################################
# ─── PARSE ARGUMENTS ──────────────────────────────────────────────────────────
################################################################################

if length(ARGS) < 2
    error("Usage: julia ppc_classic_vch.jl <model_type> <medians_csv>")
end

model_type  = ARGS[1]
medians_csv = ARGS[2]

if model_type ∉ ("2level_empiric", "2level_nominal", "3level_empiric")
    error("model_type must be one of: \"2level_empiric\", " *
          "\"2level_nominal\", \"3level_empiric\"; got: \"$model_type\"")
end

################################################################################
# ─── STIMULUS MAPPINGS ────────────────────────────────────────────────────────
# Same values as VCH_STIM_MAPPINGS in model_fitting_singleagent_forarray.jl and
# ppc_vch.jl.  Condition labels (0.0/0.25/0.5/0.75) come from :condition column
# in staged CSV; these maps convert them to signal intensities fed to the HGF.
################################################################################

const VCH_STIM_MAPPINGS = Dict(
    "2level_empiric" => Dict(
        0.0  => 0.0,
        0.25 => 0.4180444024563061,
        0.50 => 0.7115104419621175,
        0.75 => 0.8994252873563219,
    ),
    "2level_nominal" => Dict(
        0.0  => 0.0,
        0.25 => 0.25,
        0.50 => 0.50,
        0.75 => 0.75,
    ),
    "3level_empiric" => Dict(
        0.0  => 0.0,
        0.25 => 0.4180444024563061,
        0.50 => 0.7115104419621175,
        0.75 => 0.8994252873563219,
    ),
)
const STIM_MAP = VCH_STIM_MAPPINGS[model_type]

################################################################################
# ─── CONSTANTS ────────────────────────────────────────────────────────────────
################################################################################

const PROJECT    = "hppd_manuscript"
const TIMEPOINT  = "hppd"
const MODALITY   = "vch"
const N_BLOCKS   = 12
const BLOCK_SIZE = 30
const COND_LEVELS = [0.0, 0.25, 0.5, 0.75]

################################################################################
# ─── SUMMARY STATISTICS ───────────────────────────────────────────────────────
# Identical to compute_stats() in ppc_vch.jl.
################################################################################

"""
Returns (cond_rates, block_rates, block_per_cond_rates) where:
  cond_rates         — Vector{Float64} length 4 (one per COND_LEVELS)
  block_rates        — Vector{Float64} length N_BLOCKS (all conditions pooled)
  block_per_cond_rates — Dict(cond_level => Vector{Float64} length N_BLOCKS)
NaN for any block×condition cell with zero matching trials.
"""
function compute_stats(condition_labels::Vector{Float64},
                       responses::Vector{Float64})
    cond_rates = Float64[]
    for c in COND_LEVELS
        idx = findall(==(c), condition_labels)
        push!(cond_rates, isempty(idx) ? NaN : mean(responses[idx]))
    end

    block_rates = Float64[]
    block_per_cond_rates = Dict(c => Float64[] for c in COND_LEVELS)
    for b in 0:(N_BLOCKS - 1)
        block_start = b * BLOCK_SIZE + 1
        block_end   = (b + 1) * BLOCK_SIZE
        block_cond  = condition_labels[block_start:block_end]
        block_resp  = responses[block_start:block_end]
        push!(block_rates, mean(block_resp))
        for c in COND_LEVELS
            mask = block_cond .== c
            push!(block_per_cond_rates[c], any(mask) ? mean(block_resp[mask]) : NaN)
        end
    end

    return cond_rates, block_rates, block_per_cond_rates
end

################################################################################
# ─── PATHS ────────────────────────────────────────────────────────────────────
################################################################################

stim_dir   = joinpath("data_n_cmnds", PROJECT, TIMEPOINT, "$(MODALITY)_data")
out_dir    = joinpath("param_recovery", "ppc_classic", "results", model_type)
mkpath(out_dir)

hgf_string = model_type == "3level_empiric" ? "binary_3level" : "binary_2level"

################################################################################
# ─── LOAD MEDIANS ─────────────────────────────────────────────────────────────
################################################################################

medians_df = CSV.read(medians_csv, DataFrame)

println("=" ^ 60)
println("VCH PPC Classic  —  model=$model_type")
println("  $(nrow(medians_df)) participants from: $medians_csv")
println("  Output directory: $out_dir")
println("=" ^ 60)
flush(stdout)

################################################################################
# ─── MAIN LOOP ────────────────────────────────────────────────────────────────
################################################################################

n_done = 0
n_fail = 0

for mrow in eachrow(medians_df)
    record_id = Int(mrow[:record_id])
    out_file  = joinpath(out_dir, "$(record_id)_ppc_classic.csv")

    stim_file = joinpath(stim_dir, "$(record_id).csv")
    if !isfile(stim_file)
        @warn "record_id=$record_id: stimulus CSV not found at $stim_file — skipping"
        global n_fail += 1
        println("  ✗ $record_id  (stimulus CSV missing)")
        flush(stdout)
        continue
    end

    df             = CSV.read(stim_file, DataFrame)
    cond_labels    = collect(Float64, df[!, :condition])
    emp_conditions = Float64[STIM_MAP[c] for c in cond_labels]
    emp_responses  = collect(Float64, df[!, :response])

    try
        # ── Build median parameter NamedTuple ──────────────────────────────────
        ap  = Float64(mrow[:action_precision])
        ppw = Float64(mrow[:prior_posterior_weight])
        xpv = Float64(mrow[:xprob_volatility])

        if model_type == "3level_empiric"
            xvv    = Float64(mrow[:xvol_volatility])
            params = (; action_precision=ap, prior_posterior_weight=ppw,
                        xprob_volatility=xpv, xvol_volatility=xvv)
        else
            params = (; action_precision=ap, prior_posterior_weight=ppw,
                        xprob_volatility=xpv)
        end

        # ── One forward pass at median parameters ──────────────────────────────
        am    = create_agent(hgf_string)
        agent = init_agent(am)
        set_parameters!(agent, params)
        reset!(agent)

        rng           = MersenneTwister(record_id)
        sim_responses = Float64[]
        for inp in emp_conditions
            d      = agent.action_model(agent.model_attributes, inp)
            action = Float64(rand(rng, d))
            push!(sim_responses, action)
            set_actions!(agent, :action, action)
        end

        # ── Compute stats ──────────────────────────────────────────────────────
        emp_cond_rates, emp_block_rates, emp_block_per_cond =
            compute_stats(cond_labels, emp_responses)
        sim_cond_rates, sim_block_rates, sim_block_per_cond =
            compute_stats(cond_labels, sim_responses)

        # ── Build output rows ──────────────────────────────────────────────────
        rows = []
        for (source, cond_rates, block_rates, block_per_cond) in (
            ("empirical", emp_cond_rates, emp_block_rates, emp_block_per_cond),
            ("sim",       sim_cond_rates, sim_block_rates, sim_block_per_cond),
        )
            row = Dict{String,Any}(
                "record_id" => record_id,
                "source"    => source,
            )
            for (i, c) in enumerate(COND_LEVELS)
                row["det_rate_$(c)"] = cond_rates[i]
            end
            for b in 1:N_BLOCKS
                row["block_$(b)"]        = block_rates[b]
                row["block_$(b)_cond0"]  = block_per_cond[0.0][b]
                row["block_$(b)_cond25"] = block_per_cond[0.25][b]
                row["block_$(b)_cond50"] = block_per_cond[0.5][b]
                row["block_$(b)_cond75"] = block_per_cond[0.75][b]
            end
            push!(rows, row)
        end

        CSV.write(out_file, DataFrame(rows))
        global n_done += 1
        println("  ✓ $record_id")
        flush(stdout)

    catch e
        @warn "record_id=$record_id failed: $e"
        global n_fail += 1
        println("  ✗ $record_id  ($e)")
        flush(stdout)
    end
end

println("\n" * "=" ^ 60)
println("Done: $n_done processed, $n_fail failed")
println("Results in: $out_dir")
flush(stdout)
