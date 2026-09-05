### model_fitting_singleagent_forarray.jl
###
### Fit one participant's binary Hierarchical Gaussian Filter (HGF) to their
### trial-by-trial CH-task detection responses via MCMC (NUTS), then replay the
### fitted model to extract trial-by-trial belief trajectories.
###
### This is the computational core of the HGF pipeline.  It is written to be
### invoked once per participant as a single task in a SLURM/dSQ job array:
### `hgf_pipeline.py` stages one CSV per participant and emits a job file with
### one line per participant, each line calling this script with a different
### <file_index>.
###
### ── Position in the pipeline ─────────────────────────────────────────────────
###
###   hgf_pipeline.py                (local)  stage per-subject CSVs + job array
###        │
###        ▼
###   THIS SCRIPT                    (HPC)    fit HGF, save chains + states
###        │
###        ▼
###   import_hgf_results_unified.py  (local)  merge posteriors + states into dfs
###
### ── Usage ────────────────────────────────────────────────────────────────────
###
###   module load Julia/1.11.4-linux-x86_64
###   julia --project=. --threads 4 model_fitting_singleagent_forarray.jl \
###       <file_index> <modality> <skip_or_nah> <timepoint> <project> \
###       <model_type> [n_iterations]
###
###   Run `Pkg.instantiate()` once in this directory before the first job so the
###   Project.toml/Manifest.toml environment is materialised:
###       julia --project=. -e "using Pkg; Pkg.instantiate()"
###
### ── Reads ────────────────────────────────────────────────────────────────────
###
###   data_n_cmnds/<project>/<timepoint>/<modality>_data/*.csv
###     One CSV per participant, written by hgf_pipeline.py.  Required columns:
###       subj_id              participant record_id (constant within file)
###       condition            nominal contrast level {0.0, 0.25, 0.5, 0.75}
###       empirical_condition  stimulus intensity (RECOMPUTED at runtime for VCH)
###       response             binary detection response {0, 1}
###
###   The file is selected by POSITION, not by name: `file_index` is a 1-based
###   index into Glob.glob("*.csv", ...), which sorts LEXICOGRAPHICALLY.
###   "100.csv" therefore sorts before "99.csv".  Never assume index == subject.
###
### ── Writes ───────────────────────────────────────────────────────────────────
###
###   results/<model_type>/<project>/<timepoint>/<modality>/
###       <subj_id>.jls                      full MCMCChains object (serialized)
###       <subj_id>.csv                      all posterior samples, flat CSV
###       <subj_id>_medians.csv              posterior median per parameter
###       <subj_id>_rhats.csv                Gelman–Rubin R-hat per parameter
###       <subj_id>state_trajectories.csv    per-trial HGF states at the medians
###
###   Downstream QC (import_hgf_results_unified.py) sets any parameter whose
###   R-hat falls outside 0.9–1.1 to NaN in the analysis dataframe.
###
### ── Things most likely to need changing ──────────────────────────────────────
###
###   VCH_CORRECTED_MAPPING
###       The condition → stimulus-intensity lookup.  Hardcoded here
###       deliberately (see §4) and must be kept in sync with the copies in
###       hgf_pipeline.py, ppc_classic_vch.jl, and prior_recovery_vch_mcmc.jl.
###   n_iterations / n_chains / my_sampler   MCMC budget and sampler (§3).
###   priors                                  Parameter priors (§5).
###   desired_states                          Which HGF states get written (§6).

################################################################################
# ─── 1. PACKAGES AND SHARED MODEL DEFINITIONS ──────────────────────────────
# create_agent.jl defines the ActionModel: the HGF node structure, the
# belief-formation rule that blends the HGF prediction with the current
# stimulus intensity via nu, and the unit-square-sigmoid response model
# governed by beta.  Both the 2-level and 3-level architectures live there.
################################################################################

using ActionModels, HierarchicalGaussianFiltering
using CSV, Serialization
using Turing
using StatsPlots
using ProgressMeter
using DataFrames
using Glob
using MCMCChains

include("helper_functions/create_agent.jl")

################################################################################
# ─── 2. COMMAND-LINE ARGUMENTS ─────────────────────────────────────────────
# Six required positional arguments, plus an optional seventh.
#
#   file_index   1-based index into the lexicographically sorted staging
#                CSV glob.  Supplied by the dSQ job array, one per line.
#   modality     "vch" (visual) or "ach" (auditory).  Only VCH gets its
#                empirical_condition recomputed at runtime (see §7).
#   skip_or_nah  "skip"  → leave already-complete subjects untouched,
#                anything else → refit and overwrite.  "skip" makes a
#                partially-failed array safe to resubmit wholesale.
#   timepoint    Study timepoint label; "hppd" for this cross-sectional
#                manuscript.  Becomes a results path component.
#   project      Project key; "hppd_manuscript" here.  Also a path
#                component, keeping projects' results fully separate.
#   model_type   Which HGF variant to fit — selects BOTH the hierarchy
#                depth (§5) and the stimulus mapping (§4).
#   n_iterations Optional post-warmup MCMC samples per chain; default
#                1000.  Raised (e.g. 1500) for catch-up reruns of
#                subjects that failed the R-hat gate on the first pass.
################################################################################

# Read command line argument for file index
if length(ARGS) < 6
    error("Usage: julia model_fitting_singleagent_forarray.jl <file_index> <modality> <skip_or_nah> <timepoint> <project> <model_type>")
end

file_index = parse(Int, ARGS[1])
modality = ARGS[2]     # e.g., "vch" or "ach"
skip_or_nah = ARGS[3]
timepoint = ARGS[4]    # e.g., "hppd" (cross-sectional), "hyp", "acu", "sub", "pers"
project = ARGS[5]      # e.g., "hppd_manuscript", "aim2", "aim1_rpt"
model_type = ARGS[6]   # one of: "2level_empiric", "2level_nominal",
                       #          "3level_empiric", "3level_nominal"

################################################################################
# ─── 3. MCMC SETTINGS ──────────────────────────────────────────────────────
# Held constant across every participant and every model variant so that
# posterior estimates are comparable, and so that model comparison (BMS,
# BIC) is not confounded by differing sampling budgets.
#
#   use_optim = true  initialises each chain at the maximum a posteriori
#                     estimate rather than a prior draw.  This markedly
#                     shortens warmup and reduces the rate of chains that
#                     fail to converge within the iteration budget.
#   n_chains  = 4     four chains, so R-hat is well estimated.
#   NUTS with ForwardDiff — the HGF update equations are small, smooth,
#                     and forward-mode differentiable.
################################################################################

### SETTINGS ###
use_optim = true
n_iterations = length(ARGS) >= 7 ? parse(Int, ARGS[7]) : 1000
n_chains = 4
my_sampler = NUTS(adtype=AutoForwardDiff())

################################################################################
# ─── 4. VCH STIMULUS MAPPINGS  (condition → HGF input intensity) ───────────
# Each trial enters the HGF as a stimulus intensity in [0, 1] representing
# the ground-truth probability that a signal is present.  Two conventions
# exist, and the model_type argument selects between them:
#
#   VCH_CORRECTED_MAPPING  "_empiric"  — out-of-set normative hit rates
#   (no mapping)           "_nominal"  — the raw QUEST-targeted proportions
#                                        {0, 0.25, 0.5, 0.75}, read straight
#                                        from the condition column
#
# The manuscript's primary analysis uses "2level_empiric".  The remaining
# variants exist so that the choice of convention and hierarchy depth can
# be adjudicated by formal model comparison rather than assumed.
################################################################################

# ── VCH stimulus mappings (condition → empirical_condition) ───────────────────
# Hardcoded here so the HPC script is the single source of truth for which
# empirical values each model type uses.  After loading the staging CSV,
# empirical_condition is *recomputed* from the condition column at runtime —
# regardless of what the CSV column contains — preventing cross-contamination
# between pipeline runs with different mappings.
#
# 2026-06-17: VCH_CORRECTED_MAPPING updated to non-hallucinator empirical values
# for the 25 %, 50 %, and 75 % conditions (max_vh_freq == 0 AND max_ah_freq == 0,
# n=29 QC-passing; NaN treated as 0). Rationale: likelihood should reflect
# perceptual behavior of participants without hallucination histories.  See:
#   out_of_set_data/README_empirical_condition_redetermination.md
#   out_of_set_data/figures/empirical_condition_by_condition.png
#
# IMPORTANT — 0 % condition is ALWAYS 0.0, never the empirical false-alarm rate.
# At 0 % contrast there is literally no stimulus; the ground-truth detection
# probability is therefore exactly 0.0.  Using the empirical FA rate here would
# misrepresent the task structure to the HGF.
#
# Previous VCH_CORRECTED_MAPPING values (QC-filtered full sample; used through 2026-06-16):
#   0.0  => 0.0,                   # ground-truth: no stimulus → 0.0
#   0.25 => 0.44184090362893536,
#   0.5  => 0.7010130961205832,
#   0.75 => 0.8690302144249513,
const VCH_CORRECTED_MAPPING = Dict(   # non-hallucinator empirical for 25/50/75 %; 0 % fixed at 0.0
    0.0  => 0.0,                  # ground-truth: no stimulus → 0.0 (NOT the empirical FA rate)
    0.25 => 0.4180444024563061,
    0.5  => 0.7115104419621175,
    0.75 => 0.8994252873563219,
)
# Maps model_type → mapping dict (nothing = use condition column directly)
const VCH_STIM_MAPPINGS = Dict{String, Union{Nothing, Dict{Float64, Float64}}}(
    "2level_empiric"             => VCH_CORRECTED_MAPPING,
    "3level_empiric"             => VCH_CORRECTED_MAPPING,
    "2level_nominal"             => nothing,
    "3level_nominal"             => nothing,
)

# Stimulus column (used throughout): nominal models read condition directly;
# all others read empirical_condition (which is overwritten below for VCH).
stim_col = endswith(model_type, "_nominal") ? "condition" : "empirical_condition"

################################################################################
# ─── 5. MODEL SPECIFICATION: AGENT ARCHITECTURE AND PRIORS ─────────────────
# model_type determines two things at once:
#
#   Hierarchy depth
#     "2level_*"  → binary_2level: xbin ← xprob.  Three free parameters.
#     "3level_*"  → binary_3level: xbin ← xprob ← xvol.  Four free
#                   parameters; xvol tracks the volatility of xprob.
#
#   Priors
#     Truncated normals, identical across every variant so no variant is
#     advantaged in model comparison.  Centres are taken from posteriors
#     of an earlier CH-task HGF, widened (sd = 1) and truncated to the
#     regions where each parameter is defined:
#
#       action_precision        beta   — response consistency; > 0.001.
#                                        Higher beta = percepts track the
#                                        belief more deterministically.
#                                        Note the ^(-1): the stored value
#                                        0.29350739 is an inverse precision,
#                                        so the prior centre is its
#                                        reciprocal (≈ 3.41).
#       prior_posterior_weight  nu     — prior-vs-likelihood weighting; > 0.
#                                        Higher nu = perception dominated by
#                                        the top-down HGF prediction.
#       xprob_volatility        omega2 — log tonic volatility of xprob;
#                                        < -0.5.  Higher (less negative)
#                                        = faster belief updating.
#       xvol_volatility         omega3 — log tonic volatility of xvol,
#                                        3-level only; < -0.5.
#
# An unrecognised model_type is a hard error rather than a silent fallback,
# so a typo in a job file cannot quietly fit the wrong model.
################################################################################

## Create ActionModel & set priors ##
if model_type in ("2level_empiric", "2level_nominal")
    am = create_agent("binary_2level")
    priors = (;
        #β^(-1)
        action_precision = truncated(Normal(0.29350739^(-1), 1), lower = 0.001),
        #ν
        prior_posterior_weight = truncated(Normal(0.72646851, 1), lower = 0),
        #ω₂ (HGF submodel parameter, accessed via symbol :xprob_volatility)
        xprob_volatility = truncated(Normal(-5.1682685, 1), upper = -0.5),
    )

elseif model_type in ("3level_empiric", "3level_nominal")
    am = create_agent("binary_3level")
    priors = (;
        action_precision = truncated(Normal(0.29350739^(-1), 1), lower = 0.001),
        prior_posterior_weight = truncated(Normal(0.72646851, 1), lower = 0),
        xprob_volatility = truncated(Normal(-5.1682685, 1), upper = -0.5),
        #ω₃ (HGF submodel parameter, accessed via symbol :xvol_volatility)
        xvol_volatility = truncated(Normal(-6, 1), upper = -0.5),
    )

else
    error("Model type misspecified: \"$model_type\". Valid types: " *
          "\"2level_empiric\", \"2level_nominal\", " *
          "\"3level_empiric\", \"3level_nominal\"")
end

################################################################################
# ─── 6. PER-SUBJECT WORKER ─────────────────────────────────────────────────
# Runs the whole fit for one participant in two phases, each independently
# skippable and independently recoverable:
#
#   Phase 1  MCMC sampling      → chains, R-hats, posterior medians
#   Phase 2  State extraction   → per-trial belief trajectories
#
# The two-phase split matters operationally: state extraction has failed on
# its own in the past while sampling succeeded.  Because Phase 1 reloads
# medians from disk when chains already exist, such a subject can be
# recovered without paying for the MCMC a second time.
################################################################################

# Function to process a single subject
function process_subject(single_data, model_type, am, priors, n_iterations, n_chains, my_sampler, use_optim)
    # Extract subject ID from the data
    # ── 6.1  Resolve output paths and decide whether any work is needed ──────
    # Every artefact for this subject shares one path stem, which encodes the
    # model variant, project, timepoint, and modality — so runs of different
    # variants can never overwrite one another.
    subj_id = single_data[1, :subj_id]

    file_path   = joinpath("results", model_type, project, timepoint, modality, "$(subj_id)")
    chains_file = "$file_path.jls"
    states_file = "$(file_path)state_trajectories.csv"

    # Create results directory if it doesn't exist
    mkpath(dirname(file_path))

    # Skip entirely if both chains AND state trajectories already exist
    if skip_or_nah == "skip" && isfile(chains_file) && isfile(states_file)
        println("Subject $(subj_id) fully processed, skipping.")
        return nothing
    end

    param_names        = keys(priors)
    param_names_strings = [string(n) for n in param_names]

    # ── PHASE 1: MCMC sampling ─────────────────────────────────────────────────
    # Skip if chains already exist (e.g. state extraction failed last time)
    if skip_or_nah == "skip" && isfile(chains_file)
        println("Subject $(subj_id): chains already exist, loading medians for state extraction...")
        df_med_existing = CSV.read("$(file_path)_medians.csv", DataFrame)
        medians_vec     = [df_med_existing[1, Symbol(n)] for n in param_names_strings]
    else
        println("Processing subject $(subj_id)")
        println("Inputs: ",  single_data[:, stim_col], "  [column: $(stim_col)]")
        println("Actions: ", single_data[:, "response"])

        # Create model using ActionModels 0.7 single-session API
        model_fit = create_model(
            am,
            priors,
            collect(Float64, single_data[:, stim_col]),
            collect(Float64, single_data[:, "response"]),
            check_parameter_rejections=true,
            verbose=false,
        )

        try
            # Sample posterior
            if use_optim
                println("Sampling posterior with MAP initialization...")
                result = sample_posterior!(model_fit, MCMCSerial();
                    sampler     = my_sampler,
                    init_params = :MAP,
                    n_samples   = n_iterations,
                    n_chains    = n_chains,
                    progress    = true,
                )
            else
                println("Sampling posterior...")
                result = sample_posterior!(model_fit, MCMCSerial();
                    sampler     = my_sampler,
                    init_params = :sample_prior,
                    n_samples   = n_iterations,
                    n_chains    = n_chains,
                    progress    = true,
                )
            end

            # Save raw chains
            serialize(chains_file, result)
            CSV.write("$file_path.csv", DataFrame(result))
            println("Chains saved for subject $(subj_id)")

            # Rhats
            rhats    = DataFrame(rhat(result))
            df_rhats = unstack(rhats, :parameters, :rhat)
            df_rhats[!, :subj_id] .= subj_id
            CSV.write("$(file_path)_rhats.csv", df_rhats)
            println("Rhats saved for subject $(subj_id)")

            # Medians across all chains/iterations
            session_params  = get_session_parameters!(model_fit)
            session_id_sym  = Symbol(session_params.session_ids[1])
            medians_vec     = [median(session_params.value[pn][session_id_sym]) for pn in param_names]
            row             = (; subj_id=subj_id, [Symbol(n)=>m for (n,m) in zip(param_names_strings, medians_vec)]...)
            CSV.write("$(file_path)_medians.csv", DataFrame([row]))
            println("Medians saved for subject $(subj_id)")

        catch e
            println("Error during MCMC for subject $(subj_id): $e")
            rethrow(e)
        end
    end

    # ── PHASE 2: State trajectory extraction ───────────────────────────────────
    try
        println("Extracting state trajectories for subject $(subj_id)...")

        desired_states = [
            :xbin_prediction_precision,
            :xbin_prediction_mean,
            :xbin_posterior_mean,
            :xbin_value_prediction_error,
            :xprob_posterior_precision,
            :xprob_posterior_mean,
            :xprob_prediction_mean,
            :xprob_prediction_precision,
            :xprob_effective_prediction_precision,
            :xprob_value_prediction_error,
            :xprob_precision_prediction_error,
            :belief,
        ]

        # For 3-level models, also extract xvol (3rd-level volatility) states
        if startswith(model_type, "3level")
            append!(desired_states, [
                :xvol_posterior_mean,
                :xvol_posterior_precision,
                :xvol_prediction_mean,
                :xvol_prediction_precision,
                :xvol_effective_prediction_precision,
                :xvol_value_prediction_error,
                :xvol_precision_prediction_error,
            ])
        end

        # Build parameters NamedTuple from medians
        parameters_nt = NamedTuple(Symbol(n) => m for (n,m) in zip(param_names_strings, medians_vec))

        # Use the same HGF level as was fitted — 3-level models have xvol_volatility
        # in parameters_nt; passing that to a 2-level agent throws "Parameter not found"
        sim_hgf_type = startswith(model_type, "3level") ? "binary_3level" : "binary_2level"
        am_sim       = create_agent(sim_hgf_type)
        agent_sim    = init_agent(am_sim, save_history=desired_states)

        set_parameters!(agent_sim, parameters_nt)
        reset!(agent_sim)

        inputs  = single_data[:, Symbol(stim_col)]
        actions = single_data[:, :response]

        for (input, action) in zip(inputs, actions)
            agent_sim.action_model(agent_sim.model_attributes, input)
            set_actions!(agent_sim, :action, Float64(action))
            for (state_name, state_vec) in pairs(agent_sim.history)
                push!(state_vec, get_states(agent_sim, state_name))
            end
        end

        # Skip index 1 (pre-loop initial state) so all vectors align at T entries
        history_nt   = get_history(agent_sim)
        state_history = Dict(string(k) => collect(v[2:end]) for (k,v) in pairs(history_nt))
        filter!(!ismissing, state_history["belief"])

        states_df          = DataFrame(state_history)
        states_df.timestep = 0:(nrow(states_df)-1)
        states_df.ID       = fill(subj_id, nrow(states_df))
        states_df.record_id = fill(subj_id, nrow(states_df))
        states_df.subj_id  = fill(subj_id, nrow(states_df))
        states_df.modality = fill(modality, nrow(states_df))

        CSV.write(states_file, states_df)
        println("State trajectories saved for subject $(subj_id)")

    catch e
        println("Error during state extraction for subject $(subj_id): $e")
        rethrow(e)
    end

    return nothing
end


################################################################################
# ─── 7. MAIN: LOCATE AND LOAD THIS ARRAY TASK'S PARTICIPANT ────────────────────
# Everything above is definition; execution starts here.
#
# The job array passes an index rather than a subject ID, so this block
# resolves index → file → dataframe, and fails loudly if the index is out of
# range (a stale job file pointed at a smaller staging directory).
################################################################################
# Load data files
data_files = glob("*.csv", "./data_n_cmnds/$(project)/$(timepoint)/$(modality)_data")

# Check if file index is valid
if file_index < 1 || file_index > length(data_files)
    error("File index $file_index is out of range. Found $(length(data_files)) files.")
end

# Select the specific file
selected_file = data_files[file_index]
println("Processing file $file_index: $selected_file")

# Load the single dataset
single_data = CSV.read(selected_file, DataFrame)


# ── 7.1  Cross-contamination guard: recompute the stimulus mapping ──────────
# The staging CSV carries an empirical_condition column, but its value is NOT
# trusted.  It is overwritten here from the condition column using this
# model_type's mapping, so a staging file written under an older mapping can
# never silently propagate the wrong stimulus intensities into a new fit.
# This is why catch-up reruns are safe without re-staging.
# For VCH: recompute empirical_condition at runtime from the condition column
# using the hardcoded mapping for this model_type.  This overrides whatever
# value was stored in the CSV so catch-up reruns and future pipeline changes
# cannot accidentally use the wrong mapping.
if modality == "vch"
    vch_map = get(VCH_STIM_MAPPINGS, model_type, VCH_CORRECTED_MAPPING)
    if vch_map !== nothing
        single_data[!, :empirical_condition] = Float64[vch_map[c] for c in single_data[!, :condition]]
        println("Stimulus (VCH/$model_type): empirical_condition recomputed from condition column at runtime")
    end
end


# ── 7.2  Input validation ───────────────────────────────────────────────────
# Fail loudly on a malformed staging file rather than fitting the model to
# whatever columns happen to be present.
# Ensure required columns exist
required_columns = ["response", "empirical_condition", "condition", "subj_id"]
missing_columns = setdiff(required_columns, names(single_data))
if !isempty(missing_columns)
    error("Missing required columns: $missing_columns")
end


################################################################################
# ─── 8. RUN ───────────────────────────────────────────────────────────────────
################################################################################
# Process the single subject
process_subject(single_data, model_type, am, priors, n_iterations, n_chains, my_sampler, use_optim)

println("Subject processing completed!")
