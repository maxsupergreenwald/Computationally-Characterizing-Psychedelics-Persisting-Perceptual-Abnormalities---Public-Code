#!/usr/bin/env python3
"""
generate_hpc_jobs.py
====================
Unified HPC job-array generator. Combines the logic of:
  • generate_job_arrays.ipynb          (nonsp single-path regressions)
  • generate_mediation_job_arrays.ipynb (Bayesian mediation models)

All enabled sections are merged into ONE combined job-array .txt file
(COMBINED_STEM) in OUTPUT_NONSP, which is rsync'd to HPC_BASE.
Mediation-only assets (R script, helper_scripts, augmented CSV) are
rsync'd separately to HPC_MED_DIR.

Run from 03_hpc/:
    /usr/local/bin/python3.12 generate_hpc_jobs.py

The master CSV ships RAW
------------------------
Since 2026-08-22 this script performs NO normalization.  df_foranalysis_master.csv
carries raw values and every Gelman 2SD transform happens on the cluster, in
gelman_normalization.R, after the R scripts have applied their row filters and
dropped incomplete cases — so each model is centred and scaled on exactly the
rows it is fitted on.  This script still owns the MEMBERSHIP of each rule and
exports it in normalization_vars.R.  See the NORMALIZATION section below.

Outputs
-------
  data/final/nonsp_predictor_analyses/
      df_foranalysis_master.csv         ← master df (RAW) for ALL nonsp + mediation jobs
      nonsp_predictors.R                ← copied R script
      gelman_normalization.R            ← copied R normalizer
      convergence_gate.R                ← copied Rhat/ESS/divergence gate (the transform)
      convergence_gate.R                ← copied Rhat/ESS/divergence gate
      normalization_vars.R              ← generated rule membership (from master_config)
      categorical_factor_vars.R         ← generated factor list (from master_config)
      monotonic_covariates_vars.R       ← generated mo() list (from master_config)
      master_allpreds_x_hppd_caps.txt   ← Part 1 job lines (if enabled)
      master_sp_x_vch_dvs.txt           ← Part 2 job lines (if enabled)
      {COMBINED_STEM}.txt               ← all enabled jobs merged

  03_hpc/output/mediation_analyses/
      df_foranalysis_master.csv         ← same master df (RAW), copied
      hpc_mediation.R                   ← copied R script
      helper_scripts/                   ← copied helper scripts
      gelman_normalization.R            ← copied R normalizer
      convergence_gate.R                ← copied Rhat/ESS/divergence gate
      normalization_vars.R              ← generated rule membership
      categorical_factor_vars.R         ← generated factor list
      monotonic_covariates_vars.R       ← generated mo() list
      {med_stem}.txt  ...               ← individual mediation job files
      mediation_combined.txt            ← all mediation lines merged
"""


# ==============================================================================
# ████████████████████████  CONFIG — edit here only  ███████████████████████████
# ==============================================================================

# ── Base covariate model keys ─────────────────────────────────────────────────
# Each entry is resolvable by _resolve_base_covs(); canonical definitions live in
# BASE_COVARIATE_SETS in modules/master_config.py.
#
# Every entry is ALSO the HPC results directory name, so each must stay unique —
# two entries sharing a name overwrite each other's results on the cluster,
# silently.
#
# Keep this list in step with MODEL_VARIANTS in
# 04_visualizations/supplement/sensitivity_analyses.py (Supp. Fig. S4) and the
# model-type list in sensitivity_analyses_mediation.py (Supp. Fig. S5): a model
# type named there but absent here is never submitted, and shows up as a missing
# heatmap column.
BASE_MODELS = [
    # Primary manuscript model — Figures 4, 5, 6; Tables S4, S5.
    'nice_covariates_spusers',
    # Leading heatmap columns in S4/S5, beside the primary model.
    'nice_covariates',                            # full sample (adds 42 SP-naive)
    'nice_covariates_spusers_iqr',                # IQR fence on the focal predictor
    # Sensitivity columns — S4 (MODEL_VARIANTS) and S5.
    'empirical_covariates_spusers',               # strongest available caps_vision confounders
    'age_control_spusers',
    'true_univariate_spusers',
    'nice_covariates_beta_spusers',               # decision-noise (vch_beta) control
    'nice_covariates_beta_spusers_iqr',           # beta control + IQR fence
    'nice_covariates_spusers_nonan_caps',         # participants with CAPS data
    'drugs_month_spusers',                        # past-month drug classes, all kept separate
    'drugs_trimmed_month_spusers',                # past-month, grouped by receptor pharmacology
    'nice_covariates_spusers_hardware_control',   # + display class (reviewer request)
]

# ── Subsetting modifiers ──────────────────────────────────────────────────────
# R-side subsetting keywords (see 03_hpc/README.md — Row-Level Subsetting Keywords).
# The script generates ALL permutations of subsets from this list and appends
# them to each base model, e.g.:
#   ['spusers', 'iqr'] → base, base_spusers, base_iqr, base_spusers_iqr
# Supported keywords are exactly R_SIDE_SUFFIXES in modules/master_config.py:
#            'spusers', 'iqr', 'nocurrenthppd', 'nopsychosis', 'nonan_caps'.
# Outlier exclusion is not a subsetting keyword; use 'iqr' (see R_SIDE_SUFFIXES
# in master_config.py). Naming 'nooutlier' here raises at generation time.
MODIFIERS = []



# ── Variant flags ─────────────────────────────────────────────────────────────
# If True, generate _beta variants for each model (adds vch_beta as covariate).
ADD_BETA_VARIANTS = False

# If True, append _nocurrenthppd to every non-nopsychosis generated variant,
# matching the notebook's all_model_variants_nocurrenthppd logic.
ADD_NOCURRENTHPPD_VARIANTS = False

# If True, extend the VCH predictor/DV groups to include 3-level Julia HGF
# parameters (vch_comps_3lev) and MATLAB HGF parameters (vch_comp_mat).
INCLUDE_3LEVEL_VCH = False

# ── Job section toggles ────────────────────────────────────────────────────────
# All enabled sections are merged into COMBINED_STEM.
INCLUDE_HPPD_CAPS_JOBS = True    # nonsp: all predictors → HPPD/CAPS DVs (see PART1_DVS / PART1_DV_GROUPS)
INCLUDE_SP_VCH_JOBS    = False  # SP -> VCH single paths 
INCLUDE_MEDIATION_JOBS = True   # mediation


# ── Predictor / DV overrides ──────────────────────────────────────────────────
# If None, defaults apply:
#   Part 1 predictors: SP_PREDS + VCH_BEH + VCH_COMP
#   Part 1 DVs:        persist_vis_yn + caps_vision
#   Part 2 predictors: SP_PREDS
#   Part 2 DVs:        VCH_BEH + VCH_COMP
# Set to a list of iv_type_dict keys to override (resolved after data is loaded).
# PART1_PREDICTOR_GROUPS = ["sp_predictors","vch_behavior",'vch_computations']#, 'vch_comp_nominal']   # empirical + nominal params → HPPD/CAPS
# PART1_PREDICTOR_GROUPS = ["vch_comp_nominal","vch_comp_avg"]#, 'vch_comp_nominal']   # empirical + nominal params → HPPD/CAPS
# PART1_PREDICTOR_GROUPS = ["sdt_hppd"]#, 'vch_comp_nominal']   # empirical + nominal params → HPPD/CAPS
# Primary-manuscript reproduction run (2026-08-21): the four groups reported in
# the paper.  'vch_comp_nominal' (nominal-parameter HGF variants) is a supplementary
# set and is excluded here — restore the commented line below to include it.
# PART1_PREDICTOR_GROUPS = ["sp_predictors","vch_behavior",'vch_computations','sdt_hppd','vch_comp_nominal']
# 'vch_comp_nominal' = vch_nu_nominal, vch_beta_nominal, vch_omega_nominal —
# the HGF parameters refitted under the nominal likelihood, read by
# supplement/regression_results_table_nominal_sensitivity.py.
PART1_PREDICTOR_GROUPS = ["sp_predictors", "vch_behavior", "vch_computations",
                          "sdt_hppd", "vch_comp_nominal"]
PART1_DV_GROUPS        = None                                         # default: hppd_binary + caps_vision. Set to list of iv_type_dict keys to override.
PART2_PREDICTOR_GROUPS = None                                         # default: sp_predictors
PART2_DV_GROUPS        = ["vch_behavior",'vch_computations','vch_comp_nominal']    # SP → empirical + nominal params

# Direct DV list for Part 1 (column names, NOT iv_type_dict keys).
# Overrides PART1_DV_GROUPS if not None.
# Canonical default is ['hppd_binary', 'caps_vision']; set to None to use that default.
# PART1_DVS = ['hppd_binary', 'caps_vision', 'baggot_binary','caps_intensity','caps_total',
#              'caps_maximum_frequency','caps_vision_maximum_frequency','caps_intensity_maximum_frequency','lshs_total', 'baggot_total']
PART1_DVS = ['caps_vision','hppd_binary']#,'lshs_total', 'baggot_total']

# PART1_DVS = ['lshs_total']#,'lshs_total', 'baggot_total']

# caps_vision_formed and caps_vision_bottomup removed: pp_check(type="error_scatter_avg_vs_x")
# crashes for these DVs due to chain divergence, which also blocks subsequent bundled DVs via &&.

# ── Mediation analyses ─────────────────────────────────────────────────────────
# List of dicts. Each dict generates one .txt job file. Keys:
#   spvar     : SP predictor shorthand — see SP_PREDICTORS dict in setup
#               ('spage', 'lifenomic', 'avgdose', 'pc1', 'pc1ranked')
#   mediator  : VCH mediator shorthand — see MEDIATORS dict in setup
#               ('vchthreshold', 'vchrate', 'vchrate75', 'vchnu', 'vchbeta',
#                'vchomega', 'vchnu3lev', 'vchbeta3lev', 'vchnumat', 'vchbetamat')
#   dv        : Outcome shorthand — see DVS dict in setup
#               ('hppd_binary', 'persist_vis_yn', 'caps_vision', 'caps_total',
#                'hppd_true_chronicity', 'persistvis_duration', 'baggot_total', 'baggot_binary')
#   cov_types : Tuple of strings parsed by get_covs(); base types: 'nice_covariates',
#               'main', 'age_only', 'univariate', 'empiric_covariates'.
#               Append R-side suffixes: 'nice_covariates_spusers',
#               'nice_covariates_iqr', 'nice_covariates_spusers_iqr', etc.
#   stem      : (optional) output .txt stem; auto-derived as '{dv}_{spvar}_{mediator}' if omitted
#
# Default: spage/avgdose × [vchthreshold, vchrate, vchbeta, vchnu] × [persist_vis_yn, caps_vision]
# with nice_covariates variants — matches Cell 16 of generate_mediation_job_arrays.ipynb.
# Extend this list to replicate other cells from that notebook.

# ── Mediation covariate types ─────────────────────────────────────────────────
# If non-empty, all mediation specs use exactly these covariate type strings.
# If empty (default), _MED_COV_TYPES defaults to MODEL_VARIANTS at runtime —
# the same set of model variants used for the nonsp analysis.
# Example override: ['nice_covariates', 'nice_covariates_spusers']
#
# Every mediation model type any supplement figure reads: SENSITIVITY_MED_TYPES
# from 04_visualizations/supplement/sensitivity_analyses_mediation.py, plus its
# CANONICAL_MED_TYPE.  This list and that one must stay in step — a type present
# there but missing here leaves a permanently blank column in the NIE heatmap,
# with no error raised at either end.
# Note this set has no nice_covariates_beta_spusers_iqr — that variant is read by
# the nonsp figures only, and so appears in BASE_MODELS rather than here.
CUSTOM_MED_TYPES = [
    'nice_covariates_spusers',              # CANONICAL
    'nice_covariates',                      # SECOND
    'nice_covariates_spusers_iqr',          # THIRD
    'empirical_covariates_spusers',
    'age_control_spusers',
    'true_univariate_spusers',
    'nice_covariates_beta_spusers',
    'nice_covariates_spusers_nonan_caps',
    'drugs_month_spusers',                  # + past-month drug-use binaries (full)
    'drugs_trimmed_month_spusers',          # + past-month drug-use binaries (collapsed)
    'nice_covariates_spusers_hardware_control',  # + display class
]


# Primary-manuscript reproduction run (2026-08-21): the four core VCH mediators.
# The SDT mediators (dprimeoverall, criterionoverall, meanconffas) and the nominal
# HGF variants (vchnunominal, vchbetanominal) are excluded here — restore the
# commented line below to include them.
# _MED_MEDIATORS = ['dprimeoverall', 'criterionoverall', 'meanconffas','vchthreshold', 'vchrate', 'vchbeta', 'vchnu', 'vchnunominal', 'vchbetanominal']
# Full reporting run (2026-08-23): the four core VCH mediators plus the two HGF
# nominal-likelihood variants read by
# 04_visualizations/supplement/mediation_results_table_nominal_sensitivity.py
# (MEDIATORS = ["vchbetanominal", "vchnunominal"]).
# There is no 'vchomeganominal' key in the MEDIATORS dict, so vch_omega_nominal
# is a nonsp predictor only — it has no mediation counterpart to generate.
# The SDT mediators (dprimeoverall, criterionoverall, meanconffas) remain out:
# no supplement figure reads them as mediators.
_MED_MEDIATORS = ['vchthreshold', 'vchrate', 'vchbeta', 'vchnu',
                  'vchnunominal', 'vchbetanominal']

# NOTE: MEDIATION_ANALYSES is defined below END CONFIG, immediately after
# MODEL_VARIANTS is built, because _MED_COV_TYPES defaults to MODEL_VARIANTS.

# ── MCMC settings ─────────────────────────────────────────────────────────────
NONSP_ITERATIONS = 10000
NONSP_WALLTIME   = '90:00'
MED_ITER         = 10000
MED_WARMUP       = 6000
# Mediation models need 360:00, not 90:00 -- some run long, and a timed-out
# job has to be resubmitted from scratch.  Walltime is per job, not per array.
MED_WALLTIME     = '360:00'

# ── HPC paths ─────────────────────────────────────────────────────────────────
HPC_PARENT      = '/nfs/roberts/scratch/pi_arp29/msg74/aim1_baseline_final'
HPC_BASE        = f'{HPC_PARENT}/nonsp_predictor_analyses'
HPC_MED_PARENT  = f'{HPC_PARENT}/hppd_manuscript_mediation_analyses'
HPC_MED_DIR     = f'{HPC_MED_PARENT}/mediation_analyses'
HPC_MED_RESULTS = f'{HPC_MED_DIR}/results'
HPC_MED_HELPERS = f'{HPC_MED_DIR}/helper_scripts'
HPC_USER        = 'msg74'
HPC_TRANSFER    = f'{HPC_USER}@transfer-bouchet.ycrc.yale.edu'
HPC_LOGIN       = f'{HPC_USER}@bouchet.ycrc.yale.edu'

# ── Output ────────────────────────────────────────────────────────────────────
# Combined file lives in OUTPUT_NONSP (rsync'd to HPC_BASE); dSQ runs from HPC_BASE.
# Renamed 2026-08-21: the previous stem ('sdt_mediation_analyses') described an
# SDT-mediator batch that this config no longer generates.  The stem becomes the
# .txt filename passed to dsq, so change it whenever the scope of the run changes.
COMBINED_STEM = 'combined_all_analyses'

# ── Post-retrieve: compile results ────────────────────────────────────────────
# If True, after printing the RETRIEVE commands the script will ask whether to
# run compile_nonsp_results.py subset immediately.  Set to False to suppress.
PROMPT_TO_COMPILE = True

# ==============================================================================
# END CONFIG ───────────────────────────────────────────────────────────────────
# ==============================================================================


# ==============================================================================
# SETUP  [replicates generate_job_arrays.ipynb Setup Cell 1]
# ==============================================================================
import sys, os, shutil

# Set by compile_nonsp_results.py / compile_mediation_results.py, which exec this
# file purely to read its CONFIG.  When true the data load is skipped and
# execution stops at the END OF CONFIG DERIVATION marker below, so importing the
# settings can never write a file, hit the network, or prompt.
_CONFIG_ONLY = os.environ.get('HPC_JOBS_CONFIG_ONLY') == '1'
from pathlib import Path
from itertools import chain, combinations
from difflib import get_close_matches

_SCRIPT_DIR   = Path(__file__).resolve().parent          # .../03_hpc/
_PROJECT_ROOT = _SCRIPT_DIR.parent                       # .../hppd_manuscript_public/
sys.path.insert(0, str(_PROJECT_ROOT / 'modules'))

import pandas as pd
import numpy as np

# normalize_analysis_df() is deliberately NOT imported: normalization moved to R
# (2026-08-22).  BINARY_FACTOR_VARS is still needed, to pass R the same skip list.
from data_prep import (most_recent_public_df, BINARY_FACTOR_VARS)
from master_config import (
    dv_to_lab, dv_to_lab_short, iv_type_dict,
    caps_types, hppd_variables, severity_vars,
    main_predictors, all_predictors,
)
from master_config import (
    R_SIDE_SUFFIXES, BASE_COVARIATE_SETS, VARIABLE_REGISTRY,
    NEED_NON_NORMALIZED, INPLACE_NORMALIZED,
    CATEGORICAL_FACTOR_VARS, MONOTONIC_COVARIATES,
)

# Local output directories
OUTPUT_NONSP = _PROJECT_ROOT / 'data' / 'final' / 'nonsp_predictor_analyses'
OUTPUT_MED   = _SCRIPT_DIR / 'output' / 'mediation_analyses'
OUTPUT_NONSP.mkdir(parents=True, exist_ok=True)

# Augment DV label dict (matches notebook)
dv_to_lab['hppd_sx_count'] = '# of HPPD Experiences'

# Load data
# The shipped analysis dataframe is already fully prepared — every derived column
# is present in the CSV, so there is no preparation step to run here.  The R
# scripts on the cluster likewise read df_foranalysis_master.csv verbatim.
# most_recent_public_df() picks the newest data/final/df_public_*.csv and raises
# a clear FileNotFoundError if the directory is empty.
if not _CONFIG_ONLY:
    print('Loading analysis dataframe...')
    df = pd.read_csv(most_recent_public_df(_PROJECT_ROOT / 'data' / 'final'), low_memory=False)
    print(f'  df shape: {df.shape}')

# `monitor_check_operationalized_final` is a hand coding of the free-text
# `monitor_check` field and cannot be recomputed by rule, so it is not derived
# here — it ships as a column of data/final/df_public_*.csv and is read straight
# off it.  Needed by the covariate set `nice_covariates_hardware_control`.  The
# coding itself, with its five documented per-string overrides, is defined in
# 04_visualizations/supplement/hardware_keydown_check.py, which also plots it as
# a supplementary figure.
# Printed here so a changed or absent coding is visible in the run log.
if not _CONFIG_ONLY:
    print(f"  + {'monitor_check_operationalized_final'}: "
          f"{dict(df['monitor_check_operationalized_final'].value_counts(dropna=False))}")


# ==============================================================================
# COVARIATE DICT + DF_FORANALYSIS BUILD
# ==============================================================================

# ── Canonical covariate sets ──────────────────────────────────────────────────
# Single source of truth: ALL covariate sets are defined in
# modules/master_config.py :: BASE_COVARIATE_SETS (including aliases and derived
# sets such as 'univariate', 'empirical_covariates', and the beta-control
# variants).  Do NOT add entries here — edit master_config.py instead.
#
# IMPORTANT: If a covariate list includes psycheduse_life_nomic (or any variable
# matching '_nomic', 'psychedelic_age', 'PC1', 'lastuse_dose'), generate_nonsp_job_array()
# will replace that variable with the focal predictor (spvar replacement logic).
# This is intentional for model types like true_univariate and age_control.
# nice_covariates has no SP variable, so the focal predictor is simply appended.
#
# Runtime auto-generated R-side-suffix variants are appended to this dict below.
covariates_master_dict = dict(BASE_COVARIATE_SETS)


def _covs_for_variant(variant_name):
    """
    Derive the covariate list for a full variant name.

    Strips all R-side subsetting suffixes from variant_name (these are handled
    inside R, not Python), then strips any _beta suffix, then looks up the
    resulting base key in covariates_master_dict.  If _beta was present, appends
    vch_beta to the covariate list.

    All suffix stripping uses exact replacement from R_SIDE_SUFFIXES (master_config),
    so the result is always an exact covariates_master_dict key — no startswith
    fallback needed.
    """
    stem = variant_name
    for sfx in R_SIDE_SUFFIXES:
        stem = stem.replace(sfx, '')
    add_beta = '_beta' in stem and 'beta_both' not in stem
    base_key = stem.replace('_beta', '')
    if base_key not in covariates_master_dict:
        raise ValueError(
            f'Cannot resolve covariate set for variant {variant_name!r} '
            f'(stripped to {base_key!r}). '
            f'Known keys: {sorted(covariates_master_dict.keys())}'
        )
    covs = covariates_master_dict[base_key].copy()
    if add_beta:
        covs = covs + ['vch_beta']
    return covs


# ── Model variant auto-generation from CONFIG ─────────────────────────────────
def _powerset(s):
    """All subsets of s, starting from the empty set."""
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))


def build_model_variants(base_models, modifiers, add_beta, add_nocurrenthppd):
    """
    Generate model variant names from BASE_MODELS × MODIFIERS permutations.

    For each base model and each subset of modifiers (including the empty set),
    generates 'base_mod1_mod2_...' in the order modifiers appear in the list.
    If add_beta=True, also generates 'base_beta_mod1_mod2_...' variants.
    If add_nocurrenthppd=True, appends _nocurrenthppd to all non-nopsychosis variants.
    """
    variants = []
    for base in base_models:
        for mod_combo in _powerset(modifiers):
            suffix = ''.join(f'_{m}' for m in mod_combo)
            v = f'{base}{suffix}'
            variants.append(v)
            if add_beta:
                variants.append(f'{base}_beta{suffix}')

    if add_nocurrenthppd:
        nocurrenthppd = [f'{v}_nocurrenthppd' for v in variants if 'nopsychosis' not in v]
        variants.extend(nocurrenthppd)

    # Deduplicate preserving order
    seen, result = set(), []
    for v in variants:
        if v not in seen:
            seen.add(v)
            result.append(v)
    return result


MODEL_VARIANTS = build_model_variants(BASE_MODELS, MODIFIERS, ADD_BETA_VARIANTS, ADD_NOCURRENTHPPD_VARIANTS)

# ── Retired keyword guard ─────────────────────────────────────────────────────
# 'nooutlier' is not a supported keyword (see R_SIDE_SUFFIXES in
# master_config.py).  _covs_for_variant() would already fail to resolve such a
# variant, but the message would be about an unknown covariate-set key and would
# not say why.  Fail here instead, naming the reason.
_RETIRED_KEYWORDS = {
    'nooutlier': ('not supported — it thresholded |spvar| in Gelman-normalized '
                  'space, which has no meaning now that normalization happens in R '
                  'after subsetting. Use "iqr", whose fence is quantile-based and '
                  'therefore scale-free.'),
}
for _kw, _why in _RETIRED_KEYWORDS.items():
    _offenders = sorted({v for v in MODEL_VARIANTS if _kw in v}
                        | {t for t in (CUSTOM_MED_TYPES or []) if _kw in t})
    if _offenders:
        raise ValueError(
            f'Retired subsetting keyword {_kw!r} requested by: {_offenders}. {_why}'
        )

print(f'\nModel variants ({len(MODEL_VARIANTS)}):')
for v in MODEL_VARIANTS:
    print(f'  {v}')

# Resolve mediation covariate types.
# CUSTOM_MED_TYPES (set in CONFIG) overrides; empty → same variants as nonsp analysis.
_MED_COV_TYPES = tuple(CUSTOM_MED_TYPES) if CUSTOM_MED_TYPES else tuple(MODEL_VARIANTS)
print(f'Mediation covariate types ({len(_MED_COV_TYPES)}): {list(_MED_COV_TYPES)}')

# Primary-manuscript reproduction run (2026-08-21): only the two pairings the
# manuscript reports, which are also the two that compile_mediation_results.py
# mirrors.  The lshs_total and baggot_total blocks are secondary outcomes and are
# commented out below — uncomment them to regenerate those job files.
MEDIATION_ANALYSES = [
    # psychedelic_age → VCH mediator → PPA history (hppd_binary)
    *[{'spvar': 'spage',   'mediator': m, 'dv': 'hppd_binary',  'cov_types': _MED_COV_TYPES}
      for m in _MED_MEDIATORS],
    # avg_life_dose → VCH mediator → CAPS vision
    *[{'spvar': 'avgdose', 'mediator': m, 'dv': 'caps_vision',  'cov_types': _MED_COV_TYPES}
      for m in _MED_MEDIATORS],
    # # psychedelic_age → VCH mediator → LSHS total (student_t DV)
    # *[{'spvar': 'spage',   'mediator': m, 'dv': 'lshs_total',   'cov_types': _MED_COV_TYPES}
    #   for m in _MED_MEDIATORS],
    # # avg_life_dose → VCH mediator → LSHS total
    # *[{'spvar': 'avgdose', 'mediator': m, 'dv': 'lshs_total',   'cov_types': _MED_COV_TYPES}
    #   for m in _MED_MEDIATORS],
    # # psychedelic_age → VCH mediator → Baggot total (hurdle_negbinom_huvary DV)
    # *[{'spvar': 'spage',   'mediator': m, 'dv': 'baggot_total', 'cov_types': _MED_COV_TYPES}
    #   for m in _MED_MEDIATORS],
    # # avg_life_dose → VCH mediator → Baggot total
    # *[{'spvar': 'avgdose', 'mediator': m, 'dv': 'baggot_total', 'cov_types': _MED_COV_TYPES}
    #   for m in _MED_MEDIATORS],
]

# ── DV and predictor resolution ───────────────────────────────────────────────
# Hoisted above the END OF CONFIG DERIVATION marker 2026-08-31 so the compile
# scripts can read them.  Pure lookups into iv_type_dict — no data required.

# DV shorthand → column_name only.
# dv_family is derived at use-time from VARIABLE_REGISTRY[column_name]['distribution'].
# To update a DV's family, change the distribution field in master_config.py.
DVS = {
    'hppd_binary':          'hppd_binary',
    'persist_vis_yn':       'persist_vis_yn',
    'caps_vision':          'caps_vision',
    'caps_total':           'caps_total',
    'hppd_true_chronicity': 'hppd_true_chronicity',
    'persistvis_duration':  'persistvis_duration',
    'baggot_total':         'baggot_total',
    'baggot_binary':        'baggot_binary',
    'lshs_total':           'lshs_total',
}

# Predictor groups (base defaults)
SP_PREDS  = iv_type_dict['sp_predictors']
VCH_BEH   = iv_type_dict['vch_behavior']
VCH_COMP  = iv_type_dict['vch_computations']
if INCLUDE_3LEVEL_VCH:
    VCH_COMP = VCH_COMP + iv_type_dict['vch_comps_3lev'] + iv_type_dict['vch_comp_mat']

# Resolve CONFIG override groups → flat variable lists
def _resolve_groups(group_keys, default):
    if group_keys is None:
        return default
    result = []
    for key in group_keys:
        result += iv_type_dict[key]
    return list(dict.fromkeys(result))   # deduplicate preserving order

ALL_PREDS     = _resolve_groups(PART1_PREDICTOR_GROUPS, SP_PREDS + VCH_BEH + VCH_COMP)
_PART1_DV_DEFAULT = ['hppd_binary', 'caps_vision']  # canonical default — do NOT change this line
HPPD_CAPS_DVS = (
    PART1_DVS if PART1_DVS is not None
    else _resolve_groups(PART1_DV_GROUPS, _PART1_DV_DEFAULT)
)
PART2_PREDS   = _resolve_groups(PART2_PREDICTOR_GROUPS,  SP_PREDS)
VCH_AS_DVS    = _resolve_groups(PART2_DV_GROUPS,         VCH_BEH + VCH_COMP)

# ══════════════════════════════════════════════════════════════════════════════
# ██████████████████  END OF CONFIG DERIVATION  ███████████████████████████████
# ══════════════════════════════════════════════════════════════════════════════
# Everything the compile scripts need is defined above:
#   MODEL_VARIANTS  ALL_PREDS  HPPD_CAPS_DVS  PART2_PREDS  VCH_AS_DVS
#   MEDIATION_ANALYSES  _MED_COV_TYPES  DVS
# compile_nonsp_results.load_generator_config() execs this file with
# HPC_JOBS_CONFIG_ONLY=1 and stops here, so those two scripts always compile
# exactly what this file generates.  Everything BELOW writes files.
if _CONFIG_ONLY:
    raise SystemExit(0)


# Register auto-generated variants (R-side suffix combos) in covariates_master_dict
_registered = 0
for variant in MODEL_VARIANTS:
    if variant not in covariates_master_dict:
        covariates_master_dict[variant] = _covs_for_variant(variant)
        _registered += 1
print(f'Registered {_registered} new variants in covariates_master_dict. '
      f'Total: {len(covariates_master_dict)}')

# ── Build df_foranalysis  ─────────────────────────────────────────────────────
def _suggest_cols(name, cols):
    return get_close_matches(name, list(cols), n=5, cutoff=0.6)



df_foranalysis = df.copy()

# Zero-fill dose/frequency for non-SP users
for col in ['psyched_lastuse_dose', 'vasdose_bl']:
    if col in df_foranalysis.columns:
        df_foranalysis.loc[df['psychedelicuse_lifetimetot'] < 1, col] = 0

df_foranalysis.loc[df_foranalysis['psycheduse_life_nomic'] < 1,
                   ['psyched_yearsofuse', 'psychedelic_life_freq']] = 0
df_foranalysis.loc[df_foranalysis['psycheduse_life_nomic'] > 1, 'psychedelic_life_freq'] = (
    df_foranalysis.loc[df_foranalysis['psycheduse_life_nomic'] > 1, 'psycheduse_life_nomic']
    / df_foranalysis.loc[df_foranalysis['psycheduse_life_nomic'] > 1, 'psyched_yearsofuse']
)

if 'life_exposure' not in df_foranalysis.columns:
    for level, dose in zip(['micro', 'low', 'medium', 'heavy', 'vheavy'], [25, 50, 175, 400, 800]):
        df_foranalysis[f'{level}_wt'] = dose * df_foranalysis[f'psyched_percent_{level}']
    df_foranalysis['avg_life_dose'] = df_foranalysis[
        [f'{l}_wt' for l in ['micro', 'low', 'medium', 'heavy', 'vheavy']]
    ].sum(axis=1)
    df_foranalysis['life_exposure'] = df_foranalysis['avg_life_dose'] * df_foranalysis['psycheduse_life_nomic']

if len(df_foranalysis.loc[df_foranalysis['motivation'] == 0]) > 0:
    df_foranalysis.loc[df_foranalysis['motivation'] < 1, 'motivation'] = 5

# BRMS does not accept variables starting with 'sp_'
if 'sp_primary' in df_foranalysis.columns:
    df_foranalysis['psychedelic_primary'] = df_foranalysis['sp_primary']
if 'psychedelic_primary' not in df_foranalysis.columns and 'psychedelic_primary' in df.columns:
    df_foranalysis['psychedelic_primary'] = df['psychedelic_primary']

# ── Categorical variables — sourced from modules/master_config.py ─────────────
# Edit CATEGORICAL_FACTOR_VARS in master_config.py to add/remove variables.
# Do NOT add variables here directly.  This set drives:
#   (1) Exclusion from Gelman normalization below
#   (2) The low-cardinality safety check before saving the master CSV
#   (3) The auto-generated categorical_factor_vars.R file sourced by both R scripts
CATEGORICAL_COVARIATES = {x for x in CATEGORICAL_FACTOR_VARS if x in df_foranalysis.columns}
categorical_vars = CATEGORICAL_COVARIATES  # alias used downstream

# ── Ordinal-as-continuous suppressions ────────────────────────────────────────
# These variables have ≤3 observed levels but are INTENTIONALLY Gelman-normalized
# and passed to brms as continuous (ordinal regression treated linearly).
# They are excluded from the low-cardinality check below but are NOT exempt from
# normalization.  If you decide to factor any of these in R, move them to
# CATEGORICAL_COVARIATES instead and add as.factor() in the R scripts.
_ORDINAL_AS_CONTINUOUS = {x for x in [
    'highest_education_3level_final',
    'ses_employ_3level',
] if x in df_foranalysis.columns}

# ── Covariate and predictor inventories ───────────────────────────────────────
# finalcovs           : every covariate referenced by any model type
# PREDICTORS_FOR_JOBS : every predictor any job array can request.  Includes the
#                       groups named in PART1_PREDICTOR_GROUPS / PART2_PREDICTOR_GROUPS
#                       so non-default groups (e.g. sdt_hppd) are covered.
finalcovs = set().union(*covariates_master_dict.values())

_extra_pred_groups: list = []
for _pg in (PART1_PREDICTOR_GROUPS or []) + (PART2_PREDICTOR_GROUPS or []):
    if isinstance(_pg, str) and _pg in iv_type_dict:
        _extra_pred_groups += iv_type_dict[_pg]
PREDICTORS_FOR_JOBS = sorted(set(
    iv_type_dict['sp_predictors']
    + iv_type_dict['vch_behavior']
    + iv_type_dict['vch_computations_expanded']
    + iv_type_dict['vch_comp_nominal']
    + iv_type_dict['vch_comp_extended_nominal']
    + _extra_pred_groups
))
main_predictors_norm = [x for x in main_predictors if x not in NEED_NON_NORMALIZED]

# ── Columns that R will create ────────────────────────────────────────────────
# Since 2026-08-22 the master CSV ships RAW and every {col}_normalized column is
# created on the cluster by gelman_normalization.R, after subsetting.  Nothing
# below can therefore check `in df_foranalysis.columns` for a normalized name —
# it would always be absent.  _WILL_EXIST_IN_R is the authoritative prediction of
# the columns the R normalizer will produce, derived from the SAME rules
# gelman_normalize_df() applies (see the NORMALIZATION section below), and is
# what every validation in this script checks against instead.
_R_CREATED_NORMALIZED_COLS = (
    {f'{c}_normalized' for c in NEED_NON_NORMALIZED if c in df_foranalysis.columns}   # rule 1
    | {f'{c}_normalized' for c in INPLACE_NORMALIZED if c in df_foranalysis.columns}  # rule 2
    | {f'{c}_normalized' for c in PREDICTORS_FOR_JOBS}                                # rule 4
)
_WILL_EXIST_IN_R = set(df_foranalysis.columns) | _R_CREATED_NORMALIZED_COLS

# Validate covariates.  A covariate set may name a {col}_normalized column that
# R creates (e.g. a covariate set naming caps_vision_normalized), so those count.
_will_exist = _WILL_EXIST_IN_R
missing_covs = [c for c in finalcovs if c not in _will_exist]
if missing_covs:
    in_raw_df = [c for c in missing_covs if c in df.columns]
    raise ValueError(
        f'Missing covariates in df_foranalysis: {missing_covs}. '
        f'Present in raw df: {in_raw_df}. '
        f'Suggestions: { {c: _suggest_cols(c, df_foranalysis.columns) for c in missing_covs} }'
    )
# Continuous covariates must be numeric; categorical ones are exempt
for c in (finalcovs & set(df_foranalysis.columns)) - CATEGORICAL_COVARIATES:
    if not pd.api.types.is_numeric_dtype(df_foranalysis[c]):
        raise TypeError(f'Non-numeric covariate {c}: {df_foranalysis[c].dropna().astype(str).unique()[:5]}')


# ══════════════════════════════════════════════════════════════════════════════
# NORMALIZATION — MOVED TO R (2026-08-22)
# ══════════════════════════════════════════════════════════════════════════════
# This script no longer normalizes anything.  df_foranalysis_master.csv ships
# with RAW values and every Gelman 2SD transform happens on the cluster, in
# gelman_normalization.R, sourced by both nonsp_predictors.R and hpc_mediation.R.
#
# WHY.  The row-level subsetting keywords ("spusers", "iqr", ...) are applied in
# R, after the CSV is read.  Normalizing here meant every subset model was fitted
# on predictors centred and scaled on a LARGER sample than the one being fitted —
# mean was not 0 and sd was not 0.5 on the rows that actually entered the
# likelihood, so "+1 raw SD" in a reported counterfactual was +1 SD of the wrong
# sample.  R now normalizes last, after subsetting AND after dropping rows
# incomplete on any model term, so the normalization sample and the estimation
# sample are the same set of rows.
#
# WHAT THIS SCRIPT STILL OWNS.  The MEMBERSHIP of each rule.  The four lists
# below are derived here from VARIABLE_REGISTRY / BASE_COVARIATE_SETS in
# modules/master_config.py, exactly as they were when they were passed to
# normalize_analysis_df(), and are written verbatim into normalization_vars.R
# (further down, beside categorical_factor_vars.R) for R to read.  master_config
# therefore remains the single source of truth; R only supplies the arithmetic.
#
#   need_non_normalized=True → raw column preserved (its brms family requires the
#                              raw scale), normalized values in {col}_normalized
#   inplace_normalized=True  → raw column normalized, {col}_normalized aliases it
#
# Everything else in finalcovs / main_predictors is normalized in place, except
# categorical variables (R applies as.factor) and monotonic ones (mo() needs raw
# integers).  Every predictor is guaranteed a {col}_normalized column.
#
# To change how a variable is normalized, edit its VARIABLE_REGISTRY entry.
# ══════════════════════════════════════════════════════════════════════════════
NORMALIZE_IN_PLACE_VARS     = sorted(set(finalcovs) | set(main_predictors_norm))
ENSURE_NORMALIZED_COPY_VARS = list(PREDICTORS_FOR_JOBS)

# BINARY_FACTOR_VARS is a SECOND, shorter categorical list that lives in
# data_prep.py rather than master_config.py, and gelman_standardize() skipped it
# unconditionally.  It is passed through to R so the port cannot diverge from the
# Python it is verified against.  The two lists disagree (BINARY_FACTOR_VARS is a
# 7-variable subset of the 16 in CATEGORICAL_FACTOR_VARS); the difference is
# inert today because no categorical variable appears in any rule's list, but the
# duplication is a pre-existing wart worth collapsing into master_config.py.
GELMAN_SKIP_VARS = list(BINARY_FACTOR_VARS)

print('\nNormalization: SKIPPED here by design — performed in R after subsetting.')
print(f'  rule 3 normalize_in_place     : {len(NORMALIZE_IN_PLACE_VARS)} vars')
print(f'  rule 4 ensure_normalized_copy : {len(ENSURE_NORMALIZED_COPY_VARS)} vars')
print(f'  R will create {len(_R_CREATED_NORMALIZED_COLS)} _normalized columns on the cluster')
# ══════════════════════════════════════════════════════════════════════════════
# END NORMALIZATION
# ══════════════════════════════════════════════════════════════════════════════

# Log-transform PRL behavioral variables (positive skew) [matches notebook Cell 2]
for col in ['regressive_errors', 'perseverative_errors']:
    if col in df_foranalysis.columns:
        df_foranalysis[f'{col}_log'] = np.log(df_foranalysis[col] + 1)

# Ensure psycheduse_yn is present (needed for spusers R-side subsetting)
assert 'psycheduse_yn' in df_foranalysis.columns, \
    "CRITICAL: psycheduse_yn missing — spusers subsetting will fail in R!"

# hppd_binary: 0/1 alias for persist_vis_yn, added so all mediation and nonsp-predictor
# HPC results route to results/hppd_binary/ consistently.  persist_vis_yn is the original
# REDCap column name; hppd_binary is the canonical name used everywhere in this pipeline.
# Verified identical: 0 ↔ PPA(-) and 1 ↔ PPA(+) with no exceptions.
#
# NOTE FOR ANYONE READING THE SHIPPED DATA FILE
# ---------------------------------------------
# In data/final/df_public_*.csv, hppd_binary is stored as the *display* strings
# "PPA (+)" / "PPA (-)" rather than 1 / 0, because that file was written after the
# figure code's label remap had been applied.  The two encodings carry identical
# information -- "PPA (-)" is the comparator level and corresponds exactly to
# persist_vis_yn == 0 -- so no participant's outcome differs between them.
#
# The reassignment below is nevertheless UNCONDITIONAL rather than guarded by
# `if 'hppd_binary' not in df_foranalysis.columns`, because brms needs a numeric
# 0/1 (or a factor) response for the bernoulli family, and the mapping from the
# string labels back to 0/1 would otherwise be left to R's as.factor(), whose
# level ordering is locale-dependent:
#
#     LC_COLLATE=en_US.UTF-8  ->  levels: "PPA (-)", "PPA (+)"   (correct: - is 0)
#     LC_COLLATE=C            ->  levels: "PPA (+)", "PPA (-)"   (INVERTED)
#
# Compute nodes frequently run non-interactive R under the C locale, so relying on
# as.factor() would silently flip the sign of every HPPD effect on some machines
# and not others.  Deriving the column from persist_vis_yn here makes the encoding
# deterministic regardless of where the job runs.
df_foranalysis['hppd_binary'] = df_foranalysis['persist_vis_yn'].astype(int)
print("  Set hppd_binary column (0/1 alias for persist_vis_yn)")

print(f'df_foranalysis shape: {df_foranalysis.shape}')
print(f"psycheduse_yn counts:\n{df_foranalysis['psycheduse_yn'].value_counts()}")


# ==============================================================================
# JOB GENERATION FUNCTIONS
# ==============================================================================

_MODULE_PREFIX = ('module purge && module load foss/2022b && '
                  'module load R/4.4.1-foss-2022b && export R_LIBS_USER=$HOME/R/4.4')

# ── DV settings helpers ────────────────────────────────────────────────────────
# These two lists identify DVs that need pipeline-specific 'settings' strings
# that differ from a plain distribution name.  All other DVs get their family
# from VARIABLE_REGISTRY['distribution'] — do not add new entries here; add
# the variable to VARIABLE_REGISTRY in master_config.py instead.
#
# 'binomial_hierarchial'  — binomial with separate trial-count column
# 'lognormal_hierarchial' — lognormal with more MCMC iterations (pwPE variables)
# 'student_t_hierarchial' — student-t with more MCMC iterations (xprob variables)
_LOGNORMAL_HIER_VARS = {
    'ach_pwPE_negative', 'ach_pwPE_positive', 'ach_pwPE_negative_0', 'ach_pwPE_positive_0',
    'ach_pwPE_negative_75', 'ach_pwPE_positive_75', 'ach_pwPE_negative_present', 'ach_pwPE_positive_present',
    'vch_pwPE_negative', 'vch_pwPE_positive', 'vch_pwPE_negative_0', 'vch_pwPE_positive_0',
    'vch_pwPE_negative_75', 'vch_pwPE_positive_75', 'vch_pwPE_negative_present', 'vch_pwPE_positive_present',
    'pwPE_negative', 'pwPE_positive', 'pwPE_negative_0', 'pwPE_positive_0',
    'pwPE_negative_1', 'pwPE_positive_1', 'pwPE_negative_correct', 'pwPE_positive_correct',
    'pwPE_negative_incorrect', 'pwPE_positive_incorrect',
}
_STUDENT_HIER_VARS = {
    'vch_xprob_0', 'vch_xprob_75', 'vch_xprob_present',
    'ach_xprob_0', 'ach_xprob_75', 'ach_xprob_present',
}
_BINOMIAL_HIER_VARS = {'total_ch_trials', 'hits_75', 'total_vch_trials', 'vch_hits_75'}
# DVs in either hierarchial set need 'long' MCMC runs (more iterations).
_LONG_RUN_VARS = _LOGNORMAL_HIER_VARS | _STUDENT_HIER_VARS


def _dv_settings(dv):
    """
    Return the nonsp_predictors.R 'settings' string for a DV.

    For most DVs this is just VARIABLE_REGISTRY[dv]['distribution'].
    A small set of pwPE/xprob variables need pipeline-specific 'hierarchial'
    suffixes that signal more MCMC iterations to nonsp_predictors.R; these are
    captured in the module-level _*_HIER_VARS sets above rather than duplicating
    distribution information that already lives in VARIABLE_REGISTRY.

    Raises ValueError for any DV not in VARIABLE_REGISTRY — add it to
    master_config.py rather than adding a hardcoded special case here.
    """
    if dv in _BINOMIAL_HIER_VARS:
        return 'binomial_hierarchial'
    if dv in _LOGNORMAL_HIER_VARS:
        return 'lognormal_hierarchial'
    if dv in _STUDENT_HIER_VARS:
        return 'student_t_hierarchial'
    reg = VARIABLE_REGISTRY.get(dv)
    if reg is None:
        raise ValueError(
            f'DV {dv!r} is not in VARIABLE_REGISTRY (master_config.py). '
            f'Add a registry entry with the correct distribution before running.'
        )
    dist = reg.get('distribution')
    if not isinstance(dist, str):
        raise ValueError(
            f'DV {dv!r} has distribution={dist!r} in VARIABLE_REGISTRY — '
            f'must be a non-null string. Fix master_config.py.'
        )
    return dist


def generate_nonsp_job_array(
    nonsp_predictors,
    df_foranalysis,
    modeltypes,
    dvs,
    walltime=NONSP_WALLTIME,
    dvs_to_exclude=None,
    dv_dict=None,
    bundle_by=3,
    replace=True,
    iterations=NONSP_ITERATIONS,
    output_basepath=str(OUTPUT_NONSP),
    bouchet_path=HPC_BASE,
    rscript_name='nonsp_predictors.R',
    output='nonsp_predictor_analyses',
    df_path='df_foranalysis_master.csv',
    df_name='df_foranalysis_master',
):
    """
    Generate a SLURM job-array .txt for nonsp single-path regressions.

    DV brms families are read from VARIABLE_REGISTRY['distribution'] via _dv_settings().
    Covariate lists come from covariates_master_dict (keyed by model type).

    If a model type's covariate list contains a variable matching '_nomic',
    'psychedelic_age', 'PC1', or 'lastuse_dose', that variable is replaced by the
    focal predictor (spvar replacement).  For model types without such a variable
    (e.g. nice_covariates = demographics only), the focal predictor is appended.

    Predictor _normalized columns do NOT exist in df_foranalysis — the master CSV
    ships raw and R creates them after subsetting (2026-08-22).  This function
    checks each predictor's RAW column is present and that the R normalizer will
    produce the matching _normalized name, and raises clearly if either fails.
    """
    _dv_to_lab = dv_dict if dv_dict else dict(dv_to_lab)
    # binomvars are handled by _dv_settings → 'binomial_hierarchial' but still
    # need labels added to _dv_to_lab so they appear in all_dvs.
    for var, lab in zip(
        list(_BINOMIAL_HIER_VARS),
        ['CH Rate', '75% Detection Probability', 'VCH Rate', '75% Detection Probability'],
    ):
        _dv_to_lab.setdefault(var, lab)

    # DV list — filter by requested dvs, then remove exclusions
    all_dvs = list(_dv_to_lab.keys())
    if dvs:
        all_dvs = [d for d in dvs if d in all_dvs]
    if dvs_to_exclude:
        all_dvs = [d for d in all_dvs if d not in dvs_to_exclude]

    print(f'  Predictors ({len(nonsp_predictors)}): {nonsp_predictors}')
    print(f'  Model types ({len(modeltypes)}): {modeltypes}')
    print(f'  DVs ({len(all_dvs)}): {all_dvs}')

    # Resolve each predictor to its _normalized column.
    # Since 2026-08-22 those columns do NOT exist in df_foranalysis — the master CSV
    # ships raw and R creates them after subsetting.  The check is therefore against
    # _WILL_EXIST_IN_R (the predicted output of gelman_normalize_df()) plus the raw
    # column, not against df_foranalysis.columns, which would reject everything.
    # Categorical predictors are kept as-is (R factors them, never normalizes them).
    nonsp_predictors_normalized = []
    for term in nonsp_predictors:
        if term in categorical_vars:
            nonsp_predictors_normalized.append(term)
            # Ensure categorical dtype is written correctly to CSV.
            df_foranalysis[term] = df[term].astype('category')
            continue
        base_term = term.replace('_normalized', '')
        norm_term = term if term.endswith('_normalized') else f'{base_term}_normalized'
        if base_term not in df_foranalysis.columns:
            raise ValueError(
                f'Predictor {base_term!r}: raw column not found in df_foranalysis. '
                f'R cannot create {norm_term!r} without it. '
                f'Similar columns: {_suggest_cols(base_term, df_foranalysis.columns)}'
            )
        if norm_term not in _WILL_EXIST_IN_R:
            raise ValueError(
                f'Predictor {base_term!r}: R will not create {norm_term!r}. '
                f'A predictor gets a normalized column only via rule 1 '
                f'(NEED_NON_NORMALIZED), rule 2 (INPLACE_NORMALIZED) or rule 4 '
                f'(PREDICTORS_FOR_JOBS).  Add {base_term!r} to one of the predictor '
                f'groups in iv_type_dict, or to VARIABLE_REGISTRY, in master_config.py.'
            )
        nonsp_predictors_normalized.append(norm_term)

    # sex_v2 must be present for R-side binary factor conversion.
    for col in ['sex_v2']:
        if col not in df_foranalysis.columns:
            df_foranalysis = df_foranalysis.merge(df[['record_id', col]], on='record_id', how='left')

    os.makedirs(output_basepath, exist_ok=True)
    df_foranalysis.to_csv(f'{output_basepath}/{df_name}.csv', index=False)

    rscript_source = str(_SCRIPT_DIR / rscript_name)
    if os.path.exists(rscript_source):
        shutil.copy2(rscript_source, output_basepath)
        print(f'  Copied {rscript_name} to {output_basepath}')

    dvs_bundled = [all_dvs[i:i+bundle_by] for i in range(0, len(all_dvs), bundle_by)]
    commands = []

    for predictor in nonsp_predictors_normalized:
        vartype = 'factor' if predictor in categorical_vars else 'numerical'
        base_predictor = predictor.replace('_normalized', '')

        for model in modeltypes:
            orig_covs = covariates_master_dict[model].copy()
            missing_covs = [c for c in orig_covs if c not in df_foranalysis.columns]
            if missing_covs:
                in_raw = [c for c in missing_covs if c in df.columns]
                raise ValueError(
                    f'Model "{model}" missing covariates: {missing_covs}. '
                    f'In raw df: {in_raw}. '
                    f'Suggestions: { {c: _suggest_cols(c, df_foranalysis.columns) for c in missing_covs} }'
                )
            # Spvar replacement: if the covariate list contains an SP variable
            # (identified by name pattern), replace it with the focal predictor.
            # For model types with no SP variable (e.g. nice_covariates), append.
            spvar_list = [x for x in orig_covs
                          if '_nomic' in x or 'psychedelic_age' in x
                          or 'PC1' in x or 'lastuse_dose' in x]
            if spvar_list:
                spvar = spvar_list[0]
                covs = ([x for x in orig_covs if x != spvar] + [predictor]) if replace else (orig_covs + [predictor])
            else:
                covs = orig_covs + [predictor]
            # Apply mo() wrapping for ordinal monotonic covariates
            # (source of truth: MONOTONIC_COVARIATES in modules/master_config.py).
            covs_mo = [f'mo({c})' if c in MONOTONIC_COVARIATES else c for c in covs]
            sp_and_covs = ' + '.join(covs_mo)

            for dv_bundle in dvs_bundled:
                cmd = _MODULE_PREFIX
                for dv in dv_bundle:
                    if dv == predictor or dv == base_predictor:
                        continue
                    settings = _dv_settings(dv)
                    longornot = 'long' if (dv in _LONG_RUN_VARS or base_predictor in _LONG_RUN_VARS) else 'short'
                    cmd += (f' && Rscript {rscript_name} {settings} {model} {dv} {iterations}'
                            f' "{sp_and_covs}" "{predictor}" "{longornot}" "{vartype}" "{df_path}"')
                if cmd != _MODULE_PREFIX:
                    commands.append(cmd)

    txt_path = f'{output_basepath}/{output}.txt'
    with open(txt_path, 'w') as f:
        f.write(f'# to run: module load dSQ; dsq --job-file {output}.txt --mem-per-cpu 4g -t {walltime} --mail-type ALL\n')
        for cmd in commands:
            f.write(cmd + '\n')

    print(f'\n  Generated {len(commands)} job commands → {os.path.abspath(txt_path)}')
    return txt_path


def combine_job_files(file_stems, output_stem, base_path=OUTPUT_NONSP):
    """Concatenate nonsp .txt files from base_path into one combined file."""
    combined = base_path / f'{output_stem}.txt'
    total = 0
    with open(combined, 'w') as out:
        out.write(f'# Combined job array: {output_stem}\n')
        for stem in file_stems:
            src = base_path / f'{stem}.txt'
            if src.exists():
                lines = [l for l in src.read_text().splitlines() if not l.startswith('#')]
                out.write(f'### {stem} ###\n')
                out.write('\n'.join(lines) + '\n')
                total += len(lines)
            else:
                print(f'  WARNING: {src} not found — skipping')
    print(f'  Combined {len(file_stems)} files → {combined} ({total} job lines)')
    return combined


# ==============================================================================
# MEDIATION FUNCTIONS  [replicates generate_mediation_job_arrays.ipynb Cells 2-3]
# ==============================================================================

# SP predictor shorthand → (column_name, sample_default)
SP_PREDICTORS = {
    'spage':      ('psychedelic_age',          'spusers'),
    'lifenomic':  ('psycheduse_life_nomic',    'spusers'),
    'avgdose':    ('avg_life_dose',            'spusers'),
    'pc1':        ('psychedelic_use_PC1',      'full'),
    'pc1ranked':  ('psychedelic_rank_use_PC1', 'full'),
}

# Mediator shorthand → (response_col, predictor_in_dv_col)
# response_col   : raw column name used in the mediator formula (DV of the first path)
# predictor_in_dv: column used as the mediator predictor in the DV formula
#   - need_non_normalized vars: response_col=raw, predictor_in_dv=_normalized
#     (raw is non-normal so mediator formula uses raw family; _normalized is used in DV formula)
#   - inplace_normalized vars: response_col == predictor_in_dv (normalized in-place, same col)
# mediator_family is derived at use-time from VARIABLE_REGISTRY[response_col]['distribution'].
# To update a mediator's family, change the distribution field in master_config.py.
MEDIATORS = {
    'vchthreshold': ('vch_threshold', 'vch_threshold'),             # inplace_normalized
    'vchrate':      ('vch_bl_yes_0',  'vch_bl_yes_0_normalized'),   # need_non_normalized
    'vchrate75':    ('vch_bl_yes_75', 'vch_bl_yes_75_normalized'),  # need_non_normalized
    'vchnu':        ('vch_nu',        'vch_nu_normalized'),          # need_non_normalized
    'vchbeta':      ('vch_beta',      'vch_beta'),                   # inplace_normalized
    'vchomega':     ('vch_omega',     'vch_omega'),                  # inplace_normalized
    'vchnu3lev':    ('vch_nu_3lev',   'vch_nu_3lev_normalized'),    # need_non_normalized
    'vchbeta3lev':  ('vch_beta_3lev', 'vch_beta_3lev_normalized'),  # need_non_normalized
    'vchnumat':     ('vch_short_psychedelic_bl_nu',   'vch_short_psychedelic_bl_nu_normalized'),
    'vchbetamat':   ('vch_short_psychedelic_bl_beta', 'vch_short_psychedelic_bl_beta_normalized'),
    'vchnuavg':     ('vch_nu_avg',     'vch_nu_avg_normalized'),      # need_non_normalized
    'vchnunominal':   ('vch_nu_nominal',   'vch_nu_nominal_normalized'),    # need_non_normalized
    'vchbetanominal': ('vch_beta_nominal', 'vch_beta_nominal'),            # inplace_normalized
    # SDT / metacognition mediators (sdt_hppd group)
    'dprimeoverall':    ('d_prime_overall',   'd_prime_overall'),            # inplace_normalized (student_t)
    'criterionoverall': ('criterion_overall', 'criterion_overall'),          # inplace_normalized (student_t)
    'meanconffas':      ('mean_conf_fas',     'mean_conf_fas_normalized'),   # need_non_normalized (gamma)
}



def _registry_distribution(col, role='variable'):
    """
    Look up distribution for col in VARIABLE_REGISTRY.  Raises ValueError with
    a clear message if the column is missing or its distribution is unset — so
    the user knows to update master_config.py rather than chasing a silent wrong family.
    """
    reg = VARIABLE_REGISTRY.get(col)
    if reg is None:
        raise ValueError(
            f'Mediation {role} {col!r} not found in VARIABLE_REGISTRY (master_config.py). '
            f'Add a registry entry with the correct distribution field before running.'
        )
    dist = reg.get('distribution')
    if not isinstance(dist, str):
        raise ValueError(
            f'Mediation {role} {col!r} has distribution={dist!r} in VARIABLE_REGISTRY — '
            f'must be a non-null string. Fix the distribution field in master_config.py.'
        )
    return dist

def _strip_r_side(cov_type):
    stem = cov_type.lower()
    for sfx in R_SIDE_SUFFIXES:
        stem = stem.replace(sfx, '')
    return stem


def get_covs(cov_type, spvar_col):
    """
    Return covariate formula string for a mediation model.
    Strips R-side subsetting suffixes before resolving base type.
    Supports: nice_covariates, main, age_only/age_control, univariate/true_univariate,
              empiric_covariates, empirical_covariates,
              drugs_month, drugs_trimmed_month.
              Append _beta to add vch_beta as control.

    Covariate lists are sourced from BASE_COVARIATE_SETS in modules/master_config.py.
    To change a covariate set, update it there — do NOT add local lists here.

    Variables in MONOTONIC_COVARIATES (master_config.py) are automatically wrapped
    in mo() for brms monotonic effects, unless they equal spvar_col.
    """
    base = _strip_r_side(cov_type)
    add_beta = False
    if base.endswith('_beta'):
        add_beta = True
        base = base[:-5]

    if base == 'nice_covariates':
        covs = [c for c in BASE_COVARIATE_SETS['nice_covariates'] if c != spvar_col]
    elif base == 'main':
        covs = [c for c in BASE_COVARIATE_SETS['main'] if c != spvar_col]
    elif base in ('age_only', 'age_control'):
        covs = ['age_v2'] if spvar_col != 'age_v2' else []
    elif base == 'empiric_covariates':
        covs = [c for c in BASE_COVARIATE_SETS['empiric_covariates'] if c != spvar_col]
    elif base == 'empirical_covariates':
        covs = [c for c in BASE_COVARIATE_SETS[base] if c != spvar_col]
        # Wrap ordinal covariates in mo() for brms monotonic effects.
        covs = [f'mo({c})' if c in MONOTONIC_COVARIATES else c for c in covs]
    elif base in ('univariate', 'true_univariate'):
        covs = []
    # Past-month drug-use sensitivity sets.  Neither contains an SP variable, so
    # the spvar_col filter is a no-op here — kept for consistency with the
    # branches above, which rely on it.
    elif base in ('drugs_month', 'drugs_trimmed_month'):
        covs = [c for c in BASE_COVARIATE_SETS[base] if c != spvar_col]
    # Hardware (display-class) control.  `nice_covariates` + the categorical
    # monitor_check_operationalized_final.  Reached from the model type
    # `nice_covariates_spusers_hardware_control` once _strip_r_side() has
    # removed `_spusers`.  Contains no SP variable, so the spvar_col filter is a
    # no-op — kept for consistency with the branches above.
    elif base == 'nice_covariates_hardware_control':
        covs = [c for c in BASE_COVARIATE_SETS[base] if c != spvar_col]
    # 'no_outlier' / 'nooutlier' are not supported — they mapped an outlier
    # KEYWORD onto an empty covariate set, which was always a coincidence of
    # naming rather than a covariate decision.  A variant naming either now
    # raises below rather than silently fitting a univariate model.
    else:
        raise ValueError(f'Unknown cov_type: {cov_type!r} (base parsed as {base!r})')

    if add_beta and spvar_col != 'vch_beta' and 'vch_beta' not in covs:
        covs = covs + ['vch_beta']
    return ' + '.join(covs) if covs else ''


_med_generated_stems = []

_MED_MODULE_LOAD = ('module purge && module load foss/2022b && '
                    'module load R/4.4.1-foss-2022b && export R_LIBS_USER=$HOME/R/4.4')


def generate_mediation_jobs(
    spvar_short,
    mediator_short,
    dv_short,
    cov_types=('nice_covariates',),
    stem=None,
    iter_=MED_ITER,
    warmup=MED_WARMUP,
    verbose=True,
):
    """
    Generate a .txt job file for one (spvar, mediator, dv) combination.

    mediator_family is derived from VARIABLE_REGISTRY[response_col]['distribution'].
    dv_family       is derived from VARIABLE_REGISTRY[dv_col]['distribution'].
    To change a family, update the distribution field in master_config.py.
    """
    spvar_col, _              = SP_PREDICTORS[spvar_short]
    med_response, med_in_dv   = MEDIATORS[mediator_short]
    dv_col                    = DVS[dv_short]
    med_family                = _registry_distribution(med_response, role='mediator')
    dv_family                 = _registry_distribution(dv_col,       role='dv')

    if stem is None:
        stem = f'{dv_short}_{spvar_short}_{mediator_short}'

    hpc_df_path = f'{HPC_MED_DIR}/df_foranalysis_master.csv'

    job_lines, model_names = [], []

    for cov_type_raw in cov_types:
        cov_type  = cov_type_raw.lower()
        covs_str  = get_covs(cov_type, spvar_col)
        model_name = f'{dv_col}_{spvar_short}_{mediator_short}_{cov_type}'

        hpc_results = f'{HPC_MED_RESULTS}/{dv_col}/mediation_models'

        cmd = (
            f'{_MED_MODULE_LOAD} && '
            f'cd {HPC_MED_DIR} && '
            f'Rscript {HPC_MED_DIR}/hpc_mediation.R '
            f'"{spvar_col}" '
            f'"{med_response}" '
            f'"{med_in_dv}" '
            f'"{dv_col}" '
            f'"{covs_str}" '
            f'"{med_family}" '
            f'"{dv_family}" '
            f'"{model_name}" '
            f'"{hpc_df_path}" '
            f'"{hpc_results}" '
            f'"{HPC_MED_HELPERS}" '
            f'"{iter_}" '
            f'"{warmup}"'
        )
        job_lines.append(cmd)
        model_names.append(model_name)

    if not job_lines:
        print(f'  WARNING: No jobs generated for stem={stem}')
        return None

    OUTPUT_MED.mkdir(parents=True, exist_ok=True)
    txt_path = OUTPUT_MED / f'{stem}.txt'
    txt_path.write_text('\n'.join(job_lines) + '\n')
    _med_generated_stems.append(stem)

    if verbose:
        print(f'\n  {stem} ({len(job_lines)} jobs)')
        for m in model_names:
            print(f'    {m}')
        print(f'  Written: {txt_path}')

    return txt_path


def combine_mediation_jobs(stems, combined_stem):
    """Concatenate mediation .txt files from OUTPUT_MED into one file."""
    combined_path = OUTPUT_MED / f'{combined_stem}.txt'
    lines = []
    for stem in stems:
        p = OUTPUT_MED / f'{stem}.txt'
        if p.exists():
            lines.extend(l for l in p.read_text().splitlines() if l.strip())
        else:
            print(f'  WARNING: {p} not found — skipping')
    combined_path.write_text('\n'.join(lines) + '\n')
    print(f'  Combined {len(lines)} mediation jobs → {combined_path}')
    return combined_path


# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

MASTER_DF_PATH = 'df_foranalysis_master.csv'
MASTER_DF_NAME = 'df_foranalysis_master'


# ── Low-cardinality safety check ─────────────────────────────────────────────
# Any numeric column with ≤3 observed levels that is NOT explicitly declared as
# a categorical variable is almost certainly a multi-level factor being passed as
# continuous — which will silently produce wrong brms results.  Check only columns
# that are actually passed to brms (covariates + predictors), not every raw column.
print('\nLow-cardinality check on model covariates and predictors...')
_cols_to_check = (
    finalcovs                         # all covariate columns across all model types
    | set(PREDICTORS_FOR_JOBS)        # all predictor columns passed to R
) - CATEGORICAL_COVARIATES - _ORDINAL_AS_CONTINUOUS  # exempt declared categoricals + known ordinals
_low_card_problems = []
for _col in sorted(_cols_to_check):
    if _col not in df_foranalysis.columns:
        continue
    if not pd.api.types.is_numeric_dtype(df_foranalysis[_col]):
        continue  # non-numeric columns — not passed to brms as continuous
    _n_unique = df_foranalysis[_col].dropna().nunique()
    if _n_unique <= 3:
        _low_card_problems.append((_col, _n_unique))

if _low_card_problems:
    print('\n' + '!' * 70)
    print('CRITICAL: The following numeric columns have ≤3 observed levels')
    print('and are NOT in CATEGORICAL_COVARIATES.  They may be multi-level')
    print('factors being treated as continuous — brms results would be WRONG.')
    print()
    for _col, _n in _low_card_problems:
        _vals = sorted(df_foranalysis[_col].dropna().unique())
        print(f'  {_col!r}: {_n} levels → {_vals}')
    print()
    print('Fix: add the variable to CATEGORICAL_COVARIATES (and to')
    print('categorical_factor_vars in nonsp_predictors.R / hpc_mediation.R),')
    print('OR confirm it is genuinely binary/continuous and suppress this check')
    print('by adding it to CATEGORICAL_COVARIATES with a comment.')
    print('!' * 70 + '\n')
    _answer = input('Continue anyway? [y/N] ').strip().lower()
    if _answer != 'y':
        raise SystemExit('Aborted — fix categorical variable declarations before saving master CSV.')
else:
    print('  OK — no unexpected low-cardinality columns found.')

# Save nonsp master CSV
df_foranalysis.to_csv(str(OUTPUT_NONSP / MASTER_DF_PATH), index=False)
print(f'\nNonsp master CSV saved: {OUTPUT_NONSP / MASTER_DF_PATH} {df_foranalysis.shape}')

# ── Write categorical_factor_vars.R ──────────────────────────────────────────
# Auto-generated from modules/master_config.py :: CATEGORICAL_FACTOR_VARS.
# Both nonsp_predictors.R and hpc_mediation.R source this file so there is
# exactly ONE place to edit the categorical variable list.
_cat_r_lines = [
    '# Auto-generated by generate_hpc_jobs.py — DO NOT EDIT MANUALLY.',
    '# Edit CATEGORICAL_FACTOR_VARS in modules/master_config.py instead.',
    'binary_factor_vars <- c(',
]
for _v in CATEGORICAL_FACTOR_VARS:
    _cat_r_lines.append(f'  {_v!r},')
# Remove trailing comma from last entry
_cat_r_lines[-1] = _cat_r_lines[-1].rstrip(',')
_cat_r_lines.append(')')
_cat_r_content = '\n'.join(_cat_r_lines) + '\n'
for _dest_dir in [OUTPUT_NONSP, OUTPUT_MED]:
    _dest_dir.mkdir(parents=True, exist_ok=True)
    _cat_r_path = _dest_dir / 'categorical_factor_vars.R'
    _cat_r_path.write_text(_cat_r_content)
    print(f'categorical_factor_vars.R written → {_cat_r_path}')

# ── Write monotonic_covariates_vars.R ─────────────────────────────────────────
# Auto-generated from modules/master_config.py :: MONOTONIC_COVARIATES.
# Sourced by both nonsp_predictors.R and hpc_mediation.R; variables listed here
# are converted to as.ordered() before model fitting and wrapped in mo() in the
# formula string.  Edit MONOTONIC_COVARIATES in master_config.py — not here.
_mono_r_lines = [
    '# Auto-generated by generate_hpc_jobs.py — DO NOT EDIT MANUALLY.',
    '# Edit MONOTONIC_COVARIATES in modules/master_config.py instead.',
    'monotonic_covariate_vars <- c(',
]
for _v in sorted(MONOTONIC_COVARIATES):
    _mono_r_lines.append(f'  {_v!r},')
_mono_r_lines[-1] = _mono_r_lines[-1].rstrip(',')
_mono_r_lines.append(')')
_mono_r_content = '\n'.join(_mono_r_lines) + '\n'
for _dest_dir in [OUTPUT_NONSP, OUTPUT_MED]:
    _dest_dir.mkdir(parents=True, exist_ok=True)
    _mono_r_path = _dest_dir / 'monotonic_covariates_vars.R'
    _mono_r_path.write_text(_mono_r_content)
    print(f'monotonic_covariates_vars.R written → {_mono_r_path}')


# ── Write normalization_vars.R ────────────────────────────────────────────────
# Auto-generated from modules/master_config.py (VARIABLE_REGISTRY →
# NEED_NON_NORMALIZED / INPLACE_NORMALIZED, BASE_COVARIATE_SETS → finalcovs) and
# from the predictor groups configured at the top of this script.
#
# This file carries the MEMBERSHIP of each normalization rule; the arithmetic
# lives in gelman_normalization.R.  Both nonsp_predictors.R and hpc_mediation.R
# source both files, so there is exactly one place to edit each concern.
# Edit master_config.py — never this generated file, and never a variable list
# inside an R script.
def _r_char_vector(name, values, comment):
    """Render a Python list as an R character vector assignment."""
    lines = [f'# {comment}', f'{name} <- c(']
    if values:
        for v in values:
            lines.append(f'  {v!r},')
        lines[-1] = lines[-1].rstrip(',')
    lines.append(')')
    return '\n'.join(lines)


_norm_r_content = '\n'.join([
    '# Auto-generated by generate_hpc_jobs.py — DO NOT EDIT MANUALLY.',
    '# Edit VARIABLE_REGISTRY / BASE_COVARIATE_SETS in modules/master_config.py instead.',
    '#',
    '# Consumed by gelman_normalize_df() in gelman_normalization.R.  The rule numbers',
    '# match the four rules documented there and in data_prep.py :: normalize_analysis_df().',
    '',
    _r_char_vector('need_non_normalized_vars', list(NEED_NON_NORMALIZED),
                   'Rule 1 — raw column preserved (family constrains the scale); '
                   '{col}_normalized created alongside.'),
    '',
    _r_char_vector('inplace_normalized_vars', list(INPLACE_NORMALIZED),
                   'Rule 2 — raw column normalized in place; {col}_normalized aliases it.'),
    '',
    _r_char_vector('normalize_in_place_vars', NORMALIZE_IN_PLACE_VARS,
                   'Rule 3 — every covariate across every model type, plus the main '
                   'predictors; normalized in place. Categorical and monotonic '
                   'variables are excluded by gelman_normalize_df().'),
    '',
    _r_char_vector('ensure_normalized_copy_vars', ENSURE_NORMALIZED_COPY_VARS,
                   'Rule 4 — every predictor any job array can request; guaranteed a '
                   '{col}_normalized column.'),
    '',
    _r_char_vector('gelman_skip_vars', GELMAN_SKIP_VARS,
                   'Never scaled even when named by a rule. Mirrors BINARY_FACTOR_VARS '
                   'in modules/data_prep.py (a shorter list than CATEGORICAL_FACTOR_VARS '
                   '— see the NORMALIZATION section of generate_hpc_jobs.py).'),
    '',
]) + '\n'

for _dest_dir in [OUTPUT_NONSP, OUTPUT_MED]:
    _dest_dir.mkdir(parents=True, exist_ok=True)
    _norm_r_path = _dest_dir / 'normalization_vars.R'
    _norm_r_path.write_text(_norm_r_content)
    print(f'normalization_vars.R written → {_norm_r_path}')

# ── Copy the hand-written shared R files ──────────────────────────────────────
# Hand-written (not generated): it holds the transform, not the variable lists.
# Both R execution scripts source it from their own directory, so each job bundle
# needs its own copy — compute nodes only see what is rsync'd to the cluster.
# Edit 03_hpc/gelman_normalization.R, never a copy under output/.
# convergence_gate.R travels the same way: hand-written, sourced by both R
# scripts, one copy per job bundle. It holds the Rhat / ESS / divergence
# thresholds that decide whether a model may write its summary tables.
for _shared_name in ('gelman_normalization.R', 'convergence_gate.R'):
    _shared_src = _SCRIPT_DIR / _shared_name
    if not _shared_src.exists():
        raise FileNotFoundError(
            f'{_shared_src} not found. Both R scripts source it; without it '
            f'every job fails at startup.'
        )
    for _dest_dir in [OUTPUT_NONSP, OUTPUT_MED]:
        shutil.copy2(str(_shared_src), str(_dest_dir / _shared_name))
        print(f'{_shared_name} copied → {_dest_dir / _shared_name}')



print(f'\nPart 1 predictors ({len(ALL_PREDS)}): {ALL_PREDS}')
print(f'Part 1 DVs: {HPPD_CAPS_DVS}')
print(f'Part 2 predictors ({len(PART2_PREDS)}): {PART2_PREDS}')
print(f'Part 2 DVs ({len(VCH_AS_DVS)}): {VCH_AS_DVS}')

nonsp_stems = []   # track all nonsp .txt stems generated

# ── Part 1: predictors → persist_vis_yn + caps_vision ────────────────────────
if INCLUDE_HPPD_CAPS_JOBS:
    print(f'\n{"="*60}')
    print(f'PART 1: {ALL_PREDS} → {HPPD_CAPS_DVS}')
    print(f'{"="*60}')
    generate_nonsp_job_array(
        nonsp_predictors=ALL_PREDS,
        modeltypes=MODEL_VARIANTS,
        df_foranalysis=df_foranalysis.copy(),
        dvs=HPPD_CAPS_DVS,
        walltime=NONSP_WALLTIME,
        replace=True,
        iterations=NONSP_ITERATIONS,
        output='master_allpreds_x_hppd_caps',
        df_path=MASTER_DF_PATH,
        df_name=MASTER_DF_NAME,
    )
    nonsp_stems.append('master_allpreds_x_hppd_caps')

# ── Part 2: SP predictors → VCH behavior + VCH comps (as DVs) ────────────────
if INCLUDE_SP_VCH_JOBS:
    print(f'\n{"="*60}')
    print(f'PART 2: {PART2_PREDS} → VCH DVs as DVs')
    print(f'{"="*60}')
    generate_nonsp_job_array(
        nonsp_predictors=PART2_PREDS,
        modeltypes=MODEL_VARIANTS,
        df_foranalysis=df_foranalysis.copy(),
        dvs=VCH_AS_DVS,
        walltime=NONSP_WALLTIME,
        replace=True,
        iterations=NONSP_ITERATIONS,
        output='master_sp_x_vch_dvs',
        df_path=MASTER_DF_PATH,
        df_name=MASTER_DF_NAME,
    )
    nonsp_stems.append('master_sp_x_vch_dvs')

# ── Mediation jobs ────────────────────────────────────────────────────────────
med_stems = []

if INCLUDE_MEDIATION_JOBS:
    print(f'\n{"="*60}')
    print('MEDIATION: augmenting df + generating job files')
    print(f'{"="*60}')

    OUTPUT_MED.mkdir(parents=True, exist_ok=True)

    # Validate that every column the configured MEDIATION_ANALYSES will ask for
    # either exists in the raw master df or will be created by the R normalizer.
    # The raw mediator column must be in the CSV; its _normalized counterpart is
    # produced on the cluster by gelman_normalize_df(), so it is checked against
    # _WILL_EXIST_IN_R rather than df_foranalysis.columns.
    _med_col_pairs = {MEDIATORS[spec['mediator']] for spec in MEDIATION_ANALYSES}

    _required_med_raw  = {raw  for raw, _  in _med_col_pairs}
    _required_med_norm = {norm for _,  norm in _med_col_pairs}
    _required_med_dvs  = {DVS[spec['dv']] for spec in MEDIATION_ANALYSES}

    _must_be_in_csv = sorted(_required_med_raw | _required_med_dvs | {'record_id'})
    missing_med = [c for c in _must_be_in_csv if c not in df_foranalysis.columns]
    if missing_med:
        raise KeyError(f'Mediation master df missing required raw columns: {missing_med}')

    _missing_med_norm = sorted(c for c in _required_med_norm if c not in _WILL_EXIST_IN_R)
    if _missing_med_norm:
        raise KeyError(
            f'Mediation mediator column(s) {_missing_med_norm} will not be created by '
            f'the R normalizer.  A {{col}}_normalized column comes from rule 1 '
            f'(NEED_NON_NORMALIZED), rule 2 (INPLACE_NORMALIZED) or rule 4 '
            f'(PREDICTORS_FOR_JOBS) — check the mediator entry in VARIABLE_REGISTRY '
            f'and the predictor groups in master_config.py.'
        )

    # Save augmented mediation master CSV
    med_master_path = OUTPUT_MED / 'df_foranalysis_master.csv'
    df_foranalysis.to_csv(str(med_master_path), index=False)
    print(f'\nMediation master CSV saved: {med_master_path}')

    # Legacy aliases (byte-identical copies)
    # These four names predate the single-master-CSV design; every one is a
    # byte-identical copy kept only so an old job line still finds its df.  The
    # the two '_nooutlier' names are dead (the keyword is not supported)
    # but are left in place: deleting them would break a stale job file loudly at
    # read time instead of at the keyword check, which is the worse failure.
    for alias in ['df_foranalysis_full.csv', 'df_foranalysis_spusers.csv',
                  'df_foranalysis_full_nooutlier.csv', 'df_foranalysis_spusers_nooutlier.csv']:
        shutil.copy2(str(med_master_path), str(OUTPUT_MED / alias))
    print(f'  Legacy aliases copied.')

    # Copy hpc_mediation.R + helper_scripts into OUTPUT_MED.
    # Both live alongside this generator in 03_hpc/; the cluster job bundle needs
    # its own copy because compute nodes only see what is rsync'd to HPC_MED_DIR.
    _med_rscript_src = _SCRIPT_DIR / 'hpc_mediation.R'
    shutil.copy2(str(_med_rscript_src), str(OUTPUT_MED / 'hpc_mediation.R'))
    _helpers_src  = _SCRIPT_DIR / 'helper_scripts'
    _helpers_dest = OUTPUT_MED / 'helper_scripts'
    if _helpers_dest.exists():
        shutil.rmtree(str(_helpers_dest))
    shutil.copytree(str(_helpers_src), str(_helpers_dest))
    print(f'  hpc_mediation.R + helper_scripts copied into {OUTPUT_MED}')

    # Generate individual mediation job files
    for spec in MEDIATION_ANALYSES:
        stem = spec.get('stem')
        generate_mediation_jobs(
            spvar_short    = spec['spvar'],
            mediator_short = spec['mediator'],
            dv_short       = spec['dv'],
            cov_types      = spec['cov_types'],
            stem           = stem,
            verbose        = True,
        )

    # Combine all mediation jobs into one file in OUTPUT_MED
    med_combined_path = combine_mediation_jobs(_med_generated_stems, 'mediation_combined')
    med_stems = _med_generated_stems
    n_med_jobs = sum(1 for l in med_combined_path.read_text().splitlines() if l.strip())
    print(f'\n  Total mediation jobs: {n_med_jobs}')


# ==============================================================================
# COMBINE ALL INTO ONE MASTER JOB ARRAY
# ==============================================================================

print(f'\n{"="*60}')
print(f'COMBINING ALL JOBS → {COMBINED_STEM}.txt')
print(f'{"="*60}')

combined_path = OUTPUT_NONSP / f'{COMBINED_STEM}.txt'
all_lines = []

# Collect nonsp lines
for stem in nonsp_stems:
    src = OUTPUT_NONSP / f'{stem}.txt'
    if src.exists():
        lines = [l for l in src.read_text().splitlines() if l.strip() and not l.startswith('#')]
        all_lines.extend(lines)
        print(f'  + {stem}.txt  ({len(lines)} lines)')
    else:
        print(f'  WARNING: {src} not found — skipping')

# Collect mediation lines
if INCLUDE_MEDIATION_JOBS and med_stems:
    med_combined = OUTPUT_MED / 'mediation_combined.txt'
    if med_combined.exists():
        lines = [l for l in med_combined.read_text().splitlines() if l.strip()]
        all_lines.extend(lines)
        print(f'  + mediation_combined.txt  ({len(lines)} lines)')

# Determine walltime for combined submission
combined_walltime = MED_WALLTIME if INCLUDE_MEDIATION_JOBS else NONSP_WALLTIME

with open(combined_path, 'w') as f:
    f.write(f'# Combined job array: {COMBINED_STEM}\n')
    f.write(f'# to run: module load dSQ; dsq --job-file {COMBINED_STEM}.txt '
            f'--mem-per-cpu 4g -t {combined_walltime} --mail-type ALL\n')
    for line in all_lines:
        f.write(line + '\n')

print(f'\n  {len(all_lines)} total job lines → {combined_path}')


# ==============================================================================
# PRINT TRANSFER + HPC COMMANDS  (semicolon-separated for copy-paste)
# ==============================================================================

local_nonsp    = os.path.abspath(OUTPUT_NONSP)
local_med      = os.path.abspath(OUTPUT_MED)
local_results  = str(_PROJECT_ROOT / 'results')

transfer_cmds = [f'rsync -av --progress {local_nonsp}/ {HPC_TRANSFER}:{HPC_BASE}/']
if INCLUDE_MEDIATION_JOBS:
    transfer_cmds.append(f'rsync -av --progress {local_med}/ {HPC_TRANSFER}:{HPC_MED_DIR}/')

hpc_cmds = [
    f'cd {HPC_BASE}',
    'module load dSQ',
    f'dsq --job-file {COMBINED_STEM}.txt --mem-per-cpu 4g -t {combined_walltime} --mail-type ALL',
]

print(f'\n{"="*70}')
print('TRANSFER COMMANDS  (run locally — copy-paste all in one go):')
print(f'{"="*70}')
print(' ; \\\n  '.join(transfer_cmds))

print(f'\n{"="*70}')
print(f'HPC COMMANDS  (run on {HPC_LOGIN} after transfer):')
print(f'{"="*70}')
print(f'ssh {HPC_LOGIN}')
print(' ; \\\n  '.join(hpc_cmds))
print()
print(f'  ↳ dSQ will print the sbatch command to submit — run that to launch the array.')
print(f'{"="*70}')

print(f'\n{"="*70}')
print('RETRIEVE COMMANDS  (run locally after jobs finish — copy-paste all in one go):')
print(f'{"="*70}')

retrieve_cmds = [
    f'rsync -av --progress {HPC_TRANSFER}:{HPC_BASE}/results/ {local_nonsp}/results/',
]
if INCLUDE_MEDIATION_JOBS:
    # Pulls directly into results/ at the repo root — where 0X_all_figures.py reads from.
    retrieve_cmds.append(
        f'rsync -av --progress {HPC_TRANSFER}:{HPC_MED_RESULTS}/ {local_results}/'
    )
print(' ; \\\n  '.join(retrieve_cmds))
print(f'{"="*70}')

if PROMPT_TO_COMPILE:
    import subprocess

    _nonsp_compile  = _SCRIPT_DIR / 'compile_nonsp_results.py'
    _med_compile    = _SCRIPT_DIR / 'compile_mediation_results.py'

    run_nonsp = INCLUDE_HPPD_CAPS_JOBS or INCLUDE_SP_VCH_JOBS
    run_med   = INCLUDE_MEDIATION_JOBS

    if run_nonsp:
        # Derive compile filters from this script's own CONFIG — no external config file needed.
        # DVs and predictors reflect exactly what was just submitted.
        _compile_dvs = list(dict.fromkeys(
            (HPPD_CAPS_DVS if INCLUDE_HPPD_CAPS_JOBS else [])
            + (VCH_AS_DVS  if INCLUDE_SP_VCH_JOBS    else [])
        ))
        _compile_preds = list(dict.fromkeys(
            f'{p}_normalized' for p in (
                (ALL_PREDS    if INCLUDE_HPPD_CAPS_JOBS else [])
                + (PART2_PREDS if INCLUDE_SP_VCH_JOBS    else [])
            )
        ))
        print(f'\nRun compile_nonsp_results.py now to pull + compile results from Bouchet?')
        print(f'  (Requires an active SSH ControlMaster: ssh -MNf bouchet)')
        print(f'  DVs ({len(_compile_dvs)}):         {_compile_dvs}')
        print(f'  Model types ({len(MODEL_VARIANTS)}): {MODEL_VARIANTS}')
        _answer = input('  [y/N] ').strip().lower()
        if _answer == 'y':
            _compile_cmd = [
                sys.executable, str(_nonsp_compile),
                '--dvs',         *_compile_dvs,
                '--predictors',  *_compile_preds,
                '--model-types', *MODEL_VARIANTS,
                '--no-confirm',   # skip redundant DUO prompt — user already confirmed above
            ]
            subprocess.run(_compile_cmd, check=False)
        else:
            _dv_str    = ' '.join(_compile_dvs)
            _pred_str  = ' '.join(_compile_preds)
            _model_str = ' '.join(MODEL_VARIANTS)
            print(f'  Skipped. Run manually:')
            print(f'    python {_nonsp_compile.name} \\')
            print(f'      --dvs {_dv_str} \\')
            print(f'      --predictors {_pred_str} \\')
            print(f'      --model-types {_model_str}')

    if run_med:
        print(f'\nRun compile_mediation_results.py subset now to pull mediation results from Bouchet?')
        print(f'  (Reads its model list from generate_hpc_jobs.py CONFIG)')
        print(f'  (Requires an active SSH ControlMaster: ssh -MNf bouchet)')
        _answer = input('  [y/N] ').strip().lower()
        if _answer == 'y':
            subprocess.run([sys.executable, str(_med_compile), 'subset'], check=False)
        else:
            print(f'  Skipped. Run manually: python {_med_compile.name} subset')

print(f'\nDone. {len(all_lines)} jobs ready in {combined_path}')
