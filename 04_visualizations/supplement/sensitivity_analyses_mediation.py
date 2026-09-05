#!/usr/bin/env /usr/local/bin/python3.12
"""
sensitivity_analyses_mediation.py
==================================
Supplementary Figure S5 — the compound NIE sensitivity heatmap for the
mediation models.

Aggregates the fitted mediation models into one long-format CSV, then draws the
natural indirect effect across covariate specifications for both DVs. It
submits no jobs: 03_hpc/generate_hpc_jobs.py fits every model type listed in
its CUSTOM_MED_TYPES, which must stay in step with SENSITIVITY_MED_TYPES here.
A type named here but not there is never fitted, and shows up as a permanently
missing heatmap column with no error at either end.

Canonical reference column (leftmost in each heatmap):
    nice_covariates_spusers   — the primary mediation model

Stages
------
    compile    — read path/mc/summary CSVs from results/{dv}/mediation_models/,
                 pulling anything missing from Bouchet over SSH (needs an open
                 ControlMaster). Skips the pull when everything is present
                 locally. Writes compiled_sensitivity_mediation.csv.
    heatmap    — draws one compound heatmap per entry in PATHS_TO_PLOT.

PATHS_TO_PLOT is ['NIE'] — the supplementary figure. The A, B and C' path
heatmaps come off the same compiled CSV through the same plotting function;
add them back to that list to render them.

Usage (from anywhere)
---------------------
    # figure only, from the CSV already on disk
    /usr/local/bin/python3.12 .../sensitivity_analyses_mediation.py heatmap

    # re-aggregate the fitted models first, then draw
    /usr/local/bin/python3.12 .../sensitivity_analyses_mediation.py all

IMPORTANT: Keep in parallel with sensitivity_analyses.py (single-path).
Any change to covariate sets, DHARMa flagging thresholds, or heatmap style must
be replicated in both. See 04_visualizations/supplement/README.md.
"""

# ==============================================================================
# ████████████████████████  CONFIG — edit here only  ███████████████████████████
# ==============================================================================

# ── Sensitivity covariate types ───────────────────────────────────────────────
# All types must be valid keys in BASE_COVARIATE_SETS (master_config.py) after
# stripping R-side subsetting suffixes (e.g. _spusers, _nooutlier, etc.).
# Every type here must also be in CUSTOM_MED_TYPES in 03_hpc/generate_hpc_jobs.py,
# which is what fits them. A type listed here but not there is never submitted,
# and shows up as a permanently missing heatmap column with no error at either end.
SENSITIVITY_MED_TYPES = [
    'empirical_covariates_spusers',
    'nice_covariates',            # full sample; promoted to SECOND_CANONICAL slot in heatmap
    'nice_covariates_spusers_iqr', # SP users, IQR outlier removed; → THIRD_CANONICAL slot
    'age_control_spusers',
    'true_univariate_spusers',
    'nice_covariates_beta_spusers', # beta-control: nice_covariates + vch_beta, SP users only
    'nice_covariates_spusers_nonan_caps', # caps_bl_1 non-NaN subset: participants with CAPS data
    'drugs_month_spusers',          # + alc/ghb/opioids/mj/atypicals/stimulants binaries
    'drugs_trimmed_month_spusers',  # + depressants/mj/stimulants binaries
    'nice_covariates_spusers_hardware_control',  # + monitor_check_operationalized_final
]

# ── Canonical columns ─────────────────────────────────────────────────────────
# Three-level hierarchy matches sensitivity_analyses.py column ordering:
#   CANONICAL  → SP users, nice covariates      (primary results)
#   SECOND     → full sample, nice covariates   (generalizability check)
#   THIRD      → SP users, IQR outlier removed  (robustness check)
# Set any to None to omit that column.
CANONICAL_MED_TYPE        = 'nice_covariates_spusers'
SECOND_CANONICAL_MED_TYPE = 'nice_covariates'
THIRD_CANONICAL_MED_TYPE  = 'nice_covariates_spusers_iqr'

# ── Mediation analyses to run ─────────────────────────────────────────────────
# Each dict specifies one (spvar, mediator, dv) combination.
# spvar keys: 'spage', 'avgdose', 'lifenomic', 'pc1', 'pc1ranked'
# mediator keys: 'vchthreshold', 'vchrate', 'vchbeta', 'vchnu', 'vchomega'
# dv keys: 'hppd_binary', 'caps_vision', 'caps_total', 'lshs_total', 'baggot_total'
_MED_MEDIATORS = ['vchthreshold', 'vchrate', 'vchbeta', 'vchnu']

MED_ANALYSES = [
    # psychedelic_age → VCH mediator → HPPD history
    *[{'spvar': 'spage',   'mediator': m, 'dv': 'hppd_binary'} for m in _MED_MEDIATORS],
    # avg_life_dose → VCH mediator → CAPS vision
    *[{'spvar': 'avgdose', 'mediator': m, 'dv': 'caps_vision'} for m in _MED_MEDIATORS],
]

# ── Paths to visualize ────────────────────────────────────────────────────────
# A path  = SP predictor → mediator (response scale)
# B path  = mediator → DV (response scale, spvar at mean)
# C' path = direct SP → DV (mediator fixed at mean)
# NIE     = natural indirect effect (MC integration)
# Only the NIE compound heatmap is a supplementary figure. The A / B / C' path
# heatmaps render from the same compiled CSV and the same plotting function, so
# add them back here if they are ever wanted again.
PATHS_TO_PLOT = ['NIE']

# Output filename per path. Only NIE has a supplementary figure number; any
# other path added to PATHS_TO_PLOT falls back to a descriptive stem.
PATH_FILE_STEMS = {'NIE': 'supplementary_figure_s5'}

# ── Stage toggles ─────────────────────────────────────────────────────────────
# Jobs are submitted by 03_hpc/generate_hpc_jobs.py, which fits every model type
# in CUSTOM_MED_TYPES; this script only aggregates the fitted models and draws
# the figure.
RUN_COMPILE        = True
RUN_HEATMAP        = True

# ── Output stems ──────────────────────────────────────────────────────────────
MED_SENSITIVITY_STEM = 'sensitivity_analyses_mediation'   # job-array .txt stem
COMPILED_CSV_STEM    = 'compiled_sensitivity_mediation.csv'

# ── Heatmap ───────────────────────────────────────────────────────────────────
DHARMA_ALPHA    = 0.05    # p < this → DHARMa flag (red border)
RHAT_MAX        = 1.01    # Rhat > this → MCMC flag (red border)
DIV_MAX         = 0       # divergents > this → MCMC flag

# (cov_type, dv) pairs to gray out — logically inapplicable.
# nonan_caps restricts to non-NaN caps_bl_1, which is the full sample for
# any caps_* DV, making it redundant with the primary model.
INAPPLICABLE_CELLS = {
    ('nice_covariates_spusers_nonan_caps', 'caps_vision'),
}

# ── MCMC settings ─────────────────────────────────────────────────────────────
MED_ITER   = 10000
MED_WARMUP = 6000
MED_WALLTIME = '360:00'   # per job (mediation models can run ~22 min for caps_vision)

# ── HPC paths ─────────────────────────────────────────────────────────────────

# ── Y-axis display labels ─────────────────────────────────────────────────────
# Used in heatmap row labels. Fallback: spvar/mediator shorthand is used as-is.
_SPVAR_DISPLAY = {
    'spage':     'Psych. age',
    'avgdose':   'Avg. dose',
    'lifenomic': 'Life uses',
    'pc1':       'PC1',
    'pc1ranked': 'PC1 (ranked)',
}
_MEDIATOR_DISPLAY = {
    'vchthreshold': 'VCH threshold',
    'vchrate':      'VCH rate (0%)',
    'vchrate75':    'VCH rate (75%)',
    'vchbeta':      'VCH \u03b2',
    'vchnu':        'VCH \u03bd',
    'vchomega':     'VCH \u03c9',
}

# ==============================================================================
# END CONFIG ───────────────────────────────────────────────────────────────────
# ==============================================================================


# ==============================================================================
# SETUP
# ==============================================================================

import os
import shutil
import sys
import tarfile
import tempfile
import warnings
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')

# ── CLI mode ──────────────────────────────────────────────────────────────────
_MODE_ARG = sys.argv[1].lower() if len(sys.argv) > 1 else None
if _MODE_ARG and _MODE_ARG not in ('compile', 'heatmap', 'all'):
    print(f"Unknown mode {_MODE_ARG!r}. Valid: compile heatmap all")
    sys.exit(1)
if _MODE_ARG:
    RUN_COMPILE        = _MODE_ARG in ('compile', 'all')
    RUN_HEATMAP        = _MODE_ARG in ('heatmap',  'all')

# ── Paths ─────────────────────────────────────────────────────────────────────
_SCRIPT_DIR   = Path(__file__).resolve().parent    # 04_visualizations/supplement/
_PROJECT_ROOT = _SCRIPT_DIR.parent.parent          # hppd_manuscript_public/

sys.path.insert(0, str(_PROJECT_ROOT / 'modules'))

from master_config import POINT_ESTIMATE_COL, MC_EFFECT_POINT_ESTIMATE_COL, point_estimate

# Shared staging dir — same as primary mediation pipeline.
# generate_mediation_sensitivity_jobs() only adds the .txt job file here;
# hpc_mediation.R, helper_scripts/, and df_foranalysis_master.csv are assumed
# to already exist (placed by generate_hpc_jobs.py).
OUTPUT_MED_SENSITIVITY = _PROJECT_ROOT / 'data' / 'final' / 'mediation_analyses'

# Sensitivity results land in the primary results tree (model names include cov_type
# so there is no collision with primary models).
LOCAL_RESULTS_DIR = _PROJECT_ROOT / 'results'

# Compiled CSV goes in a dedicated subdir; heatmaps go under supplement/
LOCAL_COMPILE_DIR = _PROJECT_ROOT / 'results' / 'sensitivity_analyses_mediation'
LOCAL_HEATMAP_DIR = _PROJECT_ROOT / 'results' / 'supplement' / 'sensitivity_analyses_mediation'

# Primary pipeline results root (canonical models; same tree as LOCAL_RESULTS_DIR)
PRIMARY_RESULTS_DIR = _PROJECT_ROOT / 'results'

# ── Shorthand lookup tables (subset of generate_hpc_jobs.py) ─────────────────
# SP predictor shorthand → (column_name, sample_default)
SP_PREDICTORS = {
    'spage':      ('psychedelic_age',          'spusers'),
    'lifenomic':  ('psycheduse_life_nomic',    'spusers'),
    'avgdose':    ('avg_life_dose',            'spusers'),
    'pc1':        ('psychedelic_use_PC1',      'full'),
    'pc1ranked':  ('psychedelic_rank_use_PC1', 'full'),
}

# Mediator shorthand → (response_col, predictor_in_dv_col)
# response_col   : raw column in mediator formula (DV of first path)
# predictor_in_dv: column used as predictor in DV formula
MEDIATORS = {
    'vchthreshold': ('vch_threshold', 'vch_threshold'),
    'vchrate':      ('vch_bl_yes_0',  'vch_bl_yes_0_normalized'),
    'vchrate75':    ('vch_bl_yes_75', 'vch_bl_yes_75_normalized'),
    'vchbeta':      ('vch_beta',      'vch_beta'),
    'vchnu':        ('vch_nu',        'vch_nu_normalized'),
    'vchomega':     ('vch_omega',     'vch_omega'),
}

# DV shorthand → column name
DVS = {
    'hppd_binary':  'hppd_binary',
    'caps_vision':  'caps_vision',
    'caps_total':   'caps_total',
    'lshs_total':   'lshs_total',
    'baggot_total': 'baggot_total',
}


def _model_name(analysis, cov_type):
    """Canonical model name: {dv_col}_{spvar}_{mediator}_{cov_type}."""
    dv_col = DVS[analysis['dv']]
    return f"{dv_col}_{analysis['spvar']}_{analysis['mediator']}_{cov_type}"


def _all_model_info():
    """
    Return list of (dv_col, spvar_short, med_short, cov_type, model_name) for
    all MED_ANALYSES × SENSITIVITY_MED_TYPES combinations.
    """
    rows = []
    for a in MED_ANALYSES:
        for ct in SENSITIVITY_MED_TYPES:
            rows.append((
                DVS[a['dv']], a['spvar'], a['mediator'], ct, _model_name(a, ct)
            ))
    return rows


# ── DHARMa columns in mediation summary CSVs ─────────────────────────────────
# Naming differs from nonsp diagnostics: 'dharma_ks_pval' not 'dharma_uniformity_pval'.
# dharma_heteroscedasticity_* tests quantiles vs *fitted values*.
# dharma_quantiles_vs_spvar_* tests quantiles vs the *predictor* (spvar).
# dharma_quantiles_vs_mediator_* tests quantiles vs the *mediator* (DV equation only).
# The vs-spvar and vs-mediator tests can flag issues that the vs-fitted test misses,
# because residual spread may be uniform across fitted values but vary with a
# specific predictor.  All three sets must be checked for complete coverage.
_MED_DHARMA_COLS = [
    'dharma_ks_pval',
    'dharma_outlier_pval',
    'dharma_dispersion_pval',
    'dharma_zeroinflation_pval',
    'dharma_heteroscedasticity_pval',
    'dharma_heteroscedasticity_q25_pval',
    'dharma_heteroscedasticity_q50_pval',
    'dharma_heteroscedasticity_q75_pval',
    'dharma_quantiles_vs_spvar_pval',
    'dharma_quantiles_vs_spvar_q25_pval',
    'dharma_quantiles_vs_spvar_q50_pval',
    'dharma_quantiles_vs_spvar_q75_pval',
    'dharma_quantiles_vs_mediator_pval',
    'dharma_quantiles_vs_mediator_q25_pval',
    'dharma_quantiles_vs_mediator_q50_pval',
    'dharma_quantiles_vs_mediator_q75_pval',
]

# ── Color map (matches sensitivity_analyses.py) ───────────────────────────────
_CMAP = LinearSegmentedColormap.from_list(
    'sens_prob_cmap',
    ['#f7f7f7', '#fff7bc', '#fee391', '#a1dab4', '#41b6c4', '#225ea8'],
    N=100,
)

# ── HPC module load prefix ────────────────────────────────────────────────────
_MED_MODULE_LOAD = (
    'module purge && module load foss/2022b && '
    'module load R/4.4.1-foss-2022b && export R_LIBS_USER=$HOME/R/4.4'
)

print(f'\n{"="*65}')
print('sensitivity_analyses_mediation.py')
print(f'  Stages: compile={RUN_COMPILE}  heatmap={RUN_HEATMAP}')
print(f'  Sensitivity types ({len(SENSITIVITY_MED_TYPES)}): {SENSITIVITY_MED_TYPES}')
print(f'  Analyses: {len(MED_ANALYSES)} (× {len(SENSITIVITY_MED_TYPES)} types = '
      f'{len(MED_ANALYSES) * len(SENSITIVITY_MED_TYPES)} models)')
print(f'{"="*65}')


# ==============================================================================
# DATA LOAD (lazy — only when generate stage is requested)
# ==============================================================================

_df_loaded = False




# ==============================================================================
# GENERATE FUNCTIONS
# ==============================================================================

def _strip_r_side(name):
    """Strip all R_SIDE_SUFFIXES from a model type name to get the base key."""
    try:
        from master_config import R_SIDE_SUFFIXES
    except ImportError:
        # Fallback if modules not yet loaded
        _KNOWN_SUFFIXES = [
            '_spusers', '_iqr', '_nocurrenthppd', '_nopsychosis', '_nonan_caps',
            '_beta',
        ]
        stem = name
        for sfx in sorted(_KNOWN_SUFFIXES, key=len, reverse=True):
            stem = stem.replace(sfx, '')
        return stem
    stem = name
    for sfx in R_SIDE_SUFFIXES:
        stem = stem.replace(sfx, '')
    return stem










# ==============================================================================
# COMPILE FUNCTIONS
# ==============================================================================







def _read_equation_diagnostics(csv_path):
    """
    Read DHARMa p-values and MCMC diagnostics from one equation's summary CSV.

    Returns (dharma_flagged, max_rhat, min_ess, num_divergents) or
    (False, nan, nan, 0) if the file is absent / unreadable.
    DHARMa flag: any p-value in _MED_DHARMA_COLS < DHARMA_ALPHA.
    MCMC flag:   max Rhat > RHAT_MAX OR min ESS < 1000 OR num_divergents > DIV_MAX.
    """
    if not csv_path.exists():
        return False, float('nan'), float('nan'), 0
    try:
        df = pd.read_csv(csv_path)
    except Exception:
        return False, float('nan'), float('nan'), 0
    if df.empty:
        return False, float('nan'), float('nan'), 0

    # DHARMa — p-values are broadcast to all rows; read from first row
    dharma_flagged = any(
        float(df[col].iloc[0]) < DHARMA_ALPHA
        for col in _MED_DHARMA_COLS
        if col in df.columns and not pd.isna(df[col].iloc[0])
    )

    max_rhat = float(pd.to_numeric(df['Rhat'], errors='coerce').max()) \
        if 'Rhat' in df.columns else float('nan')
    min_ess = float(
        df[['Bulk_ESS', 'Tail_ESS']].apply(pd.to_numeric, errors='coerce').min().min()
    ) if 'Bulk_ESS' in df.columns else float('nan')
    num_div = int(df['num_divergents'].iloc[0]) \
        if 'num_divergents' in df.columns else 0

    return dharma_flagged, max_rhat, min_ess, num_div


def _read_model_data(model_dir, dv_col, spvar_short, med_short, med_response, cov_type, model_nm):
    """
    Read path effects and diagnostics from one local mediation model directory.

    Reads:
      path_counterfactual_summary.csv → A, B, C' path rows
      mc_mediation_summary.csv        → NIE row (posterior MEDIAN — see below)
      summary_{dv_col}_{model_nm}.csv        → DV equation DHARMa + MCMC
      summary_{med_response}_{model_nm}.csv  → mediator equation DHARMa + MCMC

    Returns list of row dicts (one per path: A path, B path, C' path, NIE),
    or [] if the critical CSVs are absent.
    """
    pcf_path = model_dir / 'path_counterfactual_summary.csv'
    mcf_path = model_dir / 'mc_mediation_summary.csv'

    if not pcf_path.exists() or not mcf_path.exists():
        return []

    try:
        pcf = pd.read_csv(pcf_path)
        mcf = pd.read_csv(mcf_path)
    except Exception as e:
        print(f'  WARNING: could not read CSVs in {model_dir}: {e}')
        return []

    # Summary CSVs use original column names (underscores preserved, unlike brms internal names)
    dv_diag_path  = model_dir / f'summary_{dv_col}_{model_nm}.csv'
    med_diag_path = model_dir / f'summary_{med_response}_{model_nm}.csv'

    (dv_flag,  dv_rhat,  dv_ess,  dv_divs)  = _read_equation_diagnostics(dv_diag_path)
    (med_flag, med_rhat, med_ess, med_divs) = _read_equation_diagnostics(med_diag_path)

    dharma_flagged_any = dv_flag or med_flag
    max_rhat    = max(x for x in [dv_rhat, med_rhat] if not pd.isna(x)) \
                  if not (pd.isna(dv_rhat) and pd.isna(med_rhat)) else float('nan')
    min_ess     = min(x for x in [dv_ess, med_ess] if not pd.isna(x)) \
                  if not (pd.isna(dv_ess) and pd.isna(med_ess)) else float('nan')
    num_divergents = max(dv_divs, med_divs)

    # MCMC flagging (also contributes to red border in heatmap)
    mcmc_flagged = (
        (not pd.isna(max_rhat) and max_rhat > RHAT_MAX)
        or (not pd.isna(min_ess) and min_ess < 1000)
        or (num_divergents > DIV_MAX)
    )

    base_row = dict(
        dv=dv_col, spvar=spvar_short, mediator=med_short, cov_type=cov_type,
        model_name=model_nm,
        dharma_flagged_dv=dv_flag, dharma_flagged_med=med_flag,
        dharma_flagged_any=dharma_flagged_any,
        mcmc_flagged=mcmc_flagged,
        max_rhat=max_rhat, min_ess=min_ess, num_divergents=num_divergents,
    )

    rows = []
    # A, B, C' paths from path_counterfactual_summary.csv
    for path_label in ['A path', 'B path', "C' path"]:
        path_rows = pcf[pcf['effect'] == path_label]
        if not path_rows.empty:
            r = path_rows.iloc[0]
            rows.append({**base_row,
                         'path': path_label,
                         # Reported point estimate = posterior mean; the compiled
                         # column is named for whichever POINT_ESTIMATE_COL is.
                         # See modules/master_config.py.
                         POINT_ESTIMATE_COL: float(
                             point_estimate(r, source=str(pcf_path))),
                         'hdi_low':   float(r['hdi_low']),
                         'hdi_high':  float(r['hdi_high']),
                         'p_above_0': float(r['p_above_0']),
                         'p_below_0': float(r['p_below_0']),
                         'p_direction': float(r['p_direction'])})

    # NIE from mc_mediation_summary.csv
    nie_rows = mcf[mcf['effect'].str.startswith('NIE')]
    if not nie_rows.empty:
        r = nie_rows.iloc[0]
        rows.append({**base_row,
                     'path': 'NIE',
                     # NIE/NDE/TE/PMed come out of hpc_mediation.R's Monte-Carlo
                     # integration over the mediator's posterior predictive, whose
                     # heavy tails make the posterior mean unusable — so these
                     # report MC_EFFECT_POINT_ESTIMATE_COL ('median'), not
                     # POINT_ESTIMATE_COL ('mean') like the A/B/C' rows above.
                     # See modules/master_config.py. The column name follows the
                     # summary, so the compiled CSV says on its face which is which.
                     # The heatmap itself renders p_direction only and is unaffected;
                     # this keeps compiled_sensitivity_mediation.csv honest and
                     # matches mediation_results_table.py.
                     MC_EFFECT_POINT_ESTIMATE_COL: float(
                         point_estimate(r, source=str(mcf_path), mc_integrated=True)),
                     'hdi_low':   float(r['hdi_low']),
                     'hdi_high':  float(r['hdi_high']),
                     'p_above_0': float(r['p_above_0']),
                     'p_below_0': float(r['p_below_0']),
                     'p_direction': float(r['p_direction'])})

    return rows


def compile_sensitivity_mediation():
    """
    Compile sensitivity mediation results into a single long-format DataFrame.

    Sensitivity model results live in the same local tree as primary mediation:
      results/{dv}/mediation_models/{model_name}/
    Model names encode the covariate type, so there is no collision with primary
    models.  Canonical models (CANONICAL_MED_TYPE) are read from the same tree.

    Steps:
    1. Report which sensitivity models are present in results/{dv}/mediation_models/.
    2. Read all sensitivity model directories.
    4. Also read canonical (CANONICAL_MED_TYPE) models from the same results tree.
    5. Assemble and save compiled CSV to results/sensitivity_analyses_mediation/.

    Returns the compiled DataFrame.
    """
    model_info = _all_model_info()  # (dv_col, spvar_short, med_short, cov_type, model_name)

    # ── Report which models are present locally ───────────────────────────────
    # Nothing is fetched here.  Pulling fitted models back from the cluster is
    # 03_hpc/compile_mediation_results.py's job; this script only reads what is
    # already in results/.  A model listed as absent means that pull has not been
    # run for it, not that anything is wrong here.
    missing_pairs = [
        (dv_col, mn) for dv_col, _, _, _, mn in model_info
        if not (LOCAL_RESULTS_DIR / dv_col / 'mediation_models' / mn
                / 'path_counterfactual_summary.csv').exists()
    ]
    if missing_pairs:
        print(f'\n{len(missing_pairs)}/{len(model_info)} sensitivity mediation models '
              f'are not in results/ and will be absent from the heatmap:')
        for dv_col, mn in missing_pairs[:10]:
            print(f'    {dv_col}/mediation_models/{mn}')
        if len(missing_pairs) > 10:
            print(f'    ... and {len(missing_pairs) - 10} more')
        print('  Pull them with 03_hpc/compile_mediation_results.py subset.')
    else:
        print(f'\nAll {len(model_info)} sensitivity models found in results/.')

    # ── Read sensitivity models ───────────────────────────────────────────────
    all_rows = []
    n_found = 0
    for dv_col, spvar_short, med_short, cov_type, mn in model_info:
        med_response = MEDIATORS[med_short][0]
        model_dir = LOCAL_RESULTS_DIR / dv_col / 'mediation_models' / mn
        rows = _read_model_data(model_dir, dv_col, spvar_short, med_short,
                                med_response, cov_type, mn)
        if rows:
            n_found += 1
        all_rows.extend(rows)

    print(f'\n  Read {n_found}/{len(model_info)} sensitivity models '
          f'({len(all_rows)} path-level rows).')

    # ── Read canonical models from primary results dir ────────────────────────
    if CANONICAL_MED_TYPE:
        n_canon = 0
        for a in MED_ANALYSES:
            dv_col    = DVS[a['dv']]
            med_short = a['mediator']
            med_response = MEDIATORS[med_short][0]
            mn = _model_name(a, CANONICAL_MED_TYPE)
            model_dir = PRIMARY_RESULTS_DIR / dv_col / 'mediation_models' / mn
            rows = _read_model_data(model_dir, dv_col, a['spvar'], med_short,
                                    med_response, CANONICAL_MED_TYPE, mn)
            if rows:
                n_canon += 1
            all_rows.extend(rows)
        print(f'  Read {n_canon}/{len(MED_ANALYSES)} canonical ({CANONICAL_MED_TYPE}) models.')

    if not all_rows:
        print('\nWARNING: No data compiled — check that models have finished on HPC.')
        return pd.DataFrame()

    compiled = pd.DataFrame(all_rows)

    LOCAL_COMPILE_DIR.mkdir(parents=True, exist_ok=True)
    out_csv = LOCAL_COMPILE_DIR / COMPILED_CSV_STEM
    compiled.to_csv(out_csv, index=False)
    print(f'\nSaved {len(compiled):,} rows → {out_csv}')
    print(f'  DVs:          {sorted(compiled["dv"].unique())}')
    print(f'  Covariate types: {sorted(compiled["cov_type"].unique())}')
    print(f'  Paths:        {sorted(compiled["path"].unique())}')
    return compiled


# ==============================================================================
# HEATMAP FUNCTIONS
# ==============================================================================

# Modifier display labels (matches sensitivity_analyses.py)
_MODIFIER_LABELS = {
    'spusers':              'SP users only',
    'nooutlier':            'no outliers',
    'nooutliers':           'no outliers',
    'nopsychosis':          'no psychosis',
    'nocurrenthppd':        'excl. current HPPD',
    'iqr':                  'IQR outlier filter',
    'beta':                 '+ vch_beta covariate',
}


def _format_cov_type_label(s):
    """
    Format a covariate type name for x-axis display.
    'empirical_covariates_spusers'      → 'empirical_covariates\n(SP users only)'
    Falls back to s unchanged if no modifiers recognised.
    """
    remaining = s
    found = []
    for mod in sorted(_MODIFIER_LABELS.keys(), key=len, reverse=True):
        suffix = f'_{mod}'
        if remaining.endswith(suffix):
            remaining = remaining[:-len(suffix)]
            found.insert(0, mod)
    if not found:
        return s.replace('_', ' ')
    descriptions = [_MODIFIER_LABELS[m] for m in found]
    return f'{remaining}\n({", ".join(descriptions)})'


def _row_label(spvar_short, med_short):
    """Y-axis label: 'Psych. age → VCH β'"""
    spvar_lbl = _SPVAR_DISPLAY.get(spvar_short, spvar_short)
    med_lbl   = _MEDIATOR_DISPLAY.get(med_short, med_short)
    return f'{spvar_lbl} \u2192 {med_lbl}'


def _is_flagged(row):
    """Return True if the row should get a red border (DHARMa or MCMC issues)."""
    if row.get('dharma_flagged_any', False):
        return True
    if row.get('mcmc_flagged', False):
        return True
    return False


def _get_cell(compiled_df, dv_col, spvar_short, med_short, cov_type, path_type):
    """
    Return (p_direction, sign, status) for one heatmap cell.
    status ∈ {'ok', 'flagged', 'missing'}
    """
    # Gray out logically inapplicable (cov_type, dv) combos
    if (cov_type, dv_col) in INAPPLICABLE_CELLS:
        return float('nan'), 'pos', 'missing'

    mask = (
        (compiled_df['dv']       == dv_col)
        & (compiled_df['spvar']    == spvar_short)
        & (compiled_df['mediator'] == med_short)
        & (compiled_df['cov_type'] == cov_type)
        & (compiled_df['path']     == path_type)
    )
    sub = compiled_df[mask]
    if sub.empty:
        return float('nan'), 'pos', 'missing'
    r = sub.iloc[0]
    p_above = float(r.get('p_above_0', 0) or 0)
    p_below = float(r.get('p_below_0', 0) or 0)
    p_dir   = float(r.get('p_direction', max(p_above, p_below)) or 0)
    sign    = 'pos' if p_above >= p_below else 'neg'
    status  = 'flagged' if _is_flagged(r) else 'ok'
    return p_dir, sign, status


def _plot_med_heatmap_panel(ax, compiled_df, analyses_for_dv, cov_types, dv_col,
                             path_type, show_xticklabels=True,
                             x_tick_fs=14, y_tick_fs=14, cell_text_fs=11):
    """
    Render one mediation heatmap panel (one DV) onto ax.

    Rows = (spvar × mediator) combos in analyses_for_dv order.
    Columns = cov_types.
    Returns the AxesImage (for shared colorbar) or None if no data.
    """
    try:
        from master_config import COVARIATE_SET_LABELS
    except ImportError:
        COVARIATE_SET_LABELS = {}

    # Build cell arrays
    row_labels = [_row_label(a['spvar'], a['mediator']) for a in analyses_for_dv]
    n_rows = len(row_labels)
    n_cols = len(cov_types)

    prob_arr = np.full((n_rows, n_cols), np.nan)
    sign_arr = [['pos'] * n_cols for _ in range(n_rows)]
    stat_arr = [['missing'] * n_cols for _ in range(n_rows)]

    any_data = False
    for ri, a in enumerate(analyses_for_dv):
        for ci, ct in enumerate(cov_types):
            p, s, st = _get_cell(compiled_df, dv_col, a['spvar'], a['mediator'], ct, path_type)
            prob_arr[ri, ci] = p
            sign_arr[ri][ci] = s
            stat_arr[ri][ci] = st
            if not np.isnan(p):
                any_data = True

    if not any_data:
        return None

    im = ax.imshow(prob_arr, cmap=_CMAP, aspect='auto', vmin=0.5, vmax=1.0)

    # Grey wash for NaN cells
    nan_mask = np.ma.masked_where(~np.isnan(prob_arr), np.ones_like(prob_arr))
    ax.imshow(nan_mask, cmap=plt.get_cmap('Greys'), aspect='auto',
              vmin=0, vmax=1, alpha=0.2)

    ax.set_xticks(range(n_cols))
    ax.set_yticks(range(n_rows))

    if show_xticklabels:
        col_labels = [
            COVARIATE_SET_LABELS.get(ct, _format_cov_type_label(ct))
            for ct in cov_types
        ]
        ax.set_xticklabels(col_labels, rotation=35, ha='right',
                           fontsize=x_tick_fs, multialignment='center')
    else:
        ax.set_xticklabels([])
        ax.tick_params(axis='x', which='both', bottom=False, top=False)

    ax.set_yticklabels(row_labels, fontsize=y_tick_fs)

    # Cell text + borders
    for ri in range(n_rows):
        for ci in range(n_cols):
            p  = prob_arr[ri, ci]
            sg = sign_arr[ri][ci]
            st = stat_arr[ri][ci]

            if not np.isnan(p):
                sign_str = '+' if sg == 'pos' else '\u2212'
                txt_col  = 'white' if p > 0.92 else 'black'
                ax.text(ci, ri, f'{sign_str}{p:.2f}',
                        ha='center', va='center',
                        fontsize=cell_text_fs, color=txt_col, fontweight='bold')

            if st == 'flagged':
                ax.add_patch(mpatches.Rectangle(
                    (ci - 0.48, ri - 0.48), 0.96, 0.96,
                    fill=False, edgecolor='red', linewidth=2.0,
                ))
            elif st == 'missing':
                ax.add_patch(mpatches.Rectangle(
                    (ci - 0.48, ri - 0.48), 0.96, 0.96,
                    fill=False, edgecolor='#e991c8', linewidth=2.0,
                ))

    # Minor grid
    ax.set_xticks(np.arange(n_cols + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(n_rows + 1) - 0.5, minor=True)
    ax.grid(which='minor', color='#cccccc', linestyle='-', linewidth=0.35)
    ax.tick_params(which='minor', length=0)

    return im


def plot_compound_path_heatmap(compiled_df, path_type, all_cov_types, out_dir):
    """
    Compound two-panel heatmap for one path type.

    Top panel: hppd_binary analyses.
    Bottom panel: caps_vision (or whichever other DV is in MED_ANALYSES).

    X-axis column order: canonical type first, then sensitivity types.
    X-tick labels shown on bottom panel only.
    Single shared y-axis label; single colorbar on the right.
    Layout mirrors plot_compound_sensitivity_heatmap() in sensitivity_analyses.py.

    Output: {out_dir}/{safe_path_type}_compound_heatmap.png
    """
    import matplotlib as mpl

    try:
        from master_config import dv_to_lab_short
    except ImportError:
        dv_to_lab_short = {}

    # Ordered list of DVs appearing in MED_ANALYSES
    dv_order = list(dict.fromkeys(DVS[a['dv']] for a in MED_ANALYSES))

    # Analyses per DV (preserving MED_ANALYSES order within each DV)
    analyses_per_dv = {dv: [] for dv in dv_order}
    for a in MED_ANALYSES:
        analyses_per_dv[DVS[a['dv']]].append(a)

    # Keep only types with any data for this path
    def _has_data(ct):
        return not compiled_df[
            (compiled_df['cov_type'] == ct)
            & (compiled_df['path'] == path_type)
        ].empty

    present_types = [ct for ct in all_cov_types if _has_data(ct)]
    if not present_types:
        print(f'  [{path_type}] No data in compiled CSV — skipping heatmap')
        return

    n_cols   = len(present_types)
    row_counts = {dv: len(analyses_per_dv[dv]) for dv in dv_order}
    n_panels = len([dv for dv in dv_order if row_counts[dv] > 0])
    if n_panels == 0:
        return

    # Figure dimensions
    cell_w  = max(2.0, min(3.2, 36 / max(n_cols, 1)))
    cell_h  = max(0.55, min(1.0, 14 / max(sum(row_counts.values()), 1)))
    fig_w   = n_cols * cell_w + 8.5
    fig_h   = sum(row_counts.values()) * cell_h + 5.0

    X_TICK_FS  = 14
    Y_TICK_FS  = 14
    TITLE_FS   = 17
    XLABEL_FS  = 15
    YLABEL_FS  = 18
    CELL_FS    = 11

    # Path type → human-readable title fragment
    _PATH_TITLES = {
        'A path':   'A path  (SP predictor \u2192 mediator)',
        'B path':   'B path  (mediator \u2192 DV)',
        "C' path":  "C\u2019 path  (direct: SP \u2192 DV)",
        'NIE':      'NIE  (natural indirect effect)',
    }
    path_title = _PATH_TITLES.get(path_type, path_type)

    height_ratios = [max(row_counts[dv], 1) for dv in dv_order]

    with mpl.rc_context({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    }):
        fig = plt.figure(figsize=(fig_w, fig_h))

        gs = GridSpec(
            n_panels, 2,
            figure=fig,
            width_ratios=[1, 0.03],
            height_ratios=height_ratios[:n_panels],
            hspace=0.10,
            wspace=0.03,
            left=0.22, right=0.93,
            bottom=0.22, top=0.92,
        )

        axes = [fig.add_subplot(gs[i, 0]) for i in range(n_panels)]
        cax  = fig.add_subplot(gs[:, 1])

        ims = []
        for i, dv_col in enumerate(dv_order[:n_panels]):
            ax        = axes[i]
            is_bottom = (i == n_panels - 1)
            dv_label  = dv_to_lab_short.get(dv_col, dv_col)

            im = _plot_med_heatmap_panel(
                ax, compiled_df, analyses_per_dv[dv_col],
                present_types, dv_col, path_type,
                show_xticklabels=is_bottom,
                x_tick_fs=X_TICK_FS,
                y_tick_fs=Y_TICK_FS,
                cell_text_fs=CELL_FS,
            )
            if im is not None:
                ims.append(im)

            ax.set_title(dv_label, fontsize=TITLE_FS, fontweight='bold', pad=6)

        # Supertitle for the path type
        fig.suptitle(path_title, fontsize=TITLE_FS + 2, fontweight='bold', y=0.975)

        # X-axis label on bottom panel
        axes[-1].set_xlabel(
            'Model type (sensitivity)', fontsize=XLABEL_FS,
            fontweight='bold', labelpad=8,
        )

        # Shared y-axis label (left margin).  Moved right from 0.03 to 0.157 on
        # 2026-09-02 to sit just outside the y-tick labels rather than at the
        # far edge of the margin.  Axes start at left=0.22 and the longest
        # y-tick label ("Avg. dose -> VCH threshold") begins at ~0.165 of
        # fig_w, so 0.157 leaves ~0.3 in of clearance.  Do not raise this past
        # ~0.16: at 0.18 the label renders on top of the tick labels.
        fig.text(
            0.157, 0.5, 'Mediation model',
            va='center', ha='center', rotation='vertical',
            fontsize=YLABEL_FS, fontweight='bold',
        )

        # Colorbar in its own GridSpec column
        if ims:
            cbar = fig.colorbar(ims[0], cax=cax)
            # No colorbar title (kept in sync with sensitivity_analyses.py):
            # the figure caption carries the meaning of the scale, and the
            # axis label crowded the tick labels once they were enlarged.
            # Ticks at 24 pt (3x the previous 8 pt).
            cbar.ax.tick_params(labelsize=24)

        out_dir.mkdir(parents=True, exist_ok=True)
        # Only the NIE heatmap is a numbered supplementary figure. The other
        # three paths keep a descriptive stem, so adding one back to
        # PATHS_TO_PLOT cannot silently claim a figure number.
        stem = PATH_FILE_STEMS.get(
            path_type,
            (path_type.replace("'", 'prime').replace(' ', '_')
                      .replace('/', '').lower() + '_compound_heatmap'),
        )
        out_path = out_dir / f'{stem}.png'
        plt.savefig(out_path, dpi=200, bbox_inches='tight')
        tiff_path = out_dir / f'{stem}.tiff'
        plt.savefig(tiff_path, dpi=200, bbox_inches='tight')
        plt.close()
        print(f'  Saved: {out_path}')
        print(f'  Saved: {tiff_path}')


# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

if __name__ == '__main__':
    # ── Stage 1: Compile ──────────────────────────────────────────────────────
    if RUN_COMPILE:
        print(f'\n{"="*65}')
        print('STAGE 1: COMPILE')
        print(f'{"="*65}')
        compiled_df = compile_sensitivity_mediation()

    # ── Stage 2: Heatmap ──────────────────────────────────────────────────────
    if RUN_HEATMAP:
        print(f'\n{"="*65}')
        print('STAGE 2: HEATMAP')
        print(f'{"="*65}')

        compiled_csv = LOCAL_COMPILE_DIR / COMPILED_CSV_STEM
        if not compiled_csv.exists():
            print(f'  ERROR: compiled CSV not found at {compiled_csv}')
            print(f'  Run the compile stage first:  '
                  f'python 04_visualizations/supplement/sensitivity_analyses_mediation.py compile')
        else:
            compiled_df = pd.read_csv(compiled_csv, low_memory=False)
            print(f'  Loaded {len(compiled_df):,} rows from {compiled_csv.name}')
            print(f'  DVs:          {sorted(compiled_df["dv"].unique())}')
            print(f'  Covariate types: {sorted(compiled_df["cov_type"].unique())}')
            print(f'  Paths:        {sorted(compiled_df["path"].unique())}')

            # Column order: three canonical slots first, then sensitivity variants.
            # Matches sensitivity_analyses.py ordering:
            #   CANONICAL | SECOND_CANONICAL | THIRD_CANONICAL | SENSITIVITY_MED_TYPES...
            present = set(compiled_df['cov_type'].values)
            canon_present  = CANONICAL_MED_TYPE is not None        and CANONICAL_MED_TYPE        in present
            second_present = SECOND_CANONICAL_MED_TYPE is not None and SECOND_CANONICAL_MED_TYPE in present
            third_present  = THIRD_CANONICAL_MED_TYPE is not None  and THIRD_CANONICAL_MED_TYPE  in present
            _all_canonical = {ct for ct in [CANONICAL_MED_TYPE,
                                            SECOND_CANONICAL_MED_TYPE,
                                            THIRD_CANONICAL_MED_TYPE] if ct is not None}
            ordered_cols = (
                ([CANONICAL_MED_TYPE]        if canon_present  else [])
                + ([SECOND_CANONICAL_MED_TYPE] if second_present else [])
                + ([THIRD_CANONICAL_MED_TYPE]  if third_present  else [])
                + [ct for ct in SENSITIVITY_MED_TYPES
                   if ct in present and ct not in _all_canonical]
            )
            all_expected = list(_all_canonical) + SENSITIVITY_MED_TYPES
            missing_cols = [ct for ct in all_expected if ct not in present]
            if missing_cols:
                print(f'  NOTE: types not yet compiled: {missing_cols}')

            # One compound heatmap per path type
            for path_type in PATHS_TO_PLOT:
                if path_type not in compiled_df['path'].values:
                    print(f'\n  [{path_type}] Not in compiled data — skipping')
                    continue
                print(f'\n  Path: {path_type}')
                plot_compound_path_heatmap(
                    compiled_df, path_type, ordered_cols, LOCAL_HEATMAP_DIR
                )

    print(f'\n{"="*65}')
    print('Done.')
    print(f'{"="*65}')
