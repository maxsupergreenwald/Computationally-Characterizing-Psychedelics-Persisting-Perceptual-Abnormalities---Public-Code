#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Results narrative pipeline.

Run:
    python results_narrative.py

Output:
    results_narrative_output_editted.txt  (OUTPUT_TXT below)

The output file is the Results section of the manuscript, byte for byte: no
banner, no separators, no markup.  Section headings are plain lines; paragraphs
are separated by a single blank line.

─────────────────────────────────────────────────────────────────────────────
Figure & Table Numbering Logic
─────────────────────────────────────────────────────────────────────────────
Main figure numbers derive from FIRST_FIGURE_NUM in Section 2; change that one
value to renumber every main figure. Supplementary references are literals in
the same block — the supplement is not ordered by the same sequence.

Main Tables
  TABLE_1 = Demographics (psychiatric & substance history merged in)

Figures 2-7  (FIRST_FIGURE_NUM = 1 by default)
    Figure 2  (+1): PPA history descriptive
    Figure 3  (+2): CAPS item descriptive
    Figure 4  (+3): SP predictors of PPA history & CAPS vision
    Figure 5  (+4): VCH behaviour, incl. its mediation diagrams
    Figure 6  (+5): HGF parameters, incl. their mediation diagrams
    Figure 7  (+6): beta / signal-detection exploratory panel

    Figures 4 and 5 share a panel scheme:
      a = PPA history, nonparametric      c = PPA history, regression
      b = CAPS vision, nonparametric      d = CAPS vision, regression

    Figure 5 mediation panels (row-wise; mirrors MED_PANEL_MODELS in
    04_visualizations/0X_all_figures.py):
      e = psychedelic_age -> vch_threshold -> PPA history
      f = avg_life_dose   -> vch_threshold -> CAPS vision
      g = psychedelic_age -> vch_bl_yes_0  -> PPA history
      h = avg_life_dose   -> vch_bl_yes_0  -> CAPS vision

    Figure 6 panels: b/h/j = omega/nu/beta group comparisons; c/i/k = the CAPS
    correlations; d & f, e & g = belief trajectories; l/m = regressions;
    n = psychedelic_age -> vch_beta -> PPA history;
    o = avg_life_dose   -> vch_beta -> CAPS vision.

    Every panel reference in the prose comes from a named constant defined at
    the bottom of the matching block in Section 4 -- never written inline.

Supplementary figures and tables
    Every supplementary reference that reaches the narrative is a named constant
    in the "FIGURES AND TABLES LISTED HERE" block below.  Renumber there.
─────────────────────────────────────────────────────────────────────────────
"""

# ================================================================
# SECTION 1: IMPORTS & SETUP
# ================================================================

import sys
import math
import warnings
from pathlib import Path
from datetime import datetime

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

# Every in-repo path below is anchored to this script's location, not to the
# process CWD, so the script runs identically no matter where it is launched
# from (repo root, an IDE "Run File" button, a cron job). Launching from
# anywhere other than 05_results_narrative/ previously failed at the import
# below with "ModuleNotFoundError: No module named 'visualization_helpers'".
_REPO_ROOT = Path(__file__).resolve().parent.parent

sys.path.insert(0, str(_REPO_ROOT / "modules"))

from visualization_helpers import *
from master_config import (
    dv_to_lab_short,
    dv_to_lab,
    iv_type_dict,
    caps_types,
    hppd_variables,
    severity_vars,
    SP_USER_COL,
    SP_USER_VALUE,
    RECRUIT_CSV,
)
# point_estimate() reads POINT_ESTIMATE_COL — the reported posterior summary —
# from a brms summary row, and raises PointEstimateColumnMissing (naming the
# file and the required refit) if that column is absent. Every reported Δ below
# goes through it; none reads a column name directly.
from master_config import point_estimate


# ================================================================
# SECTION 2: DATA LOADING & CONFIG OPTIONS
# ================================================================

# The shipped analysis dataframe is already fully prepared — every derived column
# is present in the CSV, so there is no preparation step to run.
# most_recent_public_df() picks the newest data/final/df_public_*.csv and prints
# which one it chose.
from data_prep import most_recent_public_df
_DATA_DIR    = _REPO_ROOT / "data" / "final"
_RESULTS_DIR = _REPO_ROOT / "results"
# String form, for the f-strings that build result paths by interpolation.
_RESULTS = str(_RESULTS_DIR)
df = pd.read_csv(most_recent_public_df(_DATA_DIR), low_memory=False)
df_sp = df[df[SP_USER_COL] == SP_USER_VALUE].copy()   # canonical SP-user filter
print(f"Loaded: df={df.shape}, df_sp={df_sp.shape}")

dv_to_lab["hppd_sx_count"] = "# of HPPD Experiences"
df_hppd = df[df["persist_vis_yn"] > 0].copy()

df_recruit_raw = pd.read_csv(RECRUIT_CSV).copy()

# Audit table (pre-filtering)
recruit_audit = df_recruit_raw.copy()
recruit_audit["excluded_testing_filter"] = (
    (recruit_audit["record_id"] < 203) & (recruit_audit["student_yn"] < 1)
)
recruit_audit["excluded_cutoff_filter"] = (
    recruit_audit["record_id"] >= df["record_id"].max()
)
recruit_audit["kept_initial_recruit_filters"] = ~(
    recruit_audit["excluded_testing_filter"] | recruit_audit["excluded_cutoff_filter"]
)

# Working df_recruit with filters applied
df_recruit = df_recruit_raw.copy()
df_recruit = df_recruit[~((df_recruit['record_id'] < 203) & (df_recruit['student_yn'] < 1))].copy()
df_recruit = df_recruit[df_recruit["record_id"] <= df["record_id"].max()].copy()
df_recruit.loc[df_recruit['record_id'] == 1858, ['salvage_yn', 'qc_passed']] = 0
df_recruit['raven_total'] = (
    df_recruit[[f"correct_answer{x}_v2" for x in [2, 3]]].sum(axis=1, min_count=1)
    + df_recruit["correct_answer_v2"]
    + df_recruit['raven_total_score_v2']
)
df_recruit.loc[(df_recruit["raven_total"] < 1) & (df_recruit['student_yn'] > 0), 'raven_total'] = np.nan

# Counterfactual CSV: response-scale marginal contrasts (Δ per +1 Gelman SD).
# Columns: dv, spvar, cov, median, hdi_lower_94, hdi_upper_94,
#          prob_above_0, prob_below_0, N_obs, scale, ...
# Note: hppd_binary is stored as 'hppd_binary' in this CSV (updated Jun 2026).
# mu (count/continuous) effects only — hu (hurdle) component is in
#   hu_paths_summary.csv inside each mediation model directory
#   (columns: path, var_brms_col, median, hdi_lower_94, hdi_upper_94, ...).
COUNTERFACTUAL_CSV = (
    f"{_RESULTS}/"
    "sensitivity_analyses_single_paths/existingresults_manuscript_counterfactual.csv"
)


# ── Figure & Table Number Constants ─────────────────────────────────────────
# Adjust FIRST_FIGURE_NUM to renumber all main figures in one place
# (e.g. set to 2 if a new figure is inserted before these results).
# Supplementary numbers are literals below — the supplement is not ordered by
# the same sequence, so deriving them from an offset only obscured them.

FIRST_FIGURE_NUM     = 1   # First main results figure number

# File the Results narrative is written to. Anchored to this script's directory
# like every other path here, so launching from the repo root does not scatter a
# second copy of the output there.
OUTPUT_TXT = Path(__file__).resolve().parent / "results_narrative_output_editted.txt"

# ── Effect scale ─────────────────────────────────────────────────────────────
# Every reported effect is a response-scale marginal contrast (Δ per +1 Gelman SD
# in the predictor), read from COUNTERFACTUAL_CSV for single paths and from
# path_counterfactual_summary.csv (A/B/C' paths) + mc_mediation_summary.csv (NIE)
# for mediation.  CAPS analyses use the SP-user sample (nice_covariates_spusers).

_cf_df = pd.read_csv(COUNTERFACTUAL_CSV)

# ── Mann-Whitney verbosity toggle ────────────────────────────────────────────
# True  → full string: U, p, rrb, Med(+) [IQR], Med(-) [IQR]
# False → compact string: U, p only (default)
MANN_WHITNEY_VERBOSE = False

##############################################################################
### FIGURES AND TABLES LISTED HERE!
################################################################################
# Main tables — change numbers here to renumber throughout the narrative
TABLE_1 = "Table 1"   # Demographics (psychiatric & substance history merged in)
# SP_FIG is a full label ("Supplementary Figure S6"), used bare in the text.
# PPA_FIG / CAPS_FIG and everything below are integers — the surrounding prose
# already supplies "Fig.", and panel refs are built as f"Fig. {FIG_X}a".
# PPA and CAPS descriptive figures sit between SP use and SP predictors, so
# all later figures are offset by +3 relative to FIRST_FIGURE_NUM.
SP_FIG               = "Supplementary Figure S6"   # SP use patterns
SP_FIG_SHORT         = SP_FIG.replace("Figure", "Fig.")  # second, parenthetical mention
PPA_FIG              = FIRST_FIGURE_NUM + 1             # integer, e.g. 3
CAPS_FIG             = FIRST_FIGURE_NUM + 2             # integer, e.g. 4
FIG_SP_PREDICTORS    = FIRST_FIGURE_NUM + 3             # integer, e.g. 5
FIG_VCH_BEHAVIOR     = FIRST_FIGURE_NUM + 4             # integer, e.g. 6
FIG_VCH_COMPUTATIONS = FIRST_FIGURE_NUM + 5             # integer, e.g. 7
FIG_MEDIATION        = FIRST_FIGURE_NUM + 5            # COMBINED INTO VCH_COMPUTATIONS NOW!
FIG_BETA = FIRST_FIGURE_NUM + 6                            # integer, e.g. 7


# Supplementary references, one per sentence that cites one. Renumber here.
sensitivity_analysis_heatmap_single_path = "Supplementary Figure S4"   # single-path outlier sensitivity heatmap
task_engagement_fig_threshold_error      = "Supplementary Figure S7"   # VCH task-engagement measures
supp_fig_consort                         = "Supplementary Figure S9"   # CONSORT participant flow diagram
mann_whitney_u_table                     = "Supplementary Table S2"    # all binary Mann-Whitney U tests
supp_table_mediation_sensitivity         = "Supplementary Table S3"    # full mediation sensitivity suite

# ── Reference-list citation numerals ────────────────────────────────────────
# Superscript numbers exactly as they appear against the manuscript reference
# list. Grouped citations are assembled in the narrative below ("6,7,10", or
# "6–9" with an en dash for a contiguous run), so renumbering the reference list
# only ever requires editing this dict.
CITE = {
    "izmi_2024":           "3",
    "carhart_nutt_2010":   "6",
    "baggott_2011":        "7",
    "zhou_2025":           "8",
    "muller_2022":         "9",
    "kvam_2023":          "10",
    "kessler_2005":       "56",
    "hirschfeld_2023":    "57",
    # Cited for the classical reading of a flat psychometric curve as low
    # sensitivity, in the beta/SDT exploratory section.
    "sdt_low_sensitivity": "58",
}



# ================================================================
# SECTION 3: UTILITY FUNCTIONS
# ================================================================

# ── Formatting helpers ───────────────────────────────────────────

def _fmt_num(x, digits=3):
    if pd.isna(x):
        return "nan"
    s = f"{float(x):.{digits}f}"
    rounded = float(s)
    if rounded.is_integer():
        return f"{int(rounded)}"
    return s


def _fmt_trim_num(x, digits=2):
    if pd.isna(x):
        return "nan"
    s = f"{float(x):.{digits}f}"
    return s.rstrip("0").rstrip(".") if "." in s else s


def _fmt_p(p):
    if pd.isna(p):
        return "nan"
    return "< 0.001" if float(p) < 0.001 else f"= {float(p):.3f}"


def _fmt_prob(prob):
    if pd.isna(prob):
        return "nan"
    prob = float(prob)
    if prob > 1:
        prob = prob / 100
    return f"{_fmt_trim_num(100 * prob, 1)}%"


def _num_pct_from_mask(mask, denom,semicolon=None):
    n = int(mask.sum())
    if denom == 0:
        return "n/a"
    pct = (n / denom) * 100
    if semicolon:
        return f"{n}; {_fmt_trim_num(pct, 1)}%"
    else:
        return f"{n} ({_fmt_trim_num(pct, 1)}%)"


def _pct_str(df, mask, denom_mask=None):
    if denom_mask is None:
        denom_mask = pd.Series(True, index=df.index)
    d = df.loc[denom_mask]
    return str(int(round(mask.loc[d.index].mean() * 100)))


# ── Row-selection helper ─────────────────────────────────────────

def _select_single_row(df, mask, context):
    sub = df.loc[mask].copy()
    if len(sub) == 0:
        raise ValueError(f"No rows found for {context}")
    if len(sub) > 1:
        if "N" in sub.columns:
            sub = sub.sort_values("N", ascending=False)
        print(f"Note: Multiple rows found for {context}; using first row after sorting.")
    return sub.iloc[0]


# ── Result-string builders ───────────────────────────────────────

def _mw_result_string(nonparam_df, var_name, verbose=None):
    """Mann-Whitney U string from a boxplot_grid CSV.

    When MANN_WHITNEY_VERBOSE=False (default): returns 'U = X, p Y' only.
    When MANN_WHITNEY_VERBOSE=True: also includes rrb and group medians with IQRs.
    Pass verbose=True/False to override the global toggle for a specific call.
    """
    row = _select_single_row(
        nonparam_df,
        (nonparam_df["dv"] == var_name) & (nonparam_df["test"] == "mann_whitney"),
        context=f"Mann-Whitney row for var='{var_name}'",
    )
    base = f"U = {_fmt_num(row['statistic'], 1)}, p {_fmt_p(row['p_value'])}"
    use_verbose = MANN_WHITNEY_VERBOSE if verbose is None else verbose
    if not use_verbose:
        return base
    g1_label = str(row.get("group_1_label", ""))
    g2_label = str(row.get("group_2_label", ""))
    if "(+)" in g1_label:
        plus_med, plus_q1, plus_q3 = row["group_1_median"], row["group_1_iqr_q1"], row["group_1_iqr_q3"]
        minus_med, minus_q1, minus_q3 = row["group_2_median"], row["group_2_iqr_q1"], row["group_2_iqr_q3"]
    elif "(+)" in g2_label:
        plus_med, plus_q1, plus_q3 = row["group_2_median"], row["group_2_iqr_q1"], row["group_2_iqr_q3"]
        minus_med, minus_q1, minus_q3 = row["group_1_median"], row["group_1_iqr_q1"], row["group_1_iqr_q3"]
    else:
        print(f"Warning: Could not identify (+) group for var='{var_name}'. Labels: '{g1_label}', '{g2_label}'")
        return np.nan
    return (
        f"{base}, "
        f"rrb = {_fmt_num(row['rank_biserial_correlation'], 3)}, "
        f"Med(+) = {_fmt_num(plus_med)} [{_fmt_num(plus_q1)}, {_fmt_num(plus_q3)}], "
        f"Med(-) = {_fmt_num(minus_med)} [{_fmt_num(minus_q1)}, {_fmt_num(minus_q3)}]"
    )


def _perm_mean_diff(plus_vals, minus_vals, n_boot=10000, seed=42):
    """Bootstrap CI on the difference in group means.

    Reports the observed Δ = mean(+) − mean(−) and a 94% bootstrap CI on Δ
    (each group resampled separately with replacement, n_boot iterations).

    Use this when the figure shows bootstrapped group-mean trajectories.
    """
    plus_arr  = np.asarray(plus_vals,  dtype=float)
    minus_arr = np.asarray(minus_vals, dtype=float)
    plus_arr  = plus_arr[~np.isnan(plus_arr)]
    minus_arr = minus_arr[~np.isnan(minus_arr)]
    n_plus  = len(plus_arr)
    n_minus = len(minus_arr)
    obs_diff = float(np.mean(plus_arr) - np.mean(minus_arr))

    # Bootstrap CI on Δ — resample each group separately with replacement
    rng = np.random.default_rng(seed)
    boot_plus  = plus_arr[rng.integers(0, n_plus,  size=(n_boot, n_plus))]
    boot_minus = minus_arr[rng.integers(0, n_minus, size=(n_boot, n_minus))]
    boot_diffs = boot_plus.mean(axis=1) - boot_minus.mean(axis=1)
    ci_lo = float(np.percentile(boot_diffs, 3))   # 94% CI: 3rd–97th percentile
    ci_hi = float(np.percentile(boot_diffs, 97))

    return f"Δ = {_fmt_num(obs_diff)}, 94% CI [{_fmt_num(ci_lo)}, {_fmt_num(ci_hi)}]"


def _spearman_string(nonparam_df, row_var, column_var):
    """Spearman ρ string from a correlation_grid CSV."""
    row = _select_single_row(
        nonparam_df,
        (nonparam_df["row_var"] == row_var) & (nonparam_df["column_var"] == column_var),
        context=f"Spearman row for var='{row_var}'",
    )
    return f"ρ = {_fmt_num(row['rho'])}, p {_fmt_p(row['p_value'])}"


def _partial_spearman_str(df_sub, x_col, y_col, covariate='age_v2'):
    """
    Partial Spearman correlation between x_col and y_col, controlling for
    `covariate`. Matches the method used by correlation_matrix_plot in
    visualization_helpers_parts/correlation_matrix_plot.py: all variables are
    rank-transformed first, then partial Pearson correlation is computed on
    the ranks via pingouin.partial_corr (partial Pearson on ranks = partial
    Spearman). Returns a formatted string matching _spearman_string output.
    """
    import pingouin as _pg
    from scipy.stats import rankdata as _rankdata
    all_vars = [x_col, y_col, covariate]
    d = df_sub[all_vars].dropna()
    ranked = pd.DataFrame(
        {v: _rankdata(d[v]) for v in all_vars},
        index=d.index,
    )
    pc = _pg.partial_corr(data=ranked, x=x_col, y=y_col, covar=covariate, method='pearson')
    rho  = float(pc['r'].iloc[0])
    pval = float(pc['p-val'].iloc[0])
    return f"ρ = {_fmt_num(rho)}, p {_fmt_p(pval)}"


# ── Counterfactual result helpers ────────────────────────────────

def _cf_result_string(cf_df, spvar_name, cov_name, dv_name):
    """
    Counterfactual result string: Δ = X, P(Δ?0) = Y%, 94% HDI [A, B].

    Effect = E[Y|X=mean+1SD] − E[Y|X=mean] (response-scale marginal contrast).
    dv_name: use 'hppd_binary' for the HPPD binary outcome (updated Jun 2026).
    Returns (result_str, n_obs, estimate).
    """
    mask = (
        (cf_df['dv'] == dv_name)
        & (cf_df['spvar'] == spvar_name)
        & (cf_df['cov'] == cov_name)
    )
    if not mask.any():
        available = cf_df.loc[
            (cf_df['dv'] == dv_name) & (cf_df['spvar'] == spvar_name), 'cov'
        ].tolist()
        raise ValueError(
            f"NO CF RESULT for dv='{dv_name}', spvar='{spvar_name}', cov='{cov_name}'.\n"
            f"  Available covariates for this dv×spvar: {available}"
        )
    row = _select_single_row(cf_df, mask,
                             context=f"CF result dv='{dv_name}', spvar='{spvar_name}', cov='{cov_name}'")
    est = float(point_estimate(row, source='single-path counterfactual CSV (COUNTERFACTUAL_CSV)'))
    lci = float(row['hdi_lower_94'])
    uci = float(row['hdi_upper_94'])
    direction = '>' if est >= 0 else '<'
    prob = float(row['prob_above_0']) if est >= 0 else float(row['prob_below_0'])
    n = int(row['N_obs'])
    return (
        f"Δ = {_fmt_num(est)}, P(Δ{direction}0) = {_fmt_prob(prob)}, "
        f"94% HDI [{_fmt_num(lci)}, {_fmt_num(uci)}]"
    ), n, est


def _cf_result_string_or_fallback(cf_df, spvar_name, preferred_cov, fallback_cov, dv_name):
    """
    Try preferred_cov first; if missing, fall back to fallback_cov with a warning.
    Returns (result_str, n_obs, estimate, used_cov).
    """
    try:
        s, n, e = _cf_result_string(cf_df, spvar_name, preferred_cov, dv_name)
        return s, n, e, preferred_cov
    except ValueError:
        print(
            f"WARNING: CF result not found for cov='{preferred_cov}'; "
            f"falling back to '{fallback_cov}' for dv='{dv_name}', spvar='{spvar_name}'."
        )
        s, n, e = _cf_result_string(cf_df, spvar_name, fallback_cov, dv_name)
        return f"{s} [full-sample fallback — {preferred_cov} not available]", n, e, fallback_cov


def _cf_ppa_pct_parts(cf_df, spvar_name, cov_name, dv_name):
    """
    Like _cf_result_string but scales Δ and HDI by ×100 for probability-scale outcomes.
    Use for hppd_binary DVs where marginal contrasts are in probability units (0–1).

    Returns (full_pct_str, delta_abs_pct_str, otherstats_str, n_obs, est_raw).
      full_pct_str:      "Δ = -12.3%, P(Δ<0) = 95.2%, 94% HDI [-20.1%, -4.5%]"
      delta_abs_pct_str: "12.3%"  — absolute value, for inline text ("a 12.3% lower...").
                         Trailing ".0" is trimmed (8.0 → "8%") to match the
                         manuscript; the full/HDI strings keep one decimal.
      otherstats_str:    "P(Δ<0) = 95.2%, 94% HDI [-20.1%, -4.5%]"  — for parenthetical
      n_obs:             int, sample size
      est_raw:           float, raw (un-scaled) estimate
    """
    mask = (
        (cf_df['dv'] == dv_name)
        & (cf_df['spvar'] == spvar_name)
        & (cf_df['cov'] == cov_name)
    )
    if not mask.any():
        available = cf_df.loc[
            (cf_df['dv'] == dv_name) & (cf_df['spvar'] == spvar_name), 'cov'
        ].tolist()
        raise ValueError(
            f"NO CF RESULT for dv='{dv_name}', spvar='{spvar_name}', cov='{cov_name}'.\n"
            f"  Available covariates for this dv×spvar: {available}"
        )
    row = _select_single_row(cf_df, mask,
                             context=f"CF ppa_pct dv='{dv_name}', spvar='{spvar_name}', cov='{cov_name}'")
    est = float(point_estimate(row, source='single-path counterfactual CSV (COUNTERFACTUAL_CSV)'))
    lci = float(row['hdi_lower_94'])
    uci = float(row['hdi_upper_94'])
    direction = '>' if est >= 0 else '<'
    prob = float(row['prob_above_0']) if est >= 0 else float(row['prob_below_0'])
    n = int(row['N_obs'])
    otherstats_str = (
        f"P(Δ{direction}0) = {_fmt_prob(prob)}, "
        f"94% HDI [{100 * lci:.1f}%, {100 * uci:.1f}%]"
    )
    full_str = f"Δ = {100 * est:.1f}%, {otherstats_str}"
    delta_abs_str = f"{_fmt_trim_num(abs(100 * est), 1)}%"
    return full_str, delta_abs_str, otherstats_str, n, est


def _cf_path_string(cf_path_df, effect_label):
    """
    A/B/C' path string from path_counterfactual_summary.csv.
    effect_label: 'A path', 'B path', or "C' path".
    Columns: effect, median, hdi_low, hdi_high, p_above_0, p_below_0.
    """
    row = _select_single_row(
        cf_path_df, cf_path_df['effect'] == effect_label,
        context=f"CF path effect='{effect_label}'",
    )
    est = float(point_estimate(row, source='path_counterfactual_summary.csv'))
    lci = float(row['hdi_low'])
    uci = float(row['hdi_high'])
    direction = '>' if est >= 0 else '<'
    prob = float(row['p_above_0']) if est >= 0 else float(row['p_below_0'])
    return (
        f"Δ = {_fmt_num(est)}, P(Δ{direction}0) = {_fmt_prob(prob)}, "
        f"94% HDI [{_fmt_num(lci)}, {_fmt_num(uci)}]"
    )


def _cf_nie_string(mc_df):
    """
    Indirect effect (NIE) string from mc_mediation_summary.csv.
    Selects the row whose 'effect' column starts with 'NIE'.
    Columns: effect, median, hdi_low, hdi_high, p_above_0, p_below_0.
    """
    nie_rows = mc_df[mc_df['effect'].str.startswith('NIE')]
    if len(nie_rows) == 0:
        raise ValueError("No NIE row found in mc_mediation_summary.")
    row = nie_rows.iloc[0]
    est = float(point_estimate(row, source='mc_mediation_summary.csv', mc_integrated=True))
    lci = float(row['hdi_low'])
    uci = float(row['hdi_high'])
    direction = '>' if est >= 0 else '<'
    prob = float(row['p_above_0']) if est >= 0 else float(row['p_below_0'])
    # Reported as "Δmed", not "Δ": this is the posterior MEDIAN
    # (MC_EFFECT_POINT_ESTIMATE_COL), whereas every other Δ in this file is the
    # posterior mean. The label is explicit so the two are never read as the same
    # summary. The HDI and P(Δ) are quantile-based and unaffected, so they keep
    # the plain Δ. NOTE: the manuscript-sync tool in the master repository
    # (06_submission/reviewer_comments/update_manuscript_stats.py) matches these
    # tokens with a regex — its STAT_PAT must list "Δmed = " before the bare
    # "Δ = " branch, or these statistics stop syncing to the manuscript.
    return (
        f"Δmed = {_fmt_num(est)}, P(Δ{direction}0) = {_fmt_prob(prob)}, "
        f"94% HDI [{_fmt_num(lci)}, {_fmt_num(uci)}]"
    )


# ── Section order ────────────────────────────────────────────────

# The recruitment paragraph opens the Results section and carries no heading.
LEAD_PARAGRAPH_VAR = "recruitment_results_text"

# (manuscript heading, narrative variable) in manuscript order.  The heading
# string is what the output file prints, so it must match the manuscript exactly.
RESULT_SECTION_ORDER = [
    ("Demographics, clinical, and SP use history", "clinical_demographic_results"),
    ("PPA History & Current PPAs", "ppa_history_results"),
    ("Earlier age at first SP use is associated with a greater lifetime risk of SP-associated PPAs.", "ppa_hx_sp_results"),
    ("Higher average SP doses are associated with more current PPAs", "caps_sp_results"),
    ("Lower visual detection thresholds and higher VCH Rates are associated with both PPA history and current CAPS visual symptoms.", "vch_behavior_results"),
    ("Decreased decision precision is predictive of both prior SP-associated PPA risk and current PPAs, and increased prior weighting may be associated with current PPAs.", "vch_computations_results"),
    ("Decreased decision precision, PPA history, and current PPAs are related to diminished criterion, signal-versus-noise discriminability, and more confident VCHs.", "beta_sdt_results"),
]


# ================================================================
# SECTION 4: VARIABLE CALCULATIONS
# ================================================================

# ── Recruitment ──────────────────────────────────────────────────

df_eligibility = df_recruit.copy()
df_eligibility["ineligibile_reason"] = ""

for index, row in df_eligibility[
    ~(df_eligibility["qc_passed"] == 1) & (df_eligibility["screening_survey_complete"] > 0)
].iterrows():
    df_eligibility = df_eligibility.copy()
    if row["screening_pass"] < 1:
        df_eligibility.loc[index, "ineligibile_reason"] += "Fraud-associated phone # or IP associated"
    elif pd.isna(row["geo_crit"]):
        df_eligibility.loc[index, "ineligibile_reason"] += "Fraud-associated phone # or IP associated"
    else:
        nonnegotiables = False
        if row["no_computer"] > 0:
            nonnegotiables = True
            df_eligibility.loc[index, "ineligibile_reason"] += "No computer"
        if row["english_fluency"] < 1:
            nonnegotiables = True
            df_eligibility.loc[index, "ineligibile_reason"] += "Non-English speaking"
        if not nonnegotiables:
            majorcriteria = False
            if row["age_v2"] > 65:
                majorcriteria = True
                df_eligibility.loc[index, "ineligibile_reason"] += ">65 years old"
            if row["age_v2"] < 18:
                majorcriteria = True
                df_eligibility.loc[index, "ineligibile_reason"] += "<18 years old"
            if row["cognition_screener_v2"] > 0:
                majorcriteria = True
                df_eligibility.loc[index, "ineligibile_reason"] += "Neurocognitive Impairment"
            if row["seizure_hx_v2"] > 0:
                majorcriteria = True
                df_eligibility.loc[index, "ineligibile_reason"] += "Epilepsy"
            if row["intox_screen_v2"] > 0:
                majorcriteria = True
                df_eligibility.loc[index, "ineligibile_reason"] += "Active intoxication "
            if row["raven_total_score_v2"] < 1:
                majorcriteria = True
                df_eligibility.loc[index, "ineligibile_reason"] += "Low RAVEN score"
            if not majorcriteria:
                if (
                    (row["activecannabisuse_lastuse"] < 28 and row["cannabis_frequency"] > 9)
                    or (row["activecannabisuse_lastuse"] < 14 and row["cannabis_frequency"] in [8, 9])
                    or (row["activecannabisuse_lastuse"] < 7 and row["cannabis_frequency"] == 7)
                    or row["activecannabisuse_lastuse"] < 3
                ):
                    df_eligibility.loc[index, "ineligibile_reason"] += "Recent heavy cannabis use"
                if row["atypical_since_sp"] > 0:
                    df_eligibility.loc[index, "ineligibile_reason"] += "Atypical psychedelic more recent than SP"
                if (row["psycheduse_yn"] > 1) or (row["sp_naiive"] < 1):
                    df_eligibility.loc[index, "ineligibile_reason"] += "No SP Use "

df_eligibility["ineligibile_reason"] = df_eligibility["ineligibile_reason"].str.strip()
df_eligibility = df_eligibility[df_eligibility["ineligibile_reason"] != ""].copy()
print(f"df_eligibility: {len(df_eligibility)} ineligible participants")

studentswhocompleted = (
    (df_recruit["student_yn"] == 1)
    & (df_recruit["qc_passed"] > 0)
    & (df_recruit["task_data_prltask_present"] == 1)
)

# Student task-fail records: non-SP students (student_yn==1, psycheduse_yn==2)
# who failed QC (qc_passed<1) but for whom the RA did not flag bad data
# (qc_bad_data not set). These students failed the behavioral tasks; they are
# awarded credit per institutional requirement but their data cannot be used.
# They are not genuine QC failures — there is no recoverable data-quality reason.
#
# The three filters work together to isolate exactly this group:
#   psycheduse_yn==2 — restricts to non-SP users; no SP-related QC issue can apply,
#     and also excludes students with atypical psychedelic use (psycheduse_yn==3)
#     who have recoverable QC reasons in automated checks.
#   qc_bad_data.fillna(0) < 1 — RA did not flag the data as bad quality. This is
#     the key condition that distinguishes these 4 task-fail students from record
#     2028, which matches the first two conditions but has qc_bad_data==1 because
#     the post-hoc audit found critical data inconsistencies (SP type mismatch,
#     bizarre route of administration). See consort_diagram.py for details.
RSTF_MASK = (
    (df_recruit['student_yn'] == 1)
    & (df_recruit['psycheduse_yn'] == 2)
    & (df_recruit['qc_passed'] < 1)
    & (df_recruit['raven_total_score_v2'] >= 1)
    & (df_recruit['raven_total_score_v2'].notna())
    & (df_recruit['qc_bad_data'].fillna(0) < 1)
)

# continue_date and timestamp_survey_bl are wall-clock dates and are NOT shipped.
# Both were only ever used through a threshold -- `continue_date < today` and
# `days_since_survey > 14` -- and both comparisons had saturated (734/734 and
# 455/455) well before release, so the recruitment CSV carries the frozen
# booleans instead. See 06_submission/deidentify_recruit_csv.py.
df_recruit["continue_date_passed"] = df_recruit["continue_date_passed"] == 1
df_recruit["timesincesurveystart"] = df_recruit["timesincesurveystart_gt14"]

eligibiles_nocomplete = df_recruit[
    (df_recruit["ineligibile_reason"] == "Eligibile") & (df_recruit["honesty_qc"].isna())
] if "ineligibile_reason" in df_recruit.columns else df_recruit[
    (~df_recruit["si_2_v2"].isna()) & (df_recruit["honesty_qc"].isna())
]

passed_qc_count   = len(df_recruit[(df_recruit["qc_passed"] > 0) | (df_recruit["salvage_yn"] == 1)])
missing_raven_count = len(df_recruit[
    (df_recruit["qc_passed"] > 0)
    & ((df_recruit["raven_total"] < 1) | (df_recruit["raven_total"].isna()))
])

df_rec_count = pd.DataFrame({
    "Stage": [
        "Opened Consent", "Signed Consent", "Screened", "Eligible", "Completed",
        "Failed QC", "Timed Out but Enough Data to Analyze",
        "Waiting to do Longitudinal Study", "Lost to F/U (Longitudinal Study)",
        "Timed Out (>72hr)", "Passed QC", "Missing RAVEN", "Final Dataset",
    ],
    "Participants": [np.nan] * 13,
})
df_rec_count.loc[df_rec_count["Stage"] == "Opened Consent",   "Participants"] = len(df_recruit)
df_rec_count.loc[df_rec_count["Stage"] == "Signed Consent",   "Participants"] = len(df_recruit[df_recruit["consent_baseline_complete"] > 0])
df_rec_count.loc[df_rec_count["Stage"] == "Screened",         "Participants"] = len(df_recruit[~(df_recruit["screening_survey_complete"] < 2)])
df_rec_count.loc[df_rec_count["Stage"] == "Eligible",         "Participants"] = (
    len(df_recruit[~(df_recruit["screening_survey_complete"] < 2)])
    - len(df_eligibility[df_eligibility["si_2_v2"].isna()])
)
df_rec_count.loc[df_rec_count["Stage"] == "Completed",        "Participants"] = len(df_recruit[(df_recruit["honesty_qc"].notna()) | studentswhocompleted])
df_rec_count.loc[df_rec_count["Stage"] == "Failed QC",        "Participants"] = -1 * len(df_recruit[(df_recruit["qc_passed"] < 1) & (df_recruit["honesty_qc"].notna() | (df_recruit["student_yn"] > 0)) & (df_recruit["raven_total_score_v2"] >= 1) & (df_recruit["raven_total_score_v2"].notna()) & ~RSTF_MASK])
df_rec_count.loc[df_rec_count["Stage"] == "Timed Out but Enough Data to Analyze", "Participants"] = len(df_recruit[df_recruit["salvage_yn"] == 1])
df_rec_count.loc[df_rec_count["Stage"] == "Waiting to do Longitudinal Study",     "Participants"] = -1 * len(eligibiles_nocomplete[(eligibiles_nocomplete["waiting_emailed_yn"] > 0) & (eligibiles_nocomplete["continue_date_passed"] == False)])
df_rec_count.loc[df_rec_count["Stage"] == "Lost to F/U (Longitudinal Study)",     "Participants"] = -1 * len(eligibiles_nocomplete[(eligibiles_nocomplete["waiting_emailed_yn"] > 0) & (eligibiles_nocomplete["continue_date_passed"] == True)])
df_rec_count.loc[df_rec_count["Stage"] == "Timed Out (>72hr)",                    "Participants"] = -1 * len(eligibiles_nocomplete[(eligibiles_nocomplete["waiting_emailed_yn"].isna()) & (df_recruit["timesincesurveystart"] == 1)])
df_rec_count.loc[df_rec_count["Stage"] == "Passed QC",        "Participants"] = passed_qc_count
df_rec_count.loc[df_rec_count["Stage"] == "Missing RAVEN",    "Participants"] = -1 * missing_raven_count
df_rec_count.loc[df_rec_count["Stage"] == "Final Dataset",    "Participants"] = passed_qc_count - missing_raven_count
print("Recruitment pipeline:\n", df_rec_count.to_string(index=False))

RS  = len(df_recruit[~(df_recruit['screening_survey_complete'] < 2)])
# RE: eligible participants = those who completed screening minus those the loop
# identified as truly ineligible (df_eligibility rows with si_2_v2 null = 757).
# The old definition (si_2_v2 not null = 448) undercounted by ~44 participants
# who passed every screening criterion but never answered si_2_v2 because they
# were eligible early-cohort participants or dropped out before that item.
RE  = RS - len(df_eligibility[df_eligibility["si_2_v2"].isna()])
RSTF = int(RSTF_MASK.sum())
# Students whose honesty_qc is filled in also reach the end-of-survey QC items.
# Printed so the two populations can be told apart when auditing the flow.
_rstf_in_rq = int(df_recruit[RSTF_MASK & df_recruit['honesty_qc'].notna()].shape[0])
print(f"Student task-fail records (RSTF={RSTF}, in_RQ={_rstf_in_rq}): "
      f"IDs = {sorted(df_recruit.loc[RSTF_MASK, 'record_id'].tolist())}")

df_passed = df_recruit[df_recruit['qc_passed'] > 0]
dishonest = df_passed[df_passed['honesty_qc'] <= 2]
attention   = df_passed[df_passed['attentionq_qc'] <= 2]['record_id'].tolist()
effort      = df_passed[df_passed['effort_qc'] <= 2]['record_id'].tolist()
distraction = df_passed[df_passed['distraction_qc'] <= 2]['record_id'].tolist()
failed_twoplus = list(
    set(attention) & set(effort)
    | set(attention) & set(distraction)
    | set(effort) & set(distraction)
)
if len(failed_twoplus) > 0 or len(dishonest) > 0:
    # WARNING ONLY — df is NOT mutated here. The shipped dataframe is the
    # canonical data source; any QC-based exclusion belongs upstream in the data
    # export, not here. Neither branch fires on the current data.
    print(f"WARNING: {len(dishonest)} participants failed the honesty check and "
          f"{len(failed_twoplus)} failed 2+ QC checks: {failed_twoplus}")
    print("  → df has NOT been filtered; apply the exclusion upstream in the data export.")

final_dfrec_count = int(df_rec_count.loc[df_rec_count['Stage'] == 'Final Dataset', 'Participants'].values[0])
final_df_count    = df[df['raven_total'] > 0].shape[0]
if final_dfrec_count != final_df_count:
    print(f"Warning: Final N mismatch — df_rec_count={final_dfrec_count}, df={final_df_count}")
else:
    RG = final_df_count

RSP      = len(df[df['psycheduse_yn'] == "Yes"])

# ── Stage-conditioned recruitment counts ────────────────────────────────────
# RQ / RFQ / RTS above are *marginal*: each is computed over everyone in
# df_recruit rather than within the preceding stage, so they do not subtract.
# The counts below place every participant at exactly one stage, so the chain
# closes exactly:
#
#   screened − ineligible                     = eligible
#   eligible − failed minimum study measures  = completed minimum measures
#   completed minimum measures − failed QC    = final analytic sample
#
# The quality-control group is restricted to participants found eligible at
# screening. Someone screened ineligible who then also failed QC leaves the flow
# at the screening step and is counted there. consort_diagram.py applies the
# same si_2_v2 restriction, so the supplementary table and this count describe
# the same population. consort_diagram.py recomputes all of this independently
# and asserts the closure — keep the three in sync.

_IN_DF      = df_recruit['record_id'].isin(set(df['record_id']))
_SCREENED   = ~(df_recruit['screening_survey_complete'] < 2)
_NO_RAVEN   = (df_recruit['raven_total'] < 1) | (df_recruit['raven_total'].isna())
_SALVAGED   = df_recruit['salvage_yn'] == 1
# Truly ineligible = assigned a screening reason AND not found eligible by REDCap.
# The si_2_v2 filter drops the early-cohort participants the loop mislabels as
# fraud because geo_crit did not yet exist. Same filter as consort_diagram.py.
_INELIGIBLE = df_recruit['record_id'].isin(
    set(df_eligibility.loc[df_eligibility['si_2_v2'].isna(), 'record_id'])
)
_ELIGIBLE   = _SCREENED & ~_INELIGIBLE
_QC_FAILED  = (
    (df_recruit['qc_passed'] < 1)
    & (df_recruit['honesty_qc'].notna() | (df_recruit['student_yn'] > 0))
    & (df_recruit['raven_total_score_v2'] >= 1)
    & (df_recruit['raven_total_score_v2'].notna())
    & ~RSTF_MASK
    & df_recruit['si_2_v2'].notna()      # eligible at screening
)

_LOST         = _ELIGIBLE & ~_IN_DF
_QC_ELIGIBLE  = _QC_FAILED & _ELIGIBLE
_MIN_MEASURES = _LOST & ~_QC_ELIGIBLE
# Disjoint by fixed precedence so the parts sum exactly.
_B_NO_RAVEN   = _MIN_MEASURES & _NO_RAVEN
_B_NO_RECORD  = _MIN_MEASURES & df_recruit['qc_passed'].isna() & ~_B_NO_RAVEN
_B_TASK       = _MIN_MEASURES & ~_B_NO_RAVEN & ~_B_NO_RECORD

RE_INELIGIBLE   = int(_INELIGIBLE.sum())        # 757
RE_MINMEASURES  = int(_MIN_MEASURES.sum())      # 195
RE_NORAVEN      = int(_B_NO_RAVEN.sum())        # 68 — screening battery never scored
RE_NORECORD     = int(_B_NO_RECORD.sum())       # 118 — no QC record; did not complete
RE_TASKFAIL     = int(_B_TASK.sum())            # 9
RE_AFTERMIN     = RE - RE_MINMEASURES           # 297
RE_QCFAIL       = int(_QC_ELIGIBLE.sum())       # 69 (cf. marginal RFQ = 83)
RE_SALVAGED     = int((_ELIGIBLE & _SALVAGED & _IN_DF).sum())   # 27

# The recruitment paragraph and the CONSORT diagram describe the same flow by
# two orthogonal routes. The diagram subtracts exclusions from the eligible pool
# (492 - 195 = 297 analysable, - 69 QC = 228). The paragraph instead counts the
# completers, removes the quality-control failures, and adds the salvage records
# back at the end (270 - 69 = 201, + 27 = 228), preserving the sentence structure
# of the submitted manuscript. Both close on 228; the assertions below tie them
# together so the two descriptions cannot drift apart.
_ANALYSABLE     = _ELIGIBLE & ~_MIN_MEASURES                    # 297
_SALV           = _ELIGIBLE & _SALVAGED & _IN_DF                # 27
RE_COMPLETERS   = int((_ANALYSABLE & ~_SALV).sum())             # 270
RE_RETAINED     = RE_COMPLETERS - RE_QCFAIL                     # 201

_NUMBER_WORDS = {
    1: "One", 2: "Two", 3: "Three", 4: "Four", 5: "Five", 6: "Six",
    7: "Seven", 8: "Eight", 9: "Nine", 10: "Ten", 11: "Eleven",
    12: "Twelve", 13: "Thirteen", 14: "Fourteen", 15: "Fifteen",
    16: "Sixteen", 17: "Seventeen", 18: "Eighteen", 19: "Nineteen",
    20: "Twenty",
}

_TENS_WORDS = {
    20: "Twenty", 30: "Thirty", 40: "Forty", 50: "Fifty",
    60: "Sixty", 70: "Seventy", 80: "Eighty", 90: "Ninety",
}


def spell_sentence_initial(n):
    """Spell *n* for use at the start of a sentence.

    Journal style forbids opening a sentence with a numeral. Covers 1-99, which
    is every count the narrative currently opens a sentence with. Anything
    outside that range falls through to digits rather than being silently
    reworded — if that ever fires, rephrase the sentence instead of extending
    this blindly.
    """
    n = int(n)
    if n in _NUMBER_WORDS:
        return _NUMBER_WORDS[n]
    if 21 <= n <= 99:
        tens, units = divmod(n, 10)
        return f"{_TENS_WORDS[tens * 10]}-{_NUMBER_WORDS[units].lower()}"
    return str(n)

# VCH task QC exclusions: attempted the task vs. yielded a usable threshold.
RVCH_QC_EXCLUDED    = int(df['task_data_vch_short_psychedelic_bl'].notna().sum()
                          - df['vch_threshold'].notna().sum())
# Opens a sentence in the recruitment paragraph, so it is spelled out.
VCH_QC_EXCLUDED_WORD = spell_sentence_initial(RVCH_QC_EXCLUDED)

_flow_checks = [
    ("screened - ineligible = eligible", RS - RE_INELIGIBLE, RE),
    ("eligible - failed minimum measures = completed minimum measures",
     RE - RE_MINMEASURES, RE_AFTERMIN),
    ("minimum-measure buckets sum",
     RE_NORAVEN + RE_NORECORD + RE_TASKFAIL, RE_MINMEASURES),
    ("completers + salvaged = analysable records",
     RE_COMPLETERS + RE_SALVAGED, RE_AFTERMIN),
    # Every quality-control failure sits inside the completers, so the two
    # routes stay disjoint and the salvage records are never double-subtracted.
    ("no quality-control failure is a salvage record",
     int((_QC_ELIGIBLE & _SALV).sum()), 0),
]
if "RG" in dir():
    _flow_checks.append(("completed minimum measures - QC failures = final",
                         RE_AFTERMIN - RE_QCFAIL, RG))
    _flow_checks.append(("completers - QC failures + salvaged = final",
                         RE_RETAINED + RE_SALVAGED, RG))
_flow_bad = [(lbl, got, want) for lbl, got, want in _flow_checks if got != want]
if _flow_bad:
    print("WARNING: recruitment flow no longer closes —")
    for lbl, got, want in _flow_bad:
        print(f"  {lbl}: got {got}, expected {want}")
else:
    print(f"Recruitment flow closes: {RS} → {RE} → {RE_AFTERMIN} → "
          f"{RE_AFTERMIN - RE_QCFAIL}")
    print(f"  narrative route: {RE} → {RE_COMPLETERS} completers "
          f"− {RE_QCFAIL} QC = {RE_RETAINED}, + {RE_SALVAGED} salvaged = "
          f"{RE_RETAINED + RE_SALVAGED}")


# ── Clinical & Demographic ───────────────────────────────────────

_sp_use = pd.to_numeric(df_sp["psycheduse_life_nomic"], errors="coerce")
_edu    = pd.to_numeric(df_sp["highest_education"], errors="coerce")

FOUR_YR_OR_MORE_PCT = int(round((_edu.gt(6).mean(skipna=True)) * 100))

def _corr_str(x, y):
    tmp = pd.DataFrame({"x": x, "y": y}).dropna()
    rho, p = spearmanr(tmp["x"], tmp["y"], nan_policy="omit")
    return f"ρ = {rho:.3f}, p {_fmt_p(p)}"

M_CORR                = _corr_str(df_sp['sex_v2'], _sp_use)
age_v2_CORR           = _corr_str(pd.to_numeric(df_sp["age_v2"], errors="coerce"), _sp_use)
drugs       = ['alc', 'mj', 'stimulants', 'sedatives', 'atypicals']
drugs_month = [f"{drug}_month_yn" for drug in drugs]
drugs_life  = [f"{drug}_lifetime" for drug in ['alc', 'mj']] + [f"{drug}_life_yn" for drug in ['stimulants', 'atypicals', 'sedatives']]

df['pdu_life_yn']  = (df[drugs_life] > 0).any(axis=1).astype(int)
df['pdu_month_yn'] = (df[drugs_month] > 0).any(axis=1).astype(int)

def _spearman_vs_life_nomic(var_name):
    rho, p = spearmanr(df_sp[var_name], df_sp["psycheduse_life_nomic"], nan_policy="omit")
    return f"ρ = {rho:.3f}, p {_fmt_p(p)}"

mj_lifetime_corr       = _spearman_vs_life_nomic("mj_lifetime")
atypicals_life_yn_corr = _spearman_vs_life_nomic("atypicals_life_yn")
stimulants_life_yn_corr = _spearman_vs_life_nomic("stimulants_life_yn")
atypicals_month_yn_corr = _spearman_vs_life_nomic("atypicals_month_yn")
stimulants_month_yn_corr = _spearman_vs_life_nomic("stimulants_month_yn")

print("FOUR_YR_OR_MORE_PCT:", FOUR_YR_OR_MORE_PCT)
print("M_CORR:", M_CORR, "| age_v2_CORR:", age_v2_CORR)
print("mj_lifetime_corr:", mj_lifetime_corr)
print("atypicals_life_yn_corr:", atypicals_life_yn_corr, "| stimulants_life_yn_corr:", stimulants_life_yn_corr)
print("atypicals_month_yn_corr:", atypicals_month_yn_corr, "| stimulants_month_yn_corr:", stimulants_month_yn_corr)


# ── SP Use ───────────────────────────────────────────────────────

_age_first = pd.to_numeric(df_sp["psychedelic_age"], errors="coerce")
_years_use = pd.to_numeric(df_sp["psyched_yearsofuse"], errors="coerce")
_life_uses = pd.to_numeric(df_sp["psycheduse_life_nomic"], errors="coerce")
_avg_dose  = pd.to_numeric(df_sp["avg_life_dose"], errors="coerce")
_perceived = pd.to_numeric(df_sp['perceived_benefit'], errors="coerce")

psyched_yearsofuse_PCT = _pct_str(df_sp, _years_use < 10, _years_use.notna())
PERCENT_ABOVE_40       = _pct_str(df_sp, _age_first >= 40, _age_first.notna())

_nov = df_sp["sp_experience"].astype(str).str.strip().str.lower().eq("novice")
_sp_yes = df["psycheduse_yn"].astype(str).str.strip().eq("Yes")
_age_first_all_sp = pd.to_numeric(df["psychedelic_age"], errors="coerce")
PERCENT_BELOW_25_EXP = _pct_str(df, _age_first_all_sp < 25, _sp_yes & _age_first_all_sp.notna())
PERCENT_BELOW_25_NOV = _pct_str(df_sp, _age_first < 25, _nov & _age_first.notna())
psychedelic_age_sd  = f"{df_sp["psychedelic_age"].std():.1f}"
psycheduse_life_nomic_mean = f"{_life_uses.mean():.1f}"
psycheduse_life_nomic_sd = f"{_life_uses.std():.1f}"
avg_life_dose_sd = f"{_avg_dose.std():.1f}"
avg_life_dose_mean = f"{_avg_dose.mean():.1f}"


perceived_benefit_pct_below50 = _pct_str(df_sp, _perceived < 50, _perceived.notna())
perceived_benefit_pct_50      = _pct_str(df_sp, pd.Series(np.isclose(_perceived, 50.0, atol=1e-8), index=df_sp.index), _perceived.notna())

avg_life_dose_OVER300_PCT = str(int(round((_avg_dose[_avg_dose.notna()] > 300).mean() * 100)))

# Lifetime SP use range (min and max number of uses across SP users)
sp_lifetime_uses_min = int(_life_uses.min())
sp_lifetime_uses_max = int(_life_uses.max())

# Non-serotonergic (atypical) psychedelic use: lifetime and past-month prevalence among SP users
_nonsp_psychedelic_life_yn  = pd.to_numeric(df_sp['atypicals_life_yn'],  errors='coerce') > 0
_nonsp_psychedelic_month_yn = pd.to_numeric(df_sp['atypicals_month_yn'], errors='coerce') > 0
nonsp_psychedelic_lifetime_pct   = int(round(_nonsp_psychedelic_life_yn.mean()  * 100))
nonsp_psychedelic_past_month_pct = int(round(_nonsp_psychedelic_month_yn.mean() * 100))

# Cannabis use: lifetime prevalence and past-month rate among lifetime cannabis users (SP users)
_cannabis_lifetime_users = pd.to_numeric(df_sp['mj_lifetime'],  errors='coerce') > 0
_cannabis_past_month_yn  = pd.to_numeric(df_sp['mj_month_yn'],  errors='coerce') > 0
cannabis_lifetime_pct                    = int(round(_cannabis_lifetime_users.mean() * 100))
cannabis_past_month_pct_of_lifetime_users = int(round(_cannabis_past_month_yn[_cannabis_lifetime_users].mean() * 100))

print("psyched_yearsofuse_PCT:", psyched_yearsofuse_PCT, "| PERCENT_ABOVE_40:", PERCENT_ABOVE_40)
print("PERCENT_BELOW_25_EXP:", PERCENT_BELOW_25_EXP, "| PERCENT_BELOW_25_NOV:", PERCENT_BELOW_25_NOV)
print("avg_life_dose_OVER300_PCT:", avg_life_dose_OVER300_PCT)
print("perceived_benefit_pct_below50:", perceived_benefit_pct_below50, "| pct_50:", perceived_benefit_pct_50)
print(f"sp_lifetime_uses_min={sp_lifetime_uses_min}, sp_lifetime_uses_max={sp_lifetime_uses_max}")
print(f"nonsp_psychedelic_lifetime_pct={nonsp_psychedelic_lifetime_pct}%, nonsp_psychedelic_past_month_pct={nonsp_psychedelic_past_month_pct}%")
print(f"cannabis_lifetime_pct={cannabis_lifetime_pct}%, cannabis_past_month_pct_of_lifetime_users={cannabis_past_month_pct_of_lifetime_users}%")


# ── PPA History ──────────────────────────────────────────────────

persist_vis_yn_num_PCT = _num_pct_from_mask(df["persist_vis_yn"] == 1, df["persist_vis_yn"].notna().sum())

_psychdoses = pd.to_numeric(df_hppd["persistvis_psychdoses"], errors="coerce")
_q1, _q3   = _psychdoses.quantile(0.25), _psychdoses.quantile(0.75)
persistvis_psychdoses_median    = f"median {_fmt_trim_num(_psychdoses.median(), 2)}"
persistvis_psychdoses_iqr       = f"IQR [{_fmt_trim_num(_q1, 2)},{_fmt_trim_num(_q3, 2)}]"
persistvis_psychdoses_median_iqr = f"{persistvis_psychdoses_median}; {persistvis_psychdoses_iqr}"
_n1, _d1 = int((_psychdoses == 1).sum()), int(_psychdoses.notna().sum())
_one_dose_pct = (100 * _n1 / _d1) if _d1 > 0 else float("nan")
# Opens a sentence in the narrative, so the count is spelled out.
persistvis_psychdoses_one_dose_sentence = (
    f"{spell_sentence_initial(_n1)} participants ({_fmt_trim_num(_one_dose_pct, 1)}%)"
)

# _persistvis_most_map = {
#     13: "I have never experienced any of the above visual effects",
#     1:  "Halos or auras around things",
#     2:  "Stationary things appearing to move, breathe, grow, or shrink",
#     3:  "Moving objects appear to not be moving",
#     4:  "Afterimages left behind moving objects",
#     5:  "Brighter or more intense colors",
#     6:  "Seeing patterns or textures that aren't there with eyes open",
#     7:  "Seeing patterns or textures that aren't there with eyes closed",
#     8:  "Seeing objects that aren't really there",
#     9:  "Increased intensity in oscillating or flashing lights (TV, lightbulbs, etc.)",
#     10: "Distortion, movement, or patterns in grids, gratings, or closely spaced lines",
#     11: "Noticing more things in your environment",
#     12: "Things just looked different",
# }

_persistvis_most_map = {
    13: "I have never experienced any of the above visual effects",
    1:  "Halos or auras",
    2:  "Stationary objects appear to move",
    3:  "Moving objects appear to not move",
    4:  "Afterimages for moving objects",
    5:  "Greater color intensity",
    6:   "Patterns w/eyes open",
    7:  "Patterns w/eyes closed",
    8:  "Objects that aren't really there",
    9:  "Accentuation of light oscillations",
    10:  "Distorted/moving lines/grids",
    11: "Noticing more things in environment",
    12: "Things just look different",
}

_persistvis_most_series = df_hppd["persistvis_most"].dropna()
if len(_persistvis_most_series) == 0:
    persistvis_most_mode     = "n/a"
    persistvis_most_mode_PCT = "n/a"
else:
    _mode_val = _persistvis_most_series.mode().iloc[0]
    _raw_most = _persistvis_most_map.get(_mode_val, str(_mode_val))
    persistvis_most_mode     = _raw_most[:1].lower() + _raw_most[1:]
    persistvis_most_mode_PCT = _num_pct_from_mask(_persistvis_most_series == _mode_val, len(_persistvis_most_series))

_chronicity_series = df_hppd["hppd_true_chronicity"].dropna()
_n_below_5, _d_below_5 = int((_chronicity_series < 5).sum()), int(len(_chronicity_series))
_p_below_5 = (100 * _n_below_5 / _d_below_5) if _d_below_5 > 0 else np.nan
hppd_true_chronicity_BELOW_5_PCT = f"{_n_below_5}; {_fmt_trim_num(_p_below_5, 1)}%"
hppd_true_chronicity_ABOVE_7_PCT = _num_pct_from_mask(_chronicity_series > 7, len(_chronicity_series))

_persistvis_duration_map = {
    1: "brief spurts (seconds to minutes)",
    2: "Longer periods (hours to days)",
    3: "Constant or near-constant (weeks to years)",
}
_dur_series = df_hppd["persistvis_duration"].dropna()
if len(_dur_series) == 0:
    persistvis_duration_MODE         = "n/a"
    persistvis_duration_MODE_num_PCT = "n/a"
else:
    _dur_mode_val = _dur_series.mode().iloc[0]
    _raw_dur = _persistvis_duration_map.get(_dur_mode_val, str(_dur_mode_val))
    persistvis_duration_MODE         = _raw_dur[:1].lower() + _raw_dur[1:] # lowercase first letter
    persistvis_duration_MODE_num_PCT = _num_pct_from_mask(_dur_series == _dur_mode_val, len(_dur_series))

_baggot_total_series  = df_hppd["baggot_total"].dropna()
baggot_total_mode     = _baggot_total_series.mode().iloc[0] if len(_baggot_total_series) else "n/a"
baggot_total_mode_num_PCT  = f"{int(baggot_total_mode)} ({((_baggot_total_series == baggot_total_mode).sum() / len(_baggot_total_series) * 100):.1f}%)"
baggot_total_ABOVE_4_num_PCT = _num_pct_from_mask(_baggot_total_series >= 4, len(_baggot_total_series))

_dist_series = df_hppd["persistvis_distress"].dropna()
persistvis_distress_PCT = _num_pct_from_mask(_dist_series == 2, len(_dist_series),semicolon=True)  # 2 = no distress
_tx_series = df_hppd["persistvis_txseek"].dropna()
persistvis_txseek_PCT = _num_pct_from_mask(_tx_series == 2, len(_tx_series),semicolon=True)        # 2 = no treatment seeking

# Self-reported formal HPPD diagnosis, as a share of the whole analytic sample
# (not of the PPA-history subgroup) — the denominator the manuscript reports.
hppd_diagnosis_num_PCT = _num_pct_from_mask(df["hppd_ever"] > 0, len(df))

print("persist_vis_yn_num_PCT:", persist_vis_yn_num_PCT)
print("persistvis_psychdoses_median_iqr:", persistvis_psychdoses_median_iqr)
print("baggot_total_mode_num_PCT:", baggot_total_mode_num_PCT, "| ABOVE_4:", baggot_total_ABOVE_4_num_PCT)
print("hppd_true_chronicity_BELOW_5_PCT:", hppd_true_chronicity_BELOW_5_PCT)
print("persistvis_duration_MODE:", persistvis_duration_MODE, "|", persistvis_duration_MODE_num_PCT)


# ── CAPS Descriptive ─────────────────────────────────────────────
# Denominator for all percentages in this block: SP users who have a
# non-missing caps_bl_1 value (CAPS completion indicator).
_df_caps  = df_sp[df_sp["caps_bl_1"].notna()].copy()
_n_caps   = len(_df_caps)

# Participants with at least one CAPS vision item endorsed (caps_vision > 0).
_df_caps_pos = _df_caps[_df_caps["caps_vision"] > 0].copy()
_n_caps_pos  = len(_df_caps_pos)

# CAPS vision item indices and derived column lists used throughout this block.
#   binary endorsement  → caps_bl_{x}
#   frequency           → caps_bl_{x}c  (scale 1–5: Hardly at all … All the Time)
#   distress            → caps_bl_{x}a  (scale 0=absent, 1=Not at all … 5=Very)
_CAPS_ITEMS = [4, 26, 31, 23, 19, 22]
_freq_cols  = [f"caps_bl_{x}c" for x in _CAPS_ITEMS]

# caps_pos_majority / caps_pos_pct:
# "majority" or "minority" label and % of CAPS-valid SP users who endorse
# at least one current visual PPA.
caps_pos_majority = "majority" if (_n_caps_pos / _n_caps) > 0.50 else "minority"
caps_pos_pct      = f"{100 * _n_caps_pos / _n_caps:.1f}"

# caps_ppa_pos_pct:
# Of those endorsing any current CAPS item, % who also have prior SP-induced
# PPA history (persist_vis_yn == 1).
caps_ppa_pos_pct = f"{100 * (_df_caps_pos['persist_vis_yn'] == 1).mean():.1f}"

# hppd_current_pct:
# Of those endorsing any current CAPS item, % who currently attribute PPAs
# to prior SP use (hppd_current > 0 → derived from persistvis_time > 8,
# i.e., the "I still experience these effects" response option).
hppd_current_pct = f"{100 * (_df_caps_pos['hppd_current'] > 0).mean():.1f}"

# n_capsnonnan_caps_bl_31 / pct_capsnonnan_caps_bl_31:
# Count and % of CAPS-valid SP users endorsing caps_bl_31 ("seeing things
# that others cannot") — the CAPS item closest to a formed visual hallucination.
n_capsnonnan_caps_bl_31   = int((_df_caps["caps_bl_31"] > 0).sum())
pct_capsnonnan_caps_bl_31 = f"{100 * (_df_caps['caps_bl_31'] > 0).mean():.1f}"

# caps_freq_fives:
# Of those endorsing any current CAPS item, % who rated at least one
# frequency subdimension item at 5 ("All the Time").
caps_freq_fives = f"{100 * (_df_caps_pos[_freq_cols].max(axis=1) > 4).mean():.1f}"

# pct_caps_vision_distress_0 / pct_caps_vision_distress_over3:
# Distress ratings evaluated only over items the participant actually endorsed
# (caps_bl_x == 1).  Two complementary thresholds:
#   distress_0    → max distress across endorsed items < 2 (i.e., all "Not at all")
#   distress_over3 → any endorsed item rated > 3 ("Firmly" or "Very")

def _max_distress_endorsed(row):
    """Maximum distress rating across items the participant endorsed."""
    endorsed_vals = [
        row[f"caps_bl_{x}a"]
        for x in _CAPS_ITEMS
        if row[f"caps_bl_{x}"] == 1
    ]
    return max(endorsed_vals) if endorsed_vals else np.nan

def _any_high_distress(row):
    """True if any endorsed item has distress > 3 ('Firmly' or 'Very')."""
    return any(
        row[f"caps_bl_{x}a"] > 3
        for x in _CAPS_ITEMS
        if row[f"caps_bl_{x}"] == 1
    )

_max_dist_series  = _df_caps_pos.apply(_max_distress_endorsed, axis=1)
_high_dist_series = _df_caps_pos.apply(_any_high_distress, axis=1)

pct_caps_vision_distress_0     = f"{100 * (_max_dist_series < 2).mean():.1f}"
pct_caps_vision_distress_over3 = f"{100 * _high_dist_series.mean():.1f}"

# chi_square_caps_ppa:
# Permutation chi-square test of homogeneity comparing the relative distribution
# of CAPS item endorsements between PPA history+ and PPA history− groups.
# Restricted to participants with ≥1 current CAPS item endorsed (caps_vision > 0)
# so the test reflects *which* items are endorsed (phenomenology), not whether
# items are endorsed at all.
#
# Classical chi-square is unreliable here because ~9 PPA− participants have
# current symptoms (small expected cell counts).  The permutation approach
# sidesteps that assumption: observed χ² is compared against an empirical null
# built by shuffling PPA group labels 10,000 times at the person level,
# which also preserves within-person item correlations.
# Phipson & Smyth (2010) continuity: p = (B+1) / (N_perm+1).

# Pre-compute a (n_people × 6) boolean endorsement matrix for caps_vision > 0.
_pheno_labels = _df_caps_pos["persist_vis_yn"].values
_pheno_matrix = np.array(
    [(_df_caps_pos[f"caps_bl_{x}"].values == 1) for x in _CAPS_ITEMS],
    dtype=float,
).T   # shape: (n_people, 6)

def _chi2_from_labels(labels, matrix):
    """Chi-square statistic for a 2×6 endorsement table given group labels."""
    row_neg = matrix[labels == 0].sum(axis=0)
    row_pos = matrix[labels == 1].sum(axis=0)
    table   = np.array([row_neg, row_pos])
    row_sums = table.sum(axis=1, keepdims=True)
    col_sums = table.sum(axis=0, keepdims=True)
    total    = table.sum()
    if total == 0 or (row_sums == 0).any():
        return np.nan
    expected = row_sums * col_sums / total
    mask = expected > 0
    return float(((table[mask] - expected[mask]) ** 2 / expected[mask]).sum())

_obs_chi2    = _chi2_from_labels(_pheno_labels, _pheno_matrix)
_obs_dof     = (2 - 1) * (len(_CAPS_ITEMS) - 1)   # = 5

_rng         = np.random.default_rng(seed=42)
_n_perm      = 10_000
_perm_chi2s  = np.array([
    _chi2_from_labels(_rng.permuted(_pheno_labels), _pheno_matrix)
    for _ in range(_n_perm)
])
_valid        = _perm_chi2s[~np.isnan(_perm_chi2s)]
_perm_p       = ((_valid >= _obs_chi2).sum() + 1) / (len(_valid) + 1)

chi_square_caps_ppa = (
    f"χ²({_obs_dof}) = {_obs_chi2:.2f}, "
    f"p_permutation {_fmt_p(_perm_p)} (10,000 permutations)"
)

print("caps_pos_pct:", caps_pos_pct, "| majority/minority:", caps_pos_majority)
print("caps_ppa_pos_pct:", caps_ppa_pos_pct)
print("hppd_current_pct:", hppd_current_pct)
print("n_capsnonnan_caps_bl_31:", n_capsnonnan_caps_bl_31, "| pct:", pct_capsnonnan_caps_bl_31)
print("caps_freq_fives:", caps_freq_fives)
print("pct_caps_vision_distress_0:", pct_caps_vision_distress_0,
      "| pct_caps_vision_distress_over3:", pct_caps_vision_distress_over3)
print("chi_square_caps_ppa:", chi_square_caps_ppa)


# ── SP Predictors of PPA Risk ────────────────────────────────────

_hppd_dv            = "hppd_binary"
_hppd_ivtype        = "sp_predictors"
_hppd_nonparam_csv  = (
    f"{_RESULTS}/"
    f"{_hppd_dv}/{_hppd_ivtype}/data_visualization/summary_results/boxplot_grid.csv"
)
_hppd_nonparam_df  = pd.read_csv(_hppd_nonparam_csv)

PA         = _mw_result_string(_hppd_nonparam_df, "psychedelic_age")
PB_lifetime = _mw_result_string(_hppd_nonparam_df, "psycheduse_life_nomic")
PB_dose    = _mw_result_string(_hppd_nonparam_df, "avg_life_dose")
PB         = f"lifetime SP uses: {PB_lifetime}; average SP dose: {PB_dose}"
_cf_hppd_cov = "nice_covariates_spusers"
_cf_hppd_dv  = "hppd_binary"

# Mean age at first SP use, for the "young adulthood (X ± Y years)" parenthetical.
psychedelic_age_mean = f"{df_sp['psychedelic_age'].mean():.1f}"

# All hppd_binary CFs scaled ×100 — Δ and HDI reported as percentages
PC,      PC_abs_pct, _,            _, _ = _cf_ppa_pct_parts(_cf_df, "psychedelic_age_normalized",       _cf_hppd_cov, _cf_hppd_dv)
PD,      _,          _,            _, _ = _cf_ppa_pct_parts(_cf_df, "psycheduse_life_nomic_normalized",  _cf_hppd_cov, _cf_hppd_dv)
PD_dose, _,          _,            _, _ = _cf_ppa_pct_parts(_cf_df, "avg_life_dose_normalized",          _cf_hppd_cov, _cf_hppd_dv)

# IQR-outlier sensitivity, read from the same CF CSV under
# cov='nice_covariates_spusers_iqr'. The drop-% denominator is the matching
# base-model row (_cf_hppd_cov), so both N come from the same CSV and the
# same dv x spvar — never from a separately computed sample size.
PE, _, _, _n_life_iqr,  _ = _cf_ppa_pct_parts(_cf_df, "psycheduse_life_nomic_normalized", "nice_covariates_spusers_iqr", _cf_hppd_dv)
_,  _, _, _n_life_base, _ = _cf_ppa_pct_parts(_cf_df, "psycheduse_life_nomic_normalized", _cf_hppd_cov,                  _cf_hppd_dv)
PEN = f"{100 * (_n_life_base - _n_life_iqr) / _n_life_base:.1f}"

PF, _, _, _n_dose_iqr,  _ = _cf_ppa_pct_parts(_cf_df, "avg_life_dose_normalized", "nice_covariates_spusers_iqr", _cf_hppd_dv)
_,  _, _, _n_dose_base, _ = _cf_ppa_pct_parts(_cf_df, "avg_life_dose_normalized", _cf_hppd_cov,                  _cf_hppd_dv)
PFN = f"{100 * (_n_dose_base - _n_dose_iqr) / _n_dose_base:.1f}"

print("PA:", PA)
print("PB:", PB)
print("PC:", PC, "| PD:", PD, "| PD_dose:", PD_dose)
print("PE:", PE, "| PEN:", PEN, "| PF:", PF, "| PFN:", PFN)

# Figure & table references — SP Predictors of PPA History
sp_predictors_ppa_hx_nonparametric = f"Fig. {FIG_SP_PREDICTORS}a"


# ── SP Predictors of CAPS Vision ───────────────────────────────────

_caps_dv = "caps_vision"

# Age-controlled partial Spearman correlations for the nonparametric CAPS section.
# These match the panel actually shown in the manuscript figure, which is the
# age-controlled grid (results/caps_vision/sp_predictors/data_visualization/
# correlation_grid_age_control.png), not the zero-order correlation_grid.png.
# Sample = SP users with a non-missing caps_bl_1 (n = 130); _partial_spearman_str
# rank-transforms all three variables and then partials age_v2 out via pingouin.
_df_sp_caps = df_sp[df_sp["caps_bl_1"].notna()]
CA          = _partial_spearman_str(_df_sp_caps, "psychedelic_age",       _caps_dv)
CB_lifetime = _partial_spearman_str(_df_sp_caps, "psycheduse_life_nomic", _caps_dv)
CB_dose     = _partial_spearman_str(_df_sp_caps, "avg_life_dose",         _caps_dv)

_cf_caps_cov = "nice_covariates_spusers"   # mediation/forest models are SP users only

CC,       n_age_base,  _ = _cf_result_string(_cf_df, "psychedelic_age_normalized",       _cf_caps_cov, _caps_dv)
CD,       n_dose_base, _ = _cf_result_string(_cf_df, "avg_life_dose_normalized",         _cf_caps_cov, _caps_dv)
CSE_life, n_life_base, _ = _cf_result_string(_cf_df, "psycheduse_life_nomic_normalized", _cf_caps_cov, _caps_dv)


# IQR-outlier sensitivity, from the same CF CSV. n_life_base above is the
# matching nice_covariates_spusers row, so the two N are directly comparable.
CSE_life_iqr, n_life_iqr, _ = _cf_result_string(
    _cf_df, "psycheduse_life_nomic_normalized", "nice_covariates_spusers_iqr", _caps_dv)

print("CA:", CA, "| CB_lifetime:", CB_lifetime, "| CB_dose:", CB_dose)
print("CC:", CC, "| CD:", CD, "| CSE_life:", CSE_life, "| CSE_life_iqr:", CSE_life_iqr)

# ── Spearman correlations: age confounders for caps_sp_results text ───────────
# These explain the discrepancy between nonparametric and regression associations
# with lifetime SP uses and age of first use in the CAPS vision section.

_mask_life_age = df_sp[["psycheduse_life_nomic", "age_v2"]].notna().all(axis=1)
_rho_life_age  = spearmanr(
    df_sp.loc[_mask_life_age, "psycheduse_life_nomic"],
    df_sp.loc[_mask_life_age, "age_v2"],
)
psycheduse_life_nomic_age_v2 = f"ρ = {_fmt_num(_rho_life_age.statistic)}, p {_fmt_p(_rho_life_age.pvalue)}"

_mask_agefirst_age = df_sp[["psychedelic_age", "age_v2"]].notna().all(axis=1)
_rho_agefirst_age  = spearmanr(
    df_sp.loc[_mask_agefirst_age, "psychedelic_age"],
    df_sp.loc[_mask_agefirst_age, "age_v2"],
)
psyched_agefirst_age_v2 = f"ρ = {_fmt_num(_rho_agefirst_age.statistic)}, p {_fmt_p(_rho_agefirst_age.pvalue)}"

_df_sp_caps_age = df_sp[df_sp["caps_bl_1"].notna()][["age_v2", "caps_vision"]].dropna()
_rho_age_caps   = spearmanr(_df_sp_caps_age["age_v2"], _df_sp_caps_age["caps_vision"])
age_v2_caps_corr = f"ρ = {_fmt_num(_rho_age_caps.statistic)}, p {_fmt_p(_rho_age_caps.pvalue)}"

print("psycheduse_life_nomic_age_v2:", psycheduse_life_nomic_age_v2)
print("psyched_agefirst_age_v2:", psyched_agefirst_age_v2)
print("age_v2_caps_corr:", age_v2_caps_corr)

# Figure & table references — SP Predictors of CAPS Vision
sp_predictors_caps_nonparametric = f"Fig. {FIG_SP_PREDICTORS}b"


# ── VCH Behavior ─────────────────────────────────────────────────

_hppd_vch_nonparam_csv = (
    f"{_RESULTS}/"
    f"hppd_binary/vch_behavior/data_visualization/summary_results/boxplot_grid.csv"
)
_caps_vch_nonparam_csv = (
    f"{_RESULTS}/"
    f"caps_vision/vch_behavior/data_visualization/summary_results/correlation_grid.csv"
)
_hppd_vch_nonparam_df = pd.read_csv(_hppd_vch_nonparam_csv)
_caps_vch_nonparam_df = pd.read_csv(_caps_vch_nonparam_csv)

HB_thresh_np = _mw_result_string(_hppd_vch_nonparam_df, "vch_threshold")
HB_bl0_np    = _mw_result_string(_hppd_vch_nonparam_df, "vch_bl_yes_0")
HB_bl75_np   = _mw_result_string(_hppd_vch_nonparam_df, "vch_bl_yes_75")

CB_thresh_np = _spearman_string(_caps_vch_nonparam_df, "vch_threshold",  "caps_vision")
CB_bl0_np    = _spearman_string(_caps_vch_nonparam_df, "vch_bl_yes_0",   "caps_vision")
CB_bl75_np   = _spearman_string(_caps_vch_nonparam_df, "vch_bl_yes_75",  "caps_vision")

_cf_vch_hppd_cov  = "nice_covariates_spusers"
_cf_vch_caps_cov  = "nice_covariates_spusers"
_cf_vch_hppd_dv   = "hppd_binary"

# ── hppd_binary outcomes: Δ and HDI scaled ×100 as percentages ──────────
HB_thresh_reg, HB_thresh_delta_pct, HB_thresh_otherstats, HB_thresh_n, _ = \
    _cf_ppa_pct_parts(_cf_df, "vch_threshold_normalized", _cf_vch_hppd_cov, _cf_vch_hppd_dv)
HB_bl0_reg,    HB_bl0_delta_pct,    HB_bl0_otherstats,    HB_bl0_n,    _ = \
    _cf_ppa_pct_parts(_cf_df, "vch_bl_yes_0_normalized",  _cf_vch_hppd_cov, _cf_vch_hppd_dv)
HB_bl75_reg,   HB_bl75_delta_pct,   _,                    HB_bl75_n,   _ = \
    _cf_ppa_pct_parts(_cf_df, "vch_bl_yes_75_normalized", _cf_vch_hppd_cov, _cf_vch_hppd_dv)

# Concrete VCH threshold SD for the "Lower (-{vchthreshold_sd_str}) thresholds"
# parenthetical. Model contrast: +0.5 Gelman units = +1 raw SD (NOT 2 raw SDs).
_vch_thresh_1rawsd  = df_sp["vch_threshold"].dropna().std()
vchthreshold_sd_str = f"{_vch_thresh_1rawsd:.3f}"

# ── caps_vision outcomes: raw-scale Δ split into delta and otherstats ───
# vch_threshold falls back to full-sample nice_covariates when SP-user model unavailable.
CB_thresh_reg, CB_thresh_n, _, _ = _cf_result_string_or_fallback(
    _cf_df, "vch_threshold_normalized", _cf_vch_caps_cov, "nice_covariates", "caps_vision"
)

# Detection threshold → CAPS: compute split (CB_thresh_reg_delta / CB_thresh_reg_otherstats)
_cb_thresh_mask = (
    (_cf_df['dv'] == 'caps_vision')
    & (_cf_df['spvar'] == 'vch_threshold_normalized')
    & (_cf_df['cov'] == _cf_vch_caps_cov)
)
if not _cb_thresh_mask.any():
    _cb_thresh_mask = (
        (_cf_df['dv'] == 'caps_vision')
        & (_cf_df['spvar'] == 'vch_threshold_normalized')
        & (_cf_df['cov'] == 'nice_covariates')
    )
_cb_thresh_row  = _cf_df[_cb_thresh_mask].iloc[0]
_cb_thresh_est  = float(point_estimate(_cb_thresh_row, source=str(COUNTERFACTUAL_CSV)))
_cb_thresh_lci  = float(_cb_thresh_row['hdi_lower_94'])
_cb_thresh_uci  = float(_cb_thresh_row['hdi_upper_94'])
_cb_thresh_dir  = '>' if _cb_thresh_est >= 0 else '<'
_cb_thresh_prob = float(_cb_thresh_row['prob_above_0']) if _cb_thresh_est >= 0 else float(_cb_thresh_row['prob_below_0'])
CB_thresh_reg_delta      = _fmt_num(abs(_cb_thresh_est))
CB_thresh_reg_otherstats = (
    f"P(Δ{_cb_thresh_dir}0) = {_fmt_prob(_cb_thresh_prob)}, "
    f"94% HDI [{_fmt_num(_cb_thresh_lci)}, {_fmt_num(_cb_thresh_uci)}]"
)

# VCH rate → CAPS: compute split (vchrate_caps_delta / vchrate_caps_otherstats)
_cb_bl0_mask = (
    (_cf_df['dv'] == 'caps_vision')
    & (_cf_df['spvar'] == 'vch_bl_yes_0_normalized')
    & (_cf_df['cov'] == _cf_vch_caps_cov)
)
_cb_bl0_row  = _select_single_row(_cf_df, _cb_bl0_mask, context="CF VCH rate → caps_vision")
_cb_bl0_est  = float(point_estimate(_cb_bl0_row, source=str(COUNTERFACTUAL_CSV)))
_cb_bl0_lci  = float(_cb_bl0_row['hdi_lower_94'])
_cb_bl0_uci  = float(_cb_bl0_row['hdi_upper_94'])
_cb_bl0_dir  = '>' if _cb_bl0_est >= 0 else '<'
_cb_bl0_prob = float(_cb_bl0_row['prob_above_0']) if _cb_bl0_est >= 0 else float(_cb_bl0_row['prob_below_0'])
vchrate_caps_delta     = _fmt_num(_cb_bl0_est)
vchrate_caps_otherstats = (
    f"P(Δ{_cb_bl0_dir}0) = {_fmt_prob(_cb_bl0_prob)}, "
    f"94% HDI [{_fmt_num(_cb_bl0_lci)}, {_fmt_num(_cb_bl0_uci)}]"
)
CB_bl0_reg = f"Δ = {vchrate_caps_delta}, {vchrate_caps_otherstats}"

# Hit rate → CAPS
_cb_bl75_mask = (
    (_cf_df['dv'] == 'caps_vision')
    & (_cf_df['spvar'] == 'vch_bl_yes_75_normalized')
    & (_cf_df['cov'] == _cf_vch_caps_cov)
)
_cb_bl75_row  = _select_single_row(_cf_df, _cb_bl75_mask, context="CF hit rate → caps_vision")
_cb_bl75_est  = float(point_estimate(_cb_bl75_row, source=str(COUNTERFACTUAL_CSV)))
_cb_bl75_lci  = float(_cb_bl75_row['hdi_lower_94'])
_cb_bl75_uci  = float(_cb_bl75_row['hdi_upper_94'])
_cb_bl75_dir  = '>' if _cb_bl75_est >= 0 else '<'
_cb_bl75_prob = float(_cb_bl75_row['prob_above_0']) if _cb_bl75_est >= 0 else float(_cb_bl75_row['prob_below_0'])
CB_bl75_reg_otherstats = (
    f"P(Δ{_cb_bl75_dir}0) = {_fmt_prob(_cb_bl75_prob)}, "
    f"94% HDI [{_fmt_num(_cb_bl75_lci)}, {_fmt_num(_cb_bl75_uci)}]"
)
CB_bl75_reg = f"Δ = {_fmt_num(_cb_bl75_est)}, {CB_bl75_reg_otherstats}"

# Partial Spearman: vch_threshold → caps_vision controlling for age_v2
_df_thresh_partial = (
    df_sp[df_sp["caps_bl_1"].notna()][["vch_threshold", "caps_vision", "age_v2"]].dropna()
)
_n_partial     = len(_df_thresh_partial)
_rho_xy = spearmanr(_df_thresh_partial["vch_threshold"], _df_thresh_partial["caps_vision"]).statistic
_rho_xz = spearmanr(_df_thresh_partial["vch_threshold"], _df_thresh_partial["age_v2"]).statistic
_rho_yz = spearmanr(_df_thresh_partial["caps_vision"],   _df_thresh_partial["age_v2"]).statistic
_denom_partial = np.sqrt((1 - _rho_xz ** 2) * (1 - _rho_yz ** 2))
_rho_partial   = (_rho_xy - _rho_xz * _rho_yz) / _denom_partial
_t_partial     = _rho_partial * np.sqrt(_n_partial - 3) / np.sqrt(1 - _rho_partial ** 2)
_p_partial     = 2 * stats.t.sf(abs(_t_partial), df=_n_partial - 3)
CB_thresh_np_age_control = f"ρ_partial = {_fmt_num(_rho_partial)}, p {_fmt_p(_p_partial)}"

print("HB_thresh_np:", HB_thresh_np, "| HB_bl0_np:", HB_bl0_np, "| HB_bl75_np:", HB_bl75_np)
print("HB_thresh_reg:", HB_thresh_reg, "| HB_bl0_reg:", HB_bl0_reg, "| HB_bl75_reg:", HB_bl75_reg)
print("CB_thresh_np:", CB_thresh_np, "| CB_bl0_np:", CB_bl0_np, "| CB_bl75_np:", CB_bl75_np)
print("CB_thresh_reg:", CB_thresh_reg, "| CB_bl0_reg:", CB_bl0_reg, "| CB_bl75_reg:", CB_bl75_reg)

# Figure & table references — VCH Behavior
vch_behavior_ppa_hx_nonparametric = f"Fig. {FIG_VCH_BEHAVIOR}a"
vch_behavior_ppa_hx_regression    = f"Fig. {FIG_VCH_BEHAVIOR}c"
vch_behavior_caps_nonparametric   = f"Fig. {FIG_VCH_BEHAVIOR}b"
vch_behavior_caps_regression      = f"Fig. {FIG_VCH_BEHAVIOR}d"


# ── VCH Computations ─────────────────────────────────────────────

_hppd_comp_nonparam_csv = (
    f"{_RESULTS}/"
    f"hppd_binary/vch_computations/data_visualization/summary_results/boxplot_grid.csv"
)
_caps_comp_nonparam_csv = (
    f"{_RESULTS}/"
    f"caps_vision/vch_computations/data_visualization/summary_results/correlation_grid.csv"
)
_hppd_comp_nonparam_df = pd.read_csv(_hppd_comp_nonparam_csv)
_caps_comp_nonparam_df = pd.read_csv(_caps_comp_nonparam_csv)

HC_beta_np  = _mw_result_string(_hppd_comp_nonparam_df, "vch_beta")
HC_nu_np    = _mw_result_string(_hppd_comp_nonparam_df, "vch_nu")
HC_omega_np = _mw_result_string(_hppd_comp_nonparam_df, "vch_omega")

CC_nu_np    = _spearman_string(_caps_comp_nonparam_df, "vch_nu",    "caps_vision")
CC_beta_np  = _spearman_string(_caps_comp_nonparam_df, "vch_beta",  "caps_vision")
CC_omega_np = _spearman_string(_caps_comp_nonparam_df, "vch_omega", "caps_vision")

# xprob_by_ppa / xprob_by_caps — Mann-Whitney U on per-participant MEAN xprob
# loaded directly from the primary states CSV.
#
# WHY not vch_xprob_median from df_sp:
#   vch_xprob_median is a median-of-block-medians from an earlier model fit.
#   The trajectory figures (results/*/vch_computations/trajectories/xprob.png)
#   use MEANS (trial→block mean → group mean) from the most recent states CSV.
#   To match the test to the figure, we compute per-participant mean xprob here
#   using the same CSV and the same aggregation method (mean across all trials).
#
# Groups mirror the trajectory figures:
#   hppd_binary:  PPA(-) = 'PPA (-)'  vs  PPA(+) = 'PPA (+)'
#   caps_vision:  CAPS(-) = caps_vision==0  vs  CAPS(+) = caps_vision>0
#                 (restricted to CAPS completers, i.e. caps_vision is not NaN)
# The public repo ships a single trial-level states file, vch_master_public.csv.
# It is the public counterpart of the master repo's dated
# vch_master_withstates_<date>.csv exports, so there is no date-glob or
# most-recent selection to do here — name it directly and fail loudly if absent.
# Verified identical to vch_master_withstates_06-22-2026.csv on the three
# columns read below (216 participants, 77,760 rows, no NaN xprob).
_states_csv = _DATA_DIR / 'vch_master_public.csv'
if not _states_csv.exists():
    raise FileNotFoundError(
        f"Trial-level states file not found: {_states_csv}\n"
        "It supplies record_id/block/xprob for the belief-trajectory comparisons "
        "(xprob_by_ppa, xprob_by_caps)."
    )
_states_df = pd.read_csv(_states_csv, usecols=['record_id', 'block', 'xprob'])
# Per-participant mean xprob across all trials (equivalent to mean of block means
# since all blocks have equal trial counts — 30 trials × 12 blocks = 360 per subject)
_xprob_subject_mean = (
    _states_df.groupby('record_id')['xprob'].mean()
    .reset_index()
    .rename(columns={'xprob': 'xprob_mean'})
)
print(f"[xprob_by_*] States CSV: {_states_csv.name}, "
      f"{_xprob_subject_mean['record_id'].nunique()} subjects with states data")

_xprob_base = df_sp[['record_id', 'hppd_binary', 'caps_vision']].merge(
    _xprob_subject_mean, on='record_id', how='inner'
)

_xprob_ppa_plus  = _xprob_base.loc[_xprob_base['hppd_binary'] == 'PPA (+)', 'xprob_mean']
_xprob_ppa_minus = _xprob_base.loc[_xprob_base['hppd_binary'] == 'PPA (-)', 'xprob_mean']
xprob_by_ppa     = _perm_mean_diff(_xprob_ppa_plus, _xprob_ppa_minus)

_xprob_caps_base  = _xprob_base.dropna(subset=['caps_vision'])
_xprob_caps_plus  = _xprob_caps_base.loc[_xprob_caps_base['caps_vision'] > 0,  'xprob_mean']
_xprob_caps_minus = _xprob_caps_base.loc[_xprob_caps_base['caps_vision'] == 0, 'xprob_mean']
xprob_by_caps     = _perm_mean_diff(_xprob_caps_plus, _xprob_caps_minus)

print("xprob_by_ppa:", xprob_by_ppa)
print("xprob_by_caps:", xprob_by_caps)

_cf_comp_hppd_cov = "nice_covariates_spusers"
_cf_comp_caps_cov = "nice_covariates_spusers"
_cf_comp_hppd_dv  = "hppd_binary"

HC_beta_reg,  HC_beta_n,  _ = _cf_result_string(_cf_df, "vch_beta_normalized",  _cf_comp_hppd_cov, _cf_comp_hppd_dv)
HC_nu_reg,    HC_nu_n,    _ = _cf_result_string(_cf_df, "vch_nu_normalized",    _cf_comp_hppd_cov, _cf_comp_hppd_dv)
HC_omega_reg, HC_omega_n, _ = _cf_result_string(_cf_df, "vch_omega_normalized", _cf_comp_hppd_cov, _cf_comp_hppd_dv)

CC_nu_reg,    CC_nu_n,    _ = _cf_result_string(_cf_df, "vch_nu_normalized",    _cf_comp_caps_cov, "caps_vision")
CC_beta_reg,  CC_beta_n,  _ = _cf_result_string(_cf_df, "vch_beta_normalized",  _cf_comp_caps_cov, "caps_vision")
CC_omega_reg, CC_omega_n, _ = _cf_result_string(_cf_df, "vch_omega_normalized", _cf_comp_caps_cov, "caps_vision")

# CC_beta_reg split (raw CAPS count scale)
_cc_beta_mask = (
    (_cf_df['dv'] == 'caps_vision')
    & (_cf_df['spvar'] == 'vch_beta_normalized')
    & (_cf_df['cov'] == _cf_comp_caps_cov)
)
_cc_beta_row  = _cf_df[_cc_beta_mask].iloc[0]
_cc_beta_est  = float(point_estimate(_cc_beta_row, source=str(COUNTERFACTUAL_CSV)))
_cc_beta_lci  = float(_cc_beta_row['hdi_lower_94'])
_cc_beta_uci  = float(_cc_beta_row['hdi_upper_94'])
_cc_beta_dir  = '>' if _cc_beta_est >= 0 else '<'
_cc_beta_prob = float(_cc_beta_row['prob_above_0']) if _cc_beta_est >= 0 else float(_cc_beta_row['prob_below_0'])
CC_beta_reg_delta = _fmt_num(_cc_beta_est)
CC_beta_reg_stats = (
    f"P(Δ{_cc_beta_dir}0) = {_fmt_prob(_cc_beta_prob)}, "
    f"94% HDI [{_fmt_num(_cc_beta_lci)}, {_fmt_num(_cc_beta_uci)}]"
)
CC_beta_reg = f"Δ = {CC_beta_reg_delta}, {CC_beta_reg_stats}"

# ×100 %-point overrides for HPPD binary regressions + split vars for β
HC_beta_reg, HC_beta_reg_delta, HC_beta_reg_otherstats, HC_beta_n, _ = \
    _cf_ppa_pct_parts(_cf_df, "vch_beta_normalized",  _cf_comp_hppd_cov, _cf_comp_hppd_dv)
HC_nu_reg,   _, _, HC_nu_n,    _ = \
    _cf_ppa_pct_parts(_cf_df, "vch_nu_normalized",    _cf_comp_hppd_cov, _cf_comp_hppd_dv)
HC_omega_reg, _, _, HC_omega_n, _ = \
    _cf_ppa_pct_parts(_cf_df, "vch_omega_normalized", _cf_comp_hppd_cov, _cf_comp_hppd_dv)

# CC_nu_reg split (raw CAPS count scale)
_cc_nu_mask = (
    (_cf_df['dv'] == 'caps_vision')
    & (_cf_df['spvar'] == 'vch_nu_normalized')
    & (_cf_df['cov'] == _cf_comp_caps_cov)
)
_cc_nu_row   = _cf_df[_cc_nu_mask].iloc[0]
_cc_nu_est   = float(point_estimate(_cc_nu_row, source=str(COUNTERFACTUAL_CSV)))
_cc_nu_lci   = float(_cc_nu_row['hdi_lower_94'])
_cc_nu_uci   = float(_cc_nu_row['hdi_upper_94'])
_cc_nu_dir   = '>' if _cc_nu_est >= 0 else '<'
_cc_nu_prob  = float(_cc_nu_row['prob_above_0']) if _cc_nu_est >= 0 else float(_cc_nu_row['prob_below_0'])
CC_nu_reg_delta      = _fmt_num(_cc_nu_est)
CC_nu_reg_otherstats = (
    f"P(Δ{_cc_nu_dir}0) = {_fmt_prob(_cc_nu_prob)}, "
    f"94% HDI [{_fmt_num(_cc_nu_lci)}, {_fmt_num(_cc_nu_uci)}]"
)

# IQR-outlier sensitivity, from the same CF CSV. Denominator is the base
# model row already selected above (_cc_nu_row), so both N describe the same
# dv x spvar and differ only by the IQR exclusion.
nu_caps_iqr_results, _n_nu_iqr, _ = _cf_result_string(
    _cf_df, "vch_nu_normalized", "nice_covariates_spusers_iqr", "caps_vision")
_n_nu_base = int(_cc_nu_row['N_obs'])
nu_caps_iqr_drop = f"{100 * (_n_nu_base - _n_nu_iqr) / _n_nu_base:.1f}"

# IQR + beta-control sensitivity: nu -> caps_vision, with vch_beta as additional covariate.
# cov name is nice_covariates_beta_spusers_iqr (beta before spusers).
# Read from COUNTERFACTUAL_CSV. The type is in BASE_MODELS in
# 03_hpc/generate_hpc_jobs.py, so the array produces it natively; if this raises
# "NO CF RESULT ... cov='nice_covariates_beta_spusers_iqr'", the model type was
# dropped from BASE_MODELS or the array has not been re-pulled.
nu_caps_iqr_beta_results, _, _ = _cf_result_string(_cf_df, "vch_nu_normalized", "nice_covariates_beta_spusers_iqr", "caps_vision")

print("HC_beta_np:", HC_beta_np, "| HC_nu_np:", HC_nu_np, "| HC_omega_np:", HC_omega_np)
print("HC_beta_reg:", HC_beta_reg, "| HC_nu_reg:", HC_nu_reg, "| HC_omega_reg:", HC_omega_reg)
print("CC_nu_np:", CC_nu_np, "| CC_nu_reg:", CC_nu_reg)
print("CC_beta_np:", CC_beta_np, "| CC_beta_reg:", CC_beta_reg)
print("CC_omega_np:", CC_omega_np, "| CC_omega_reg:", CC_omega_reg)

# Figure & table references — VCH Computations
vch_computations_ppa_hx_nonparametric_ppa_beta = f"Fig. {FIG_VCH_COMPUTATIONS}j"
vch_computations_ppa_hx_nonparametric_ppa_nu  = f"Fig. {FIG_VCH_COMPUTATIONS}h"
vch_computations_ppa_hx_nonparametric_ppa_omega = f"Fig. {FIG_VCH_COMPUTATIONS}b"
vch_computations_ppa_hx_regression    = f"Fig. {FIG_MEDIATION}l" #A
vch_computations_caps_regression      = f"Fig. {FIG_MEDIATION}m" #B
vch_computations_caps_nonparametric_nu = f"Fig. {FIG_VCH_COMPUTATIONS}i"
vch_computations_caps_nonparametric_beta = f"Fig. {FIG_VCH_COMPUTATIONS}k"
vch_computations_trajectory_ppa_hx = f"Fig. {FIG_VCH_COMPUTATIONS}d & f"
vch_computations_trajectory_caps   = f"Fig. {FIG_VCH_COMPUTATIONS}e & g"
# ω → CAPS spans two panels of the same figure: the nonparametric panel (c)
# and the regression panel (m).
vch_computations_caps_omega_panels = f"Fig. {FIG_VCH_COMPUTATIONS}c & m"


# ── Beta SDT Correlates ──────────────────────────────────────────────────────
# _corr: Spearman correlations between vch_beta and iv_type_dict["sdt_hppd"]
#   variables (criterion_overall, d_prime_overall, mean_conf_fas), read from
#   the vch_beta × sdt_hppd correlation grid CSV.
# _reg:  counterfactual regression results (SDT variable → hppd_binary /
#   caps_vision) from existingresults_manuscript_counterfactual.csv, using
#   nice_covariates_spusers.  HPPD binary uses _cf_ppa_pct_parts (×100 scale);
#   CAPS vision uses _cf_result_string (response-count scale).
# Variable naming rule: beta_{sdt_var}_corr / {outcome}_beta_{sdt_var}_reg,
#   where sdt_var is the row_var in the correlation grid.

_beta_sdt_corr_csv = (
    f"{_RESULTS}/"
    "vch_beta/sdt_hppd/data_visualizations/summary_results/correlation_grid.csv"
)
_beta_sdt_corr_df = pd.read_csv(_beta_sdt_corr_csv)

beta_d_prime_overall_corr   = _spearman_string(_beta_sdt_corr_df, "d_prime_overall",   "vch_beta")
beta_criterion_overall_corr = _spearman_string(_beta_sdt_corr_df, "criterion_overall", "vch_beta")
beta_mean_conf_fas_corr     = _spearman_string(_beta_sdt_corr_df, "mean_conf_fas",     "vch_beta")

# Partial Spearman: vch_threshold vs. vch_hit_rate / vch_beta / d_prime_overall,
# controlling for age_v2. SP users only; _partial_spearman_str handles pairwise dropna.
vch_hit_rate_vch_threshold_corr_age_v2_control    = _partial_spearman_str(df_sp, "vch_threshold", "vch_hit_rate")
vch_beta_vch_threshold_corr_age_v2_control        = _partial_spearman_str(df_sp, "vch_threshold", "vch_beta")
d_prime_overall_vch_threshold_corr_age_v2_control = _partial_spearman_str(df_sp, "vch_threshold", "d_prime_overall")

_sdt_hppd_cov = "nice_covariates_spusers"
_sdt_caps_cov = "nice_covariates_spusers"

# HPPD binary — probability scale (×100); full_pct_str from _cf_ppa_pct_parts
ppa_beta_d_prime_overall_reg, _, _, _, _   = _cf_ppa_pct_parts(
    _cf_df, "d_prime_overall_normalized",   _sdt_hppd_cov, "hppd_binary"
)
ppa_beta_criterion_overall_reg, _, _, _, _ = _cf_ppa_pct_parts(
    _cf_df, "criterion_overall_normalized", _sdt_hppd_cov, "hppd_binary"
)
ppa_beta_mean_conf_fas_reg, _, _, _, _     = _cf_ppa_pct_parts(
    _cf_df, "mean_conf_fas_normalized",     _sdt_hppd_cov, "hppd_binary"
)

# CAPS vision — response-count scale
caps_beta_d_prime_overall_reg, _, _   = _cf_result_string(
    _cf_df, "d_prime_overall_normalized",   _sdt_caps_cov, "caps_vision"
)
caps_beta_criterion_overall_reg, _, _ = _cf_result_string(
    _cf_df, "criterion_overall_normalized", _sdt_caps_cov, "caps_vision"
)
caps_beta_mean_conf_fas_reg, _, _     = _cf_result_string(
    _cf_df, "mean_conf_fas_normalized",     _sdt_caps_cov, "caps_vision"
)

print("beta_d_prime_overall_corr:", beta_d_prime_overall_corr)
print("beta_criterion_overall_corr:", beta_criterion_overall_corr)
print("beta_mean_conf_fas_corr:", beta_mean_conf_fas_corr)
print("vch_hit_rate_vch_threshold_corr_age_v2_control:", vch_hit_rate_vch_threshold_corr_age_v2_control)
print("vch_beta_vch_threshold_corr_age_v2_control:", vch_beta_vch_threshold_corr_age_v2_control)
print("d_prime_overall_vch_threshold_corr_age_v2_control:", d_prime_overall_vch_threshold_corr_age_v2_control)
print("ppa_beta_d_prime_overall_reg:", ppa_beta_d_prime_overall_reg)
print("ppa_beta_criterion_overall_reg:", ppa_beta_criterion_overall_reg)
print("ppa_beta_mean_conf_fas_reg:", ppa_beta_mean_conf_fas_reg)
print("caps_beta_d_prime_overall_reg:", caps_beta_d_prime_overall_reg)
print("caps_beta_criterion_overall_reg:", caps_beta_criterion_overall_reg)
print("caps_beta_mean_conf_fas_reg:", caps_beta_mean_conf_fas_reg)


# ── Mediation Models ─────────────────────────────────────────────

_med_base = _RESULTS

_med_cov = "nice_covariates_spusers"

def _load_cf_med(subdir):
    """Load (cf_path_df, mc_df) for a mediation model directory."""
    d = Path(f"{_med_base}/{subdir}")
    return (
        pd.read_csv(d / "path_counterfactual_summary.csv"),
        pd.read_csv(d / "mc_mediation_summary.csv"),
    )

_hppd_thresh_cf,  _hppd_thresh_mc  = _load_cf_med(f"hppd_binary/mediation_models/hppd_binary_spage_vchthreshold_{_med_cov}")
_hppd_vchrate_cf, _hppd_vchrate_mc = _load_cf_med(f"hppd_binary/mediation_models/hppd_binary_spage_vchrate_{_med_cov}")
_hppd_beta_cf,    _hppd_beta_mc    = _load_cf_med(f"hppd_binary/mediation_models/hppd_binary_spage_vchbeta_{_med_cov}")
_hppd_vchnu_cf,   _hppd_vchnu_mc   = _load_cf_med(f"hppd_binary/mediation_models/hppd_binary_spage_vchnu_{_med_cov}")
_caps_vchrate_cf, _caps_vchrate_mc = _load_cf_med(f"caps_vision/mediation_models/caps_vision_avgdose_vchrate_{_med_cov}")
_caps_vchnu_cf,   _caps_vchnu_mc   = _load_cf_med(f"caps_vision/mediation_models/caps_vision_avgdose_vchnu_{_med_cov}")
_caps_thresh_cf,  _caps_thresh_mc  = _load_cf_med(f"caps_vision/mediation_models/caps_vision_avgdose_vchthreshold_{_med_cov}")

# Indirect effects (NIE) — from MC integration summary
MHA = _cf_nie_string(_hppd_thresh_mc)   # psychedelic_age → vch_threshold → hppd_binary
MHV = _cf_nie_string(_hppd_vchrate_mc)  # psychedelic_age → vch_bl_yes_0  → hppd_binary
MHB = _cf_nie_string(_hppd_beta_mc)     # psychedelic_age → vch_beta       → hppd_binary
M_nuhppd = _cf_nie_string(_hppd_vchnu_mc)  # psychedelic_age → vch_nu    → hppd_binary
MDV = _cf_nie_string(_caps_vchrate_mc)  # avg_life_dose   → vch_bl_yes_0  → caps_vision
MDN = _cf_nie_string(_caps_vchnu_mc)    # avg_life_dose   → vch_nu         → caps_vision
MDT = _cf_nie_string(_caps_thresh_mc)   # avg_life_dose   → vch_threshold  → caps_vision
M_nucaps = MDN                          # alias: prior weighting mediation of dose → caps_vision

# ── Split NIE strings (delta / other-stats), for inline narrative use ───
# Every one of these is the posterior MEDIAN of the MC-integrated NIE, so every
# one is labelled "Δmed" — the same label _cf_nie_string() applies. Nothing in
# the narrative reports an indirect effect as a bare "Δ".
# MHA_delta / MHA_otherstats: HPPD binary NIE via threshold (×100 for prob scale)
_nie_ht          = _hppd_thresh_mc[_hppd_thresh_mc['effect'].str.startswith('NIE')].iloc[0]
_nie_ht_est      = float(point_estimate(_nie_ht, source='mc_mediation_summary.csv', mc_integrated=True))
_nie_ht_lci      = float(_nie_ht['hdi_low'])
_nie_ht_uci      = float(_nie_ht['hdi_high'])
_nie_ht_dir      = '>' if _nie_ht_est >= 0 else '<'
_nie_ht_prob     = float(_nie_ht['p_above_0']) if _nie_ht_est >= 0 else float(_nie_ht['p_below_0'])
MHA_delta      = f"{100 * _nie_ht_est:.1f}%"
MHA_otherstats = (
    f"P(Δ{_nie_ht_dir}0) = {_fmt_prob(_nie_ht_prob)}, "
    f"94% HDI [{100 * _nie_ht_lci:.1f}%, {100 * _nie_ht_uci:.1f}%]"
)
MHA_allstats   = f"Δmed = {MHA_delta}, {MHA_otherstats}"

# MDV_delta / MDV_otherstats: CAPS NIE via vchrate (raw count scale)
_nie_cv          = _caps_vchrate_mc[_caps_vchrate_mc['effect'].str.startswith('NIE')].iloc[0]
_nie_cv_est      = float(point_estimate(_nie_cv, source='mc_mediation_summary.csv', mc_integrated=True))
_nie_cv_lci      = float(_nie_cv['hdi_low'])
_nie_cv_uci      = float(_nie_cv['hdi_high'])
_nie_cv_dir      = '>' if _nie_cv_est >= 0 else '<'
_nie_cv_prob     = float(_nie_cv['p_above_0']) if _nie_cv_est >= 0 else float(_nie_cv['p_below_0'])
MDV_delta      = _fmt_num(_nie_cv_est)
MDV_otherstats = (
    f"P(Δ{_nie_cv_dir}0) = {_fmt_prob(_nie_cv_prob)}, "
    f"94% HDI [{_fmt_num(_nie_cv_lci)}, {_fmt_num(_nie_cv_uci)}]"
)
MDV_allstats   = f"Δmed = {MDV_delta}, {MDV_otherstats}"

# MDT_delta / MDT_otherstats: CAPS NIE via threshold (raw count scale)
_nie_ct          = _caps_thresh_mc[_caps_thresh_mc['effect'].str.startswith('NIE')].iloc[0]
_nie_ct_est      = float(point_estimate(_nie_ct, source='mc_mediation_summary.csv', mc_integrated=True))
_nie_ct_lci      = float(_nie_ct['hdi_low'])
_nie_ct_uci      = float(_nie_ct['hdi_high'])
_nie_ct_dir      = '>' if _nie_ct_est >= 0 else '<'
_nie_ct_prob     = float(_nie_ct['p_above_0']) if _nie_ct_est >= 0 else float(_nie_ct['p_below_0'])
MDT_delta      = _fmt_num(_nie_ct_est)
MDT_otherstats = (
    f"P(Δ{_nie_ct_dir}0) = {_fmt_prob(_nie_ct_prob)}, "
    f"94% HDI [{_fmt_num(_nie_ct_lci)}, {_fmt_num(_nie_ct_uci)}]"
)
MDT_allstats   = f"Δmed = {MDT_delta}, {MDT_otherstats}"

# ×100 %-point overrides for HPPD binary NIE strings (probability scale)
MHA = (
    f"Δmed = {100 * _nie_ht_est:.1f}%, P(Δ{_nie_ht_dir}0) = {_fmt_prob(_nie_ht_prob)}, "
    f"94% HDI [{100 * _nie_ht_lci:.1f}%, {100 * _nie_ht_uci:.1f}%]"
)
_nie_hv      = _hppd_vchrate_mc[_hppd_vchrate_mc['effect'].str.startswith('NIE')].iloc[0]
_nie_hv_est  = float(point_estimate(_nie_hv, source='mc_mediation_summary.csv', mc_integrated=True))
_nie_hv_lci  = float(_nie_hv['hdi_low'])
_nie_hv_uci  = float(_nie_hv['hdi_high'])
_nie_hv_dir  = '>' if _nie_hv_est >= 0 else '<'
_nie_hv_prob = float(_nie_hv['p_above_0']) if _nie_hv_est >= 0 else float(_nie_hv['p_below_0'])
MHV = (
    f"Δmed = {100 * _nie_hv_est:.1f}%, P(Δ{_nie_hv_dir}0) = {_fmt_prob(_nie_hv_prob)}, "
    f"94% HDI [{100 * _nie_hv_lci:.1f}%, {100 * _nie_hv_uci:.1f}%]"
)
MHV_delta      = f"{abs(100 * _nie_hv_est):.1f}%"
MHV_otherstats = (
    f"P(Δ{_nie_hv_dir}0) = {_fmt_prob(_nie_hv_prob)}, "
    f"94% HDI [{100 * _nie_hv_lci:.1f}%, {100 * _nie_hv_uci:.1f}%]"
)
MHV_allstats   = f"Δmed = {MHV_delta}, {MHV_otherstats}"
_nie_hb      = _hppd_beta_mc[_hppd_beta_mc['effect'].str.startswith('NIE')].iloc[0]
_nie_hb_est  = float(point_estimate(_nie_hb, source='mc_mediation_summary.csv', mc_integrated=True))
_nie_hb_lci  = float(_nie_hb['hdi_low'])
_nie_hb_uci  = float(_nie_hb['hdi_high'])
_nie_hb_dir  = '>' if _nie_hb_est >= 0 else '<'
_nie_hb_prob = float(_nie_hb['p_above_0']) if _nie_hb_est >= 0 else float(_nie_hb['p_below_0'])
MHB = (
    f"Δmed = {100 * _nie_hb_est:.1f}%, P(Δ{_nie_hb_dir}0) = {_fmt_prob(_nie_hb_prob)}, "
    f"94% HDI [{100 * _nie_hb_lci:.1f}%, {100 * _nie_hb_uci:.1f}%]"
)
_nie_hn      = _hppd_vchnu_mc[_hppd_vchnu_mc['effect'].str.startswith('NIE')].iloc[0]
_nie_hn_est  = float(point_estimate(_nie_hn, source='mc_mediation_summary.csv', mc_integrated=True))
_nie_hn_lci  = float(_nie_hn['hdi_low'])
_nie_hn_uci  = float(_nie_hn['hdi_high'])
_nie_hn_dir  = '>' if _nie_hn_est >= 0 else '<'
_nie_hn_prob = float(_nie_hn['p_above_0']) if _nie_hn_est >= 0 else float(_nie_hn['p_below_0'])
M_nuhppd = (
    f"Δmed = {100 * _nie_hn_est:.1f}%, P(Δ{_nie_hn_dir}0) = {_fmt_prob(_nie_hn_prob)}, "
    f"94% HDI [{100 * _nie_hn_lci:.1f}%, {100 * _nie_hn_uci:.1f}%]"
)

# A-paths (predictor → mediator) — from path_counterfactual_summary
AT = _cf_path_string(_hppd_thresh_cf,  'A path')   # psychedelic_age → vch_threshold
AB = _cf_path_string(_hppd_beta_cf,    'A path')   # psychedelic_age → vch_beta
DN = _cf_path_string(_caps_vchnu_cf,   'A path')   # avg_life_dose   → vch_nu
# MHB_apath: A-path narrative string for psychedelic_age → vch_beta in the
# hppd_binary_spage_vchbeta_nice_covariates_spusers model.  Used in
# vch_computations_results to justify why β does not mediate PPA risk from
# younger age of first use ("age of first use was not associated with β").
MHB_apath = AB

_caps_beta_cf, _caps_beta_mc = _load_cf_med(f"caps_vision/mediation_models/caps_vision_avgdose_vchbeta_{_med_cov}")

# M_caps_delta / M_caps_otherstats: CAPS NIE via vch_beta (raw count scale)
_nie_cb          = _caps_beta_mc[_caps_beta_mc['effect'].str.startswith('NIE')].iloc[0]
_nie_cb_est      = float(point_estimate(_nie_cb, source='mc_mediation_summary.csv', mc_integrated=True))
_nie_cb_lci      = float(_nie_cb['hdi_low'])
_nie_cb_uci      = float(_nie_cb['hdi_high'])
_nie_cb_dir      = '>' if _nie_cb_est >= 0 else '<'
_nie_cb_prob     = float(_nie_cb['p_above_0']) if _nie_cb_est >= 0 else float(_nie_cb['p_below_0'])
M_caps_delta      = _fmt_num(_nie_cb_est)
M_caps_otherstats = (
    f"P(Δ{_nie_cb_dir}0) = {_fmt_prob(_nie_cb_prob)}, "
    f"94% HDI [{_fmt_num(_nie_cb_lci)}, {_fmt_num(_nie_cb_uci)}]"
)
M_caps_allstats   = f"Δmed = {M_caps_delta}, {M_caps_otherstats}"

print("MHA:", MHA, "| MHV:", MHV, "| MHB:", MHB)
print("MDV:", MDV, "| MDN:", MDN, "| MDT:", MDT)
print("AT:", AT, "| DN:", DN)

# Figure & table references — Mediation Models
# Column 1 (A, C, E): predictor = psychedelic_age, outcome = PPA History
# Column 2 (B, D, F): predictor = avg_life_dose,   outcome = CAPS Vision
# Row 1 (A/B): mediator = vch_threshold | Row 2 (C/D): mediator = vch_bl_yes_0
# Row 3 (E):   mediator = vch_beta      | Row 3 (F):   mediator = vch_nu
mediation_ppa_hx_threshold_panel = f"Fig. {FIG_VCH_BEHAVIOR}e"
mediation_ppa_hx_vchrate_panel   = f"Fig. {FIG_VCH_BEHAVIOR}g"
mediation_ppa_hx_beta_panel      = f"Fig. {FIG_MEDIATION}n"
mediation_dose_beta_panel        = f"Fig. {FIG_MEDIATION}o"
mediation_caps_threshold_panel   = f"Fig. {FIG_VCH_BEHAVIOR}f"
mediation_caps_vchrate_panel     = f"Fig. {FIG_VCH_BEHAVIOR}h"

psychometric_curve_panel = f"Fig. {FIG_BETA}a"
beta_sdt_correlates_panel = f"Fig. {FIG_BETA}b"
sdt_ppa_hx_panel = f"Fig. {FIG_BETA}c"
sdt_caps_panel = f"Fig. {FIG_BETA}d"

# ================================================================
# SECTION 5: NARRATIVE TEXT BLOCKS
# ================================================================
#
# Each block below is one manuscript section, reproduced verbatim apart from the
# f-string placeholders.  Two conventions are load-bearing and must be preserved:
#
#   * The literal begins and ends with a newline.  write_results_txt() emits the
#     blocks without stripping, so those newlines are what separate the heading
#     from the first paragraph and one section from the next.  beta_sdt_results
#     opens with a blank line because the manuscript puts one there.
#   * Trailing spaces at the end of some paragraphs are inherited from the
#     manuscript and are intentional.  Do not let an editor strip them, or the
#     output will stop matching the submitted text.

recruitment_results_text = f"""
Of the {RS} participants who completed screening, {RE} were found eligible (see {supp_fig_consort}), and {RE_COMPLETERS} completed all study measures. An additional {RE_SALVAGED} participants completed at least the full first questionnaire and passed all QC checks and so were retained for final analysis pertaining to the completed data. Final analysis included {RG} participants ({RSP} SP users; {TABLE_1}) passing all QC checks. {VCH_QC_EXCLUDED_WORD} participants’ VCH data were excluded from analysis for failure of QC measures described above.
"""

clinical_demographic_results = f"""
The sample skewed young, male, white, United-States (US)-residing, and well-educated, with more than {FOUR_YR_OR_MORE_PCT}% completing a four-year degree or greater ({TABLE_1}). There were virtually no differences in demographic profile between participants endorsing a history of SP-associated PPAs or current, CAPS-assessed PPAs. Those with a SP-associated PPA history were less likely to be multiracial, and those with current PPAs were younger (consistent with other reports),{CITE['izmi_2024']} more likely to be European, and less likely to be US-residing. Psychiatric illness rates were similar to the general population.{CITE['kessler_2005']} As expected, polydrug use was considerable; all participants had used another psychoactive drug and the majority had used one in the past month. {'Most' if nonsp_psychedelic_lifetime_pct > 50 else f'{nonsp_psychedelic_lifetime_pct}%'} participants reported prior non-serotonergic psychedelic use, though fewer than {int(math.ceil(nonsp_psychedelic_past_month_pct / 10) * 10)}% had used in the past month. 

SP use patterns are summarized in {SP_FIG}. Total lifetime uses varied significantly ({psycheduse_life_nomic_mean} ± {psycheduse_life_nomic_sd}) and were positively skewed—from {sp_lifetime_uses_min} to {sp_lifetime_uses_max}. The weighted average subjective SP dose used was more consistent ({avg_life_dose_mean} μg ± {avg_life_dose_sd} μg), with the vast majority of uses concentrated at the highest end of doses administered in contemporary in-laboratory research (the equivalent of ~200μg of LSD).{CITE['hirschfeld_2023']} More than a quarter of participants reported doses substantially higher than this cutoff. First SP use generally occurred in young adulthood ({psychedelic_age_mean} ± {psychedelic_age_sd} years), with {PERCENT_BELOW_25_EXP}% of participants reporting initiation prior to age 25 and only {PERCENT_ABOVE_40}% at age 40 or older. The plurality of participants reported using SPs for both therapeutic and recreational reasons, whereas the remainder were split relatively evenly between primarily recreational and primarily therapeutic use (see {SP_FIG_SHORT}). 
"""

ppa_history_results = f"""
Consistent with prior findings,{CITE['carhart_nutt_2010']},{CITE['baggott_2011']},{CITE['kvam_2023']} SP-associated PPAs were common, endorsed by {persist_vis_yn_num_PCT} of the full sample (Fig. {PPA_FIG}). Also consistent,{CITE['izmi_2024']},{CITE['zhou_2025']},{CITE['muller_2022']} PPAs typically occurred after few SP uses ({persistvis_psychdoses_median_iqr}). {persistvis_psychdoses_one_dose_sentence} reported effects after just one dose. Participants generally endorsed multiple PPA phenomena; the modal symptom count was {baggot_total_mode_num_PCT}, and {baggot_total_ABOVE_4_num_PCT} endorsed 4 or more symptoms. The most endorsed PPA was "{persistvis_most_mode}", reported by {persistvis_most_mode_PCT}. Also consistent with other studies, most ({hppd_true_chronicity_BELOW_5_PCT}) reported that symptoms did not persist beyond 1 week, whereas {hppd_true_chronicity_ABOVE_7_PCT} reported durations longer than 1 year. Finally, as reported before,{CITE['carhart_nutt_2010']}–{CITE['muller_2022']} the vast majority of participants reported not experiencing distress ({persistvis_distress_PCT}) nor seeking treatment ({persistvis_txseek_PCT}) and only {hppd_diagnosis_num_PCT} reported a formal HPPD diagnosis.

The {caps_pos_majority} ({caps_pos_pct}%) of SP users with CAPS data endorsed one or more current PPAs—primarily those endorsing prior SP-associated PPAs ({caps_ppa_pos_pct}%). Consistent with SP-associated PPAs, the most endorsed CAPS item was "lights seeming brighter or colors seeming more intense", and this effect was typically rated as non-distressing and non-distracting, but also one of the most frequent symptoms (Fig. {CAPS_FIG}). Frequency was episodic, with only {caps_freq_fives}% endorsing any symptoms "all the time", and most symptoms experienced "not often" to "sometimes". Distress was rarer, with {pct_caps_vision_distress_0}% of participants reporting one or more PPAs endorsing all of their experiences as "not at all" distressing, while only {pct_caps_vision_distress_over3}% endorsed an experience as "firmly" or "very" distressing.
"""

ppa_hx_sp_results = f"""
We first examined which patterns of SP use were associated with PPA history among SP users. Compared with those who denied a history of PPAs after SP use, participants who endorsed PPAs (PPA+) had a younger age at first SP use ({sp_predictors_ppa_hx_nonparametric}; {PA}; see {mann_whitney_u_table} for additional statistics from all binary Mann-Whitney U tests). By contrast, there was less certain evidence that PPA+ individuals have had greater lifetime SP uses ({PB_lifetime}) and no credible sign of greater average doses used ({PB_dose}). We further interrogated these associations with Bayesian regressions that incorporate relevant covariates (age, sex, IQ, and mental illness history) and magnitude of both variables in order to estimate how much PPA probability cross-sectionally varies with SP use patterns. We found that younger (- {psychedelic_age_sd} years) age of first SP use was associated with a {PC_abs_pct} higher estimated probability of lifetime PPAs ({PC}). Greater lifetime SP uses ({psycheduse_life_nomic_sd} uses) were associated with greater estimated PPA probability ({PD}) but this association diminished considerably when excluding outliers ({PE}; {PEN}% of participants dropped; {sensitivity_analysis_heatmap_single_path.replace("Figure", "Fig.")}). Higher average SP dose used (+{avg_life_dose_sd} LSD μg equivalents) was not credibly associated with PPA probability ({PD_dose}); though this changed after outlier exclusion ({PF}; {PFN}% of participants dropped). In summary, only younger age of first SP use was robustly associated with small but meaningfully greater probability of prior SP-associated PPAs, and so subsequent mediation analyses examining past PPA risk focus on this predictor.
"""

caps_sp_results = f"""
We next asked how these SP patterns varied with current (past-month) visual PPAs ({sp_predictors_caps_nonparametric}). We controlled for age because age correlated with all SP use variables and strongly anticorrelated with CAPS endorsements ({TABLE_1}). Average SP dose ({CB_dose}) strongly correlated with greater CAPS vision items endorsed, whereas age of first SP use ({CA}) and lifetime SP uses ({CB_lifetime}) had 50% weaker correlations. Regression models confirmed that a higher average dose was associated with a greater estimated number of CAPS vision items endorsed ({CD}). Age of first SP use showed no association with CAPS vision ({CC}). Lifetime SP uses were associated with higher estimates of CAPS vision items ({CSE_life}), but this association was entirely outlier-driven ({CSE_life_iqr}). Subsequent mediation analyses of current PPA symptoms focus on higher average SP doses used, as this showed the strongest cross-sectional association with current PPA symptoms.
"""

vch_behavior_results = f"""
We next investigated potential behavioral correlates of SP-associated PPA history and current PPAs. Those with a history of SP-associated PPAs had lower detection thresholds (i.e., reported detection at lower visual contrasts; {HB_thresh_np}; {vch_behavior_ppa_hx_nonparametric}), associated with {HB_thresh_delta_pct} greater PPA risk ({HB_thresh_otherstats}; {vch_behavior_ppa_hx_regression}). PPA (+) participants also had greater VCH rates ({HB_bl0_np}), associated with an estimated {HB_bl0_delta_pct} greater PPA risk ({HB_bl0_otherstats}). Veridical detection in the easiest-to-detect 75% condition (hit-rate) was not related to PPA history ({HB_bl75_np}; {HB_bl75_reg}).

CAPS items showed a weaker correlation with threshold ({CB_thresh_np}) that strengthened after controlling for age ({CB_thresh_np_age_control}; {vch_behavior_caps_nonparametric}). Lower (-{vchthreshold_sd_str}) thresholds were associated with a {CB_thresh_reg_delta} more CAPS vision endorsements ({CB_thresh_reg_otherstats}; {vch_behavior_caps_regression}). Higher VCH rates ({CB_bl0_np}) correlated with current CAPS visual items endorsed and was associated with {vchrate_caps_delta} more CAPS visual items ({vchrate_caps_otherstats}). Hit rate was not correlated with CAPS items ({CB_bl75_np}), though it exhibited uncertain associations with estimated CAPS endorsements ({CB_bl75_reg}). 

We then asked if greater visual sensitivity (lower thresholds) and higher VCH rates statistically mediated the SP use patterns previously associated with past and current PPAs—i.e. could covariance between SP use and behavior and behavior and PPAs alone plausibly explain the association between SP use and PPAs via this indirect, modeled pathway. We found plausible evidence of such an indirect path between younger age of first use and higher PPA history probability through lower visual threshold ({MHA_allstats}; {mediation_ppa_hx_threshold_panel}), but not VCH rate ({MHV_allstats}; {mediation_ppa_hx_vchrate_panel}). By contrast, SP dose and current PPAs’ association was not convincingly explainable through lower threshold ({MDT_allstats}; {mediation_caps_threshold_panel}), but it was through VCH rate ({MDV_allstats}; {mediation_caps_vchrate_panel}).
"""

vch_computations_results = f"""
Finally, to determine which latent states drove these behavioral patterns, we fit HGF models to VCH behavioral data. We then examined whether HGF-derived parameter estimates varied with both PPA measures.

Those with a history of SP-associated PPAs exhibited no clear difference in estimated prior weighting, ν ({HC_nu_np}; {HC_nu_reg}; {vch_computations_ppa_hx_nonparametric_ppa_nu}). By contrast, PPA (+) participants did have more stochasticity in responding, reflected in lower β values ({HC_beta_np}; {vch_computations_ppa_hx_nonparametric_ppa_beta}). Lower β was associated with a {HC_beta_reg_delta} greater probability of prior SP-associated PPAs ({HC_beta_reg_otherstats}; {vch_computations_ppa_hx_regression}). There were no associations between PPA risk and ω ({HC_omega_np}; {HC_omega_reg}; {vch_computations_ppa_hx_nonparametric_ppa_omega}) nor belief trajectories ({xprob_by_ppa}; {vch_computations_trajectory_ppa_hx}).

In the case of current PPAs, β again negatively correlated with endorsements ({CC_beta_np}; {vch_computations_caps_nonparametric_beta}), and was associated with a {CC_beta_reg_delta} fewer current PPAs ({CC_beta_reg_stats}; {vch_computations_caps_regression}). We also saw a very weak association with prior weighting, ν ({CC_nu_np}). Prior weighting was associated with {CC_nu_reg_delta} more visual PPAs with modest confidence ({CC_nu_reg_otherstats}; {vch_computations_caps_nonparametric_nu}), but this association significantly diminished with outlier exclusion ({nu_caps_iqr_results}; {nu_caps_iqr_drop}% of participants dropped) and remained uncertain even after controlling for β ({nu_caps_iqr_beta_results}; {sensitivity_analysis_heatmap_single_path.replace("Figure", "Fig.")}). Current PPAs were also not associated with ω ({CC_omega_np}; {CC_omega_reg}; {vch_computations_caps_omega_panels}) nor belief trajectories ({xprob_by_caps}; {vch_computations_trajectory_caps}).

β did not statistically mediate the association between PPA risk and younger age of first SP use ({MHB}; {mediation_ppa_hx_beta_panel}), apparently because age of first use was not associated with β ({MHB_apath}). However, the relationship between average SP dose and greater CAPS visual item endorsements was plausibly explainable through dose’s association with low β ({M_caps_allstats}; {mediation_dose_beta_panel}). Unsurprisingly, ν did not statistically mediate associations between age of first use and PPA history ({M_nuhppd}) nor average dose and current PPAs ({M_nucaps}; {supp_table_mediation_sensitivity}). 
"""

beta_sdt_results = f"""

From the results above, low decision precision (β) appeared to be the strongest explanatory parameter. However, its interpretation is ambiguous because β was not designed after any particular theoretical process—any behavior not well explained using perceptual model parameters could influence it. Hence, low decision precision could reflect poor task engagement, impaired sensory processing, and/or altered metacognitive processing. To help determine which of these interpretations is best supported by the data, we conducted a series of exploratory analyses. 

We first considered whether lower β scores indicated noisy responses due to poor task engagement. We examined correlations between β and six separate measures of task engagement, including repeated responses, response randomness, and response time-based measures (see Methods and {task_engagement_fig_threshold_error.replace("Figure", "Fig.")}) and found no evidence that poor task engagement related to low decision precision. 

We next examined participants’ overall task performance, splitting participants by the median values of the two PPA-linked HGF parameters, β and ν ({psychometric_curve_panel}). While both were associated with higher VCH rates, lower hit rates in the target-present portion of the task drove lower β estimates, and those with the greatest VCH rates also had lower hit rates, resulting in flatter psychometric curves, classically corresponding to low sensitivity.{CITE['sdt_low_sensitivity']} Consistent with this, we found direct correlations ({beta_sdt_correlates_panel}) between β and both dʹ ({beta_d_prime_overall_corr}) and decision criterion ({beta_criterion_overall_corr})—suggesting poor sensory fidelity alongside a low threshold for target perception. Given associations between lower threshold and SP-related PPAs above, we also considered whether diminished hit rates associated with low β were due simply to the task being systematically more challenging for those with lower threshold. After controlling for age, threshold exhibited a very uncertain correlation with hit rate ({vch_hit_rate_vch_threshold_corr_age_v2_control}) and no correlation with β ({vch_beta_vch_threshold_corr_age_v2_control}). Threshold more strongly correlated with dʹ ({d_prime_overall_vch_threshold_corr_age_v2_control}). These results suggest that threshold alone does not explain beta’s association with poor sensory sensitivity. 

Finally, we considered whether there was evidence of metacognitive abnormalities by inspecting confidence ratings: because β is outside the HGF’s perceptual model, low decision precision could reflect lower confidence in the posterior derived from the perceptual model. We found the opposite: lower β correlated with greater confidence in VCHs ({beta_mean_conf_fas_corr}; {beta_sdt_correlates_panel}).

Having found evidence that lower β is associated with poor sensory discrimination coincident with a liberal criterion for signal detection and reduced metacognitive accuracy, we asked whether these behaviors were also linked to PPAs. Greater PPA history probability ({sdt_ppa_hx_panel}) was associated with lower dʹ ({ppa_beta_d_prime_overall_reg}), lower decision criterion ({ppa_beta_criterion_overall_reg}), and higher VCH confidence ({ppa_beta_mean_conf_fas_reg}), while greater current PPA count ({sdt_caps_panel}) was confidently associated with lower dʹ ({caps_beta_d_prime_overall_reg}) and decision criterion ({caps_beta_criterion_overall_reg}). It was less confidently associated with VCH confidence ({caps_beta_mean_conf_fas_reg}).
"""


# ================================================================
# SECTION 6: OUTPUT
# ================================================================

def build_results_text():
    """Assemble the full Results narrative in manuscript order.

    The recruitment paragraph opens the section and carries no heading; every
    other block is preceded by its manuscript heading.  Blocks are concatenated
    verbatim (no strip) so that the newlines and trailing spaces written into
    the literals above are what the output file actually contains.

    A missing or empty block raises rather than being skipped: an absent section
    means an upstream result failed to load, and silently dropping it would hide
    that.
    """
    parts = [globals()[LEAD_PARAGRAPH_VAR].lstrip("\n")]
    for label, var_name in RESULT_SECTION_ORDER:
        text = globals().get(var_name)
        if not isinstance(text, str) or not text.strip():
            raise ValueError(
                f"Narrative block '{var_name}' (section '{label}') is missing or empty. "
                "Every section must be present — check the result CSVs it reads."
            )
        parts.append("\n" + label + text)
    return "".join(parts)


def write_results_txt(output_path=OUTPUT_TXT):
    output_path = Path(output_path)
    text = build_results_text()
    output_path.write_text(text, encoding="utf-8")
    print(f"Wrote narrative to: {output_path.resolve()}")
    return output_path


if __name__ == "__main__":
    print()
    print(build_results_text())
    write_results_txt()
