#!/usr/bin/env /usr/local/bin/python3.12
"""
sensitivity_analyses.py
=======================
Supplementary Figure S4 — the compound single-path sensitivity heatmap.

Reads one file and draws one figure. It fits nothing, submits nothing, and
opens no connection to the cluster: every model type it renders is fitted by
03_hpc/generate_hpc_jobs.py and compiled by 03_hpc/compile_nonsp_results.py.

Input
-----
    results/sensitivity_analyses_single_paths/
        existingresults_manuscript_counterfactual.csv     (HEATMAP_SOURCE_CSV)

One CSV, deliberately. Every column — the three canonical model types and the
nine sensitivity variants alike — is read from it, so a heatmap cannot mix rows
from two different cluster runs. Cells show the posterior probability of
direction from the response-scale marginal contrast, not the log-scale
coefficient.

Output
------
    results/supplement/sensitivity_analyses/supplementary_figure_s4.png
    results/supplement/sensitivity_analyses/supplementary_figure_s4.tiff

hppd_binary on top, caps_vision below. Both DVs must be present in the source
CSV or the script raises — a one-DV figure is not the figure the manuscript
cites.

Usage (from anywhere)
---------------------
    /usr/local/bin/python3.12 04_visualizations/supplement/sensitivity_analyses.py

Two ways this can go wrong quietly, both now made loud:
  - A model type in MODEL_VARIANTS with no rows in the source CSV is dropped
    with a warning, so it shows up as a missing column rather than a blank one.
    Keep MODEL_VARIANTS in step with BASE_MODELS in generate_hpc_jobs.py.
  - A missing DHARMa column raises rather than leaving cells silently
    unflagged. All eleven are written by compile_nonsp_results.py; if one is
    absent the compile is stale and re-pulling is the fix.

Parallel-maintenance: this script and sensitivity_analyses_mediation.py share
their display conventions (colour map, cell format, flagging logic, compound
layout). Change one, change the other.
"""

# ==============================================================================
# ████████████████████████  CONFIG — edit here only  ███████████████████████████
# ==============================================================================

# ── DVs and predictor groups ──────────────────────────────────────────────────
PART1_DVS              = ['hppd_binary', 'caps_vision']
PART1_PREDICTOR_GROUPS = ['sp_predictors', 'vch_behavior', 'vch_computations']

# ── Diagnostic flagging ───────────────────────────────────────────────────────
# A cell is flagged when any DHARMa test for that model has p below this.
DHARMA_ALPHA       = 0.05

# ── Heatmap input: single source of truth ─────────────────────────────────────
# Changed 2026-08-23 (Max).  The heatmap previously stitched two compiled CSVs
# together: the three canonical columns came from the manuscript CSV and the six
# sensitivity columns from existingresults_sensitivity_counterfactual.csv, which
# only exists if this script's own compile stage has been run against its own
# separate cluster array.
#
# The 2026-08-23 full reporting run (generate_hpc_jobs.py, array 23252250) fits
# all nine model types in one array, so every column the heatmap needs is now in
# the single CSV that compile_nonsp_results.py writes.  Reading one file removes
# the risk of a heatmap mixing rows from two different cluster runs.
#
# Filename is relative to LOCAL_COMPILE_DIR (resolved at runtime below).
HEATMAP_SOURCE_CSV = 'existingresults_manuscript_counterfactual.csv'

# Column order in the heatmap: these three lead, then MODEL_VARIANTS in order.
# All are read from HEATMAP_SOURCE_CSV.
CANONICAL_MODEL_TYPE  = 'nice_covariates_spusers'
SECOND_CANONICAL_TYPE = 'nice_covariates'             # full sample; column immediately right of canonical
THIRD_CANONICAL_TYPE  = 'nice_covariates_spusers_iqr' # SP users, IQR outlier exclusion; column after nice_covariates
BULK_ESS_MIN       = 1000
NUM_DIVERGENTS_STRICT_MIN = 1

# (model_type, dv) pairs to gray out in the heatmap because they are
# logically inapplicable.  nonan_caps restricts to participants with
# non-NaN caps_bl_1 — which is the full sample for any caps_* DV, making
# that column redundant (identical to the primary model).
INAPPLICABLE_CELLS = {
    ('nice_covariates_spusers_nonan_caps', 'caps_vision'),
}

# ==============================================================================
# END CONFIG ───────────────────────────────────────────────────────────────────
# ==============================================================================


# ==============================================================================
# SETUP
# ==============================================================================

import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import warnings
from itertools import chain, combinations
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.image as mpimg
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')

# ── Paths ─────────────────────────────────────────────────────────────────────
_SCRIPT_DIR   = Path(__file__).resolve().parent   # 04_visualizations/supplement/
_PROJECT_ROOT = _SCRIPT_DIR.parent.parent         # hppd_manuscript_public/

sys.path.insert(0, str(_PROJECT_ROOT / 'modules'))

# Local results directories
LOCAL_COMPILE_DIR = _PROJECT_ROOT / 'results' / 'sensitivity_analyses_single_paths'
LOCAL_HEATMAP_DIR = _PROJECT_ROOT / 'results' / 'supplement' / 'sensitivity_analyses'
# ── Heatmap columns ───────────────────────────────────────────────────────────
# The sensitivity model types this figure renders as columns, in display order,
# after the three canonical types defined above.  Every one is fitted on the
# cluster by 03_hpc/generate_hpc_jobs.py and compiled into HEATMAP_SOURCE_CSV;
# this script only reads that CSV.
#
# This list and BASE_MODELS in 03_hpc/generate_hpc_jobs.py must stay in step. A
# type named here but absent from the CSV is dropped with a warning below rather
# than drawn as an empty column, so a silent gap shows up as a missing column,
# not a blank one.
#
#   beta variants        + vch_beta as a covariate (decision-noise control)
#   nonan_caps           restricted to participants with non-NaN caps_bl_1
#   drugs_month          + alc/ghb/opioids/mj/atypicals/stimulants past-month flags
#   drugs_trimmed_month  + depressants/mj/stimulants (sparse classes collapsed)
#   hardware_control     + monitor_check_operationalized_final (3-level display class)
MODEL_VARIANTS = [
    'empirical_covariates_spusers',
    'age_control_spusers',
    'true_univariate_spusers',
    'nice_covariates_beta_spusers',
    'nice_covariates_beta_spusers_iqr',
    'nice_covariates_spusers_nonan_caps',
    'drugs_month_spusers',
    'drugs_trimmed_month_spusers',
    'nice_covariates_spusers_hardware_control',
]

print(f'Sensitivity model types rendered as columns ({len(MODEL_VARIANTS)}):')
for v in MODEL_VARIANTS:
    print(f'  {v}')












# ==============================================================================
# COMPILE FUNCTIONS
# (adapted from compile_nonsp_results.py with per-quantile DHARMa columns added)
# ==============================================================================

# Full set of DHARMa columns to extract from each diagnostics CSV.
# Includes per-quantile p-values (q25/q50/q75) added in the June 2026 fix
# to nonsp_predictors.R — these are NOT in compile_nonsp_results.py's DHARMA_COLS.
# Columns not present in a given CSV are safely stored as NaN.
_DHARMA_COLS_COMPILE = [
    'dharma_uniformity_pval',
    'dharma_dispersion_pval',
    'dharma_outlier_pval',
    'dharma_quantiles_pval',              # combined BH-adjusted vs fitted values
    'dharma_quantiles_q25_pval',          # individual quantile line (q25)
    'dharma_quantiles_q50_pval',          # individual quantile line (q50)
    'dharma_quantiles_q75_pval',          # individual quantile line (q75)
    'dharma_heteroscedasticity_pval',     # combined BH-adjusted vs spvar
    'dharma_heteroscedasticity_q25_pval', # individual quantile line (q25) vs spvar
    'dharma_heteroscedasticity_q50_pval',
    'dharma_heteroscedasticity_q75_pval',
]

_RHAT_MIN, _RHAT_MAX = 0.9, 1.1














# ==============================================================================
# HEATMAP FUNCTIONS
# (adapted from nonsp_predictor_heatmap_summary.py; uses counterfactual CSV)
# ==============================================================================

# All DHARMa columns checked for p < DHARMA_ALPHA → red 'flagged' box.
# Includes per-quantile columns not present in older compiled CSVs.
_DHARMA_COLS_HEATMAP = _DHARMA_COLS_COMPILE  # same full set

_CMAP = LinearSegmentedColormap.from_list(
    'sens_prob_cmap',
    ['#f7f7f7', '#fff7bc', '#fee391', '#a1dab4', '#41b6c4', '#225ea8'],
    N=100,
)


def _heatmap_label(predictor):
    """Short display label for a predictor."""
    from master_config import dv_to_lab_short
    norm = f'{predictor}_normalized' if not predictor.endswith('_normalized') else predictor
    return dv_to_lab_short.get(norm, dv_to_lab_short.get(predictor, predictor.replace('_', ' ').title()))


# Human-readable descriptions for R-side subsetting keywords used in model type names.
_MODIFIER_LABELS = {
    'spusers':              'SP users only',
    'nooutlier':            'no outliers',
    'nooutliers':           'no outliers',
    'nopsychosis':          'no psychosis',
    'nocurrenthppd':        'excl. current HPPD',
    'iqr':                  'IQR outlier filter',
    'beta':                 '+ vch_beta covariate',
}


def _format_model_type_label(s):
    """
    Full covariate-set name with subset modifiers on a second line in parentheses.
    'empirical_covariates_spusers'      → 'empirical_covariates\n(SP users only)'
    'age_control_spusers'               → 'age_control\n(SP users only)'
    Unrecognised modifiers are passed through unchanged.
    """
    remaining = s
    found = []
    # Strip modifiers longest-first to avoid partial-match collisions
    for mod in sorted(_MODIFIER_LABELS.keys(), key=len, reverse=True):
        suffix = f'_{mod}'
        if remaining.endswith(suffix):
            remaining = remaining[:-len(suffix)]
            found.insert(0, mod)
    if not found:
        return s
    descriptions = [_MODIFIER_LABELS[m] for m in found]
    return f'{remaining}\n({", ".join(descriptions)})'


def _diagnostic_status_cf(cf_row):
    """
    Return 'flagged', 'missing', or 'ok' for one counterfactual row.

    Checks all _DHARMA_COLS_HEATMAP p-values (any p < DHARMA_ALPHA → flagged),
    num_divergents > 0, and (Bulk_ESS < BULK_ESS_MIN AND num_divergents > 1).

    The counterfactual CSV has one row per (spvar, cov, dv) so we check the
    scalar values directly rather than aggregating across coefficient rows.
    """
    # Sampler quality
    if 'Bulk_ESS' in cf_row.index and 'num_divergents' in cf_row.index:
        ess = pd.to_numeric(cf_row.get('Bulk_ESS'), errors='coerce')
        div = pd.to_numeric(cf_row.get('num_divergents'), errors='coerce')
        if pd.notna(ess) and pd.notna(div):
            if float(ess) < BULK_ESS_MIN and float(div) > NUM_DIVERGENTS_STRICT_MIN:
                return 'flagged'

    found_any = False
    for col in _DHARMA_COLS_HEATMAP:
        if col not in cf_row.index:
            continue
        val = pd.to_numeric(cf_row[col], errors='coerce')
        if pd.isna(val):
            continue
        found_any = True
        if float(val) < DHARMA_ALPHA:
            return 'flagged'

    if not found_any:
        return 'missing'

    # Also flag non-zero divergents when diagnostics are present
    if 'num_divergents' in cf_row.index:
        div = pd.to_numeric(cf_row.get('num_divergents'), errors='coerce')
        if pd.notna(div) and float(div) > 0:
            return 'flagged'

    return 'ok'


def _build_heatmap_data_cf(cf_df, predictors, model_types, dv):
    """
    Build (prob_pivot, dir_pivot, status_pivot) for one DV from the
    counterfactual CSV.

    Rows = predictors (in order), Columns = model_types.
    Probability is P(direction) = max(prob_above_0, prob_below_0).
    """
    sub = cf_df[cf_df['dv'] == dv].copy()

    prob_data, dir_data, status_data = {}, {}, {}

    for mt in model_types:
        mt_sub = sub[sub['cov'] == mt]
        prob_data[mt] = dir_data[mt] = status_data[mt] = {}
        prob_data[mt]   = {}
        dir_data[mt]    = {}
        status_data[mt] = {}

        for pred in predictors:
            # Skip logically inapplicable (model_type, dv) combos → gray cell
            if (mt, dv) in INAPPLICABLE_CELLS:
                continue

            spvar_name = f'{pred}_normalized' if not pred.endswith('_normalized') else pred
            cell_rows  = mt_sub[mt_sub['spvar'] == spvar_name]

            if cell_rows.empty:
                continue

            row = cell_rows.iloc[0]
            prob_above = pd.to_numeric(row.get('prob_above_0', 0), errors='coerce') or 0.0
            prob_below = pd.to_numeric(row.get('prob_below_0', 0), errors='coerce') or 0.0

            prob_data[mt][pred]   = max(prob_above, prob_below)
            dir_data[mt][pred]    = 'pos' if prob_above >= prob_below else 'neg'
            status_data[mt][pred] = _diagnostic_status_cf(row)

    prob_pivot   = pd.DataFrame(prob_data,   index=predictors)
    dir_pivot    = pd.DataFrame(dir_data,    index=predictors)
    status_pivot = pd.DataFrame(status_data, index=predictors)
    return prob_pivot, dir_pivot, status_pivot


def _plot_heatmap_ax(ax, cf_df, all_preds, sp_preds, vch_beh,
                     cols_with_data, dv, vch_comp=None,
                     show_xticklabels=True,
                     x_tick_fs=8, y_tick_fs=8.5, cell_text_fs=6.5):
    """
    Render a sensitivity heatmap for one DV onto an existing Axes.

    Called by both plot_sensitivity_heatmap (standalone figure) and
    plot_compound_sensitivity_heatmap (two-panel figure).

    Parameters
    ----------
    ax               : matplotlib Axes to draw on
    cf_df            : counterfactual DataFrame (all DVs / model types)
    all_preds        : full ordered predictor list
                       (sp + vch_beh + vch_comp + sdt_hppd)
    sp_preds         : SP predictor names (for boundary lines)
    vch_beh          : VCH behavior predictor names (for boundary lines)
    vch_comp         : VCH computational predictor names (for the boundary that
                       separates them from the SDT block).  Optional: pass None
                       to reproduce the pre-2026-08-31 two-boundary layout.
    cols_with_data   : ordered list of model type column keys to display
    dv               : DV name string (used to filter cf_df)
    show_xticklabels : if False, suppress x-tick text (e.g. top panel)
    x_tick_fs        : font size for x-tick labels
    y_tick_fs        : font size for y-tick labels

    Returns
    -------
    im               : AxesImage from imshow (needed for shared colorbar), or None
    rows_with_data   : Index of predictor rows that have any data
    """
    from master_config import COVARIATE_SET_LABELS

    prob_pivot, dir_pivot, status_pivot = _build_heatmap_data_cf(
        cf_df, all_preds, cols_with_data, dv
    )

    rows_with_data = prob_pivot.index[prob_pivot.notna().any(axis=1)]
    if rows_with_data.empty:
        return None, rows_with_data

    prob_pivot   = prob_pivot.loc[rows_with_data]
    dir_pivot    = dir_pivot.loc[rows_with_data]
    status_pivot = status_pivot.loc[rows_with_data]

    n_preds, n_models = prob_pivot.shape

    sp_end   = sum(1 for p in rows_with_data if p in sp_preds)
    vchb_end = sp_end + sum(1 for p in rows_with_data if p in vch_beh)
    # Fourth block (SDT: criterion / d' / mean confidence on FAs) added
    # 2026-08-31, so a boundary is also needed after the VCH computational rows.
    _cut = [sp_end, vchb_end]
    if vch_comp:
        _cut.append(vchb_end + sum(1 for p in rows_with_data if p in vch_comp))
    boundaries = [b - 0.5 for b in _cut if 0 < b < n_preds]

    im = ax.imshow(
        prob_pivot.values.astype(float),
        cmap=_CMAP, aspect='auto', vmin=0.5, vmax=1.0,
    )
    nan_mask = np.ma.masked_where(
        ~np.isnan(prob_pivot.values.astype(float)),
        np.ones_like(prob_pivot.values.astype(float)),
    )
    ax.imshow(nan_mask, cmap=plt.get_cmap('Greys'), aspect='auto', vmin=0, vmax=1, alpha=0.2)

    ax.set_xticks(range(n_models))
    ax.set_yticks(range(n_preds))

    if show_xticklabels:
        ax.set_xticklabels(
            [COVARIATE_SET_LABELS.get(m, m) for m in cols_with_data],
            rotation=35, ha='right', fontsize=x_tick_fs, multialignment='center',
        )
    else:
        ax.set_xticklabels([])
        ax.tick_params(axis='x', which='both', bottom=False, top=False)

    ax.set_yticklabels([_heatmap_label(p) for p in rows_with_data], fontsize=y_tick_fs)

    for boundary in boundaries:
        ax.axhline(boundary, color='#555555', linewidth=1.4, linestyle='-', zorder=5)

    for ri, pred in enumerate(rows_with_data):
        for ci, mt in enumerate(cols_with_data):
            prob   = prob_pivot.at[pred, mt]
            direc  = dir_pivot.at[pred, mt]   if pred in dir_pivot.index   else 'pos'
            status = status_pivot.at[pred, mt] if pred in status_pivot.index else 'missing'

            if not pd.isna(prob):
                sign    = '+' if direc == 'pos' else '−'
                txt_col = 'white' if prob > 0.92 else 'black'
                ax.text(ci, ri, f'{sign}{prob:.2f}', ha='center', va='center',
                        fontsize=cell_text_fs, color=txt_col, fontweight='bold')

            if status == 'flagged':
                ax.add_patch(mpatches.Rectangle(
                    (ci - 0.48, ri - 0.48), 0.96, 0.96,
                    fill=False, edgecolor='red', linewidth=2.0,
                ))
            elif status == 'missing':
                ax.add_patch(mpatches.Rectangle(
                    (ci - 0.48, ri - 0.48), 0.96, 0.96,
                    fill=False, edgecolor='#e991c8', linewidth=2.0,
                ))

    ax.set_xticks(np.arange(n_models + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(n_preds   + 1) - 0.5, minor=True)
    ax.grid(which='minor', color='#cccccc', linestyle='-', linewidth=0.35)
    ax.tick_params(which='minor', length=0)

    return im, rows_with_data


def plot_compound_sensitivity_heatmap(cf_df, model_types, out_dir):
    """
    Compound two-panel sensitivity heatmap: hppd_binary (top) + caps_vision (bottom).

    Both panels share the same x-axis column order.  X-tick labels and the x-axis
    label are shown on the bottom panel only.  A single y-axis label is placed
    between the two panels using fig.text() (2.5× the single-panel font size).
    Y-tick labels appear on both panels.  All tick label font sizes are doubled
    relative to the single-panel defaults.  Arial font is used throughout.

    Layout uses GridSpec with an explicit narrow colorbar column so the colorbar
    never overlaps the heatmap cells.  The left margin is wide enough that the
    shared y-axis label does not collide with the (doubled-size) y-tick labels.

    Output: {out_dir}/supplementary_figure_s4.png
    """
    import matplotlib as mpl
    from matplotlib.gridspec import GridSpec
    from master_config import iv_type_dict, dv_to_lab_short

    dvs       = ['hppd_binary', 'caps_vision']
    sp_preds  = iv_type_dict['sp_predictors']
    vch_beh   = iv_type_dict['vch_behavior']
    vch_comp  = iv_type_dict['vch_computations']
    # SDT / metacognition predictors (criterion_overall, d_prime_overall,
    # mean_conf_fas) added to the y-axis 2026-08-31.  Their jobs come from
    # 03_hpc/generate_hpc_jobs.py (PART1_PREDICTOR_GROUPS includes 'sdt_hppd'),
    # so this script's own generation stage is deliberately NOT extended --
    # adding the group to PART1_PREDICTOR_GROUPS here would emit a duplicate
    # copy of every SDT job.  Heatmap rows only.
    sdt_preds = iv_type_dict['sdt_hppd']
    all_preds = sp_preds + vch_beh + vch_comp + sdt_preds

    # Canonical types first, then remaining sensitivity types — same order for both panels.
    _canonical = [CANONICAL_MODEL_TYPE, SECOND_CANONICAL_TYPE, THIRD_CANONICAL_TYPE]
    _has_any   = lambda mt: any(
        not cf_df[(cf_df['dv'] == dv) & (cf_df['cov'] == mt)].empty for dv in dvs
    )
    cols_with_data = (
        [mt for mt in _canonical if _has_any(mt)]
        + [mt for mt in model_types if mt not in _canonical and _has_any(mt)]
    )
    if not cols_with_data:
        print('  [compound] no model columns with data — skipping')
        return

    n_models = len(cols_with_data)

    # Per-panel row counts (for height_ratios)
    row_counts = {}
    for dv in dvs:
        pp, _, _ = _build_heatmap_data_cf(cf_df, all_preds, cols_with_data, dv)
        row_counts[dv] = len(pp.index[pp.notna().any(axis=1)])

    cell_w = max(2.2, min(3.2, 36 / max(n_models, 1)))
    cell_h = max(0.45, min(0.75, 16 / max(sum(row_counts.values()), 1)))
    # Extra width: left margin (for y-label + doubled tick labels) + right (colorbar + label)
    fig_w  = n_models * cell_w + 9.0
    fig_h  = (row_counts[dvs[0]] + row_counts[dvs[1]]) * cell_h + 6.0

    # Font sizes
    X_TICK_FS = 16     # doubled from single-panel default of 8
    Y_TICK_FS = 17     # doubled from single-panel default of 8.5
    XLABEL_FS = 20     # doubled from 10
    YLABEL_FS = 25     # 2.5× single-panel default of 10
    TITLE_FS  = 26     # doubled from 13

    n0, n1 = row_counts[dvs[0]], row_counts[dvs[1]]

    with mpl.rc_context({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    }):
        fig = plt.figure(figsize=(fig_w, fig_h))

        # GridSpec: 2 heatmap rows × [heatmap column | narrow colorbar column].
        # Explicit left/right margins ensure the shared y-label and colorbar label
        # both have dedicated space and cannot overlap the heatmap cells.
        #
        # left=0.20  → 20 % of fig_w for y-label + doubled y-tick labels
        # right=0.93 → 7 % of fig_w for the colorbar column + its rotated label
        gs = GridSpec(
            2, 2,
            figure=fig,
            width_ratios=[1, 0.03],   # colorbar column ≈ 3 % of content width
            height_ratios=[n0, n1],
            hspace=0.08,
            wspace=0.03,
            left=0.20,
            right=0.93,
            bottom=0.22,
            top=0.95,
        )

        axes = [fig.add_subplot(gs[i, 0]) for i in range(2)]
        cax  = fig.add_subplot(gs[:, 1])

        ims = []
        for i, dv in enumerate(dvs):
            ax        = axes[i]
            is_bottom = (i == len(dvs) - 1)
            im, _     = _plot_heatmap_ax(
                ax, cf_df, all_preds, sp_preds, vch_beh,
                cols_with_data, dv, vch_comp=vch_comp,
                show_xticklabels=is_bottom,
                x_tick_fs=X_TICK_FS,
                y_tick_fs=Y_TICK_FS,
                cell_text_fs=13.0,   # doubled from single-panel default of 6.5
            )
            if im is not None:
                ims.append(im)
            ax.set_title(
                dv_to_lab_short.get(dv, dv),
                fontsize=TITLE_FS, fontweight='bold', pad=8,
            )

        # X-axis label on the bottom panel only
        axes[-1].set_xlabel(
            'Model type (sensitivity)', fontsize=XLABEL_FS,
            fontweight='bold', labelpad=10,
        )

        # Single y-axis label centred between both panels.
        # With left=0.20, the heatmap axes begin at 20 % of fig_w from the left.
        # x=0.13 sits just outside the y-tick labels (which extend no further
        # left than ~14 % of fig_w).  Moved right from 0.03 on 2026-09-02 to
        # close the wide gap the old position left between label and ticks.
        fig.text(
            0.13, 0.5, 'Predictor',
            va='center', ha='center', rotation='vertical',
            fontsize=YLABEL_FS, fontweight='bold',
        )

        # Colorbar in its own dedicated GridSpec column — no space stolen from heatmap.
        if ims:
            cbar = fig.colorbar(ims[0], cax=cax)
            # No colorbar title: the compound figure's caption carries the
            # meaning of the scale, and the axis label crowded the tick labels
            # once they were enlarged.  Ticks at 24 pt (3x the previous 8 pt)
            # so the scale stays readable at the figure's print width.
            cbar.ax.tick_params(labelsize=24)

        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / 'supplementary_figure_s4.png'
        plt.savefig(out_path, dpi=200, bbox_inches='tight')
        tiff_path = out_dir / 'supplementary_figure_s4.tiff'
        plt.savefig(tiff_path, dpi=200, bbox_inches='tight')
        plt.close()
        print(f'  Saved compound heatmap: {out_path}')
        print(f'  Saved compound heatmap: {tiff_path}')




# ==============================================================================
# PDF COMPILATION FUNCTIONS
# (adapted from compile_single_path_diagnostic_pdfs.py)
# ==============================================================================

_PDF_DPI = 100

# Keys used internally to refer to each figure type
_INDIVIDUAL_DIAG_SUFFIXES = [
    'dharma_diagnostics',       # DHARMa uniformity/dispersion/outlier overview
    'linearity_check',          # residuals vs. fitted (linearity)
    'quantiles_check',          # quantile residuals
    'heteroscedasticity_check', # heteroscedasticity vs. predictor
    'error_scatter_vs_x',       # error scatter vs. predictor
    'error_scatter_vs_y',       # error scatter vs. fitted values
    'pp_check',                 # posterior predictive check
    'convergence',              # MCMC convergence (trace/rank)
]

# Subdirectory (under {pred_normalized}/{model_type}/) and filename for each figure type.
# The 6 DHARMa plots live under results/diagnostics/; pp_check and convergence have
# their own top-level subdirectories.
_DIAG_SUBDIR = {
    'dharma_diagnostics':       'results/diagnostics',
    'linearity_check':          'results/diagnostics',
    'quantiles_check':          'results/diagnostics',
    'heteroscedasticity_check': 'results/diagnostics',
    'error_scatter_vs_x':       'results/diagnostics',
    'error_scatter_vs_y':       'results/diagnostics',
    'pp_check':                 'pp_checks',
    'convergence':              'convergence_tests',
}
_DIAG_FILENAME = {
    'dharma_diagnostics':       '{dv}_dharma_diagnostics.png',
    'linearity_check':          '{dv}_linearity_check.png',
    'quantiles_check':          '{dv}_quantiles_check.png',
    'heteroscedasticity_check': '{dv}_heteroscedasticity_check.png',
    'error_scatter_vs_x':       '{dv}_error_scatter_vs_x.png',
    'error_scatter_vs_y':       '{dv}_error_scatter_vs_y.png',
    'pp_check':                 '{dv}_posterior_check_plot.png',
    'convergence':              '{dv}_convergence_plot.png',
}

# 2 rows × 4 cols: DHARMa panels on row 1, error scatter + model checks on row 2
_PANEL_LAYOUT = [
    ['dharma_diagnostics', 'linearity_check',    'quantiles_check',   'heteroscedasticity_check'],
    ['error_scatter_vs_x', 'error_scatter_vs_y', 'pp_check',          'convergence'],
]

_PANEL_TITLES = {
    'dharma_diagnostics':       'DHARMa diagnostics',
    'linearity_check':          'Linearity (vs fitted)',
    'quantiles_check':          'Quantile residuals',
    'heteroscedasticity_check': 'Heteroscedasticity',
    'error_scatter_vs_x':       'Error scatter vs predictor',
    'error_scatter_vs_y':       'Error scatter vs fitted',
    'pp_check':                 'Posterior predictive check',
    'convergence':              'MCMC convergence (trace/rank)',
}
















# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

print(f'\n{"="*65}')
print(f'sensitivity_analyses.py — Supplementary Figure: compound sensitivity heatmap')
print(f'{"="*65}')

if True:
    # Single source of truth: every heatmap column — the three canonical model
    # types and the six sensitivity variants alike — is read from this one CSV,
    # written by 03_hpc/compile_nonsp_results.py after generate_hpc_jobs.py pulls
    # the cluster results.  See HEATMAP_SOURCE_CSV in CONFIG for why.
    cf_csv_path = LOCAL_COMPILE_DIR / HEATMAP_SOURCE_CSV
    if not cf_csv_path.exists():
        print(f'  ERROR: heatmap source CSV not found at {cf_csv_path}')
        print(f'  It is produced by the compile step of 03_hpc/generate_hpc_jobs.py')
        print(f'  (which calls 03_hpc/compile_nonsp_results.py).  Pull the HPC')
        print(f'  results first, then re-run this stage.')
    else:
        cf_df = pd.read_csv(cf_csv_path, low_memory=False)
        print(f'  Loaded {len(cf_df):,} rows from {cf_csv_path.name}')

        # The DHARMa columns the heatmap flags on are written into this CSV by
        # 03_hpc/compile_nonsp_results.py.  Stop rather than render a heatmap
        # with silently unflagged cells if a stale compile is missing any.
        missing_dharma = [c for c in _DHARMA_COLS_HEATMAP if c not in cf_df.columns]
        if missing_dharma:
            raise KeyError(
                f'{cf_csv_path.name} is missing {len(missing_dharma)} DHARMa '
                f'columns the heatmap flags on: {missing_dharma}. Re-pull the '
                'cluster results with 03_hpc/compile_nonsp_results.py rather '
                'than plotting unflagged cells.'
            )

        types_present = set(cf_df['cov'].unique())
        print(f'  DVs present:         {sorted(cf_df["dv"].unique())}')
        print(f'  Model types present: {sorted(types_present)}')

        # Column order: the three canonical types lead, then MODEL_VARIANTS.
        # Every column must be present in the source CSV — a type that is absent
        # is dropped with a warning rather than rendered as an empty column.
        _wanted = [CANONICAL_MODEL_TYPE, SECOND_CANONICAL_TYPE, THIRD_CANONICAL_TYPE] + list(MODEL_VARIANTS)
        _wanted = list(dict.fromkeys(_wanted))   # dedupe, preserve order
        model_types_in_csv = [mt for mt in _wanted if mt in types_present]
        missing_types      = [mt for mt in _wanted if mt not in types_present]
        if missing_types:
            print(f'  WARNING: these model types are not in {cf_csv_path.name} and will be '
                  f'omitted from the heatmap: {missing_types}')
            print(f'           Check that BASE_MODELS in 03_hpc/generate_hpc_jobs.py covers '
                  f'them and that their cluster tasks completed.')

        # The compound figure (hppd_binary on top, caps_vision below) is the
        # supplementary figure; both DVs must be present or it is not the
        # figure the manuscript cites.
        missing_dvs = [dv for dv in PART1_DVS if dv not in cf_df['dv'].unique()]
        if missing_dvs:
            raise ValueError(
                f'{cf_csv_path.name} has no rows for {missing_dvs}. The compound '
                'heatmap needs both DVs; re-pull the cluster results.'
            )
        plot_compound_sensitivity_heatmap(cf_df, model_types_in_csv, LOCAL_HEATMAP_DIR)

