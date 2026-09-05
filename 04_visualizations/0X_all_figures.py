##############################################################################
### MASTER RAPID FIGURE ASSEMBLY FOR MANUSCRIPT
# Regenerates all primary manuscript figures from a single entry point.
# Change MODEL_TYPE (and optionally the section toggles) at the top to
# reproduce all figures with a different covariate set or output format.
#
# FONT: ALL figures MUST use Arial.  Global rcParams are set after the
# matplotlib import (search "GLOBAL FONT").  Any code that calls
# sns.set_style() must re-enforce Arial immediately afterward.
# PIL-assembled figures (3, 6) load Arial Bold.ttf explicitly.
#
# Can be run from any directory — the script changes to its own directory
# at startup so all relative paths (../results/, figures/, linked_figures/)
# resolve correctly regardless of invocation location.
#   python 04_visualizations/0X_all_figures.py        # from project root
#   python 0X_all_figures.py                         # from 04_visualizations/
#   python ../04_visualizations/0X_all_figures.py     # from 03_hpc/
##############################################################################
import os as _os
_os.chdir(_os.path.dirname(_os.path.abspath(__file__)))


##############################################################################
### ═══════════════════════════════════════════════════════════════════════ ###
### ██████████████████████████   C O N F I G   ██████████████████████████ ###
### ═══════════════════════════════════════════════════════════════════════ ###
##############################################################################

# ---------------------------------------------------------------------------
# Model type — the covariate string used in directory and file names.
# Change this to regenerate all figures with a different covariate set.
#   Examples: "empiric_covariates", "nice_covariates"
# ---------------------------------------------------------------------------
MODEL_TYPE = "nice_covariates"

# ---------------------------------------------------------------------------
# Section toggles — set to False to skip a section entirely.
# ---------------------------------------------------------------------------
RUN_HPPD_CAPS_FIGS     = True   # boxplots, correlation grids, forest plots
RUN_TABLES             = True   # demographics, clinical, SP publication tables
RUN_DESCRIPTIVE_FIGS   = True   # SP distribution violin/bar summary figures
RUN_MEDIATION_DIAGRAMS = True   # batch-create missing mediation diagrams
RUN_FIGURE_ASSEMBLY    = True   # assemble multi-panel manuscript figures
RUN_DIAGNOSTICS        = True   # DHARMa diagnostic review for manuscript mediation panels (Figs 2+4)
RUN_DIAGNOSTIC_COMPILATION = True  # generate diagnostic_compilation.png per mediation model (calls R, slow)
RUN_SINGLE_PATH_DIAGNOSTICS = True # assemble single-path nonsp PDFs from HPC_RESULTS_MIRROR (no R needed)
RUN_LINKED_FIGURES     = True   # copy assembled figures to linked_figures/

# ---------------------------------------------------------------------------
# Categorical tests in Table 1
#
# Chi-square is invalid when any expected cell count is < 5.  For 2x2 tables
# (binary indicator rows such as "Any psychiatric diagnosis: Yes") Fisher's
# exact test is substituted automatically and recorded as
# statistic_type = 'fisher_exact' in the sidecar CSV.
#
# Multi-level categoricals (sex, race, education) cannot be: scipy provides no
# general RxC Fisher exact.  For those, chi-square runs on the subset of levels
# meeting the expected-count threshold, and 'NA' appears only when fewer than
# two levels survive.
# ---------------------------------------------------------------------------

# Forest plot model names (must match covariate strings in existingresults_manuscript.csv).
# HPPD uses the SP-user sample, so append _spusers. Add more entries here if additional
# covariate sets are present in the CSV (e.g. "univariate_spusers", "age_control_spusers").
FOREST_MODELS_HPPD = [
    f"{MODEL_TYPE}_spusers" if 'spusers' not in MODEL_TYPE else MODEL_TYPE,
]

# Canonical HPPD model-type suffix used wherever HPPD model directories are read.
HPPD_MODEL_TYPE = FOREST_MODELS_HPPD[0]

# ---------------------------------------------------------------------------
# Per-figure mediation panel configuration.
#
# FIG2_MED_*      — drives panels e + f (row 1, vch_threshold) on figure 2 (vch_behavior).
# FIG2_MED_ROW2_* — drives panels g + h (row 2, vch_bl_yes_0) on figure 2 (vch_behavior).
# FIG4_MED_*      — drives panels c + d on the vch_computations figure (figure 4).
#
# Figure 2 layout:
#   Rows a–d: boxplots (a, b) and forest plots (c, d)
#   Row e–f: vch_threshold mediator row
#     e: spage  → vchthreshold → hppd_binary
#     f: avgdose → vchthreshold → caps_vision
#   Row g–h: vch_bl_yes_0 mediator row
#     g: spage  → vchrate → hppd_binary
#     h: avgdose → vchrate → caps_vision
#
# spvar shorthand  →  full SP variable name used in diagram filenames:
#   "spage"    → "psychedelic_age"
#   "avgdose"  → "avg_life_dose"
#   "lifenomic"→ "psycheduse_life_nomic"
# ---------------------------------------------------------------------------
_SPVDR_TO_DIAGRAMNAME = {
    "spage":     "psychedelic_age",
    "avgdose":   "avg_life_dose",
    "lifenomic": "psycheduse_life_nomic",
}

# Figure 2 (vch_behavior): mediation panels e + f — row 1 (vch_threshold mediator)
FIG2_MED_HPPD_SPVDR     = "spage"         # hppd_binary_{spvar}_{mediator}_{model_type}
FIG2_MED_HPPD_MEDIATOR  = "vchthreshold"
FIG2_MED_CAPS_SPVDR     = "avgdose"        # caps_vision_{spvar}_{mediator}_{model_type}
FIG2_MED_CAPS_MEDIATOR  = "vchthreshold"

# Figure 2 (vch_behavior): mediation panels g + h — row 2 (vch_bl_yes_0 mediator)
FIG2_MED_ROW2_HPPD_SPVDR    = "spage"
FIG2_MED_ROW2_HPPD_MEDIATOR = "vchrate"
FIG2_MED_ROW2_CAPS_SPVDR    = "avgdose"
FIG2_MED_ROW2_CAPS_MEDIATOR = "vchrate"

# Figure 3 continued (panels l–o, appended below HGF composite)
FIG4_MED_HPPD_SPVDR     = "spage"
FIG4_MED_HPPD_MEDIATOR  = "vchbeta"
FIG4_MED_CAPS_SPVDR     = "avgdose"
FIG4_MED_CAPS_MEDIATOR  = "vchbeta"

# IV type groups for which we run boxplot/correlation + forest figures
IVTYPES = ["sp_predictors", "vch_behavior", "vch_computations", "sdt_hppd"]

# ---------------------------------------------------------------------------
# Paths  (relative to 04_visualizations/ — the expected working directory)
# ---------------------------------------------------------------------------
RESULTS_BASE          = "../results"
COUNTERFACTUAL_CSV    = f"{RESULTS_BASE}/sensitivity_analyses_single_paths/existingresults_manuscript_counterfactual.csv"
TABLES_OUT_DIR        = f"{RESULTS_BASE}/descriptive/tables"
FIGURES_OUT_DIR       = "figures/caps_vision_hppd_combined_4_figs"

# ---------------------------------------------------------------------------
# Output format — applies to ALL figures: assembled multi-panel figures AND
# individual panel figures (boxplots, correlation grids, forest plots,
# mediation diagrams). PNG copies are always kept alongside so assembly works.
#   Options: "png", "jpeg", "tiff", "svg", "pdf"
# ---------------------------------------------------------------------------
FIGURE_FORMAT = "png"
FIGURE_DPI    = 600

# Formats for journal-submission copies saved to results/final_figures/.
# Every format listed is written for every figure, so results/final_figures/
# carries a complete set in each.  "svg" is the vector master Scientific Reports
# asks for; "tiff" is the raster fallback some journals require; "png" is for
# quick viewing and for pasting into drafts.  Also accepts "pdf", "jpeg".
FINAL_FIGURE_FORMATS = ("svg", "tiff", "png")

# ---------------------------------------------------------------------------
# Effect scale — counterfactual marginal contrasts
#
# Every reported effect is a response-scale marginal contrast:
#   E[Y | X = mean + 1 SD] - E[Y | X = mean], computed by posterior_epred().
#     - bernoulli (hppd_binary):      risk difference (probability scale)
#     - hurdle_negbinom / negbinom:   expected count difference
#     - student_t / gaussian:         raw score difference
# These are reported in preference to standardized log-scale coefficients
# because they are interpretable in the units of the outcome.  Forest plots
# therefore load COUNTERFACTUAL_CSV and call counterfactual_forest_plot(), with
# the reference line at 0 (not 1) and no exp() or signed-reciprocal transform.
# The x-axis label comes from cfg["counterfactual_x_label"] in DV_CONFIGS.
#
# Sample: mediation and forest-plot models are fitted in SP users only
# (df_sp, n ~ 186, model nice_covariates_spusers).  The mediation question is
# about variation in exposure among the exposed, so SP-naive participants carry
# no information about it.
#
# Figure layout that follows from this:
#   Figure 2 (vch_behavior), row 3 adds panels e + f:
#     e = hppd_binary FIG2_MED_HPPD_SPVDR->FIG2_MED_HPPD_MEDIATOR diagram
#     f = caps_vision FIG2_MED_CAPS_SPVDR->FIG2_MED_CAPS_MEDIATOR diagram
#   Figure 3 continued (panels l-o, appended below the HGF composite):
#     l = hppd_binary vch_computations forest plot
#     m = caps_vision vch_computations forest plot
#     n = hppd_binary FIG4_MED_HPPD_SPVDR->FIG4_MED_HPPD_MEDIATOR diagram
#     o = caps_vision FIG4_MED_CAPS_SPVDR->FIG4_MED_CAPS_MEDIATOR diagram
#
# Column mapping applied to COUNTERFACTUAL_CSV before plotting:
#   mean     -> Estimate   (POINT_ESTIMATE_COL, modules/master_config.py)
#   spvar    -> var
#   cov      -> covariates
#   settings -> model
# ---------------------------------------------------------------------------

# Forest plot interval type (applies to standard mode only; counterfactual
# always uses hdi_lower_94 / hdi_upper_94).
# "hdi_94"  → use 94% Highest Density Interval
# "ci_95"   → use 95% equal-tailed credible interval
# 94% HDI is the reported interval throughout the manuscript.
_FOREST_LCI_COL = "hdi_lower_94"
_FOREST_UCI_COL = "hdi_upper_94"

# ---------------------------------------------------------------------------
# Standardized manuscript font sizes.
#
# TARGET_* values are the desired font sizes IN THE FINAL ASSEMBLED FIGURE
# (in points, at journal double-column width = 183 mm ≈ 7.2").  All source-panel font
# sizes are back-computed from these targets so that regardless of a panel's
# native figsize, text appears at a uniform size after assembly scaling.
#
# The math:  scale = row_display_height / source_height
#            source_pt = target_pt / scale
# See compute_source_fontsize() in modules/figure_assembly.py.
# ---------------------------------------------------------------------------
import sys as _sys
_sys.path.insert(0, '../modules')
_sys.path.insert(0, 'supplement/diagnostics')  # diagnostic PDF compilation scripts
from figure_assembly import compute_source_fontsize
TARGET_AXIS_LABEL  = 9.2    # y-axis labels, x-axis labels (+15%)
TARGET_TICK_LABEL  = 7.0    # tick labels on both axes
TARGET_SIG_MARKER  = 24.0   # significance asterisks (*, **, ***)
TARGET_XLABEL      = 8.6    # x-axis labels (+15%, matched to forest xlabel)

# ---------------------------------------------------------------------------
# Journal sizing — Scientific Reports (Nature)
# Change these two constants when targeting a different journal.
# ---------------------------------------------------------------------------
JOURNAL_SINGLE_COL_MM = 89    # single-column figure width
JOURNAL_DOUBLE_COL_MM = 183   # double-column (full-page) figure width
_MM_PER_INCH = 25.4

JOURNAL_SINGLE_COL_IN = JOURNAL_SINGLE_COL_MM / _MM_PER_INCH  # ≈ 3.504"
JOURNAL_DOUBLE_COL_IN = JOURNAL_DOUBLE_COL_MM / _MM_PER_INCH  # ≈ 7.205"

# Assembly geometry shared across all figures
_FIG_W  = JOURNAL_DOUBLE_COL_IN   # all figures use double-column width
_GAP    = 0.05   # horizontal gap between panels (inches)

# Panel dimensions for compute_source_fontsize.
# Both ASPECT RATIOS and HEIGHTS should reflect the *actual* saved panel
# dimensions so the font-size back-calculation exactly cancels the
# assembly scaling.  bbox_inches='tight' changes saved dimensions from the
# figsize — e.g. the forest plot shrinks 10.5% vertically (4.8" → 4.30")
# while the boxplot shrinks only 0.6% (15.0" → 14.90").  Using figsize
# heights causes forest fonts to render ~10% larger than boxplot fonts
# at the same target size.  Reading actual heights from disk eliminates this.
from PIL import Image as _PIL_Image_early
import glob as _glob

def _actual_panel_dims(pattern, fallback_ar, fallback_h):
    """Read AR and actual height (inches) of an existing panel from disk.

    Returns (aspect_ratio, height_inches).  Falls back to
    (fallback_ar, fallback_h) when no matching file exists (e.g. first run).
    """
    matches = _glob.glob(pattern)
    if matches:
        with _PIL_Image_early.open(matches[0]) as _img:
            w_px, h_px = _img.size
            dpi_info = _img.info.get('dpi')
            if dpi_info:
                actual_h = h_px / dpi_info[1]
            else:
                actual_h = fallback_h
        return w_px / h_px, actual_h
    return fallback_ar, fallback_h

BOXPLOT_PANEL_HEIGHT = 5  # height per row — used by multipanel_boxplot_grid

_BOXPLOT_AR, _BOXPLOT_H = _actual_panel_dims(
    '../results/hppd_binary/sp_predictors/data_visualization/boxplot_grid.png',
    6.0 / 15.0, 15.0)
_CORR_AR, _CORR_H = _actual_panel_dims(
    '../results/caps_vision/sp_predictors/data_visualization/correlation_grid_age_control.png',
    8.0 / 9.0, 9.0)
_FOREST_AR, _FOREST_H = _actual_panel_dims(
    '../results/hppd_binary/sp_predictors/forest_plots/*.png',
    9.5 / 4.8, 4.8)
_MED_AR = 10.0 / 6.0  # mediation diagrams: no existing panel to read
_MED_H  = 6.0          # figsize height for mediation diagrams (bbox_tight effect is small)

# ── Pre-computed source font sizes ────────────────────────────────────────
# Each FONT_* dict holds the source font sizes for one panel-type ×
# row-context combination.  Values are floats (rounded to int at call sites).
#
# Both aspect ratios AND heights are read from saved panels on disk, so
# the back-calculation exactly cancels the assembly scaling (including
# the bbox_inches='tight' distortion).
#
# _DATAVIZ_BOOST: empirical scaling factor for the data-viz row.  The
# assembly's bbox_inches='tight' crop interacts with the tall data-viz
# row differently than the short forest row, causing data-viz text to
# render ~25% smaller than the math predicts.  This factor compensates.
_DATAVIZ_XLABEL_BOOST = 1.3

# Row context: 2-panel data-viz row (boxplot + correlation), Figs 1/2/5
_ROW_DATAVIZ_ARS = [_BOXPLOT_AR, _CORR_AR]

FONT_BOXPLOT_3VAR = {
    'ylabel':  compute_source_fontsize(TARGET_AXIS_LABEL,                       _BOXPLOT_H, _FIG_W, _GAP, _ROW_DATAVIZ_ARS),
    'xlabel':  compute_source_fontsize(TARGET_XLABEL * _DATAVIZ_XLABEL_BOOST,   _BOXPLOT_H, _FIG_W, _GAP, _ROW_DATAVIZ_ARS),
    'ytick':   compute_source_fontsize(TARGET_TICK_LABEL,                       _BOXPLOT_H, _FIG_W, _GAP, _ROW_DATAVIZ_ARS),
    'xtick':   compute_source_fontsize(TARGET_TICK_LABEL,                       _BOXPLOT_H, _FIG_W, _GAP, _ROW_DATAVIZ_ARS),
    'sig':     compute_source_fontsize(TARGET_SIG_MARKER,                       _BOXPLOT_H, _FIG_W, _GAP, _ROW_DATAVIZ_ARS),
}

FONT_CORR_3X1 = {
    'ylab':    compute_source_fontsize(TARGET_AXIS_LABEL,                       _CORR_H, _FIG_W, _GAP, _ROW_DATAVIZ_ARS),
    'xlab':    compute_source_fontsize(TARGET_XLABEL * _DATAVIZ_XLABEL_BOOST,   _CORR_H, _FIG_W, _GAP, _ROW_DATAVIZ_ARS),
    'tick':    compute_source_fontsize(TARGET_TICK_LABEL,                        _CORR_H, _FIG_W, _GAP, _ROW_DATAVIZ_ARS),
    'sig':     compute_source_fontsize(TARGET_SIG_MARKER,                        _CORR_H, _FIG_W, _GAP, _ROW_DATAVIZ_ARS),
}

# Row context: 2-panel forest row, Figs 1/2/4/5
_ROW_FOREST_ARS = [_FOREST_AR, _FOREST_AR]

FONT_FOREST = {
    'label':   compute_source_fontsize(TARGET_AXIS_LABEL, _FOREST_H, _FIG_W, _GAP, _ROW_FOREST_ARS),
    'tick':    compute_source_fontsize(TARGET_TICK_LABEL, _FOREST_H, _FIG_W, _GAP, _ROW_FOREST_ARS),
    'xlabel':  compute_source_fontsize(TARGET_XLABEL,     _FOREST_H, _FIG_W, _GAP, _ROW_FOREST_ARS),
    'sig':     compute_source_fontsize(TARGET_SIG_MARKER, _FOREST_H, _FIG_W, _GAP, _ROW_FOREST_ARS),
}

# Row context: 2-panel mediation row, Figs 2/4
_ROW_MED_ARS = [_MED_AR, _MED_AR]

FONT_MEDIATION = {
    'box':     compute_source_fontsize(TARGET_AXIS_LABEL, _MED_H, _FIG_W, _GAP, _ROW_MED_ARS),
    'stat':    compute_source_fontsize(TARGET_TICK_LABEL, _MED_H, _FIG_W, _GAP, _ROW_MED_ARS),
}

# ── Figure 3 (PIL assembly, dpi=200) ──────────────────────────────────────
# Right-side panels are generated at dpi=200 and pasted at 1300 px wide
# (scale factor 1.0 since 6.5" × 200 = 1300).  The composite is saved with
# dpi metadata = FIGURE_DPI (600).  Source fonts rendered at 200 DPI appear
# at 200/600 = 1/3 of their physical size when printed at 600 DPI.
# Multiply target by output_dpi/source_dpi = 3 to compensate.
_FIG3_SRC_DPI = 200
_FIG3_DPI_RATIO = FIGURE_DPI / _FIG3_SRC_DPI    # 3.0

_FIG3_BOXPLOT_H  = 5.0   # per-parameter boxplots: figsize=(6.5, 5)
_FIG3_CORR_H     = 4.5   # per-parameter correlations: figsize=(6.5, 4.5)
_FIG3_TRAJ_H     = 4.5   # trajectories: figsize=(6.5, 4.5)

FONT_FIG3_BOXPLOT = {
    'ylabel':  TARGET_AXIS_LABEL * _FIG3_DPI_RATIO,
    'xlabel':  TARGET_XLABEL     * _FIG3_DPI_RATIO * _DATAVIZ_XLABEL_BOOST,
    'ytick':   TARGET_TICK_LABEL * _FIG3_DPI_RATIO,
    'xtick':   TARGET_TICK_LABEL * _FIG3_DPI_RATIO,
    'sig':     TARGET_SIG_MARKER * _FIG3_DPI_RATIO,
}

FONT_FIG3_CORR = {
    'ylab':    TARGET_AXIS_LABEL * _FIG3_DPI_RATIO,
    'xlab':    TARGET_XLABEL     * _FIG3_DPI_RATIO * _DATAVIZ_XLABEL_BOOST,
    'tick':    TARGET_TICK_LABEL * _FIG3_DPI_RATIO,
    'sig':     TARGET_SIG_MARKER * _FIG3_DPI_RATIO,
}

FONT_FIG3_TRAJECTORY = {
    'label':   TARGET_AXIS_LABEL * _FIG3_DPI_RATIO,
    'tick':    TARGET_TICK_LABEL * _FIG3_DPI_RATIO,
    'legend':  TARGET_TICK_LABEL * _FIG3_DPI_RATIO,
}

# ── Figure 6 (2×2 quadrant: detection curves + correlation | forest plots)
# Row 1: panel a = detection curves (7, 5), panel b = VCH×SDT correlation (~8, 9)
# Row 2: panel c = HPPD forest, panel d = CAPS forest
_FIG6_DETECTION_AR = 7.0 / 5.0                  # 1.4
_FIG6_CORR_SDT_AR  = _CORR_AR                   # 8/9 ≈ 0.889
_FIG6_ROW1_ARS = [_FIG6_DETECTION_AR, _FIG6_CORR_SDT_AR]

# Detection curves panel a: figsize=(7, 5)
_FIG6_DET_H = 5.0
FONT_FIG6_DETECTION = {
    'xlabel':  compute_source_fontsize(TARGET_XLABEL,     _FIG6_DET_H, _FIG_W, _GAP, _FIG6_ROW1_ARS),
    'ylabel':  compute_source_fontsize(TARGET_AXIS_LABEL, _FIG6_DET_H, _FIG_W, _GAP, _FIG6_ROW1_ARS),
    'tick':    compute_source_fontsize(TARGET_TICK_LABEL, _FIG6_DET_H, _FIG_W, _GAP, _FIG6_ROW1_ARS),
    'legend':  compute_source_fontsize(TARGET_TICK_LABEL, _FIG6_DET_H, _FIG_W, _GAP, _FIG6_ROW1_ARS),
}

# VCH×SDT correlation panel b: figsize=(8, 9)
FONT_FIG6_CORR_SDT = {
    'ylab':    compute_source_fontsize(TARGET_AXIS_LABEL, _CORR_H, _FIG_W, _GAP, _FIG6_ROW1_ARS),
    'xlab':    compute_source_fontsize(TARGET_XLABEL,     _CORR_H, _FIG_W, _GAP, _FIG6_ROW1_ARS),
    'tick':    compute_source_fontsize(TARGET_TICK_LABEL, _CORR_H, _FIG_W, _GAP, _FIG6_ROW1_ARS),
    'sig':     compute_source_fontsize(TARGET_SIG_MARKER, _CORR_H, _FIG_W, _GAP, _FIG6_ROW1_ARS),
}

# Figure 6 Row 2 = same forest panels as Figs 1/2 → use FONT_FOREST.

# ── Panel label settings for assembled manuscript figures ─────────────────
# FIGURE_LABEL_X < 0 places labels to the LEFT of each panel (outside image).
# FIGURE_LABEL_MARGIN reserves that fraction of figure width as a left gutter.
FIGURE_LABEL_FONTSIZE = 20     # pt in the assembly figure's coordinate system
FIGURE_LABEL_X        = -0.04  # axes coords; negative = outside panel to the left
FIGURE_LABEL_MARGIN   = 0.05   # fraction of figure width reserved left of panels
# PIL panel labels: convert FIGURE_LABEL_FONTSIZE to pixels at the output DPI
PIL_PANEL_LABEL_PX    = round(FIGURE_LABEL_FONTSIZE / 72 * FIGURE_DPI)

# Pre-built keyword dicts for correlation_matrix_plot and boxplot calls.
# Defined at module level so every conditional block can use them.
_corr_font_kw = dict(
    sig_fontsize=round(FONT_CORR_3X1['sig']),
    xlab_fontsize=round(FONT_CORR_3X1['xlab']),
    ylab_fontsize=round(FONT_CORR_3X1['ylab']),
    tick_size=round(FONT_CORR_3X1['tick']),
)


# ---------------------------------------------------------------------------
# Display labels for HPPD history groups (used across figures and tables)
# ---------------------------------------------------------------------------
HPPD_TERM = "PPA"
HPPD_NEG_LABEL = f"{HPPD_TERM} (-)"
HPPD_POS_LABEL = f"{HPPD_TERM} (+)"
HPPD_STATUS_LABEL = f"{HPPD_TERM} History"
HPPD_RISK_LABEL = f"{HPPD_TERM} Risk"

# ---------------------------------------------------------------------------
# DV configurations for the HPPD & CAPS figures section.
#
# Each entry drives one full pass through IVTYPES: data-visualization figure
# (boxplot or correlation grid) + standardized forest plots.
#
# Keys:
#   dv              — column name; also used as the dv filter in the forest CSV
#                     and as the subdirectory name under results/.  Must match
#                     an entry in outcome_bundles in modules/master_config.py.
#   figure_type     — "boxplot"     → multipanel_boxplot_grid (binary outcomes)
#                     "correlation" → correlation_matrix_plot (continuous outcomes)
#   results_subdir  — folder name under RESULTS_BASE (e.g. "hppd_binary")
#   data            — which DataFrame to pass to the plotting helper:
#                       "sp_plot" → df_sp (SP users, psycheduse_yn == "Yes")
#                       "full"    → df         (all QC-passing subjects)
#   forest_models   — list of covariate strings; each produces one forest-plot
#                     PNG.  Must match the "covariates" column in
#                     existingresults_manuscript.csv for this dv.
#   counterfactual_x_label — forest-plot x-axis label.
#                     Describes the response-scale effect unit for this DV.
#
# Boxplot-only keys:
#   group_var           — column used to split boxes (must be in data)
#   group_order         — display order of group labels
#
# To add a new DV:
#   1. Add its entry to outcome_bundles in modules/master_config.py.
#   2. Append a dict here.
#   3. Confirm its dv rows exist in existingresults_manuscript.csv.
# ---------------------------------------------------------------------------
DV_CONFIGS = [
    # ── Binary outcomes: boxplot grid, split by group_var ──────────────────
    {
        "dv":                      "hppd_binary",
        "figure_type":             "boxplot",
        "results_subdir":          "hppd_binary",
        "data":                    "sp_plot",          # SP users (psycheduse_yn == "Yes")
        "forest_models":           FOREST_MODELS_HPPD, # nice_covariates_spusers (+ any extras)
        "group_var":               "hppd_binary",
        "group_order":             [HPPD_NEG_LABEL, HPPD_POS_LABEL],
        # Forest x-axis: risk difference (probability scale, bernoulli family)
        "counterfactual_x_label":  "Δ in Prior SP-related PPA Risk",
    },
    # ── Continuous outcomes: correlation grid (IV rows × DV column) ─────────
    {
        "dv":                      "caps_vision",
        "figure_type":             "correlation",
        "results_subdir":          "caps_vision",
        # SP users only (n≈186), matching the mediation analysis sample.
        "data":                    "sp_plot",
        # Forest model follows the same sample: spusers.
        "forest_models":           FOREST_MODELS_HPPD,
        # Forest x-axis: expected count difference (hurdle_negbinom_huvary family)
        "counterfactual_x_label":  "Δ in CAPS Visual Items Endorsed",
    },
]

### ══════════════════════   END CONFIG   ═════════════════════════════════ ###


##############################################################################
### SETUP
##############################################################################
import sys, os
sys.path.insert(0, '../modules')

import shutil, textwrap, warnings
from collections import OrderedDict
from pathlib import Path

import matplotlib
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as sci_stats
import seaborn as sns

warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# GLOBAL FONT: All manuscript figures MUST use Arial.
# Set here, before any plotting.  sns.set_style() resets font.family to
# 'sans-serif', so this is re-enforced after every such call.
# ---------------------------------------------------------------------------
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

### Import helper functions
from visualization_helpers import (
    # Figure helpers
    correlation_matrix_plot,
    counterfactual_forest_plot,
    create_mediation_diagram,
    multipanel_boxplot_grid,
    compute_state_stats,
    plot_state_trajectories,
    # Table helpers
    generate_publication_table_thickonly,
    generate_combined_split_table_thickonly,
    # Utilities
    get_field_label_dict,
    load_most_recent_csv,
)

### Import variable dicts
from master_config import (
    binary_palette,
    caps_vision_palette,
    caps_types,
    dv_to_lab_short,
    electric_blue_palette,
    hppd_variables,
    iv_type_dict,
    severity_vars,
)



### Import assembly helper (compute_source_fontsize already imported above for font constants)
from figure_assembly import assemble_manuscript_figure, assemble_manuscript_figure_svg

### Data loading
# ---------------------------------------------------------------------------
# This repository ships the analysis dataframe already prepared: every derived
# column (hppd_binary, subtle, baggot_total, the SDT/metacognition block, the
# recalculated VCH hit rates, corrected avg_life_dose, ...) is present in the CSV
# below.  There is therefore no data-preparation step to run — read the file and
# work from it.
#
# `df_sp` is the SP-user subsample used by most figures.  `psycheduse_yn == "Yes"`
# is the project's canonical SP filter — not `!= "No"`, and not
# `psycheduse_life_nomic > 0`, both of which give different Ns.
# ---------------------------------------------------------------------------
from master_config import SP_USER_COL, SP_USER_VALUE
from master_config import POINT_ESTIMATE_COL, point_estimate
from data_prep import most_recent_public_df

DATA_DIR  = Path(__file__).resolve().parent.parent / "data" / "final"
DF_PUBLIC = most_recent_public_df(DATA_DIR)

df = pd.read_csv(DF_PUBLIC, low_memory=False)
df_sp = df[df[SP_USER_COL] == SP_USER_VALUE].copy()
print(f"Loaded {DF_PUBLIC.name}: df={df.shape}, df_sp={df_sp.shape}")

# Keep helper-level label dictionaries synchronized with script-level overrides.
try:
    multipanel_boxplot_grid.__globals__["dv_to_lab_short"] = dv_to_lab_short
except Exception:
    pass
try:
    create_mediation_diagram.__globals__["dv_to_lab_short"] = dv_to_lab_short
except Exception:
    pass

# Keep hppd_binary display labels synchronized with the PPA/PVD toggle.
_hppd_binary_remap = {
    0: HPPD_NEG_LABEL,
    1: HPPD_POS_LABEL,
    "0": HPPD_NEG_LABEL,
    "1": HPPD_POS_LABEL,
    "PPA (-)": HPPD_NEG_LABEL,
    "PPA (+)": HPPD_POS_LABEL,
    "PVD (-)": HPPD_NEG_LABEL,
    "PVD (+)": HPPD_POS_LABEL,
}
for _frame in (df, df_sp):
    if "hppd_binary" in _frame.columns:
        _frame["hppd_binary"] = _frame["hppd_binary"].replace(_hppd_binary_remap)

dv_to_lab_short["hppd_sx_count"] = "# of HPPD Experiences"


# ---------------------------------------------------------------------------
# Utility: save a copy in FIGURE_FORMAT alongside an already-saved PNG.
# Called after every constituent plot save so all panel figures exist as both
# PNG (needed for figure assembly) and FIGURE_FORMAT.
# No-op when FIGURE_FORMAT == "png".
# ---------------------------------------------------------------------------
from PIL import Image as _PIL_Image
_PIL_Image.MAX_IMAGE_PIXELS = None  # suppress DecompressionBombError for high-DPI PNGs

def _also_save_as(png_path_str):
    """Copy png_path to the same stem with FIGURE_FORMAT extension."""
    if FIGURE_FORMAT == "png":
        return
    p = Path(str(png_path_str))
    if not p.exists():
        p = Path(str(png_path_str) + '.png')   # savepath passed without extension
    if not p.exists():
        return
    fmt = FIGURE_FORMAT.lower().lstrip('.')
    dest = p.with_suffix('.' + fmt)
    with _PIL_Image.open(p) as _img:
        if fmt in {"jpg", "jpeg"}:
            # JPEG does not support alpha channels; convert first when needed.
            if _img.mode not in {"RGB", "L"}:
                _img = _img.convert("RGB")
            _img.save(str(dest), format="JPEG", quality=95, subsampling=0, optimize=True, dpi=(FIGURE_DPI, FIGURE_DPI))
        else:
            _img.save(str(dest), dpi=(FIGURE_DPI, FIGURE_DPI))

print("Setup complete.")




##############################################################################
### HPPD & CAPS FIGURES
# Driven by DV_CONFIGS (defined in the CONFIG section above).
# For each DV × IV-type group combination:
#   boxplot DVs:     multipanel_boxplot_grid + optional nocurrenthppd variant
#   correlation DVs: correlation_matrix_plot
#   both:            standardized forest plots
#
# To add a new DV: edit DV_CONFIGS in the CONFIG section — do not touch here.
##############################################################################
if RUN_HPPD_CAPS_FIGS:

    # Load forest-plot results once; reused across all DV × ivtype combinations.
    # Any persist_vis_yn rows in older compiled CSVs are aliased → hppd_binary
    # so they match the results_subdir. New runs already output hppd_binary.
    results_forest = pd.read_csv(COUNTERFACTUAL_CSV)
    # Rename counterfactual CSV columns to the convention expected by both
    # forest plot functions: Estimate, var, covariates, model.
    # The point-estimate source is POINT_ESTIMATE_COL ('mean') — the same
    # summary of the posterior that brms' own "Estimate" is in the
    # single-path coefficient CSVs, so forest plots built from either
    # source now show the same quantity. point_estimate() is called purely
    # for its presence check: rename() would silently no-op on a CSV that
    # predates the mean columns, leaving a missing 'Estimate' to surface
    # much later as an unrelated-looking error.
    point_estimate(results_forest, source=str(COUNTERFACTUAL_CSV))
    results_forest = results_forest.rename(columns={
        POINT_ESTIMATE_COL: 'Estimate',
        'spvar':    'var',
        'cov':      'covariates',
        'settings': 'model',
    })

    results_forest = results_forest[results_forest['var'] != 'Intercept'].copy()
    results_forest['prob'] = results_forest[['prob_above_0', 'prob_below_0']].max(axis=1)
    results_forest.loc[results_forest['dv'] == 'persist_vis_yn', 'dv'] = 'hppd_binary'
    # A (dv, var, covariates, contrast_label) key must identify exactly one
    # model.  Duplicates mean the compiled CSV holds two results for the same
    # specification -- usually a persist_vis_yn row and an hppd_binary row that
    # survived the rename above, or two compile runs merged without dedup.
    # Silently keeping one would pick an estimate arbitrarily, so stop instead.
    _key = ['dv', 'var', 'covariates', 'contrast_label']
    _dups = results_forest[results_forest.duplicated(subset=_key, keep=False)]
    if not _dups.empty:
        _msg = [
            f"Duplicate forest-plot rows in {COUNTERFACTUAL_CSV}",
            f"  {len(_dups)} rows across {_dups.groupby(_key, dropna=False).ngroups} "
            f"duplicated key(s); key = {_key}",
        ]
        for _k, _grp in _dups.groupby(_key, dropna=False):
            _msg.append(f"\n  key={dict(zip(_key, _k))}")
            _msg.append(f"  row indices: {list(_grp.index)}")
            for _idx, _row in _grp.iterrows():
                _msg.append(f"    [{_idx}] " + ", ".join(
                    f"{_c}={_row[_c]!r}" for _c in results_forest.columns))
        raise ValueError("\n".join(_msg))

    # Shared keyword arguments for standard forest plots
    _forest_kwargs = dict(
        binary_palette=binary_palette,
        lci_col=_FOREST_LCI_COL,
        uci_col=_FOREST_UCI_COL,
        prob_col="prob",
        prob_threshold=0.9,
        alpha_high=1.0,
        alpha_low=0.10,
        estimate_alpha_boost=0.15,
        figsize=(9.5, 4.8),
        label_fontsize=FONT_FOREST['label'],
        tick_fontsize=FONT_FOREST['tick'],
        xlabel_fontsize=FONT_FOREST['xlabel'],
        significance_fontsize=FONT_FOREST['sig'],
        y_label_map=dv_to_lab_short,
        show_plot=False,
    )
    # Shared kwargs for counterfactual forest plots — same visual style but
    # excludes params specific to the standardized forest plot (lci_col,
    # uci_col, prob_col, estimate_alpha_boost are handled internally).
    _forest_kwargs_cf = {
        k: v for k, v in _forest_kwargs.items()
        if k not in {'lci_col', 'uci_col', 'prob_col', 'estimate_alpha_boost'}
    }
    # significance_y_offset controlled by default in counterfactual_forest_plot.py (currently 0.02)

    # Map config "data" strings → live DataFrames (resolved after data loading).
    _DATA_MAP = {
        "full":    df,          # all QC-passing subjects (n ≈ 228)
        "sp_plot": df_sp,  # SP users, psycheduse_yn == "Yes" (n ≈ 186)
        "sp":      df_sp,       # SP users, no special labels
    }

    for cfg in DV_CONFIGS:
        dv           = cfg["dv"]
        subdir       = cfg["results_subdir"]
        figure_type  = cfg["figure_type"]
        data         = _DATA_MAP[cfg["data"]]
        forest_models = cfg["forest_models"]

        print(f"\n--- {dv} figures ---")

        for ivtype in IVTYPES:
            print(f"  {ivtype}")

            # ── Data visualization ─────────────────────────────────────────
            if figure_type == "boxplot":
                panel_specs = [
                    {'dv': var, 'group_var': cfg["group_var"],
                     'order': cfg["group_order"],
                     'palette': binary_palette,
                     'ylabel': dv_to_lab_short[var]}
                    for var in iv_type_dict[ivtype]
                ]
                _boxplot_kwargs = dict(
                    panel_specs=panel_specs,
                    nrows=len(panel_specs), ncols=1,
                    panel_width=6, panel_height=BOXPLOT_PANEL_HEIGHT,
                    show_stats=False,
                    ylabel_fontsize=round(FONT_BOXPLOT_3VAR['ylabel']),
                    xlabel_fontsize=round(FONT_BOXPLOT_3VAR['xlabel']),
                    ytick_fontsize=round(FONT_BOXPLOT_3VAR['ytick']),
                    xtick_fontsize=round(FONT_BOXPLOT_3VAR['xtick']),
                    ytick_font_multiplier=None,   # disabled; using explicit overrides
                    xtick_font_multiplier=None,
                    sig_marker_fontsize=round(FONT_BOXPLOT_3VAR['sig']),
                )
                multipanel_boxplot_grid(
                    data=data, **_boxplot_kwargs,
                    savepath=f"{RESULTS_BASE}/{subdir}/{ivtype}/data_visualization/boxplot_grid",
                )
                _also_save_as(f"{RESULTS_BASE}/{subdir}/{ivtype}/data_visualization/boxplot_grid.png")
            elif figure_type == "correlation":
                _corr_path = f"{RESULTS_BASE}/{subdir}/{ivtype}/data_visualization/correlation_grid.png"
                correlation_matrix_plot(
                    row_cols=iv_type_dict[ivtype],
                    column_cols=[dv],
                    dataframe=data,
                    dv_to_label_dict=dv_to_lab_short,
                    add_x_jitter=0.1,
                    savepath=_corr_path,
                    **_corr_font_kw,
                )
                _also_save_as(_corr_path)

                # ── Age-controlled variants (sp_predictors ivtype only) ────
                # Partial Spearman controlling for age_v2; shows rank-residual
                # scatter and significance from the partial correlation.
                if ivtype == "sp_predictors":
                    _corr_path_ac = (
                        f"{RESULTS_BASE}/{subdir}/{ivtype}/data_visualization/"
                        f"correlation_grid_age_control.png"
                    )
                    correlation_matrix_plot(
                        row_cols=iv_type_dict[ivtype],
                        column_cols=[dv],
                        dataframe=data,
                        dv_to_label_dict=dv_to_lab_short,
                        savepath=_corr_path_ac,
                        control_vars=['age_v2'],
                        **_corr_font_kw,
                    )
                    _also_save_as(_corr_path_ac)

            # ── Per-parameter data viz (vch_computations only) ──────────────
            # Saves one boxplot_grid or correlation_grid per HGF parameter so
            # the figure_3 assembly can place each parameter
            # at the correct level in the HGF diagram layout.
            # Files named: {param}_boxplot_grid.png / {param}_correlation_grid.png
            if ivtype == "vch_computations":
                _param_dir = f"{RESULTS_BASE}/{subdir}/{ivtype}/data_visualization"
                for _pvar in iv_type_dict["vch_computations"]:
                    if figure_type == "boxplot":
                        _pspec = [{'dv': _pvar, 'group_var': cfg["group_var"],
                                   'order': cfg["group_order"],
                                   'palette': binary_palette,
                                   'ylabel': dv_to_lab_short[_pvar]}]
                        # panel_width=6.5" × 200 DPI = 1300px — matches correlation
                        # and trajectory source widths so PIL composite is 1:1 (no
                        # upscaling).  Font sizes from FONT_FIG3_BOXPLOT account for
                        # the 200→600 DPI metadata difference in the composite.
                        multipanel_boxplot_grid(
                            data=data,
                            panel_specs=_pspec,
                            nrows=1, ncols=1,
                            panel_width=6.5, panel_height=BOXPLOT_PANEL_HEIGHT,
                            show_stats=False,
                            ylabel_fontsize=round(FONT_FIG3_BOXPLOT['ylabel']),
                            xlabel_fontsize=round(FONT_FIG3_BOXPLOT['xlabel']),
                            ytick_fontsize=round(FONT_FIG3_BOXPLOT['ytick']),
                            xtick_fontsize=round(FONT_FIG3_BOXPLOT['xtick']),
                            ytick_font_multiplier=None,
                            xtick_font_multiplier=None,
                            sig_marker_fontsize=round(FONT_FIG3_BOXPLOT['sig']),
                            xlabel=(_pvar == 'vch_beta'),
                            dpi=200,
                            savepath=f"{_param_dir}/{_pvar}_boxplot_grid",
                        )
                        _also_save_as(f"{_param_dir}/{_pvar}_boxplot_grid.png")
                    elif figure_type == "correlation":
                        _pcorr = f"{_param_dir}/{_pvar}_correlation_grid.png"
                        # figsize=(6.5, 4.5) at dpi=200 → ~1300×900px, matching
                        # trajectory dimensions so the PIL composite is 1:1.
                        # Font sizes from FONT_FIG3_CORR account for the
                        # 200→600 DPI metadata difference in the composite.
                        # X-axis label only on vch_beta (bottom row in figure_3).
                        # Bottom-row panels (show_xlabel=True) use extra
                        # pad_inches so the large xlabel is not clipped.
                        _is_bottom = (_pvar == 'vch_beta')
                        correlation_matrix_plot(
                            row_cols=[_pvar],
                            column_cols=[dv],
                            dataframe=data,
                            dv_to_label_dict=dv_to_lab_short,
                            add_x_jitter=0.1,
                            savepath=_pcorr,
                            scatter_size=150,
                            line_lw=3,
                            xlab_fontsize=round(FONT_FIG3_CORR['xlab']),
                            ylab_fontsize=round(FONT_FIG3_CORR['ylab']),
                            tick_size=round(FONT_FIG3_CORR['tick']),
                            sig_fontsize=round(FONT_FIG3_CORR['sig']),
                            show_xlabel=_is_bottom,
                            figsize=(6.5, 4.5),
                            dpi=200,
                            pad_inches=0.5 if _is_bottom else 0.1,
                        )
                        _also_save_as(_pcorr)

            # ── Forest plots ───────────────────────────────────────────────
            # Guard: only generate if results CSV has matching rows for this
            # dv × covariates × iv_type combination.  Prevents empty-plot errors
            # when HPC results are absent for a given iv_type (e.g. sdt_hppd
            # predictors that were not included in the original HPC sweep).
            # Same guard is applied uniformly across all IVTYPES — including
            # vch_behavior and sdt_hppd — via the same code path.
            forest_outdir = Path(f"{RESULTS_BASE}/{subdir}/{ivtype}/forest_plots")
            forest_outdir.mkdir(parents=True, exist_ok=True)
            for model in forest_models:
                _rf_rows = results_forest[
                    (results_forest['dv'] == dv) &
                    (results_forest['covariates'] == model) &
                    (results_forest['var'].str.replace("_normalized", "", regex=False)
                                         .isin(iv_type_dict[ivtype]))
                ]
                if _rf_rows.empty:
                    print(f"    [skip] no forest results for {dv}/{ivtype}/{model} — skipping")
                    continue
                counterfactual_forest_plot(
                    results_df=results_forest,
                    dv=dv, covariates=model,
                    x_label=cfg["counterfactual_x_label"],
                    predictors=iv_type_dict[ivtype],
                    savepath=str(forest_outdir / model),
                    **_forest_kwargs_cf,
                )
                _also_save_as(forest_outdir / f"{model}.png")

    print("HPPD & CAPS figures complete.")


##############################################################################
### VCH COMPUTATIONS — STATE TRAJECTORIES
# Generates one block-trajectory figure per state variable per DV.
# Saved to: results/{dv}/vch_computations/trajectories/{var}.png
#
# These modular per-variable PNGs feed directly into the
# figure_3 assembly (HGF diagram layout).  One file per variable makes it
# easy to add new state variables later without disturbing existing files.
#
# Governed by RUN_HPPD_CAPS_FIGS since trajectories are part of the
# vch_computations data-visualization family.
##############################################################################
if RUN_HPPD_CAPS_FIGS:
    print("\n--- VCH state trajectories ---")
    _traj_vars = ['xprob', 'xbin_pred']

    # Pre-compute shared y-axis limits so d↔e and f↔g panels align.
    # We aggregate stats for both DVs first, then take the union of the
    # data ranges (mean ± CI) across both, and pass the shared limits
    # into each plot_state_trajectories call.
    _all_traj_stats = {}
    for _traj_dv in ['hppd_binary', 'caps_vision']:
        _all_traj_stats[_traj_dv], _ = compute_state_stats(
            dv=_traj_dv, df_bl=df, data_dir='../data/final',
            variables=_traj_vars,
        )

    _shared_ylims = {}
    for _var in _traj_vars:
        _vals = []
        for _traj_dv in ['hppd_binary', 'caps_vision']:
            for _centers, _los, _his in _all_traj_stats[_traj_dv].get(_var, {}).values():
                for _arr in (_centers, _los, _his):
                    _vals.extend(_arr[~np.isnan(_arr)].tolist())
        if _vals:
            _ymin, _ymax = min(_vals), max(_vals)
            _pad = (_ymax - _ymin) * 0.08
            _shared_ylims[_var] = (_ymin - _pad, _ymax + _pad)

    for _traj_dv in ['hppd_binary', 'caps_vision']:
        print(f"  {_traj_dv}")
        plot_state_trajectories(
            dv=_traj_dv,
            df_bl=df,
            data_dir='../data/final',
            variables=_traj_vars,
            sp_only=False,
            save_dir=f"{RESULTS_BASE}/{_traj_dv}/vch_computations/trajectories",
            show_plot=False,
            dpi=200,
            # Font sizes account for 200→600 DPI metadata difference in
            # the Figure 3 PIL composite.
            label_fontsize=round(FONT_FIG3_TRAJECTORY['label']),
            tick_fontsize=round(FONT_FIG3_TRAJECTORY['tick']),
            legend_fontsize=round(FONT_FIG3_TRAJECTORY['legend']),
            # X-axis label only on xbin_pred (bottom trajectory row in figure_3)
            xlabel_variables=['xbin_pred'],
            # Shared y-limits so d↔e and f↔g panels use the same scale
            ylim=_shared_ylims,
        )
    print("State trajectories complete.")


##############################################################################
### VCH BETA × SDT_HPPD CORRELATION (figure 6 panel c)
# Generates a correlation grid with vch_beta on the x-axis and each
# iv_type_dict["sdt_hppd"] variable on the y-axis.  Points are colored by
# caps_vision using the canonical caps_vision_palette (discrete 0–6).
#
# Statistics: standard Spearman ρ per panel on non-missing (vch_beta, y, caps_vision)
# triples.  Significance markers: *** p<0.001, ** p<0.01, * p<0.05, ~ p<0.10.
#
# Outputs:
#   results/vch_beta/sdt_hppd/data_visualizations/correlation_grid.png
#   results/vch_beta/sdt_hppd/data_visualizations/summary_results/correlation_grid.csv
##############################################################################
if RUN_HPPD_CAPS_FIGS:
    print("\n--- vch_beta × sdt_hppd correlation (figure 6 panel c) ---")
    _vch_beta_sdt_path = (
        f"{RESULTS_BASE}/vch_beta/sdt_hppd/data_visualizations/correlation_grid.png"
    )
    correlation_matrix_plot(
        row_cols=iv_type_dict["sdt_hppd"],   # criterion_overall, d_prime_overall, mean_conf_fas
        column_cols=["vch_beta"],
        dataframe=df,                          # full QC-passing sample (n ≈ 228)
        palette=caps_vision_palette,
        color_var="caps_vision",               # color points by caps_vision (0–6)
        dv_to_label_dict=dv_to_lab_short,
        add_x_jitter=0.1,
        savepath=_vch_beta_sdt_path,
        sig_fontsize=round(FONT_FIG6_CORR_SDT['sig']),
        xlab_fontsize=round(FONT_FIG6_CORR_SDT['xlab']),
        ylab_fontsize=round(FONT_FIG6_CORR_SDT['ylab']),
        tick_size=round(FONT_FIG6_CORR_SDT['tick']),
        scatter_size=50,
        dpi=FIGURE_DPI,
        ylab_y=0.3,    # shift y-labels down to avoid top-label clipping in assembly
    )
    _also_save_as(_vch_beta_sdt_path)
    print(f"  Saved: {_vch_beta_sdt_path}")
    print("vch_beta × sdt_hppd correlation complete.")


##############################################################################
### TABLES AND DESCRIPTIVE FIGURES
##############################################################################
if RUN_TABLES or RUN_DESCRIPTIVE_FIGS:

    # ── Preprocessing shared by both tables and descriptive figures ────────

    df.loc[df["us_loc_v2"] == 1, "outside_us_v2"] = 0
    df["location_summary"] = "Other"
    df.loc[df["outside_us_v2"].isin([0, 32]), "location_summary"] = "United States"
    df.loc[df["outside_us_v2"].isin([64, 82, 154, 159, 127, 136, 17, 45, 9, 74, 60, 137, 139, 180]), "location_summary"] = "Europe"
    df.loc[df["outside_us_v2"].isin([122, 9]), "location_summary"] = "Oceania"
    df.loc[df["outside_us_v2"].isin([134, 110]), "location_summary"] = "Latin America"

    for spcol in ["psychedelic_primary", "sp_type_recent"]:
        df[spcol] = df[spcol].replace({"other": 6})
        df[spcol] = df[spcol].astype(float).astype("Int64")

    # ── Category label mappings ────────────────────────────────────────────
    race_mapping = {
        1: "American Indian/Alaska Native", 2: "Asian",
        3: "Native Hawaiian or Other Pacific Islander", 4: "Black or African American",
        5: "White", 6: "More than one race", 7: "Unknown/I prefer not to say", 8: "Latino/a",
    }
    sex_mapping = {1: "Male", 2: "Female", 3: "Other"}
    education_mapping = {
        1: "Less than 8th grade", 2: "Some high school", 3: "Highschool diploma or GED",
        4: "Some college", 5: "2-year degree", 6: "College student",
        7: "4-year bachelor's degree", 8: "Master's degree",
        9: "Doctoral or 4 year professional degree",
    }
    location_mapping = {
        "United States": "United States", "Europe": "Europe",
        "Oceania": "Oceania", "Latin America": "Latin America",
        "Other": "Other",
    }
    motivation_mapping = {
        1: "Recreational", 2: "Therapeutic",
        3: "Equally Recreational and Therapeutic",
        # 4: "Not Using" excluded — treated as missing data, not a motivation category
        5: "Spiritual", 6: "Equally Recreational, Therapeutic, and Spiritual",
    }
    sptype_dict = {
        1: "Psilocybin", 2: "LSD", 3: "Mescaline",
        4: "DMT", 5: "5-MeO-DMT", 6: "Other",
    }
    lastuse_dict = {0: "1 week", 1: "1 Month", 2: "6 Months", 3: "1 Year", 4: "> 1 Year"}
    binary_dict = {1: "Yes"}
    howtheyfoundusdict = get_field_label_dict("howtheyfoundus")

    # ── Sample subsets ─────────────────────────────────────────────────────
    df_wraven = df[(df["raven_total"] > 0) & (df["qc_passed"] > 0)].copy()

    month_any_candidates    = ["alc_month_yn", "ghb_month_yn", "opioids_month_yn",
                                "sedatives_month_yn", "mj_month_yn", "atypicals_month_yn", "stimulants_month_yn"]
    lifetime_any_candidates = ["alc_lifetime", "ghb_lifetime", "opioids_lifetime",
                                "sedatives_life_yn", "mj_lifetime", "atypicals_life_yn", "stimulants_life_yn"]
    month_any_cols    = [c for c in month_any_candidates    if c in df_wraven.columns]
    lifetime_any_cols = [c for c in lifetime_any_candidates if c in df_wraven.columns]
    if not month_any_cols:
        raise KeyError("No valid past-month substance columns found to compute any_substance_month_yn")
    if not lifetime_any_cols:
        raise KeyError("No valid lifetime substance columns found to compute any_substance_lifetime_yn")

    df_wraven["any_substance_month_yn"] = (
        df_wraven[month_any_cols].apply(pd.to_numeric, errors="coerce").fillna(0).gt(0).any(axis=1)
    ).astype(int)
    df_wraven["any_substance_lifetime_yn"] = (
        df_wraven[lifetime_any_cols].apply(pd.to_numeric, errors="coerce").fillna(0).gt(0).any(axis=1)
    ).astype(int)

    df_spusers = df_wraven[df_wraven["psycheduse_yn"] == "Yes"].copy()

    # Drop unused levels of 'motivation' (codes with 0 observations in the SP-user sample)
    _motivation_observed = set(df_spusers["motivation"].dropna().unique())
    motivation_mapping = {k: v for k, v in motivation_mapping.items() if k in _motivation_observed}

    df_spusers["sp_lastuse_category"] = df_spusers['psych_dayslastuse_nomicro'].apply(
        lambda x: 0 if x <= 7 else (1 if x <= 31 else (2 if x <= 180 else (3 if x <= 365 else 4)))
    )

    hppd_split_wraven = OrderedDict({
        HPPD_NEG_LABEL: df_wraven.index.isin(df_spusers[df_spusers["persist_vis_yn"] == 0].index),
        HPPD_POS_LABEL: df_wraven.index.isin(df_spusers[df_spusers["persist_vis_yn"] == 1].index),
    })
    hppd_split_spusers = OrderedDict({
        HPPD_NEG_LABEL: df_spusers["persist_vis_yn"] == 0,
        HPPD_POS_LABEL: df_spusers["persist_vis_yn"] == 1,
    })

    # ── Shared table section dicts ─────────────────────────────────────────
    DEMOGRAPHICS_SECTIONS = OrderedDict({
        "Age": "age_v2",
        "Sex": {"sex_v2": sex_mapping},
        "Race": {"race_v2": race_mapping},
        "Location": {"location_summary": location_mapping},
        "Education": {"highest_education": education_mapping},
        "RAVEN's Progressive Matrices Score (9 Item)": "raven_total",
    })

    CLINICAL_SECTIONS = OrderedDict({
        "Psychiatric History": np.nan,
        "Any Diagnosis": {"mental_illness2_v2": binary_dict},
        "Psychotic Spectrum": {"psych_spectrum_v2": binary_dict},
        "Schizophrenia Spectrum": {"schizophrenia_spectrum": binary_dict},
        "Unipolar Mood Disorder": {"mood_disorder": binary_dict},
        "Bipolar Disorder (no psychosis)": {"mental_illness_dx_current_9": binary_dict},
        "Substance Use Disorder": {"addiction": binary_dict},
        "Anxiety Disorder": {"anxiety_disorder": binary_dict},
        "Obsessive-Compulsive Disorder (OCD)": {"mental_illness_dx_current_30": binary_dict},
        "Trauma-Related Disorder (PTSD/c-PTSD)": {"trauma_disorder": binary_dict},
        "Eating Disorder": {"eating_disorder": binary_dict},
        "Personality Disorder": {"personality_disorder": binary_dict},
        "Autism Spectrum Disorder": {"mental_illness_dx_current_6": binary_dict},
        "Attention Deficit Hyperactivity Disorder": {"mental_illness_dx_current_5": binary_dict},
        "Sleep Disorder": {"sleep_disorder": binary_dict},
        "Psychiatric Medications": np.nan,
        "Any Medication": {"medication_current_v2": binary_dict},
        "Antipsychotic": {"antipsychotic": binary_dict},
        "Antidepressant": {"antidepressants": binary_dict},
        "Stimulant": {"simulant_medication": binary_dict},
        "Benzodiazepine": {"benzos": binary_dict},
        "Anxiolytic (Non-Benzodiazepine)": {"nonbenzo_anxiolytics": binary_dict},
        "Sedative (Non-Benzodiazepine)": {"nonbenzo_sedatives": binary_dict},
        "Opioid Antagonist": {"opioid_antagonists": binary_dict},
        "Other Substance Use (Past Month)": np.nan,
        "Any Substance Use (Past Month)": {"any_substance_month_yn": binary_dict},
        "Alcohol Use (Past Month)": {"alc_month_yn": binary_dict},
        "Sedative-Hypnotic Use (Past Month)": {"ghb_month_yn": binary_dict},
        "Opioid Use (Past Month)": {"opioids_month_yn": binary_dict},
        "Cannabis Use (Past Month)": {"mj_month_yn": binary_dict},
        "Atypical Psychedelics Use (Past Month)": {"atypicals_month_yn": binary_dict},
        "Stimulants Use (Past Month)": {"stimulants_month_yn": binary_dict},
        "Other Substance Use History (Lifetime)": np.nan,
        "Any Substance Use (Lifetime)": {"any_substance_lifetime_yn": binary_dict},
        "Alcohol Use (Lifetime)": {"alc_lifetime": binary_dict},
        "Sedative-Hypnotic Use (Lifetime)": {"ghb_lifetime": binary_dict},
        "Opioid Use (Lifetime)": {"opioids_lifetime": binary_dict},
        "Cannabis Use (Lifetime)": {"mj_lifetime": binary_dict},
        "Atypical Psychedelics Use (Lifetime)": {"atypicals_life_yn": binary_dict},
        "Stimulants Use (Lifetime)": {"stimulants_life_yn": binary_dict},
    })

    SP_SECTIONS = OrderedDict({
        "Lifetime Psychedelic Uses (Macrodoses)": "psycheduse_life_nomic",
        "Average SP dose used (LSD ug equivalents)": "avg_life_dose",
        "Age of first Psychedelic Use": "psychedelic_age",
        "Years Of Psychedelic Use": "psyched_yearsofuse",
        "Most Recent Psychedelic Use": {"sp_lastuse_category": lastuse_dict},
        "Most Recent Psychedelic Used": {"sp_type_recent": sptype_dict},
        "Preferred Psychedelic": {"psychedelic_primary": sptype_dict},
        "Reason for Psychedelic Use": {"motivation": motivation_mapping},
        "Perceived Benefit of SP use (VAS)": "perceived_benefit",
        f"History of SP-induced {HPPD_TERM}s": {"persist_vis_yn": {1: "Yes"}},
    })

if RUN_TABLES:
    print("\n--- Publication tables ---")
    Path(TABLES_OUT_DIR).mkdir(parents=True, exist_ok=True)

    generate_publication_table_thickonly(
        dataframe=df_wraven, columns=None, sections=CLINICAL_SECTIONS,
        table_caption="Clinical and Other Substance Use Characteristics", table_label="tab:clinical",
        savepath=f"{TABLES_OUT_DIR}/clinical_table",
        # Not a manuscript table. It exists because results_narrative.py reads
        # clinical_table.csv for the largest single-diagnosis and single-medication
        # N, so the CSV is the only format needed.
        save_formats=("csv",), include_test_stats_csv=False,
    )

    # SP split: history row becomes risk row (it is the split variable)
    _sp_history_key = f"History of SP-induced {HPPD_TERM}s"
    sp_sections_split = OrderedDict(
        {k: v for k, v in SP_SECTIONS.items() if k != _sp_history_key}
        | {HPPD_RISK_LABEL: {"persist_vis_yn": {1: "Yes"}}}
    )
    # ── CAPS-valid subsample ──────────────────────────────────────────────────
    # Participants with null caps_vision cannot be assigned to CAPS(-)/CAPS(+)
    # groups.  Dropped only here; the left-hand HPPD group keeps the full n=186.
    _n_caps_before = len(df_spusers)
    df_spusers_caps_valid = df_spusers[df_spusers["caps_vision"].notna()].copy()
    _n_caps_after  = len(df_spusers_caps_valid)
    print(f"  CAPS-valid subsample: {_n_caps_after} of {_n_caps_before} SP users retained "
          f"({_n_caps_before - _n_caps_after} dropped — null caps_vision).")
    caps_only_split = OrderedDict({
        "CAPS(-)": df_spusers_caps_valid["caps_vision"] == 0,
        "CAPS(+)": df_spusers_caps_valid["caps_vision"] > 0,
    })

    # ── Combined HPPD+CAPS split tables (always generated when RUN_TABLES=True) ─
    # Left group : PPA(-) | PPA(+) | Total | P-value  (df_spusers, n=186)
    # Right group: CAPS(-) | CAPS(+) | Total | P-value  (df_spusers_caps_valid, n≈130)
    # A thick vertical line separates the two groups in the PNG.
    # SP table left group uses sp_sections_split (HPPD_RISK_LABEL row); right
    # group keeps persist_vis_yn to show PPA prevalence within each CAPS group.
    print("  Generating combined HPPD+CAPS split tables...")
    generate_combined_split_table_thickonly(
        dataframe=df_spusers,
        columns_left=hppd_split_spusers,
        dataframe_right=df_spusers_caps_valid,
        columns_right=caps_only_split,
        sections=DEMOGRAPHICS_SECTIONS,
        table_caption=f"Demographic Characteristics by {HPPD_STATUS_LABEL} and by CAPS Vision Status",
        table_label="tab:demographics_hppd_caps_split",
        savepath=f"{TABLES_OUT_DIR}/demographics_table_hppd_split_caps_split",
        save_formats=("png", "csv", "docx"),
        include_total=True,
        continuous_test_left="mannwhitneyu",
        continuous_test_right="mannwhitneyu",
    )
    generate_combined_split_table_thickonly(
        dataframe=df_spusers,
        columns_left=hppd_split_spusers,
        dataframe_right=df_spusers_caps_valid,
        columns_right=caps_only_split,
        sections=CLINICAL_SECTIONS,
        table_caption=f"Clinical and Other Substance Use Characteristics by {HPPD_STATUS_LABEL} and by CAPS Vision Status",
        table_label="tab:clinical_hppd_caps_split",
        savepath=f"{TABLES_OUT_DIR}/clinical_table_hppd_split_caps_split",
        save_formats=("png", "csv", "docx"),
        include_total=True,
        continuous_test_left="mannwhitneyu",
        continuous_test_right="mannwhitneyu",
    )
    # SP combined table: left uses sp_sections_split (persist_vis_yn → HPPD_RISK_LABEL
    # row since it is the left split axis); right uses sp_sections_split too so the row
    # label is consistent, and shows PPA prevalence within each CAPS group.
    generate_combined_split_table_thickonly(
        dataframe=df_spusers,
        columns_left=hppd_split_spusers,
        dataframe_right=df_spusers_caps_valid,
        columns_right=caps_only_split,
        sections=sp_sections_split,
        table_caption=f"Psychedelic Use Characteristics by {HPPD_STATUS_LABEL} and by CAPS Vision Status",
        table_label="tab:sp_hppd_caps_split",
        savepath=f"{TABLES_OUT_DIR}/sp_table_hppd_split_caps_split",
        save_formats=("png", "csv", "docx"),
        include_total=True,
        continuous_test_left="mannwhitneyu",
        continuous_test_right="mannwhitneyu",
    )
    print("  Combined HPPD+CAPS split tables done.")

    # ── PPA History Distributions Figure ─────────────────────────────────────
    _ppa_color = '#6ab0c5'

    _baggot_cols_ppa = [
        "grids", "moving", "different", "oscillating", "halos",
        "still", "trails", "colors", "objects",
        "pattern_open", "pattern_closed", "attention",
    ]
    _baggott_map = {
        "halos": "Halos or auras",
        "moving": "Stationary objects appear to move",
        "still": "Moving objects appear to not move",
        "trails": "Afterimages for moving objects",
        "colors": "Greater color intensity",
        "pattern_open": "Patterns w/eyes open",
        "pattern_closed": "Patterns w/eyes closed",
        "objects": "Objects that aren't really there",
        "oscillating": "Accentuation of light oscillations",
        "grids": "Distorted/moving lines/grids",
        "attention": "Noticing more things in environment",
        "different": "Things just look different",
    }
    _baggot_xlabels_ppa = [_baggott_map[col] for col in _baggot_cols_ppa]
    _chronicity_map_ppa = {
        2: "< 1 day",
        3: "1 - 3 days",
        4: "3 days – 1 week",
        5: "1 wk – 1 month",
        6: "1 – 6 months",
        7: "6 mos – 1 year",
        8: "> 1 year",
    }
    _duration_map_ppa = {
        1: "Brief spurts\n(sec–min)",
        2: "Longer periods\n(hrs–days)",
        3: "Constant/\nnear-constant",
    }
    # Reordered to match _baggot_cols_ppa order so panels C and F share aligned xticks.
    # Key 13 ("never experienced") removed — always 0 observations.
    _most_map_ppa = {
        10: "Distorted/moving lines/grids",
        2:  "Stationary objects appear to move",
        12: "Things just look different",
        9:  "Accentuation of light oscillations",
        1:  "Halos or auras",
        3:  "Moving objects appear to not move",
        4:  "Afterimages for moving objects",
        5:  "Greater color intensity",
        8:  "Objects that aren't really there",
        6:  "Patterns w/eyes open",
        7:  "Patterns w/eyes closed",
        11: "Noticing more things in environment",
    }
    
    def plot_ppa_history_distributions(df_ppa, save_path):
        """Half-violin + barplot summary figure for PPA/HPPD history characteristics."""

        def _pviolin(ax, data, color, alpha=0.78):
            """Half-violin scaled so y-axis = expected count per unit x."""
            data = np.asarray(data, float); data = data[np.isfinite(data)]
            n = len(data)
            if n < 5: return None, None, None
            bw = max(0.4, sci_stats.gaussian_kde(data).factor)
            kde = sci_stats.gaussian_kde(data, bw_method=bw)
            lo, hi = np.percentile(data, [0.5, 99.5])
            pad = (hi - lo) * 0.12
            x_range = (lo - pad, hi + pad)
            xs = np.linspace(*x_range, 400)
            ys_raw = kde(xs)
            # Scale by n so y-axis represents count density
            ys = ys_raw * n
            ax.fill_between(xs, np.zeros_like(ys), ys, color=color, alpha=alpha, zorder=2)
            return kde, n, x_range

        def _pqmarks(ax, data, kde, n_total, color, fs=13.0):
            data = np.asarray(data, float); data = data[np.isfinite(data)]
            q25, q50, q75 = np.percentile(data, [25, 50, 75])
            r, g, b, *_ = mcolors.to_rgba(color)
            dark = (r * 0.6, g * 0.6, b * 0.6)
            peak = float(kde(np.linspace(q25, q75, 200)).max()) * n_total
            for q, ls, lw, is_med in [
                (q25, "--", 1.0, False), (q50, "-", 1.8, True), (q75, "--", 1.0, False)
            ]:
                h = float(kde(np.array([q]))[0]) * n_total
                ax.plot([q, q], [0, h], color=dark, ls=ls, lw=lw, alpha=0.9, zorder=5)
                y_txt = h + peak * (0.25 if is_med else 0.10)
                ax.text(q, y_txt, f"{round(q)}", ha="center", va="bottom",
                        fontsize=fs, color="black",
                        fontweight="bold" if is_med else "normal", clip_on=False)

        def _pviolin_style(ax, x_range, label, peak_count=None):
            ax.set_xlim(*x_range)
            if peak_count is not None:
                # Show a few integer y-ticks
                _step = max(1, round(peak_count / 4))
                _yticks = np.arange(0, peak_count * 1.4 + _step, _step)
                ax.set_yticks([int(t) for t in _yticks if t <= peak_count * 1.3])
                ax.tick_params(axis='y', labelsize=9)
                ax.set_ylabel("Count", fontsize=11, labelpad=4)
            else:
                ax.set_yticks([])
            ax.spines["left"].set_visible(False)
            ax.set_xlabel(label, fontsize=18, labelpad=4)
            ax.spines["bottom"].set_visible(False)
            ax.axhline(0, color="#999999", lw=0.9, zorder=1)
            ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
            ax.tick_params(length=0)
            ax.xaxis.grid(True, alpha=0.18, zorder=0); ax.set_axisbelow(True)

        def _pbarplot(ax, series, mapping, color, fs=9.5, rotate=False):
            codes = list(mapping.keys()); labels = list(mapping.values())
            series = pd.to_numeric(series, errors="coerce").dropna()
            counts = np.array([int((series == c).sum()) for c in codes])
            total = counts.sum() or 1; top = counts.max() or 1
            xs = np.arange(len(codes))
            bars = ax.bar(xs, counts, color=color, alpha=0.82, width=0.62,
                          edgecolor="white", linewidth=0.7, zorder=3)
            for bar, cnt in zip(bars, counts):
                if cnt:
                    xc = bar.get_x() + bar.get_width() / 2
                    y_pct = bar.get_height() + top * 0.04
                    ax.text(xc, y_pct, f"({cnt/total*100:.0f}%)",
                            ha="center", va="bottom", fontsize=fs * 1.1, color="#333",
                            clip_on=False, fontweight="bold")
                    ax.text(xc, y_pct + top * 0.07, f"{cnt}",
                            ha="center", va="bottom", fontsize=fs * 1.1, color="#333",
                            clip_on=False)
            ax.set_xticks(xs)
            if rotate:
                ax.set_xticklabels(labels, ha="right", rotation=40, fontsize=fs * 1.05)
            else:
                ax.set_xticklabels(labels, ha="center", fontsize=fs * 1.1)
            ax.set_ylabel("Count", fontsize=fs)
            ax.set_ylim(0, top * 1.55)
            ax.yaxis.grid(True, alpha=0.22, zorder=0); ax.set_axisbelow(True)
            ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False); ax.spines["bottom"].set_visible(False)
            ax.tick_params(length=0)

        sns.set_style("white")
        # Re-enforce Arial after sns.set_style (which resets font.family)
        matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
        fig = plt.figure(figsize=(18.5, 13.0))
        gs = gridspec.GridSpec(2, 3, figure=fig, left=0.07, right=0.97,
                               top=0.95, bottom=0.13, wspace=0.42, hspace=0.70,
                               width_ratios=[1.92, 1, 1.92])
        # Column-major (top-down) ordering: read down each column before moving right.
        # Physical order (row-scan): (r0c0)=a, (r0c1)=c, (r0c2)=e, (r1c0)=b, (r1c1)=d, (r1c2)=f
        panel_iter = iter("acebdf")
        lbl_kw = dict(fontsize=36, fontweight="bold", fontfamily="Arial",
                      va="bottom", ha="left", clip_on=False)

        # A (row 0, col 0) — persistvis_psychdoses violin
        ax = fig.add_subplot(gs[0, 0])
        ax.text(-0.12, 1.0, f"{next(panel_iter)}", transform=ax.transAxes, **lbl_kw)
        d = pd.to_numeric(df_ppa["persistvis_psychdoses"], errors="coerce").dropna()
        kde, n_total, xr = _pviolin(ax, d, _ppa_color)
        if xr:
            _pqmarks(ax, d, kde, n_total, _ppa_color)
            _peak_a = float(kde(np.linspace(*xr, 400)).max()) * n_total
            _pviolin_style(ax, xr, "SP Doses Before\nFirst PPA Onset",
                           peak_count=_peak_a)
            ax.set_ylim(-_peak_a * 0.08, _peak_a * 1.5)

        # B (row 0, col 1) — hppd_true_chronicity barplot  [old C]
        ax = fig.add_subplot(gs[0, 1])
        ax.text(-0.12, 1.0, f"{next(panel_iter)}", transform=ax.transAxes, **lbl_kw)
        _pbarplot(ax, df_ppa["hppd_true_chronicity"], _chronicity_map_ppa, _ppa_color, rotate=True)
        ax.set_xlabel("Total Duration of PPA Experience", fontsize=18, labelpad=6, color="#333")

        # C (row 0, col 2) — baggot symptom prevalence countplot
        ax_c = fig.add_subplot(gs[0, 2])
        ax_c.text(-0.12, 1.0, f"{next(panel_iter)}", transform=ax_c.transAxes, **lbl_kw)
        n_ppa = len(df_ppa)
        counts_c = np.array([
            int((pd.to_numeric(df_ppa[col], errors="coerce") > 0).sum())
            for col in _baggot_cols_ppa
        ])
        top_c = counts_c.max() or 1
        xs_c = np.arange(len(_baggot_cols_ppa))
        bars_c = ax_c.bar(xs_c, counts_c, color=_ppa_color, alpha=0.82, width=0.62,
                           edgecolor="white", linewidth=0.7, zorder=3)
        for bar, cnt in zip(bars_c, counts_c):
            xc = bar.get_x() + bar.get_width() / 2
            y_pct = bar.get_height() + top_c * 0.04
            ax_c.text(xc, y_pct, f"({cnt/n_ppa*100:.0f}%)",
                      ha="center", va="bottom", fontsize=9.5 * 1.1, color="#333",
                      clip_on=False, fontweight="bold")
            ax_c.text(xc, y_pct + top_c * 0.07, f"{cnt}",
                      ha="center", va="bottom", fontsize=9.5 * 1.1, color="#333",
                      clip_on=False)
        ax_c.set_xticks(xs_c)
        ax_c.set_xticklabels(_baggot_xlabels_ppa, ha="right", rotation=40, fontsize=9.5 * 1.05)
        ax_c.set_ylabel("Count", fontsize=9.5)
        ax_c.set_ylim(0, top_c * 1.55)
        ax_c.yaxis.grid(True, alpha=0.22, zorder=0); ax_c.set_axisbelow(True)
        ax_c.spines["top"].set_visible(False); ax_c.spines["right"].set_visible(False)
        ax_c.spines["left"].set_visible(False); ax_c.spines["bottom"].set_visible(False)
        ax_c.tick_params(length=0)
        ax_c.set_xlabel("PPA Symptom Endorsed", fontsize=18, labelpad=6, color="#333")

        # D (row 1, col 0) — baggot_total violin
        ax = fig.add_subplot(gs[1, 0])
        ax.text(-0.12, 1.0, f"{next(panel_iter)}", transform=ax.transAxes, **lbl_kw)
        d = pd.to_numeric(df_ppa["baggot_total"], errors="coerce").dropna()
        kde, n_total, xr = _pviolin(ax, d, _ppa_color)
        if xr:
            _pqmarks(ax, d, kde, n_total, _ppa_color)
            _peak_b = float(kde(np.linspace(*xr, 400)).max()) * n_total
            _pviolin_style(ax, xr, "# PPA Symptoms\nEndorsed",
                           peak_count=_peak_b)
            ax.set_ylim(-_peak_b * 0.08, _peak_b * 1.5)

        # E (row 1, col 1) — persistvis_duration barplot
        ax = fig.add_subplot(gs[1, 1])
        ax.text(-0.12, 1.0, f"{next(panel_iter)}", transform=ax.transAxes, **lbl_kw)
        _pbarplot(ax, df_ppa["persistvis_duration"], _duration_map_ppa, _ppa_color)
        ax.set_xlabel("PPA Timing Pattern", fontsize=18, labelpad=6, color="#333")

        # F (row 1, col 2) — persistvis_most barplot
        ax = fig.add_subplot(gs[1, 2])
        ax.text(-0.12, 1.0, f"{next(panel_iter)}", transform=ax.transAxes, **lbl_kw)
        _pbarplot(ax, df_ppa["persistvis_most"], _most_map_ppa, _ppa_color, rotate=True)
        ax.set_xlabel("Most Vivid/Intense PPA", fontsize=18, labelpad=6, color="#333")

        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=180, bbox_inches="tight")
        fig.savefig(str(Path(save_path).with_suffix('.svg')), format='svg', bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {save_path}")

    df_ppa_hist = df_spusers[df_spusers["persist_vis_yn"] == 1].copy()
    plot_ppa_history_distributions(
        df_ppa_hist,
        Path(TABLES_OUT_DIR).parent / "ppa_history_distributions.png",
    )

    print(f"Tables saved to: {TABLES_OUT_DIR}")


##############################################################################
### DESCRIPTIVE DISTRIBUTION FIGURES (SP use characteristics)
# Half-violin + countplot summary figures for the SP use table.
# Unsplit version (all SP users) and HPPD-split version (NEG vs POS).
##############################################################################
if RUN_DESCRIPTIVE_FIGS:
    print("\n--- Descriptive distribution figures ---")

    sns.set_style("white")
    # Re-enforce Arial after sns.set_style (which resets font.family)
    plt.rcParams.update({
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 11,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    PPA_NEG  = binary_palette[0]
    PPA_POS  = binary_palette[1]
    COL_NEUT = '#6ab0c5'

    FIG_DIR = Path(TABLES_OUT_DIR).parent
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # (column, y_label, kind, mapping_or_None)  kind: 'cont' | 'cat'
    SECTIONS = [
        ("psycheduse_life_nomic", "Lifetime SP Uses\n(Macrodoses)",       "cont", None),
        ("avg_life_dose",         "Avg. Dose Used\n(LSD μg eq.)",         "cont", None),
        ("psychedelic_age",       "Age at First\nPsychedelic Use",        "cont", None),
        ("psyched_yearsofuse",    "Years of\nPsychedelic Use",            "cont", None),
        ("perceived_benefit",     "Perceived Benefit\nof SP Use (VAS)",   "cont", None),
        ("sp_lastuse_category",   "Most Recent\nPsychedelic Use",         "cat",  lastuse_dict),
        ("sp_type_recent",        "Most Recent\nSP Used",                 "cat",  sptype_dict),
        ("psychedelic_primary",   "Preferred\nPsychedelic",               "cat",  sptype_dict),
        ("motivation",            "Reason for\nSP Use",                   "cat",  motivation_mapping),
        ("persist_vis_yn",        f"SP-induced {HPPD_TERM}s",             "cat",  {0: HPPD_NEG_LABEL, 1: HPPD_POS_LABEL}),
    ]
    # For HPPD-split: persist_vis_yn IS the grouping variable — omit from the plot
    SECTIONS_SPLIT = [s for s in SECTIONS if s[0] != "persist_vis_yn"]

    # ── Low-level drawing helpers ──────────────────────────────────────────

    def _clean(ax):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(length=0)

    def _darken(color, factor=0.60):
        r, g, b, *_ = mcolors.to_rgba(color)
        return (r * factor, g * factor, b * factor)

    def _half_violin(ax, data, color, flip=False, alpha=0.78, x_range=None,
                      proportion=False):
        """Draw a half-violin (KDE fill).

        When *proportion* is True the KDE is scaled so that heights represent
        the proportion of participants per unit-x (i.e. the PDF, which
        integrates to 1 over the data range).  This lets the y-axis carry
        meaningful tick labels ("Proportion of Participants").

        Returns *(kde, kde_max, x_range)*.  *kde_max* is the peak of the raw
        KDE when proportion=False, or the peak of the proportion-scaled KDE
        when proportion=True (needed by _quartile_marks).
        """
        data = np.asarray(data, float)
        data = data[np.isfinite(data)]
        if len(data) < 5:
            return None, None, x_range
        bw  = max(0.4, sci_stats.gaussian_kde(data).factor)
        kde = sci_stats.gaussian_kde(data, bw_method=bw)
        if x_range is None:
            lo, hi = np.percentile(data, [0.5, 99.5])
            pad = (hi - lo) * 0.10
            x_range = (lo - pad, hi + pad)
        xs     = np.linspace(*x_range, 400)
        ys_raw = kde(xs)
        if proportion:
            # Scale so area under curve ≈ 1 (kde is already a PDF);
            # display raw density as-is — y values = proportion per unit x.
            ys = ys_raw
            kde_max_out = float(ys_raw.max())
        else:
            ys = ys_raw / float(ys_raw.max()) * 0.46
            kde_max_out = float(ys_raw.max())
        y_lo   = -ys if flip else np.zeros_like(ys)
        y_hi   = np.zeros_like(ys) if flip else ys
        ax.fill_between(xs, y_lo, y_hi, color=color, alpha=alpha, zorder=2)
        return kde, kde_max_out, x_range

    def _quartile_marks(ax, data, kde=None, kde_max=None, flip=False, color="#1a1a1a", fs=14.0,
                         proportion=False):
        data = np.asarray(data, float)
        data = data[np.isfinite(data)]
        q25, q50, q75 = np.percentile(data, [25, 50, 75])
        sign = -1 if flip else +1
        dark = _darken(color)
        for q, ls, lw, is_median in [(q25, "--", 1.0, False), (q50, "-", 1.8, True), (q75, "--", 1.0, False)]:
            if proportion:
                h = float(kde(np.array([q]))[0]) if (kde is not None) else 0
            else:
                h = float(kde(np.array([q]))[0]) / kde_max * 0.46 if (kde is not None and kde_max) else 0.46
            y_top = sign * h
            ax.plot([q, q], [0, y_top], color=dark, ls=ls, lw=lw, alpha=0.9, zorder=5)
            # Text offset: scale to data range for proportion mode
            if proportion:
                y_offset = sign * kde_max * 0.08 * (6 if is_median else 1)
            else:
                y_offset = sign * 0.07 * (6 if is_median else 1)
            y_txt = y_top + y_offset
            ax.text(q, y_txt, f"{round(q)}", ha="center",
                    va="bottom" if sign > 0 else "top",
                    fontsize=fs, color="black", fontweight="bold" if is_median else "normal",
                    clip_on=False)

    def _style_violin_ax(ax, x_range, label, show_zero_line=True,
                          ylabel=None, yticks=None):
        ax.set_xlim(*x_range)
        if yticks is not None:
            ax.set_yticks(yticks)
            ax.tick_params(axis='y', labelsize=9)
        else:
            ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        if ylabel:
            ax.set_ylabel(ylabel, fontsize=18, labelpad=4)
        ax.set_xlabel(label, fontsize=18, labelpad=4)
        if show_zero_line:
            ax.axhline(0, color="#999999", lw=0.9, zorder=1)
        _clean(ax)
        ax.xaxis.grid(True, alpha=0.18, zorder=0)
        ax.set_axisbelow(True)

    def _wrap(s, width=18):
        return "\n".join(textwrap.wrap(str(s), width))

    def _countplot_single(ax, series, mapping, color=COL_NEUT, fs=9.0):
        cats   = list(mapping.values())
        counts = series.map(mapping).dropna().value_counts().reindex(cats, fill_value=0)
        total  = counts.sum()
        xs     = np.arange(len(cats))
        bars   = ax.bar(xs, counts.values, color=color, alpha=0.82, width=0.60,
                        edgecolor="white", linewidth=0.7, zorder=3)
        top = counts.max() if counts.max() > 0 else 1
        for bar, cnt in zip(bars, counts.values):
            if cnt:
                x = bar.get_x() + bar.get_width() / 2
                y = bar.get_height() + top * 0.028
                fs_plot = (fs - 1.5) * 1.5
                # Trailing \n anchors the empty second line to y (va="bottom"),
                # which pushes the count one line above — matching the original layout.
                ax.text(x, y, f"{cnt}\n",
                        ha="center", va="bottom", fontsize=fs_plot, color="#333")
                ax.text(x, y, f"({cnt/total*100:.0f}%)",
                        ha="center", va="bottom", fontsize=fs_plot, color="#333",
                        fontweight="bold")
        ax.set_xticks(xs)
        ax.set_xticklabels([_wrap(c) for c in cats], rotation=0, ha="center",
                            fontsize=(fs - 1.5) * 1.5)
        ax.set_ylabel("Count", fontsize=fs * 2)
        ax.yaxis.grid(True, alpha=0.22, zorder=0)
        ax.set_axisbelow(True)
        _clean(ax)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)

    def _countplot_split(ax, series_neg, series_pos, mapping, fs=9.0):
        cats  = list(mapping.values())
        c_neg = series_neg.map(mapping).dropna().value_counts().reindex(cats, fill_value=0)
        c_pos = series_pos.map(mapping).dropna().value_counts().reindex(cats, fill_value=0)
        t_neg, t_pos = int(c_neg.sum()), int(c_pos.sum())
        xs, w = np.arange(len(cats)), 0.38
        bars_neg = ax.bar(xs - w/2, c_neg.values, w, color=PPA_NEG, alpha=0.86,
                           edgecolor="white", lw=0.6, label=HPPD_NEG_LABEL, zorder=3)
        bars_pos = ax.bar(xs + w/2, c_pos.values, w, color=PPA_POS, alpha=0.86,
                           edgecolor="white", lw=0.6, label=HPPD_POS_LABEL, zorder=3)
        top = max(c_neg.max(), c_pos.max(), 1)
        for bar, cnt, tot in (list(zip(bars_neg, c_neg.values, [t_neg]*len(cats))) +
                               list(zip(bars_pos, c_pos.values, [t_pos]*len(cats)))):
            if cnt:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + top * 0.018,
                        f"{cnt}\n({cnt/tot*100:.0f}%)",
                        ha="center", va="bottom", fontsize=(fs - 2) * 1.5, color="#333")
        ax.set_xticks(xs)
        ax.set_xticklabels([_wrap(c) for c in cats], rotation=0, ha="center",
                            fontsize=(fs - 2) * 1.5)
        ax.set_ylabel("Count", fontsize=fs)
        ax.yaxis.grid(True, alpha=0.22, zorder=0)
        ax.set_axisbelow(True)
        _clean(ax)

    def _row_heights(sections):
        """Variable row heights: continuous shorter, categorical taller."""
        return [1.9 if k == "cont" else max(2.4, min(3.8, (len(m) if m else 2) * 0.6 + 1.2))
                for _, _, k, m in sections]

    # ── Figure builders ────────────────────────────────────────────────────

    def plot_sp_distributions(df, sections, save_path):
        """Two-column layout: left = continuous half-violins, right = countplots.

        Panel labels use column-major (top-to-bottom) ordering: left column is
        labelled a, b, c, ... and right column continues from where left stopped.
        """
        cont = [(c, l, k, m) for c, l, k, m in sections if k == "cont"]
        cat  = [(c, l, k, m) for c, l, k, m in sections if k == "cat"]
        fig  = plt.figure(figsize=(13.5, 12.0))
        outer    = gridspec.GridSpec(1, 2, figure=fig, left=0.10, right=0.97,
                                     top=0.97, bottom=0.03, wspace=0.32)
        gs_left  = gridspec.GridSpecFromSubplotSpec(len(cont), 1, subplot_spec=outer[0],
                                                    height_ratios=[1]*len(cont), hspace=0.60)
        gs_right = gridspec.GridSpecFromSubplotSpec(len(cat), 1, subplot_spec=outer[1],
                                                    height_ratios=_row_heights(cat), hspace=1.0)
        for i, (col, label, _, _m) in enumerate(cont):
            ax = fig.add_subplot(gs_left[i])
            # Column-major label: left col gets a, b, c, d, e (chr('a' + i))
            ax.text(-0.12, 1.0, f"{chr(ord('a') + i)}", transform=ax.transAxes,
                    fontsize=36, fontweight="bold", fontfamily="Arial",
                    va="bottom", ha="left", clip_on=False)
            kde, kde_max, x_range = _half_violin(ax, df[col].dropna(), COL_NEUT,
                                                  proportion=True)
            if x_range:
                _quartile_marks(ax, df[col].dropna(), kde=kde, kde_max=kde_max,
                                color=COL_NEUT, proportion=True)
                # Exactly 4 y-ticks (including 0) evenly spaced up to ~peak
                _peak = kde_max
                _tick_max = _peak * 1.1
                # Round tick_max up to a clean number
                _mag = 10 ** np.floor(np.log10(_tick_max)) if _tick_max > 0 else 0.01
                _tick_max_nice = np.ceil(_tick_max / _mag) * _mag
                _yticks = np.linspace(0, _tick_max_nice, 4)
                _style_violin_ax(ax, x_range, label, show_zero_line=False,
                                 ylabel="Sample\nProportion",
                                 yticks=[round(t, 3) for t in _yticks])
                ax.set_ylim(-_peak * 0.08, _peak * 1.5)
                ax.spines["bottom"].set_visible(False)
        n_cont = len(cont)
        for j, (col, label, _, mapping) in enumerate(cat):
            ax = fig.add_subplot(gs_right[j])
            # Column-major label: right col continues from left (chr('a' + n_cont + j))
            ax.text(-0.12, 1.0, f"{chr(ord('a') + n_cont + j)}", transform=ax.transAxes,
                    fontsize=36, fontweight="bold", fontfamily="Arial",
                    va="bottom", ha="left", clip_on=False)
            _countplot_single(ax, df[col].dropna(), mapping)
            ax.set_xlabel(label.replace("\n", " "), fontsize=18, labelpad=8, color="#333")
        fig.savefig(save_path, dpi=180, bbox_inches="tight")
        fig.savefig(str(Path(save_path).with_suffix('.svg')), format='svg', bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {save_path}")

    def plot_sp_distributions_split(df, df_neg, df_pos, sections, save_path):
        """Two-column layout: left = mirrored violins, right = grouped countplots."""
        cont = [(c, l, k, m) for c, l, k, m in sections if k == "cont"]
        cat  = [(c, l, k, m) for c, l, k, m in sections if k == "cat"]
        fig  = plt.figure(figsize=(13.0, 13.0))
        outer    = gridspec.GridSpec(1, 2, figure=fig, left=0.08, right=0.97,
                                     top=0.96, bottom=0.03, wspace=0.32)
        gs_left  = gridspec.GridSpecFromSubplotSpec(len(cont), 1, subplot_spec=outer[0],
                                                    height_ratios=[1]*len(cont), hspace=0.65)
        gs_right = gridspec.GridSpecFromSubplotSpec(len(cat), 1, subplot_spec=outer[1],
                                                    height_ratios=_row_heights(cat), hspace=0.70)
        fig.legend(handles=[mpatches.Patch(facecolor=PPA_NEG, label=HPPD_NEG_LABEL, alpha=0.86),
                     mpatches.Patch(facecolor=PPA_POS, label=HPPD_POS_LABEL, alpha=0.86)],
                   loc="upper right", fontsize=10, frameon=False, ncol=2,
                   bbox_to_anchor=(0.97, 0.995))
        for i, (col, label, _, _m) in enumerate(cont):
            ax  = fig.add_subplot(gs_left[i])
            s_n = df_neg[col].dropna()
            s_p = df_pos[col].dropna()
            ax.text(-0.12, 1.0, f"{chr(ord('a') + i)}", transform=ax.transAxes,
                    fontsize=36, fontweight="bold", fontfamily="Arial",
                    va="bottom", ha="left", clip_on=False)
            all_d = np.asarray(pd.concat([s_n, s_p]), float)
            all_d = all_d[np.isfinite(all_d)]
            lo, hi  = np.percentile(all_d, [0.5, 99.5])
            x_range = (lo - (hi - lo) * 0.10, hi + (hi - lo) * 0.10)
            kde_n, kdm_n, _ = _half_violin(ax, s_n, PPA_NEG, flip=False, x_range=x_range)
            kde_p, kdm_p, _ = _half_violin(ax, s_p, PPA_POS, flip=True,  x_range=x_range)
            _quartile_marks(ax, s_n, kde=kde_n, kde_max=kdm_n, flip=False, color=PPA_NEG)
            _quartile_marks(ax, s_p, kde=kde_p, kde_max=kdm_p, flip=True,  color=PPA_POS)
            _style_violin_ax(ax, x_range, label, show_zero_line=True)
            ax.set_ylim(-1.15, 1.15)
            ax.text(-0.01, 0.75, HPPD_NEG_LABEL, transform=ax.transAxes, ha="right", va="center",
                    fontsize=8.5, color=PPA_NEG, fontweight="bold")
            ax.text(-0.01, 0.25, HPPD_POS_LABEL, transform=ax.transAxes, ha="right", va="center",
                    fontsize=8.5, color=PPA_POS, fontweight="bold")
        for j, (col, label, _, mapping) in enumerate(cat):
            ax = fig.add_subplot(gs_right[j])
            ax.text(-0.12, 1.0, f"{chr(ord('a') + len(cont) + j)}", transform=ax.transAxes,
                    fontsize=36, fontweight="bold", fontfamily="Arial",
                    va="bottom", ha="left", clip_on=False)
            _countplot_split(ax, df_neg[col].dropna(), df_pos[col].dropna(), mapping)
            ax.set_xlabel(label.replace("\n", " "), fontsize=18, labelpad=8, color="#333")
        fig.savefig(save_path, dpi=180, bbox_inches="tight")
        fig.savefig(str(Path(save_path).with_suffix('.svg')), format='svg', bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {save_path}")

    # ── Generate figures ───────────────────────────────────────────────────
    df_ppa_neg = df_spusers[df_spusers["persist_vis_yn"] == 0]
    df_ppa_pos = df_spusers[df_spusers["persist_vis_yn"] == 1]

    plot_sp_distributions(
        df=df_spusers, sections=SECTIONS,
        save_path=FIG_DIR / "sp_table_distributions.png",
    )
    _also_save_as(FIG_DIR / "sp_table_distributions.png")
    plot_sp_distributions_split(
        df=df_spusers, df_neg=df_ppa_neg, df_pos=df_ppa_pos,
        sections=SECTIONS_SPLIT,
        save_path=FIG_DIR / "sp_table_hppd_split_distributions.png",
    )
    _also_save_as(FIG_DIR / "sp_table_hppd_split_distributions.png")

    # CAPS vision item distributions split by HPPD history
    from caps_item_distributions_hppd_split import (
        make_caps_item_distributions_hppd_split,
    )
    make_caps_item_distributions_hppd_split(
        df_sp=df_sp,
        binary_palette=binary_palette,
        hppd_neg_label=HPPD_NEG_LABEL,
        hppd_pos_label=HPPD_POS_LABEL,
        save_path=FIG_DIR / "caps_item_distributions_hppd_split.png",
    )

    print(f"Descriptive figures saved to: {FIG_DIR}")


##############################################################################
### MEDIATION DIAGRAMS
# Batch-creates mediation diagrams for model directories whose names contain
# MODEL_TYPE. Only diagrams that match the model type are (re)generated so
# that sensitivity-analysis variants are not accidentally overwritten.
##############################################################################
if RUN_MEDIATION_DIAGRAMS:
    print("\n--- Mediation diagrams ---")

    def _batch_create_mediation_diagrams(mediation_models_dir, model_type_filter=None,
                                         skip_existing=False):
        """
        Create mediation diagrams for all model directories under mediation_models_dir.

        model_type_filter : if set, only process directories whose name contains this string.
        skip_existing     : if True, skip directories where *_AUTO.png already exists.
        """
        model_dirs = sorted([p for p in Path(mediation_models_dir).iterdir() if p.is_dir()])
        if model_type_filter:
            model_dirs = [p for p in model_dirs if model_type_filter.lower() in p.name.lower()]
            print(f"  Filtered to {len(model_dirs)} dirs matching '{model_type_filter}'")
        else:
            print(f"  Processing all {len(model_dirs)} dirs (no filter)")

        created, skipped_exist, skipped_missing, failures = 0, 0, 0, 0
        for model_dir in model_dirs:
            model_name = model_dir.name
            med_files  = sorted(model_dir.glob(f'mediation_results_*_{model_name}.csv'))
            if not med_files:
                med_files = sorted(model_dir.glob('mediation_results_*.csv'))
            if not med_files:
                skipped_missing += 1
                continue
            for med_file in med_files:
                try:
                    med_results = pd.read_csv(med_file)
                    if med_results.empty or 'effect' not in med_results.columns:
                        raise ValueError('CSV missing expected mediation columns')
                    predictor = (med_results['predictor'].iloc[0]
                                 if 'predictor' in med_results.columns else 'predictor')
                    savepath  = model_dir / f'mediation_diagram_{predictor}_AUTO'
                    if skip_existing and Path(str(savepath) + '.png').exists():
                        skipped_exist += 1
                        continue
                    fig, ax = create_mediation_diagram(
                        med_results,
                        dv_label=None,
                        dv_to_lab_dict=dv_to_lab_short,
                        box_fontsize=FONT_MEDIATION['box'],
                        stat_fontsize=FONT_MEDIATION['stat'],
                        savepath=str(savepath),
                    )
                    plt.close(fig)
                    _also_save_as(str(savepath) + '.png')
                    created += 1
                except Exception as exc:
                    failures += 1
                    print(f'    FAIL: {med_file.name} -> {exc}')

        print(f"  Created: {created} | Skipped (exists): {skipped_exist} | "
              f"Skipped (no CSV): {skipped_missing} | Failures: {failures}")

        # Ensure TIFF copies exist for any already-created PNGs we skipped above
        if FIGURE_FORMAT != "png":
            for model_dir in model_dirs:
                for png in model_dir.glob("mediation_diagram_*_AUTO.png"):
                    _also_save_as(png)

    def _batch_create_mediation_diagrams_counterfactual(mediation_models_dir,
                                                         model_type_filter=None,
                                                         skip_existing=False,
                                                         nie_delta_as_pct=False):
        """Create mediation diagrams using response-scale counterfactual A/B/C' paths.

        Reads path_counterfactual_summary.csv (response-scale; 94% HDI; A path / B path /
        C' path effect rows) alongside mediation_results_*.csv (for predictor / mediator /
        dv labels and indirect/prop_mediated).  Falls back gracefully when
        path_counterfactual_summary.csv is absent (e.g. cumulative/ordinal models that
        emit __SKIP__).

        Column mapping from path_counterfactual_summary.csv → create_mediation_diagram:
          effect     : "A path" → "a_path", "B path" → "b_path", "C' path" → "c_prime_direct"
          mean       → estimate   (POINT_ESTIMATE_COL, modules/master_config.py)
          hdi_low    → lower_95  (shown as "94% HDI" via ci_label param)
          hdi_high   → upper_95
          p_above_0  → prob_above_0   (already matches)
          p_below_0  → prob_below_0   (already matches)

        Saves as mediation_diagram_{predictor}_COUNTERFACTUAL.png (does NOT overwrite
        the standard _AUTO.png, so both versions coexist).
        """
        _EFFECT_MAP = {"A path": "a_path", "B path": "b_path", "C' path": "c_prime_direct"}

        model_dirs = sorted([p for p in Path(mediation_models_dir).iterdir() if p.is_dir()])
        if model_type_filter:
            model_dirs = [p for p in model_dirs
                          if model_type_filter.lower() in p.name.lower()]
            print(f"  Filtered to {len(model_dirs)} dirs matching '{model_type_filter}'")
        else:
            print(f"  Processing all {len(model_dirs)} dirs (no filter)")

        created, skipped_exist, skipped_missing, skipped_no_cf, failures = 0, 0, 0, 0, 0
        for model_dir in model_dirs:
            model_name = model_dir.name
            med_files  = sorted(model_dir.glob(f'mediation_results_*_{model_name}.csv'))
            if not med_files:
                med_files = sorted(model_dir.glob('mediation_results_*.csv'))
            if not med_files:
                skipped_missing += 1
                continue

            cf_csv = model_dir / 'path_counterfactual_summary.csv'
            if not cf_csv.exists():
                # Graceful skip: model family not supported for counterfactual paths
                skipped_no_cf += 1
                continue

            for med_file in med_files:
                try:
                    med_results = pd.read_csv(med_file)
                    if med_results.empty or 'effect' not in med_results.columns:
                        raise ValueError('mediation_results CSV missing expected columns')
                    cf_results = pd.read_csv(cf_csv)
                    if cf_results.empty or 'effect' not in cf_results.columns:
                        raise ValueError('path_counterfactual_summary.csv missing expected columns')

                    predictor = (med_results['predictor'].iloc[0]
                                 if 'predictor' in med_results.columns else 'predictor')
                    savepath = model_dir / f'mediation_diagram_{predictor}_COUNTERFACTUAL'
                    if skip_existing and Path(str(savepath) + '.png').exists():
                        skipped_exist += 1
                        continue

                    # ── Build combined DataFrame ───────────────────────────────────
                    # Map effect labels and rename columns to what create_mediation_diagram expects
                    cf = cf_results.copy()
                    cf['effect'] = cf['effect'].map(_EFFECT_MAP)
                    cf = cf.dropna(subset=['effect'])
                    # Add metadata columns from mediation_results
                    cf['predictor'] = med_results['predictor'].iloc[0]
                    cf['mediator']  = (med_results['mediator'].iloc[0]
                                       if 'mediator' in med_results.columns else '')
                    cf['dv']        = (med_results['dv'].iloc[0]
                                       if 'dv' in med_results.columns else '')
                    # Rename to create_mediation_diagram expected column names.
                    # p_above_0/p_below_0 are the cf CSV names; the diagram function
                    # expects prob_above_0/prob_below_0 for significance asterisks.
                    # POINT_ESTIMATE_COL ('mean') is the reported point
                    # estimate; see modules/master_config.py. The explicit
                    # presence check keeps a pre-mean CSV from silently
                    # producing a diagram with no estimates on it.
                    point_estimate(cf, source=str(model_dir / 'path_counterfactual_summary.csv'))
                    cf = cf.rename(columns={
                        POINT_ESTIMATE_COL: 'estimate',
                        'hdi_low':   'lower_95',
                        'hdi_high':  'upper_95',
                        'p_above_0': 'prob_above_0',
                        'p_below_0': 'prob_below_0',
                    })
                    # Append indirect/prop_mediated rows from standard mediation results
                    # (these quantities are from the MC integration, not counterfactual)
                    extra_rows = med_results[
                        med_results['effect'].isin(['indirect_ab', 'prop_mediated'])
                    ].copy()
                    combined = pd.concat([cf, extra_rows], ignore_index=True)

                    # Read NIE p_direction from mc_mediation_summary.csv — this is the
                    # MC-integration-based pd for the indirect effect and is the
                    # authoritative value (path-level prob_above/below_0 in
                    # mediation_results_*.csv differs because it uses the path posterior
                    # directly rather than the MC product distribution).
                    nie_pd    = None
                    nie_delta = None
                    mc_summary_csv = model_dir / 'mc_mediation_summary.csv'
                    if mc_summary_csv.exists():
                        mc_df = pd.read_csv(mc_summary_csv)
                        nie_rows = mc_df[mc_df['effect'].str.startswith('NIE')]
                        if not nie_rows.empty:
                            if 'p_direction' in nie_rows.columns:
                                nie_pd = float(nie_rows.iloc[0]['p_direction'])
                            nie_delta = float(point_estimate(
                                nie_rows.iloc[0], source=str(mc_summary_csv),
                                mc_integrated=True))

                    fig, ax = create_mediation_diagram(
                        combined,
                        dv_label=None,
                        dv_to_lab_dict=dv_to_lab_short,
                        # No HDI — estimate + significance asterisk only (matches _AUTO style).
                        # box_text_width=14 wraps "Avg. Dose Used" / "(LSD µg eq)" correctly.
                        box_text_width=16,
                        box_fontsize=FONT_MEDIATION['box'],
                        stat_fontsize=FONT_MEDIATION['stat'],
                        show_indirect_pd=True,
                        indirect_pd_override=nie_pd,
                        nie_delta=nie_delta,
                        nie_delta_as_pct=nie_delta_as_pct,
                        # Counterfactual diagrams show marginal response-scale changes (Δ),
                        # not standardised coefficients (Β).
                        coef_label="Δ",
                        savepath=str(savepath),
                    )
                    plt.close(fig)
                    _also_save_as(str(savepath) + '.png')
                    created += 1
                except Exception as exc:
                    failures += 1
                    print(f'    FAIL: {model_dir.name} -> {exc}')

        print(f"  Created: {created} | Skipped (exists): {skipped_exist} | "
              f"Skipped (no med CSV): {skipped_missing} | "
              f"Skipped (no CF CSV): {skipped_no_cf} | Failures: {failures}")

        if FIGURE_FORMAT != "png":
            for model_dir in model_dirs:
                for png in model_dir.glob("mediation_diagram_*_COUNTERFACTUAL.png"):
                    _also_save_as(png)

    _batch_create_mediation_diagrams_counterfactual(
        f"{RESULTS_BASE}/hppd_binary/mediation_models",
        model_type_filter=HPPD_MODEL_TYPE,
        skip_existing=False,
        nie_delta_as_pct=True,   # hppd_binary NIE is on probability scale (0–1)
    )
    _batch_create_mediation_diagrams_counterfactual(
        f"{RESULTS_BASE}/caps_vision/mediation_models",
        model_type_filter=MODEL_TYPE,
        skip_existing=False,
        nie_delta_as_pct=False,  # caps_vision NIE is on count scale
    )
    print("Mediation diagrams complete.")



##############################################################################
### FINAL FIGURE ASSEMBLY
# Assembles multi-panel manuscript figures from individual panel images.
# All paths are derived from MODEL_TYPE so only one variable needs to change.
##############################################################################
if RUN_FIGURE_ASSEMBLY:
    print("\n--- Figure assembly ---")

    BASE = str(Path('..').resolve())
    Path(FIGURES_OUT_DIR).mkdir(parents=True, exist_ok=True)

    # Model names used in forest-plot and mediation-diagram panel paths.
    # caps_vision uses the SP-user model so the
    # forest-plot path resolves to nice_covariates_spusers.png instead of
    # nice_covariates.png — matching the DV_CONFIG sample change above.
    _hppd_model = HPPD_MODEL_TYPE
    _caps_model  = HPPD_MODEL_TYPE

    def _assemble(config):
        """Save figure in FIGURE_FORMAT; also save PNG if format differs (for linked_figures).
        Also saves a true-vector SVG via svgutils assembly."""
        path = assemble_manuscript_figure(config, base_dir=BASE)
        print(f"  Saved: {path}")
        if FIGURE_FORMAT != "png":
            png_cfg = {**config, "output_path": config["output_path"].rsplit(".", 1)[0] + ".png"}
            assemble_manuscript_figure(png_cfg, base_dir=BASE)
            print(f"  Saved PNG copy: {png_cfg['output_path']}")
        # True-vector SVG assembly (parallel to raster pipeline)
        try:
            svg_cfg = {**config, "output_path": config["output_path"].rsplit(".", 1)[0] + ".svg"}
            assemble_manuscript_figure_svg(svg_cfg, base_dir=BASE)
        except Exception as e:
            print(f"  WARNING: vector SVG assembly failed: {e}")

    # ── Figures 1–3: SP predictors, VCH behavior, VCH computations ────────
    # Row 1 = data visualization; Row 2 = standardized forest plots
    # Shared label settings applied to all assembled figures
    _label_cfg = {
        "label_fontsize":    FIGURE_LABEL_FONTSIZE,
        "label_color":       "black",
        "label_x":           FIGURE_LABEL_X,
        "label_left_margin": FIGURE_LABEL_MARGIN,
    }

    # ── Helper: Kafadar-style HGF composite for figure_3 ──────────────────
    # figure_3 reproduces the layout of
    # figures/kafadar_et_al_compfig.png:
    #   Left half : hgf_alone.png — spans full figure height
    #   Right half : 5-row × 2-col grid (hppd_binary | caps_vision):
    #     Row 0: vch_omega  (boxplot_grid | correlation_grid)
    #     Row 1: xprob trajectory
    #     Row 2: xbin_pred trajectory
    #     Row 3: vch_nu     (boxplot_grid | correlation_grid)
    #     Row 4: vch_beta   (boxplot_grid | correlation_grid)
    # The HGF diagram maps top → bottom to HGF parameter hierarchy top → bottom,
    # making the panel layout directly interpretable alongside the model schematic.
    def _assemble_figure3_hgf(
        hgf_path,
        right_rows,
        output_path,
        dpi=600,
        col_panel_width_px=800,
        col_gap_px=280,
        row_gap_px=12,
        panel_gap_px=270,
        bg_color=(255, 255, 255),
        panel_label_size=PIL_PANEL_LABEL_PX,
        bottom_pad_px=40,
    ):
        """
        Composite figure_3 using PIL.

        Parameters
        ----------
        hgf_path : str | Path
            Path to hgf_alone.png.
        right_rows : list[tuple[str, str]]
            Five (hppd_path, caps_path) tuples, top-to-bottom.
        output_path : str | Path
            Destination PNG.
        dpi : int
            DPI metadata written into the output PNG.
        col_panel_width_px : int
            Target pixel width for each of the two right-side columns.
            Images are downsampled / upsampled to this width; aspect ratio
            is preserved so heights vary per panel.
        col_gap_px : int
            Horizontal gap between the HGF (panel a) and the first data column.
            Must be wide enough to hold panel labels without overlapping the
            matplotlib y-axis text that lives inside each panel image's left
            margin (~280 px accommodates 240 px Arial Bold labels).
        row_gap_px : int
            Vertical gap between right-panel rows.
        panel_gap_px : int
            Horizontal gap between the first and second data columns.
            Same sizing constraint as col_gap_px (~270 px).
        panel_label_size : int
            Font size (px) for the a–k panel labels drawn in the gaps.
        """
        from PIL import ImageDraw as _IDraw, ImageFont as _IFont

        hgf_path = Path(hgf_path)
        if not hgf_path.exists():
            raise FileNotFoundError(f"HGF image not found: {hgf_path}")

        # Scale each right panel to col_panel_width_px wide; record row heights
        scaled_rows = []
        for h_path_str, c_path_str in right_rows:
            h_path, c_path = Path(h_path_str), Path(c_path_str)
            if not h_path.exists():
                raise FileNotFoundError(f"Right panel not found: {h_path}")
            if not c_path.exists():
                raise FileNotFoundError(f"Right panel not found: {c_path}")
            h_img = _PIL_Image.open(h_path).convert('RGB')
            c_img = _PIL_Image.open(c_path).convert('RGB')
            h_h = max(1, round(col_panel_width_px * h_img.size[1] / h_img.size[0]))
            c_h = max(1, round(col_panel_width_px * c_img.size[1] / c_img.size[0]))
            h_img = h_img.resize((col_panel_width_px, h_h), _PIL_Image.LANCZOS)
            c_img = c_img.resize((col_panel_width_px, c_h), _PIL_Image.LANCZOS)
            row_h = max(h_h, c_h)
            scaled_rows.append((h_img, c_img, row_h))

        n_rows   = len(scaled_rows)
        right_h  = sum(r[2] for r in scaled_rows) + (n_rows - 1) * row_gap_px
        right_w  = col_panel_width_px * 2 + panel_gap_px

        # Scale HGF to exactly match right_h (preserving aspect ratio for width)
        hgf_src  = _PIL_Image.open(hgf_path).convert('RGB')
        hgf_w    = max(1, round(right_h * hgf_src.size[0] / hgf_src.size[1]))
        hgf_img  = hgf_src.resize((hgf_w, right_h), _PIL_Image.LANCZOS)

        # Build canvas (bottom_pad_px prevents x-axis labels on the last
        # row from being clipped at the canvas edge)
        canvas_h = right_h + bottom_pad_px
        total_w  = hgf_w + col_gap_px + right_w
        canvas   = _PIL_Image.new('RGB', (total_w, canvas_h), bg_color)
        canvas.paste(hgf_img, (0, 0))

        # Paste right panels row by row, vertically centring the shorter panel.
        # Track label positions in the whitespace GAPS between columns so labels
        # never overlap the matplotlib y-axis text inside each panel image.
        # Layout: a = HGF, b–k = right panels left→right, top→bottom.
        _panel_positions = [('a', 8, 8)]   # HGF label at its top-left corner
        _panel_letter = ord('b')
        x0 = hgf_w + col_gap_px
        y  = 0
        for h_img, c_img, row_h in scaled_rows:
            # Bottom-align so x-axis labels sit at the same y-position
            h_y_offset = row_h - h_img.size[1]
            c_y_offset = row_h - c_img.size[1]
            canvas.paste(h_img, (x0, y + h_y_offset))
            canvas.paste(c_img, (x0 + col_panel_width_px + panel_gap_px, y + c_y_offset))
            # Left-column label: 10 px into col_gap (between HGF and first data column)
            # Right-column label: 10 px into panel_gap (between the two data columns)
            # Both at y+8 — identical height between columns within each row.
            _panel_positions.append((chr(_panel_letter),
                                     hgf_w + 10,
                                     y + 8))
            _panel_positions.append((chr(_panel_letter + 1),
                                     x0 + col_panel_width_px + 10,
                                     y + 8))
            _panel_letter += 2
            y += row_h + row_gap_px

        # ── Panel labels (a–k) ────────────────────────────────────────────────
        # Scale label size so it appears at FIGURE_LABEL_FONTSIZE (20pt) when
        # the figure is printed at 7" wide.  The canvas may be wider than 7"
        # (Figure 3 is ~10.7"), so labels must be proportionally larger.
        _label_px = round(FIGURE_LABEL_FONTSIZE * total_w / (_FIG_W * 72))
        _arial_candidates = [
            '/Library/Fonts/Arial Bold.ttf',
            '/System/Library/Fonts/Supplemental/Arial Bold.ttf',
            '/Library/Fonts/Arial.ttf',
            '/System/Library/Fonts/Supplemental/Arial.ttf',
            '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf',
        ]
        _font = None
        for _fp in _arial_candidates:
            if Path(_fp).exists():
                try:
                    _font = _IFont.truetype(_fp, size=_label_px)
                    break
                except Exception:
                    pass
        if _font is None:
            _font = _IFont.load_default()

        _draw = _IDraw.Draw(canvas)
        # Use baseline anchor ('ls') so letters with descenders (g, j)
        # align by their baseline rather than top, matching letters like
        # b, d, f that have no descender.
        _ascent = _font.getmetrics()[0]
        for _letter, _px, _py in _panel_positions:
            _draw.text((_px, _py + _ascent), _letter, fill=(0, 0, 0),
                       font=_font, anchor='ls', stroke_width=0)

        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)

        ext = out.suffix.lower().lstrip('.')
        if ext == 'svg':
            # PIL cannot write SVG; wrap the raster composite in an SVG
            # container with physical dimensions derived from DPI.
            from figure_assembly import raster_to_svg as _r2s
            # Save PNG first (needed for linked_figures and as SVG source)
            _png_tmp = out.with_suffix('.png')
            canvas.save(str(_png_tmp), dpi=(dpi, dpi))
            _r2s(str(_png_tmp), str(out), dpi=dpi)
        else:
            canvas.save(str(out), dpi=(dpi, dpi))

        print(f"  Saved: {out}")

        # ── True-vector SVG assembly (parallel to PIL raster) ──────────────
        try:
            import svgutils.transform as _svgt
            from figure_assembly import _parse_svg_dimension, _get_svg_aspect_ratio

            _svg_out = out.with_suffix('.svg')
            _svg_elements = []

            # Scale factor: convert PIL pixel units to SVG points.
            # Target figure width in points (matching assembled PNG at given DPI).
            _total_w_px = canvas.size[0]
            _total_h_px = canvas.size[1]
            _pts_per_px = 72.0 / dpi  # e.g. 72/600 = 0.12

            # HGF diagram: use SVG source if available, otherwise raster fallback
            _hgf_svg_path = Path(hgf_path).with_suffix('.svg')
            if _hgf_svg_path.exists():
                _hgf_svg = _svgt.fromfile(str(_hgf_svg_path))
                _hgf_svg_root = _hgf_svg.getroot()
                # Parse dimensions: try width/height attrs first, fall back to viewBox
                _hgf_svg_w = _parse_svg_dimension(_hgf_svg.width) if _hgf_svg.width else None
                _hgf_svg_h = _parse_svg_dimension(_hgf_svg.height) if _hgf_svg.height else None
                if _hgf_svg_w is None or _hgf_svg_h is None:
                    import re as _re
                    _vb_match = _re.search(r'viewBox=["\'](\S+)\s+(\S+)\s+(\S+)\s+(\S+)["\']',
                                           open(str(_hgf_svg_path)).read(2000))
                    if _vb_match:
                        _hgf_svg_w = float(_vb_match.group(3))
                        _hgf_svg_h = float(_vb_match.group(4))
                if _hgf_svg_w and _hgf_svg_h:
                    # Scale to match the raster layout: target width = hgf_w px → pts
                    _hgf_scale = (hgf_w * _pts_per_px) / _hgf_svg_w
                    _hgf_svg_root.scale(_hgf_scale)
                    _hgf_svg_root.moveto(0, 0)
                    _svg_elements.append(_hgf_svg_root)
                else:
                    raise ValueError(f"Cannot determine dimensions of {_hgf_svg_path}")
            else:
                _hgf_img_elem = _svgt.ImageElement(str(hgf_path), hgf_w * _pts_per_px, right_h * _pts_per_px)
                _hgf_img_elem.moveto(0, 0)
                _svg_elements.append(_hgf_img_elem)

            # Right-side panels: use SVG versions where available
            _svg_letter = ord('b')
            _svg_y = 0.0
            _x0_pt = (hgf_w + col_gap_px) * _pts_per_px
            _col_w_pt = col_panel_width_px * _pts_per_px
            _panel_gap_pt = panel_gap_px * _pts_per_px
            _row_gap_pt_f3 = row_gap_px * _pts_per_px

            for h_img, c_img, row_h_px in scaled_rows:
                _row_h_pt = row_h_px * _pts_per_px
                for col_idx, (img_obj, orig_path_str) in enumerate(zip(
                    [h_img, c_img],
                    [right_rows[scaled_rows.index((h_img, c_img, row_h_px))][0],
                     right_rows[scaled_rows.index((h_img, c_img, row_h_px))][1]]
                )):
                    _x_pt = _x0_pt + col_idx * (_col_w_pt + _panel_gap_pt)
                    _svg_panel_path = Path(orig_path_str).with_suffix('.svg')
                    _img_h_pt = img_obj.size[1] * _pts_per_px
                    _y_offset_pt = _row_h_pt - _img_h_pt  # bottom-align

                    if _svg_panel_path.exists():
                        _sp = _svgt.fromfile(str(_svg_panel_path))
                        _sp_root = _sp.getroot()
                        _sp_w = _parse_svg_dimension(_sp.width)
                        _sp_h = _parse_svg_dimension(_sp.height)
                        _sp_scale = _col_w_pt / _sp_w
                        _sp_root.scale(_sp_scale)
                        _sp_root.moveto(_x_pt, _svg_y + _y_offset_pt)
                        _svg_elements.append(_sp_root)
                    else:
                        _img_elem = _svgt.ImageElement(orig_path_str, _col_w_pt, _img_h_pt)
                        _img_elem.moveto(_x_pt, _svg_y + _y_offset_pt)
                        _svg_elements.append(_img_elem)

                    # Panel label
                    _label_x_pt = (_x0_pt - col_gap_px * _pts_per_px + 10 * _pts_per_px) if col_idx == 0 \
                        else (_x0_pt + _col_w_pt + 10 * _pts_per_px)
                    _lbl = _svgt.TextElement(
                        _label_x_pt, _svg_y + 8 * _pts_per_px + _label_px * _pts_per_px,
                        chr(_svg_letter),
                        size=_label_px * _pts_per_px, weight='bold', font='Arial',
                    )
                    _svg_elements.append(_lbl)
                    _svg_letter += 1

                _svg_y += _row_h_pt + _row_gap_pt_f3

            # HGF label 'a'
            _lbl_a = _svgt.TextElement(
                8 * _pts_per_px, 8 * _pts_per_px + _label_px * _pts_per_px,
                'a', size=_label_px * _pts_per_px, weight='bold', font='Arial',
            )
            _svg_elements.append(_lbl_a)

            _fig3_svg = _svgt.SVGFigure(
                f'{_total_w_px * _pts_per_px}pt',
                f'{_total_h_px * _pts_per_px}pt',
            )
            _fig3_svg.append(_svg_elements)
            _fig3_svg.save(str(_svg_out))

            # svgutils.SVGFigure.save() does not write width/height attrs;
            # inject them so scale_svg_to_dpi can read the dimensions.
            from lxml import etree as _etree_f3
            _tree_f3 = _etree_f3.parse(str(_svg_out))
            _root_f3 = _tree_f3.getroot()
            _w_pt_f3 = _total_w_px * _pts_per_px
            _h_pt_f3 = _total_h_px * _pts_per_px
            _root_f3.set('width', f'{_w_pt_f3:.4f}pt')
            _root_f3.set('height', f'{_h_pt_f3:.4f}pt')
            _root_f3.set('viewBox', f'0 0 {_w_pt_f3:.4f} {_h_pt_f3:.4f}')
            _tree_f3.write(str(_svg_out), xml_declaration=True, standalone='yes', encoding='ASCII')

            # Set viewBox and DPI-scaled display dimensions (the final
            # figures section will override with target_width_mm later).
            from figure_assembly import scale_svg_to_dpi as _scale_f3
            _scale_f3(str(_svg_out), dpi=dpi)

            print(f"  Saved vector SVG: {_svg_out}")
        except Exception as _e:
            print(f"  WARNING: figure 3 vector SVG assembly failed: {_e}")

        return str(out)

    # ── Standard loop (figures 1–2 always; figure 3 only when not USE_CF) ──
    # Includes sdt_hppd (same handling as
    # vch_behavior): row a-b = data viz, row c-d = forest plots (conditional
    # on file existence), rows e-h = mediation diagrams (vch_behavior only).
    _standard_groups = [g for g in IVTYPES if g != 'vch_computations']
    for group in _standard_groups:
        # Forest plot row — only included if both PNG files already exist on
        # disk.  Mirrors the results-availability guard in RUN_HPPD_CAPS_FIGS:
        # if forest plots were skipped during generation (no matching CSV rows),
        # the assembly row is dropped gracefully rather than raising FileNotFoundError.
        _fp_hppd = Path(BASE) / f"results/hppd_binary/{group}/forest_plots/{_hppd_model}.png"
        _fp_caps  = Path(BASE) / f"results/caps_vision/{group}/forest_plots/{_caps_model}.png"
        _forest_row = [{"panels": [
            {"path": f"results/hppd_binary/{group}/forest_plots/{_hppd_model}.png", "label": "c"},
            {"path": f"results/caps_vision/{group}/forest_plots/{_caps_model}.png", "label": "d"},
        ]}] if (_fp_hppd.exists() and _fp_caps.exists()) else []
        if not _forest_row:
            print(f"  [assembly] forest plots not found for {group} — omitting row c/d")
        _assemble({
            "title": "",
            "output_path": f"{FIGURES_OUT_DIR}/{group}_hppd_vs_caps.{FIGURE_FORMAT}",
            "fig_width_inches": _FIG_W,
            "dpi": FIGURE_DPI,
            "gap_inches": 0.05,
            "row_gap_inches": 0.10,
            **_label_cfg,
            "rows": [
                {"panels": [
                    {"path": f"results/hppd_binary/{group}/data_visualization/boxplot_grid.png",     "label": "a"},
                    # sp_predictors uses the age-controlled variant (partial Spearman | age_v2);
                    # other groups (vch_behavior, sdt_hppd, …) use the standard correlation grid.
                    {"path": (
                         f"results/caps_vision/{group}/data_visualization/"
                         f"{'correlation_grid_age_control' if group == 'sp_predictors' else 'correlation_grid'}.png"
                     ), "label": "b",
                     "label_x": FIGURE_LABEL_X + 0.01},
                ]},
                *_forest_row,
                # ── Panels e–h: mediation diagrams (vch_behavior only, USE_CF) ──
                # Driven by FIG2_MED_* and FIG2_MED_ROW2_* config variables at the top.
                # Only appended for vch_behavior; omitted for sp_predictors.
                # Row 1 (e, f): vch_threshold mediator
                # Row 2 (g, h): vch_bl_yes_0 (vchrate) mediator
                *([
                    {
                        "panels": [
                            {
                                "path": (
                                    f"results/hppd_binary/mediation_models/"
                                    f"hppd_binary_{FIG2_MED_HPPD_SPVDR}_{FIG2_MED_HPPD_MEDIATOR}_{_hppd_model}/"
                                    f"mediation_diagram_{_SPVDR_TO_DIAGRAMNAME[FIG2_MED_HPPD_SPVDR]}_COUNTERFACTUAL.png"
                                ),
                                "label": "e",
                            },
                            {
                                "path": (
                                    f"results/caps_vision/mediation_models/"
                                    f"caps_vision_{FIG2_MED_CAPS_SPVDR}_{FIG2_MED_CAPS_MEDIATOR}_{_caps_model}/"
                                    f"mediation_diagram_{_SPVDR_TO_DIAGRAMNAME[FIG2_MED_CAPS_SPVDR]}_COUNTERFACTUAL.png"
                                ),
                                "label": "f",
                            },
                        ],
                    },
                    {
                        "panels": [
                            {
                                "path": (
                                    f"results/hppd_binary/mediation_models/"
                                    f"hppd_binary_{FIG2_MED_ROW2_HPPD_SPVDR}_{FIG2_MED_ROW2_HPPD_MEDIATOR}_{_hppd_model}/"
                                    f"mediation_diagram_{_SPVDR_TO_DIAGRAMNAME[FIG2_MED_ROW2_HPPD_SPVDR]}_COUNTERFACTUAL.png"
                                ),
                                "label": "g",
                            },
                            {
                                "path": (
                                    f"results/caps_vision/mediation_models/"
                                    f"caps_vision_{FIG2_MED_ROW2_CAPS_SPVDR}_{FIG2_MED_ROW2_CAPS_MEDIATOR}_{_caps_model}/"
                                    f"mediation_diagram_{_SPVDR_TO_DIAGRAMNAME[FIG2_MED_ROW2_CAPS_SPVDR]}_COUNTERFACTUAL.png"
                                ),
                                "label": "h",
                            },
                        ],
                    },
                ] if group == "vch_behavior" else []),
            ],
        })

    # ── Figure 3: HGF composite ───────────────────────────────────────────
    _cf3_rows = [
        (
            f"results/hppd_binary/vch_computations/data_visualization/vch_omega_boxplot_grid.png",
            f"results/caps_vision/vch_computations/data_visualization/vch_omega_correlation_grid.png",
        ),
        (
            f"results/hppd_binary/vch_computations/trajectories/xprob.png",
            f"results/caps_vision/vch_computations/trajectories/xprob.png",
        ),
        (
            f"results/hppd_binary/vch_computations/trajectories/xbin_pred.png",
            f"results/caps_vision/vch_computations/trajectories/xbin_pred.png",
        ),
        (
            f"results/hppd_binary/vch_computations/data_visualization/vch_nu_boxplot_grid.png",
            f"results/caps_vision/vch_computations/data_visualization/vch_nu_correlation_grid.png",
        ),
        (
            f"results/hppd_binary/vch_computations/data_visualization/vch_beta_boxplot_grid.png",
            f"results/caps_vision/vch_computations/data_visualization/vch_beta_correlation_grid.png",
        ),
    ]
    # Resolve paths relative to BASE (project root) to match other panels
    _cf3_rows_abs = [
        (str(Path(BASE) / h), str(Path(BASE) / c))
        for h, c in _cf3_rows
    ]
    _fig3_out = str(Path(BASE) / FIGURES_OUT_DIR / f"vch_computations_hppd_vs_caps.{FIGURE_FORMAT}")
    _assemble_figure3_hgf(
        hgf_path=str(Path(BASE) / "figures" / "hgf_alone.png"),
        right_rows=_cf3_rows_abs,
        output_path=_fig3_out,
        dpi=FIGURE_DPI,
        # 1300px = 6.5" × 200 DPI — matches correlation/trajectory source
        # widths so no PIL upscaling occurs and label sizes stay consistent.
        col_panel_width_px=1300,
        # col_gap=280, panel_gap=270: wide enough to hold 240px panel labels
        # in the pure-white gap without overlapping matplotlib y-axis text.
        col_gap_px=280,
        row_gap_px=15,
        panel_gap_px=270,
        panel_label_size=PIL_PANEL_LABEL_PX,
    )

    # ── Figure 4 panels (forest + mediation, appended below figure 3) ────────
    # Labels continue from figure 3's a–k: l, m, n, o.
    # Rendered at the same pixel width as figure 3 so they combine without
    # any post-hoc scaling (which would distort panel-label sizes).
    _cf4_hppd_med_dir = (
        f"results/hppd_binary/mediation_models/"
        f"hppd_binary_{FIG4_MED_HPPD_SPVDR}_{FIG4_MED_HPPD_MEDIATOR}_{_hppd_model}"
    )
    _cf4_caps_med_dir = (
        f"results/caps_vision/mediation_models/"
        f"caps_vision_{FIG4_MED_CAPS_SPVDR}_{FIG4_MED_CAPS_MEDIATOR}_{_caps_model}"
    )
    # Compute fig_width_inches so figure 4 renders at exactly figure 3's
    # pixel width.  This avoids any post-hoc PIL resize (which would also
    # scale the panel labels).
    from PIL import Image as _CombineImage
    _fig3_png = Path(BASE) / FIGURES_OUT_DIR / f"vch_computations_hppd_vs_caps.png"
    _fig4_width_inches = _FIG_W  # fallback
    if _fig3_png.exists():
        _fig3_w_px = _CombineImage.open(_fig3_png).size[0]
        _fig4_width_inches = _fig3_w_px / FIGURE_DPI

    _fig4_tmp_path = f"{FIGURES_OUT_DIR}/mediation_hppd_vs_caps.{FIGURE_FORMAT}"
    _assemble({
        "title": "",
        "output_path": _fig4_tmp_path,
        "fig_width_inches": _fig4_width_inches,
        "dpi": FIGURE_DPI,
        "gap_inches": 0.05,
        "row_gap_inches": 0.10,
        **_label_cfg,
        # Raise panel labels 2% above the panel top edge so they clear the
        # top y-axis label of the forest plots (panels l, m).
        "label_y": 1.015,
        "rows": [
            {"panels": [
                {"path": f"results/hppd_binary/vch_computations/forest_plots/{_hppd_model}.png", "label": "l"},
                {"path": f"results/caps_vision/vch_computations/forest_plots/{_caps_model}.png", "label": "m",
                 "label_x": FIGURE_LABEL_X + 0.035},
            ]},
            {"panels": [
                {
                    "path": (
                        f"{_cf4_hppd_med_dir}/"
                        f"mediation_diagram_{_SPVDR_TO_DIAGRAMNAME[FIG4_MED_HPPD_SPVDR]}_COUNTERFACTUAL.png"
                    ),
                    "label": "n",
                },
                {
                    "path": (
                        f"{_cf4_caps_med_dir}/"
                        f"mediation_diagram_{_SPVDR_TO_DIAGRAMNAME[FIG4_MED_CAPS_SPVDR]}_COUNTERFACTUAL.png"
                    ),
                    "label": "o",
                    "label_x": FIGURE_LABEL_X + 0.025,
                },
            ]},
        ],
    })

    # ── Combine figure 3 + figure 4 into one image ────────────────────────
    # Figure 4 was rendered at figure 3's pixel width, so no scaling needed
    # — just stack vertically with figure 4 centered.
    _fig4_png = Path(BASE) / FIGURES_OUT_DIR / f"mediation_hppd_vs_caps.png"
    if _fig3_png.exists() and _fig4_png.exists():
        _img3 = _CombineImage.open(_fig3_png).convert('RGB')
        _img4 = _CombineImage.open(_fig4_png).convert('RGB')
        _combined_w = max(_img3.width, _img4.width)
        _vert_gap = 40  # pixels between fig3 and fig4
        _combined_h = _img3.height + _vert_gap + _img4.height
        _combined = _CombineImage.new('RGB', (_combined_w, _combined_h), (255, 255, 255))
        _combined.paste(_img3, ((_combined_w - _img3.width) // 2, 0))
        _combined.paste(_img4, ((_combined_w - _img4.width) // 2, _img3.height + _vert_gap))
        _combined.save(str(_fig3_png), dpi=(FIGURE_DPI, FIGURE_DPI))
        print(f"  Combined figure 3+4 saved: {_fig3_png}")

        # ── Combined SVG (figure 3 + figure 4) ────────────────────────────
        # Stack the two SVGs vertically; figure 4 SVG already matches
        # figure 3's width (rendered at same inches), so no scaling needed.
        try:
            import svgutils.transform as _svgt_c
            import re as _re_svg
            from figure_assembly import _parse_svg_dimension, scale_svg_to_dpi

            def _svg_dim(svg_obj, attr, png_fallback=None):
                """Get SVG dimension in points, preferring viewBox (stable units).

                width/height attrs may have been scaled by scale_svg_to_dpi
                to display-pixels, so viewBox is the reliable source of the
                coordinate-system size in points.
                """
                # Prefer viewBox (always in the original point coordinate system)
                root = svg_obj.root
                vb = root.get('viewBox')
                if vb:
                    parts = vb.split()
                    return float(parts[2]) if attr == 'width' else float(parts[3])
                # Fall back to width/height attrs
                val = getattr(svg_obj, attr, None)
                if val is not None and str(val) != 'None':
                    return _parse_svg_dimension(val)
                # Fall back to PNG pixel dimensions → points
                if png_fallback and png_fallback.exists():
                    _fb = _CombineImage.open(png_fallback)
                    idx = 0 if attr == 'width' else 1
                    return _fb.size[idx] * (72.0 / FIGURE_DPI)
                return None

            _fig3_svg_path = _fig3_png.with_suffix('.svg')
            _fig4_svg_path = _fig4_png.with_suffix('.svg')
            if _fig3_svg_path.exists() and _fig4_svg_path.exists():
                # Combine SVGs using pure lxml (svgutils silently drops
                # content when re-serializing complex SVGs).
                from lxml import etree as _etree_comb
                import copy as _copy_comb
                _tree3 = _etree_comb.parse(str(_fig3_svg_path))
                _tree4 = _etree_comb.parse(str(_fig4_svg_path))
                _root3 = _tree3.getroot()
                _root4 = _tree4.getroot()
                _svgns = 'http://www.w3.org/2000/svg'

                # Extract viewBox dimensions (points)
                def _vb_dims(root, png_fallback):
                    vb = root.get('viewBox')
                    if vb:
                        parts = vb.split()
                        return float(parts[2]), float(parts[3])
                    # Fallback: derive from PNG
                    _img = _CombineImage.open(png_fallback)
                    return (_img.size[0] * 72.0 / FIGURE_DPI,
                            _img.size[1] * 72.0 / FIGURE_DPI)

                _w3, _h3 = _vb_dims(_root3, _fig3_png)
                _w4, _h4 = _vb_dims(_root4, _fig4_png)
                _gap_pt = _vert_gap * (72.0 / FIGURE_DPI)
                _total_svg_w = max(_w3, _w4)
                _total_svg_h = _h3 + _gap_pt + _h4

                # Build new SVG root
                _new_root = _etree_comb.Element(
                    '{%s}svg' % _svgns,
                    nsmap={None: _svgns,
                           'xlink': 'http://www.w3.org/1999/xlink'})
                _new_root.set('version', '1.1')
                _new_root.set('width', f'{_total_svg_w:.4f}pt')
                _new_root.set('height', f'{_total_svg_h:.4f}pt')
                _new_root.set('viewBox',
                              f'0 0 {_total_svg_w:.4f} {_total_svg_h:.4f}')

                # Fig 3 group (centered)
                _g3 = _etree_comb.SubElement(_new_root, '{%s}g' % _svgns)
                _dx3 = (_total_svg_w - _w3) / 2
                _g3.set('transform', f'translate({_dx3:.4f},0)')
                for _child in _root3:
                    _g3.append(_copy_comb.deepcopy(_child))

                # Fig 4 group (centered, below fig 3)
                _g4 = _etree_comb.SubElement(_new_root, '{%s}g' % _svgns)
                _dx4 = (_total_svg_w - _w4) / 2
                _g4.set('transform',
                        f'translate({_dx4:.4f},{_h3 + _gap_pt:.4f})')
                for _child in _root4:
                    _g4.append(_copy_comb.deepcopy(_child))

                _etree_comb.ElementTree(_new_root).write(
                    str(_fig3_svg_path), xml_declaration=True,
                    standalone='yes', encoding='ASCII')
                scale_svg_to_dpi(str(_fig3_svg_path), dpi=FIGURE_DPI)
                print(f"  Combined SVG saved: {_fig3_svg_path}")
        except Exception as _svg_e:
            print(f"  WARNING: combined SVG assembly failed: {_svg_e}")
    else:
        if not _fig3_png.exists():
            print(f"  WARNING: cannot combine — missing {_fig3_png}")
        if not _fig4_png.exists():
            print(f"  WARNING: cannot combine — missing {_fig4_png}")

    print("Figure assembly complete.")


##############################################################################
### FIGURE 6 — VCH BETA / SDT / DETECTION CURVES
#
# Layout: 2×2 quadrant (assembled via assemble_manuscript_figure).
#   Row 1: a = detection curves, b = VCH×SDT correlation grid
#   Row 2: c = HPPD sdt_hppd forest plot, d = CAPS sdt_hppd forest plot
#
# Panel sources:
#   a  Detection curves by β/ν median-split quadrant
#      (replicated from the nu/beta quadrant exploration;
#       03_nu_beta_quadrants.py; x-label → "QUEST-based detection probability";
#       title removed; uses the df loaded at the top of this script)
#   b  results/vch_beta/sdt_hppd/data_visualizations/correlation_grid.png
#      (generated above by correlation_matrix_plot; color = caps_vision)
#   c  results/hppd_binary/sdt_hppd/forest_plots/nice_covariates_spusers.png
#   d  results/caps_vision/sdt_hppd/forest_plots/nice_covariates_spusers.png
##############################################################################
if RUN_FIGURE_ASSEMBLY:
    # BASE may already be set by the earlier RUN_FIGURE_ASSEMBLY block; redefine
    # here so this block works independently if the earlier one was skipped.
    BASE = str(Path('..').resolve())
    Path(FIGURES_OUT_DIR).mkdir(parents=True, exist_ok=True)

    # ── Figure 6 panel b: detection curves by β/ν median-split quadrant ──────
    # Uses the df loaded at the top of this script; does not pin to a
    # separate CSV.

    print("\n--- Figure 6 panel b: detection curves ---")

    _F6_CONTRAST_LEVELS = [0, 25, 50, 75]
    _F6_HIT_RATE_COLS   = ['vch_bl_yes_0', 'vch_bl_yes_25', 'vch_bl_yes_50', 'vch_bl_yes_75']
    _F6_QUAD_LABELS = [
        r'$\beta^-\nu^-$',
        r'$\beta^-\nu^+$',
        r'$\beta^+\nu^-$',
        r'$\beta^+\nu^+$',
    ]
    # Four evenly spaced colors from the electric_blue_palette (indices 0, 2, 4, 6)
    _F6_PALETTE_4 = [
        electric_blue_palette[0],
        electric_blue_palette[2],
        electric_blue_palette[4],
        electric_blue_palette[6],
    ]
    _F6_N_BOOT    = 2000
    _F6_BOOT_SEED = 42

    def _f6_bootstrap_ci(values, n_boot=_F6_N_BOOT, ci=95, seed=_F6_BOOT_SEED):
        """Nonparametric bootstrap CI for the mean. Returns (lower, upper)."""
        rng = np.random.default_rng(seed)
        boot_means = np.array([
            rng.choice(values, size=len(values), replace=True).mean()
            for _ in range(n_boot)
        ])
        lo = np.percentile(boot_means, (100 - ci) / 2)
        hi = np.percentile(boot_means, 100 - (100 - ci) / 2)
        return lo, hi

    def _f6_curve_stats(subdf):
        rows = []
        for contrast, col in zip(_F6_CONTRAST_LEVELS, _F6_HIT_RATE_COLS):
            vals = subdf[col].dropna().values
            if len(vals) == 0:
                rows.append({'contrast': contrast, 'mean': np.nan,
                             'ci_lo': np.nan, 'ci_hi': np.nan})
                continue
            mean = vals.mean()
            lo, hi = _f6_bootstrap_ci(vals)
            rows.append({'contrast': contrast, 'mean': mean, 'ci_lo': lo, 'ci_hi': hi})
        return pd.DataFrame(rows)

    # Compute quadrant thresholds from all participants with valid vch_beta + vch_nu
    _f6_required = ['vch_beta', 'vch_nu'] + _F6_HIT_RATE_COLS
    _f6_missing = [c for c in _f6_required if c not in df.columns]
    if _f6_missing:
        print(f"  WARNING: columns missing from df — skipping figure 6 panel b: {_f6_missing}")
        _f6_panel_b_path = None
    else:
        _f6_valid = df.dropna(subset=['vch_beta', 'vch_nu'])
        _f6_beta_median = _f6_valid['vch_beta'].median()
        _f6_nu_median   = _f6_valid['vch_nu'].median()
        print(f"  Quadrant medians: vch_beta={_f6_beta_median:.4f}, vch_nu={_f6_nu_median:.4f}")

        def _f6_assign_quad(row):
            b, n = row['vch_beta'], row['vch_nu']
            if pd.isna(b) or pd.isna(n):
                return np.nan
            beta_hi = b >= _f6_beta_median
            nu_hi   = n >= _f6_nu_median
            if not beta_hi and not nu_hi:
                return _F6_QUAD_LABELS[0]
            elif not beta_hi and nu_hi:
                return _F6_QUAD_LABELS[1]
            elif beta_hi and not nu_hi:
                return _F6_QUAD_LABELS[2]
            else:
                return _F6_QUAD_LABELS[3]

        _f6_df = df.copy()
        _f6_df['_f6_quadrant'] = _f6_df.apply(_f6_assign_quad, axis=1)

        # Rank quadrants by mean vch_bl_yes_0 to assign colors from dark → bright
        _f6_df_valid = _f6_df[_f6_df['_f6_quadrant'].notna()]
        _f6_quad_means = {
            lbl: _f6_df_valid[_f6_df_valid['_f6_quadrant'] == lbl]['vch_bl_yes_0'].mean()
            for lbl in _F6_QUAD_LABELS
        }
        _f6_quad_order  = sorted(_F6_QUAD_LABELS, key=lambda l: _f6_quad_means[l])
        _f6_quad_colors = list(_F6_PALETTE_4)

        print(f"  Quadrant order (lowest → highest vch_bl_yes_0 mean): {_f6_quad_order}")
        for lbl in _f6_quad_order:
            n_q = (_f6_df['_f6_quadrant'] == lbl).sum()
            print(f"    {lbl}: n={n_q}, mean vch_bl_yes_0={_f6_quad_means[lbl]:.4f}")

        # Build detection-curves figure for panel b
        _f6_curves_df = _f6_df[_f6_df['_f6_quadrant'].notna()].dropna(
            subset=[_F6_HIT_RATE_COLS[0]], how='all'
        )
        _f6_groups = {lbl: _f6_curves_df[_f6_curves_df['_f6_quadrant'] == lbl]
                      for lbl in _F6_QUAD_LABELS}

        _f6_fig, _f6_ax = plt.subplots(figsize=(7, 5))
        # Arial enforced by global rcParams at the top of this script

        for lbl, color in zip(_f6_quad_order, _f6_quad_colors):
            _f6_sub = _f6_groups[lbl]
            _f6_n   = len(_f6_sub)
            _f6_stats = _f6_curve_stats(_f6_sub)

            # Thin individual participant lines
            for _, _row in _f6_sub.iterrows():
                _vals = [_row.get(col, np.nan) for col in _F6_HIT_RATE_COLS]
                if any(pd.isna(v) for v in _vals):
                    continue
                _f6_ax.plot(_F6_CONTRAST_LEVELS, _vals,
                            color=color, alpha=0.07, lw=0.6)

            # 95% bootstrap CI band
            _f6_ax.fill_between(
                _f6_stats['contrast'], _f6_stats['ci_lo'], _f6_stats['ci_hi'],
                color=color, alpha=0.20,
            )
            # Group mean line
            _f6_ax.plot(
                _f6_stats['contrast'], _f6_stats['mean'],
                color=color, lw=2.5, marker='o', markersize=6,
                label=f'{lbl}  (n={_f6_n})',
            )

        _f6_ax.set_xlim(-5, 80)
        _f6_ax.set_ylim(-0.02, 1.05)
        _f6_ax.set_xticks(_F6_CONTRAST_LEVELS)
        _f6_ax.set_xticklabels([f'{c}%\nDetection\nProbability' for c in _F6_CONTRAST_LEVELS],
                               fontsize=FONT_FIG6_DETECTION['tick'])
        _f6_ax.set_xlabel('QUEST-Derived Stimulus Intensity\n(% Detection Probability)',
                          fontsize=FONT_FIG6_DETECTION['xlabel'])
        _f6_ax.set_ylabel("Proportion 'yes' responses",
                          fontsize=FONT_FIG6_DETECTION['ylabel'])
        # No title (removed per figure 6 panel b spec)
        _f6_ax.tick_params(axis='y', labelsize=FONT_FIG6_DETECTION['tick'])
        _f6_ax.legend(fontsize=FONT_FIG6_DETECTION['legend'], loc='upper left', framealpha=0.85)
        for spine in _f6_ax.spines.values():
            spine.set_visible(False)

        _f6_fig.tight_layout()
        _f6_panel_b_dir = Path('figures')
        _f6_panel_b_dir.mkdir(parents=True, exist_ok=True)
        _f6_panel_b_path = str(_f6_panel_b_dir / 'figure6b_detection_curves.png')
        _f6_fig.savefig(_f6_panel_b_path, dpi=FIGURE_DPI, bbox_inches='tight')
        _f6_fig.savefig(str(Path(_f6_panel_b_path).with_suffix('.svg')), format='svg', bbox_inches='tight')
        plt.close(_f6_fig)
        # Embed DPI metadata (matplotlib + bbox_inches='tight' drops pHYs chunk)
        _PIL_Image.open(_f6_panel_b_path).save(_f6_panel_b_path, dpi=(FIGURE_DPI, FIGURE_DPI))
        print(f"  Saved panel b: {_f6_panel_b_path}")

    # ── Figure 6 assembly ─────────────────────────────────────────────────────
    # 2×2 quadrant via assemble_manuscript_figure (replaces old PIL 3+2 layout).
    # Row 1: a = detection curves, b = VCH×SDT correlation grid
    # Row 2: c = HPPD forest, d = CAPS forest

    _f6_panel_paths = {
        'a': str(Path(_f6_panel_b_path).resolve()) if _f6_panel_b_path else None,
        'b': str(Path(BASE) / 'results/vch_beta/sdt_hppd/data_visualizations/correlation_grid.png'),
        'c': str(Path(BASE) / f'results/hppd_binary/sdt_hppd/forest_plots/{HPPD_MODEL_TYPE}.png'),
        'd': str(Path(BASE) / f'results/caps_vision/sdt_hppd/forest_plots/{HPPD_MODEL_TYPE}.png'),
    }

    _f6_missing_panels = [k for k, p in _f6_panel_paths.items()
                          if p is None or not Path(p).exists()]
    if _f6_missing_panels:
        print(f"  WARNING: figure 6 panels missing — skipping assembly: {_f6_missing_panels}")
        for k in _f6_missing_panels:
            print(f"    panel {k}: {_f6_panel_paths[k]}")
    else:
        _assemble({
            "title": "",
            "output_path": f"{FIGURES_OUT_DIR}/figure_6.{FIGURE_FORMAT}",
            "fig_width_inches": _FIG_W,
            "dpi": FIGURE_DPI,
            "gap_inches": 0.05,
            "row_gap_inches": 0.10,
            **_label_cfg,
            "rows": [
                {"panels": [
                    {"path": _f6_panel_paths['a'], "label": "a"},
                    {"path": _f6_panel_paths['b'], "label": "b",
                     "label_x": FIGURE_LABEL_X + 0.01},
                ]},
                {"panels": [
                    {"path": _f6_panel_paths['c'], "label": "c"},
                    {"path": _f6_panel_paths['d'], "label": "d",
                     "label_x": FIGURE_LABEL_X + 0.035},
                ]},
            ],
        })

    print("Figure 6 assembly complete.")


##############################################################################
### DHARMA DIAGNOSTIC REVIEW (manuscript mediation panels — Figs 2 + 4)
# Checks every constituent equation in each mediation model panel and flags
# any DHARMa tests with p < ALPHA. Prints the path to the comprehensive PNG
# for borderline cases that warrant visual inspection.
#
# Panels checked (driven by FIG2_MED_* and FIG4_MED_* config vars at top):
#   Fig 2 e: hppd_binary  FIG2_MED_HPPD_SPVDR → FIG2_MED_HPPD_MEDIATOR
#   Fig 2 f: caps_vision  FIG2_MED_CAPS_SPVDR → FIG2_MED_CAPS_MEDIATOR
#   Fig 4 c: hppd_binary  FIG4_MED_HPPD_SPVDR → FIG4_MED_HPPD_MEDIATOR
#   Fig 4 d: caps_vision  FIG4_MED_CAPS_SPVDR → FIG4_MED_CAPS_MEDIATOR
#
# NOTE: dharma_heteroscedasticity_pval = testQuantiles(sim_residuals, predictor = spvar)
#       = BH-adjusted combined p across quantiles, vs. the SP predictor variable.
#       dharma_heteroscedasticity_q25/q50/q75_pval = per-quantile pvals from the same call.
# NOTE: testQuantiles uses an internal permutation test — borderline failures
#       can be stochastic. Always inspect the PNG for borderline cases.
##############################################################################
if RUN_DIAGNOSTICS:
    print("\n--- DHARMa diagnostics ---")

    ALPHA = 0.05

    MEDIATOR_LABELS = {
        "vch_threshold": "vch_threshold (75% Detection Threshold, Student-t)",
        "vch_bl_yes_0":  "vch_bl_yes_0 (VCH Rate at 0%, Zero-Inflated Beta)",
        "vch_bl_yes_75": "vch_bl_yes_75 (VCH Rate at 75%, Beta)",
        "vch_nu":        "vch_nu (Prior Weighting ν, Gamma)",
        "vch_beta":      "vch_beta (Decision Noise β, Student-t)",
        "vch_omega":     "vch_omega (Learning Rate ω, Student-t)",
    }
    DV_LABELS = {
        "hppd_binary":    "hppd_binary (HPPD History, Bernoulli)",
        "caps_vision":    "caps_vision (CAPS Vision, Hurdle Neg-Binom [hu-vary])",
        "caps_total":     "caps_total (CAPS Total, Zero-Inflated Negbinom)",
    }
    DHARMA_TESTS = {
        "dharma_ks_pval":                  "Uniformity / KS (testUniformity)              ",
        "dharma_dispersion_pval":          "Dispersion (testDispersion)                   ",
        "dharma_outlier_pval":             "Outlier (testOutliers)                         ",
        "dharma_zeroinflation_pval":       "Zero-inflation (testZeroInflation)             ",
        # testQuantiles(sim_residuals, predictor = fit$data[[spvar]]) — combined BH-adjusted p
        "dharma_heteroscedasticity_pval":  "Quantile combined vs spvar (testQuantiles)    ",
        # Per-quantile pvals from same testQuantiles call ($pvals[1:3]); added June 2026.
        "dharma_heteroscedasticity_q25_pval": "  Quantile q25 vs spvar                    ",
        "dharma_heteroscedasticity_q50_pval": "  Quantile q50 vs spvar                    ",
        "dharma_heteroscedasticity_q75_pval": "  Quantile q75 vs spvar                    ",
    }

    # Caps model type mirrors the assembly logic: spusers.
    _diag_caps_model = HPPD_MODEL_TYPE

    MED_PANEL_MODELS = [
        # ── Figure 2: panels e + f (row 1 — vch_threshold mediator) ──────────
        dict(
            path=(f"results/hppd_binary/mediation_models/"
                  f"hppd_binary_{FIG2_MED_HPPD_SPVDR}_{FIG2_MED_HPPD_MEDIATOR}_{HPPD_MODEL_TYPE}"),
            panel="Fig2-e",
            fig_label=(f"Fig2e  hppd_binary | {FIG2_MED_HPPD_SPVDR} → {FIG2_MED_HPPD_MEDIATOR}"
                       f"  [{HPPD_MODEL_TYPE}]"),
        ),
        dict(
            path=(f"results/caps_vision/mediation_models/"
                  f"caps_vision_{FIG2_MED_CAPS_SPVDR}_{FIG2_MED_CAPS_MEDIATOR}_{_diag_caps_model}"),
            panel="Fig2-f",
            fig_label=(f"Fig2f  caps_vision | {FIG2_MED_CAPS_SPVDR} → {FIG2_MED_CAPS_MEDIATOR}"
                       f"  [{_diag_caps_model}]"),
        ),
        # ── Figure 2: panels g + h (row 2 — vch_bl_yes_0 / vchrate mediator) ─
        dict(
            path=(f"results/hppd_binary/mediation_models/"
                  f"hppd_binary_{FIG2_MED_ROW2_HPPD_SPVDR}_{FIG2_MED_ROW2_HPPD_MEDIATOR}_{HPPD_MODEL_TYPE}"),
            panel="Fig2-g",
            fig_label=(f"Fig2g  hppd_binary | {FIG2_MED_ROW2_HPPD_SPVDR} → {FIG2_MED_ROW2_HPPD_MEDIATOR}"
                       f"  [{HPPD_MODEL_TYPE}]"),
        ),
        dict(
            path=(f"results/caps_vision/mediation_models/"
                  f"caps_vision_{FIG2_MED_ROW2_CAPS_SPVDR}_{FIG2_MED_ROW2_CAPS_MEDIATOR}_{_diag_caps_model}"),
            panel="Fig2-h",
            fig_label=(f"Fig2h  caps_vision | {FIG2_MED_ROW2_CAPS_SPVDR} → {FIG2_MED_ROW2_CAPS_MEDIATOR}"
                       f"  [{_diag_caps_model}]"),
        ),
        # ── Figure 4: panels c + d ────────────────────────────────────────────
        dict(
            path=(f"results/hppd_binary/mediation_models/"
                  f"hppd_binary_{FIG4_MED_HPPD_SPVDR}_{FIG4_MED_HPPD_MEDIATOR}_{HPPD_MODEL_TYPE}"),
            panel="Fig4-c",
            fig_label=(f"Fig4c  hppd_binary | {FIG4_MED_HPPD_SPVDR} → {FIG4_MED_HPPD_MEDIATOR}"
                       f"  [{HPPD_MODEL_TYPE}]"),
        ),
        dict(
            path=(f"results/caps_vision/mediation_models/"
                  f"caps_vision_{FIG4_MED_CAPS_SPVDR}_{FIG4_MED_CAPS_MEDIATOR}_{_diag_caps_model}"),
            panel="Fig4-d",
            fig_label=(f"Fig4d  caps_vision | {FIG4_MED_CAPS_SPVDR} → {FIG4_MED_CAPS_MEDIATOR}"
                       f"  [{_diag_caps_model}]"),
        ),
    ]
    # Legacy alias so the loop below doesn't need renaming.
    FIG4_MODELS = MED_PANEL_MODELS

    def _label_for_variable(varname):
        base = varname.removesuffix("_normalized")
        for lookup, eq_type in [(DV_LABELS, "DV equation       "), (MEDIATOR_LABELS, "Mediator equation ")]:
            if base in lookup:
                return f"{eq_type}: {lookup[base]}"
            if varname in lookup:
                return f"{eq_type}: {lookup[varname]}"
            for key, label in lookup.items():
                if base == key.replace("_", "") or varname == key.replace("_", ""):
                    return f"{eq_type}: {label}"
        return f"Equation          : {varname}"

    def _comprehensive_png(model_path, varname):
        stripped  = varname.removesuffix("_normalized").replace("_", "")
        candidate = model_path / f"dharma_comprehensive_{stripped}_{model_path.name}.png"
        return candidate if candidate.exists() else None

    def check_model(model_dict):
        path  = Path(BASE) / model_dict["path"]
        panel = model_dict["panel"]
        label = model_dict["fig_label"]
        print(f"\n{'='*72}")
        print(f"  Panel ({panel.upper()})  {label}")
        print(f"  Directory: {path.name}")
        print(f"{'='*72}")
        if not path.exists():
            print("  WARNING  Directory not found — model not yet run or results not pulled.")
            return
        summary_csvs = sorted(path.glob("summary_*.csv"))
        if not summary_csvs:
            print("  WARNING  No summary_*.csv files found.")
            return
        any_fail = False
        for csv_path in summary_csvs:
            stem     = csv_path.stem
            model_id = path.name
            prefix   = "summary_"
            if stem.startswith(prefix) and stem.endswith(model_id):
                varname = stem[len(prefix):len(stem) - len(model_id)].rstrip("_")
            else:
                varname = stem[len(prefix):]
            print(f"\n  {_label_for_variable(varname)}")
            try:
                df_csv = pd.read_csv(csv_path)
            except Exception as e:
                print(f"    WARNING  Could not read CSV: {e}")
                continue
            row = df_csv.iloc[0]
            n   = int(row["N"])              if "N"              in row.index else "?"
            nd  = int(row["num_divergents"]) if "num_divergents" in row.index else "?"
            print(f"    N = {n},  num_divergents = {nd}")
            found_any = False
            for col, test_label in DHARMA_TESTS.items():
                if col not in df_csv.columns:
                    continue
                pval = row[col]
                if pd.isna(pval):
                    print(f"    -  {test_label}: p = NA")
                    continue
                found_any = True
                pval = float(pval)
                fail = pval < ALPHA
                if fail:
                    any_fail = True
                print(f"    {'FAIL' if fail else 'OK  '}  {test_label}: p = {pval:.4f}")
            if not found_any:
                print("    (no DHARMa columns found)")
            png = _comprehensive_png(path, varname)
            if png:
                print(f"    -> Comprehensive PNG: {png}")
            else:
                fallback = sorted(path.glob(f"dharma_comprehensive_*_{path.name}.png"))
                eq_s = varname.removesuffix("_normalized").replace("_", "")
                for f in ([x for x in fallback if x.name.startswith(f"dharma_comprehensive_{eq_s}_")] or fallback):
                    print(f"    -> Comprehensive PNG: {f}")
        print(f"\n  {'WARNING  One or more tests FAILED' if any_fail else 'OK  All tests passed'} (alpha = {ALPHA})")

    BASE = str(Path('..').resolve())
    print(f"Checking {len(FIG4_MODELS)} manuscript mediation panels (Figs 2+4) | alpha = {ALPHA}")
    for model in FIG4_MODELS:
        check_model(model)
    print(f"\n{'='*72}")
    print("Diagnostics complete. Borderline testQuantiles failures may be stochastic.")

    # ── Diagnostic compilation PDFs (batch, filtered to HPPD_MODEL_TYPE) ────────
    # Generates diagnostic_compilation.png for every model whose directory name
    # ends with HPPD_MODEL_TYPE, then assembles per-DV PDFs in figures/.
    # R generation runs in parallel (4 workers); already-cached models are skipped.
    if RUN_DIAGNOSTIC_COMPILATION:
        from compile_mediation_diagnostic_pdfs import compile_diagnostic_pdfs
        print(f"\n--- Diagnostic compilation PDFs (model_type={HPPD_MODEL_TYPE}) ---")
        compile_diagnostic_pdfs(model_type=HPPD_MODEL_TYPE)

    # ── Single-path nonsp diagnostic PDFs (reads from HPC_RESULTS_MIRROR) ────────
    # Step 1: retrieve *_diagnostic_compilation.png files from HPC via targeted
    #         tarball (requires active ControlMaster: ssh -MNf bouchet).
    #         User is prompted to confirm DUO — decline to skip retrieval and
    #         go straight to PDF assembly with whatever is already in the mirror.
    # Step 2: assemble figures/single_path_diagnostics_{dv}.pdf.
    if RUN_SINGLE_PATH_DIAGNOSTICS:
        import subprocess as _subprocess
        _diag_generator = Path(__file__).parent.parent / '03_hpc' / 'generate_nonsp_diagnostic_jobs.py'
        print(f"\n--- Single-path diagnostic retrieval (model_type={HPPD_MODEL_TYPE}) ---")
        print("  (Requires ssh -MNf bouchet; decline the DUO prompt to skip and use cached PNGs.)")
        _subprocess.run([sys.executable, str(_diag_generator), '--retrieve'], check=False)

        from compile_single_path_diagnostic_pdfs import compile_single_path_pdfs
        print(f"\n--- Single-path diagnostic PDFs (model_type={HPPD_MODEL_TYPE}) ---")
        compile_single_path_pdfs(model_type=HPPD_MODEL_TYPE)


##############################################################################
### LINKED FIGURES — copy assembled figures to linked_figures/ for Google Doc
##############################################################################
if RUN_LINKED_FIGURES:
    print("\n--- Updating linked figures ---")
    from PIL import Image as _PIL_Image

    LINKED_FIGURES_DIR = Path("linked_figures")

    def _ensure_size_under_mb(src_path, max_mb=1):
        max_bytes = max_mb * 1024 * 1024
        if not src_path.exists() or src_path.stat().st_size <= max_bytes:
            return True
        scale    = (max_bytes / src_path.stat().st_size) ** 0.5
        img      = _PIL_Image.open(src_path)
        new_w    = max(1, int(img.width * scale))
        new_h    = max(1, int(img.height * scale))
        img      = img.resize((new_w, new_h), _PIL_Image.LANCZOS)
        img.save(src_path, optimize=True, compress_levels=9)
        while src_path.stat().st_size > max_bytes and min(new_w, new_h) > 50:
            new_w = max(1, int(new_w * 0.9))
            new_h = max(1, int(new_h * 0.9))
            img   = img.resize((new_w, new_h), _PIL_Image.LANCZOS)
            img.save(src_path)
        return src_path.stat().st_size <= max_bytes

    def UpdateLinkedFigures(fig_map, output_dir=None):
        out = Path(output_dir) if output_dir else LINKED_FIGURES_DIR
        out.mkdir(parents=True, exist_ok=True)
        for stem, src in fig_map.items():
            src_path = Path(src)
            if src_path.is_file():
                dest = out / f"{stem}.png"
                shutil.copyfile(src_path, dest)
                _ensure_size_under_mb(dest, max_mb=1)
                # Embed DPI so Google Docs / Preview display at correct size
                _img = _PIL_Image.open(dest)
                _img.save(str(dest), dpi=(FIGURE_DPI, FIGURE_DPI))
                print(f"  Linked: {dest}")
            else:
                print(f"  Missing: {src_path}")
        return out

    BASE = str(Path('..').resolve())
    # ── Manuscript outputs: {published filename stem} -> {source image} ──────
    # The key IS the output filename, so the manuscript's numbering lives in one
    # place. Renumbering a figure means editing a key here and nothing else.
    #
    # Manuscript Figure 1 (the VCH task schematic) is hand-drawn in Illustrator
    # and has no generating script, so it is not listed; it is placed in
    # results/final_figures/figure_1.svg by hand.
    #
    # NOTE the one confusing name: the SOURCE of Figure 7 is an intermediate
    # this script assembles as figures/.../figure_6.png. That internal name is
    # unrelated to the published figure_6 (vch_computations) below.
    manuscript_outputs = {
        "figure_2": f"{RESULTS_BASE}/descriptive/ppa_history_distributions.png",
        "figure_3": f"{RESULTS_BASE}/descriptive/caps_item_distributions_hppd_split.png",
        "figure_4": str(Path(BASE) / FIGURES_OUT_DIR / "sp_predictors_hppd_vs_caps.png"),
        "figure_5": str(Path(BASE) / FIGURES_OUT_DIR / "vch_behavior_hppd_vs_caps.png"),
        # Figure 6 has the mediation panels (l-o) appended below the HGF composite.
        "figure_6": str(Path(BASE) / FIGURES_OUT_DIR / "vch_computations_hppd_vs_caps.png"),
        # Figure 7: VCH beta / SDT / detection curves (2x2 quadrant)
        #   Row 1 (a, b): detection curves, VCH x SDT correlation grid
        #   Row 2 (c, d): HPPD sdt_hppd forest plot, CAPS sdt_hppd forest plot
        "figure_7": str(Path(BASE) / FIGURES_OUT_DIR / "figure_6.png"),
        "table_1":  f"{RESULTS_BASE}/descriptive/tables/clinical_table_hppd_split_caps_split.png",
        # Routed to results/supplement/ rather than final_figures/ — see
        # _SUPPLEMENT_FIGS below.
        "supplementary_figure_s6": f"{RESULTS_BASE}/descriptive/sp_table_distributions.png",
    }
    UpdateLinkedFigures(manuscript_outputs)
    print("Linked figures updated.")

    # ── Final figures (journal-ready, full-resolution) ───────────────────
    # Saved to results/final_figures/ in every FINAL_FIGURE_FORMATS entry.
    # These are NOT downscaled (unlike linked_figures/ PNGs which are capped
    # at 1 MB for Google Docs).
    from figure_assembly import raster_to_svg as _raster_to_svg
    FINAL_FIGURES_DIR = Path(f"{RESULTS_BASE}/final_figures")
    FINAL_FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    # SP use characteristics is a SUPPLEMENTARY figure, not a main-text one, so
    # it is written to results/supplement/ alongside the other supplementary
    # figures rather than to results/final_figures/.  Everything else in
    # linked_figures goes to final_figures.  Keyed on the linked_figures key.
    SUPPLEMENT_FIGURES_DIR = Path(f"{RESULTS_BASE}/supplement")
    _SUPPLEMENT_FIGS = {"supplementary_figure_s6"}
    # Journal SVG artboard widths: 89mm single-column, 180mm double-column
    _SVG_WIDTH_MM = 89       # default: single-column
    _SVG_DOUBLE_FIGS = {"figure_6"}   # vch_computations is double-column (180mm)
    print(f"\n--- Saving final figures ({', '.join(f.upper() for f in FINAL_FIGURE_FORMATS)}) ---")
    for stem, src in manuscript_outputs.items():
        src_path = Path(src)
        if not src_path.is_file():
            print(f"  Missing: {src_path}")
            continue
        _out_dir = (SUPPLEMENT_FIGURES_DIR if stem in _SUPPLEMENT_FIGS
                    else FINAL_FIGURES_DIR)
        _out_dir.mkdir(parents=True, exist_ok=True)
        for _fmt in FINAL_FIGURE_FORMATS:
            dest = _out_dir / f"{stem}.{_fmt}"
            if _fmt == "svg":
                # Prefer true-vector SVG (assembled by svgutils) over raster-wrapped
                vector_svg = src_path.with_suffix('.svg')
                if vector_svg.exists():
                    shutil.copyfile(vector_svg, dest)
                else:
                    _raster_to_svg(str(src_path), str(dest), dpi=FIGURE_DPI)
                # Set artboard to journal column width (aspect ratio preserved)
                from figure_assembly import scale_svg_to_dpi
                _w_mm = 180 if stem in _SVG_DOUBLE_FIGS else _SVG_WIDTH_MM
                scale_svg_to_dpi(dest, target_width_mm=_w_mm, top_margin_mm=5)
            elif _fmt in ("tiff", "tif"):
                img = _PIL_Image.open(src_path).convert('RGB')
                img.save(str(dest), format='TIFF', dpi=(FIGURE_DPI, FIGURE_DPI),
                         compression='tiff_lzw')
            else:
                img = _PIL_Image.open(src_path).convert('RGB')
                img.save(str(dest), dpi=(FIGURE_DPI, FIGURE_DPI))
            print(f"  Saved: {dest}  ({dest.stat().st_size / 1024 / 1024:.1f} MB)")
    print(f"Final figures ({', '.join(f.upper() for f in FINAL_FIGURE_FORMATS)}) saved.")

print("\nAll done.")
