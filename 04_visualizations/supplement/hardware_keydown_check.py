#!/usr/bin/env /usr/local/bin/python3.12
"""
hardware_keydown_check.py
=========================
Reviewer-response figure: does the participant's DISPLAY covary with VCH task
performance?

Scientific purpose
------------------
A reviewer asked whether display hardware, uncontrolled in an unsupervised
online task, could account for the visual-task measures we report.  The task
captured hardware as a REDCap FREE-TEXT field (`monitor_check`), which cannot
be tested directly, so it is hand-coded into a 3-level DISPLAY CLASS,
`monitor_check_operationalized_final`:

    Mac                   integrated Apple displays, which are consistently
                          LED-backlit IPS LCD panels
    External Monitor      a display separate from the machine driving it, and
                          therefore the class with the most variable viewing
                          distance, panel type and setup
    Windows/Other Laptop  everything else

This is the only grouping of the free text that corresponds to a physical
property plausibly able to produce a display difference, and it is the one
coding with enough participants in every cell to test.  It is also the coding
carried into the models as the hardware-control covariate (model type
`nice_covariates_spusers_hardware_control`), so the figure and the covariate
describe the same variable.

Analysis rule
-------------
Kruskal-Wallis across the three display classes, per DV.  Pairwise Mann-Whitney
post-hoc tests are computed for every pair and written to the summary CSV with
raw and Bonferroni-adjusted p-values (adjusted within panel).  Post-hoc
brackets are drawn on the figure only when the omnibus test is significant.
Whether a display difference is carried into the reported analyses is decided
on that omnibus test and its post-hoc contrasts.

    row a   d_prime_overall  x  monitor_check_operationalized_final
    row b   vch_threshold    x  monitor_check_operationalized_final

Both rows share the grouping variable, so only panel b carries x-tick labels
and the x-axis title (xaxis_bottom_row_only=True).  Panel letters are drawn by
figure_assembly._add_panel_label -- the helper every assembled manuscript
figure uses -- so they cannot drift from the rest of the paper's.

Fonts: Arial (declared as a sans-serif stack) at the manuscript target point
sizes from 0X_all_figures.py; see the FONT_* constants in CONFIG.  The figure is
drawn at its final printed size (half a 183 mm page), so those targets are
applied directly rather than back-computed through compute_source_fontsize().

The display class
-----------------
`monitor_check_operationalized_final` is a hand coding of the free-text
`monitor_check` field. It ships as a column of `data/final/df_public_*.csv` and
is read from there -- this script does not derive it, and neither does
`03_hpc/generate_hpc_jobs.py`, which reads the same column for the
`nice_covariates_hardware_control` covariate set. Figure and covariate are
therefore the same variable by construction, not by two matching derivations.

Three levels: Mac, Windows/Other Laptop, External Monitor. See the CONFIG block
for the rule that assigns them.

Known data caveats (NOT corrected here -- the shipped coding is the analyst's
ground truth and this script does not alter it)
------------------------------------------------------------------------------
1. NON-RESPONSE IS NOT A HARDWARE CATEGORY.  15 people left the hardware
   question blank and are dropped (DROP_NON_RESPONDERS).  Non-response is
   almost perfectly confounded with VCH exclusion: 14 of the 15 already have
   missing d'/threshold, so only ONE non-responder is in the analysis sample.
   Dropping changes per-panel N from 194 to 193.  People who skipped the
   hardware question are essentially the people who did not produce usable task
   data -- consistent with partial-completion dropout, not with anything about
   hardware.  Set DROP_NON_RESPONDERS = False to test the coding as it stands.

   The exclusion is keyed on the RAW free-text field being blank, not on the
   "No Response" label, so it does not depend on how non-response was coded.

2. SINGLE-PARTICIPANT SENSITIVITY.  With only one non-responder in the analysis
   sample, results at the margin can turn on that participant.  Per-panel N and
   per-group n are written to the summary CSV and printed above every box, so
   the caption can state the cell sizes honestly.

Sample
------
228 rows in the coding CSV.  Two exclusions, in order:
  - non-responders on the hardware field (caveat 1);
  - listwise deletion on the DV (the 34 VCH-excluded participants; both DVs
    share exactly the same missing rows -- asserted at runtime).
The two overlap almost completely, so every panel ends at N = 193.  No SP-user
filter is applied: the reviewer's question is about the measurement device,
which is independent of drug exposure.

Reads
-----
  data/final/df_public_*.csv  (most recent, via most_recent_public_df)
      Columns used: monitor_check, monitor_check_operationalized_final,
      d_prime_overall, vch_threshold.

  modules/master_config.py — VARIABLE_REGISTRY (y-axis labels),
      build_linear_palette + ELECTRIC_BLUE_SPEC (house palette).

Outputs (results/supplement/hardware_keydown_check/)
----------------------------------------------------
  supplementary_figure_s8.png    manuscript drafting / Docs preview
  supplementary_figure_s8.tiff   journal submission, LZW, same DPI
  supplementary_figure_s8.svg    true-vector figure assembly
  summary_results/supplementary_figure_s8.csv                  this scheme
  summary_results/hardware_operationalization_all_schemes.csv  every active
                                                               panel, stacked

Each panel is annotated with its omnibus statistic, p, effect size, N, AND the
per-group n above every box (show_group_n), so the collapsibility of thin
levels is readable off the figure itself.

How to run
----------
  cd hppd_manuscript_public
  /usr/local/bin/python3.12 04_visualizations/supplement/hardware_keydown_check.py
"""

import re
import sys
import warnings
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal, mannwhitneyu

warnings.filterwarnings('ignore')

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT  = SCRIPT_DIR.parents[1]          # supplement/ -> 04_visualizations/ -> root

sys.path.insert(0, str(REPO_ROOT / 'modules'))
from master_config import (                                    # noqa: E402
    VARIABLE_REGISTRY, ELECTRIC_BLUE_SPEC, build_linear_palette,
)
# The canonical panel-label renderer used by every assembled manuscript figure
# (bold, Arial, black, upper-left of the axes).  Imported rather than
# reimplemented so these labels cannot drift from the rest of the paper's.
from figure_assembly import _add_panel_label                    # noqa: E402
from data_prep import most_recent_public_df                     # noqa: E402

# ══════════════════════════════════════════════════════════════════════════════
# CONFIG
# ══════════════════════════════════════════════════════════════════════════════

DATA_DIR = REPO_ROOT / 'data' / 'final'
OUT_DIR  = REPO_ROOT / 'results' / 'supplement' / 'hardware_keydown_check'

# ── The display class ─────────────────────────────────────────────────────────
# Read straight off the analysis dataframe. It is a hand coding of the free-text
# `monitor_check` field and is NOT derived here.
#
# Three levels, chosen so every cell is large enough to test:
#   Mac                   integrated Apple displays, iMac desktops included
#   Windows/Other Laptop  everything else with a built-in display
#   External Monitor      any display separate from the machine driving it,
#                         including a TV used as one
#
# One rule decides the level: it describes THE DISPLAY, not the machine. A
# response naming a separate external display codes as External Monitor whatever
# drives it. Non-responders are null, not a fourth level.
#
# The same column is the covariate behind the HPC model type
# `nice_covariates_spusers_hardware_control`, so the figure and the sensitivity
# analysis describe exactly the same variable. To change the coding, change it
# where the dataframe is built and re-export; every consumer reads the column.
GROUP_COL = 'monitor_check_operationalized_final'
RAW_FIELD = 'monitor_check'

# Row layout: row a = d', row b = threshold (DV_COLS order IS the row order).
DV_COLS = ['d_prime_overall', 'vch_threshold']

# x-axis title.
GROUP_LABELS = {GROUP_COL: 'Display Hardware Reported'}

# Exclude participants who left the hardware question blank: non-response is not
# a hardware category, so these people are dropped rather than tested as a level.
# Keyed on the RAW free-text field being blank, so it does not depend on how
# non-response was coded.
DROP_NON_RESPONDERS = True

# Draw pairwise post-hoc brackets only for panels with at most this many
# groups.  Above it the brackets overlap into illegibility; the tests are still
# written to the summary CSV.
PAIRWISE_LINE_MAX_GROUPS = 5

DPI = 300

# ── Publication font rules ───────────────────────────────────────────────────
# Arial everywhere, per 04_visualizations/README.md ("Font requirement: Arial").
# The family is declared as a sans-serif STACK rather than font.family='Arial'
# so a machine without Arial degrades predictably instead of silently falling
# back to DejaVu.  This matches 0X_all_figures.py exactly.  Nothing in this
# script calls sns.set_style(), which is the one documented thing that resets
# font.family (see the README gotcha) -- if one is ever added, re-set these two
# lines immediately afterwards.
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

# Point sizes are the manuscript targets defined at the top of
# 0X_all_figures.py: TARGET_AXIS_LABEL (9.2), TARGET_TICK_LABEL (7.0),
# TARGET_SIG_MARKER (24.0), FIGURE_LABEL_FONTSIZE (20).  Those are *apparent*
# sizes in the final 183 mm assembled figure.  This figure is drawn at its
# final printed size and is never rescaled by modules/figure_assembly.py, so
# the targets apply directly -- no compute_source_fontsize() back-calculation.
FONT_AXIS_LABEL  = 9.2     # y-axis label, x-axis label
FONT_TICK_LABEL  = 7.0     # tick labels on both axes
FONT_SIG_MARKER  = 24.0    # post-hoc asterisks (drawn at 0.5x, as before)
FONT_ANNOT       = 7.0     # omnibus stat block: small-text tier, = tick size
FONT_GROUP_N     = 7.0     # the "n = ..." row: same tier
FONT_PANEL_LABEL = 20.0    # bold a/b, = FIGURE_LABEL_FONTSIZE in 0X_all_figures

# ── Journal sizing -- Scientific Reports (Nature) ────────────────────────────
# Mirrors JOURNAL_DOUBLE_COL_MM in 0X_all_figures.py (183 mm = full page /
# double column).  This figure is drawn at HALF that width (analyst
# instruction 2026-08-31): the two panels share one grouping variable and are
# stacked, so a full-width figure would waste most of the page.
JOURNAL_DOUBLE_COL_MM = 183
_MM_PER_INCH = 25.4
FIG_WIDTH_IN = (JOURNAL_DOUBLE_COL_MM / 2) / _MM_PER_INCH   # 91.5 mm ~ 3.602"

# Height of ONE panel row, in inches.  2 rows x 2.6" = 5.2" tall at 3.602"
# wide.  Panel a carries no x-tick labels or x-label, so its plotting box is
# taller than panel b's at the same row height.
PANEL_HEIGHT_IN = 2.6


# ══════════════════════════════════════════════════════════════════════════════
# THE FIGURE FUNCTION
# ══════════════════════════════════════════════════════════════════════════════

def hardware_boxplot_grid(
    data,
    panel_specs,
    nrows=2,
    ncols=2,
    panel_height=PANEL_HEIGHT_IN,
    fig_width=None,
    min_panel_width=7.0,
    max_panel_width=13.0,
    width_per_group=0.85,
    strip_alpha=0.45,
    # Point-unit geometry (marker sizes, line widths) does NOT scale with the
    # figure, so the values inherited from the 14"-wide exploratory layout bury
    # the boxes at journal width.  Halved for FIG_WIDTH_IN; the wide schemes
    # want the originals back (strip 4, line 3, marker 14).
    strip_size=2.5,
    box_alpha=0.2,
    box_width=0.45,
    line_width=1.5,
    mean_marker='o',
    mean_marker_size=7,
    ylabel_fontsize=FONT_AXIS_LABEL,
    xlabel_fontsize=FONT_AXIS_LABEL,
    annot_fontsize=FONT_ANNOT,
    sig_marker_fontsize=FONT_SIG_MARKER,
    show_outliers=False,
    showcaps=False,
    pairwise_line_max_groups=PAIRWISE_LINE_MAX_GROUPS,
    pairwise_line_height_increment=0.11,
    annotate_stats=True,
    show_group_n=True,
    group_n_fontsize=FONT_GROUP_N,
    tick_fontsize=FONT_TICK_LABEL,
    panel_labels=None,
    panel_label_fontsize=FONT_PANEL_LABEL,
    xaxis_bottom_row_only=False,
    bbox_tight=False,
    savepath=None,
    dpi=DPI,
):
    """
    A x B grid of boxplots with a nonparametric group test per panel.

    Adapted from `modules/visualization_helpers_parts/multipanel_boxplot_grid.py`
    and visually identical to it (translucent box, per-group strip plot, mean
    marker, manually drawn quartile rules, no spines, no ticks).  It differs in
    five ways that the hardware figures require:

    1. EVERY PANEL KEEPS ITS OWN Y-AXIS AND Y-LABEL.  The parent function hides
       the y-axis for every column but the first, which assumes all panels in a
       row share a DV.  Here the DV varies BY COLUMN (d' vs. threshold), so
       hiding column 1's axis would render it unreadable.
    2. EACH PANEL CAN KEEP ITS OWN X-LABEL.  The parent labels only the bottom
       row, which assumes a shared grouping variable.  When the grouping
       variable varies BY ROW (monitor vs. headphones) every row needs its own
       label, so that is the default here; pass xaxis_bottom_row_only=True for
       the stacked single-grouping layout, where repeating the x-axis on every
       row wastes vertical space the journal width cannot spare.
    3. X-TICK LABELS ROTATE AND THE FIGURE WIDENS WITH GROUP COUNT.  Scheme 1
       has up to 17 levels with long names; the parent's horizontal, space-
       wrapped labels collide well before that.
    4. EFFECT SIZES AND POST-HOC TESTS ALWAYS REACH THE CSV.  epsilon^2 for
       Kruskal-Wallis, rank-biserial for Mann-Whitney, and every pairwise
       comparison with a Bonferroni-adjusted p-value -- drawn on the figure
       only when the panel is small enough to stay legible.
    5. IT WRITES .png, .tiff AND .svg (the parent writes .png and .svg only).
       The .tiff is the supplement's journal-submission requirement.

    Parameters
    ----------
    data : DataFrame
        Source frame.  Never modified; each panel takes its own listwise-clean
        copy.
    panel_specs : list of dict, laid out row-major (index // ncols = row)
        'dv'          : DV column name (required)
        'group_var'   : grouping column name (required)
        'order'       : explicit level order; default = descending group n
        'ylabel'      : y-axis label; default = the DV column name
        'xlabel'      : x-axis label; default = the group column name
        'force_test'  : 'kruskal' | 'mannwhitney' | None (default: kruskal for
                        3+ groups, mannwhitney for exactly 2)
    fig_width : float or None
        TOTAL figure width in inches.  When given it overrides the group-count-
        driven width below, which is what pins the figure to a journal column
        measure (FIG_WIDTH_IN).  None restores the auto-widening behaviour
        needed by the wide exploratory schemes.
    panel_labels : list of str or None
        Bold panel letters ('a', 'b', ...), row-major, one per panel spec.
        Rendered by figure_assembly._add_panel_label so they are identical to
        the labels on every assembled manuscript figure.
    xaxis_bottom_row_only : bool
        True suppresses x-tick labels and the x-label on every row but the
        last occupied one.  Only correct when all rows share one grouping
        variable (see difference 2 above).
    bbox_tight : bool
        False (default) saves at exactly `figsize`, so a figure built to a
        journal measure comes out at that measure.  True restores
        bbox_inches='tight', which crops to the artists and therefore does NOT
        preserve the requested width.
    savepath : str or Path
        Base path.  Writes '<savepath>.png', '.tiff', '.svg' and
        '<dir>/summary_results/<base>.csv'.

    Returns
    -------
    (results, summary_df)
        results    : dict keyed '<dv>_<group_var>' with the full per-panel stats
        summary_df : tidy long-format DataFrame, one row per omnibus test and
                     one per pairwise comparison
    """
    n_panels = len(panel_specs)
    if n_panels == 0:
        raise ValueError('panel_specs must contain at least one panel specification')
    if n_panels > nrows * ncols:
        raise ValueError(
            f'panel_specs has {n_panels} panels but grid only has {nrows * ncols} positions'
        )

    if fig_width is not None:
        # Fixed journal measure: the caller has decided the printed width, so
        # the group count must not be allowed to change it.
        panel_width = fig_width / ncols
    else:
        # Figure width is driven by the widest panel: enough horizontal room for
        # the largest group count anywhere in the grid, so all panels stay
        # aligned.
        max_groups = 0
        for spec in panel_specs:
            clean = data[[spec['dv'], spec['group_var']]].dropna()
            max_groups = max(max_groups, clean[spec['group_var']].nunique())
        # Clamped: scheme 1 has 17 levels, and an uncapped per-group width makes
        # a 32-inch-wide figure that no supplement page can use.  The cap trades
        # a little x-tick breathing room for a usable aspect ratio.
        panel_width = min(max_panel_width,
                          max(min_panel_width, width_per_group * max_groups))

    # Index of the last row that actually holds a panel (the grid may be only
    # partially filled), used by xaxis_bottom_row_only.
    last_occupied_row = (len(panel_specs) - 1) // ncols

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(panel_width * ncols, panel_height * nrows),
        squeeze=False,
    )

    results = {}

    for panel_idx, spec in enumerate(panel_specs):
        row, col = divmod(panel_idx, ncols)
        ax = axes[row, col]

        dv         = spec['dv']
        group_var  = spec['group_var']
        ylabel     = spec.get('ylabel', dv)
        xlabel     = spec.get('xlabel', group_var)
        force_test = spec.get('force_test', None)

        # Listwise deletion on the DV and the grouping column only.  Nothing is
        # imputed, recoded or reordered -- the group labels are used verbatim.
        clean = data[[dv, group_var]].dropna(subset=[dv, group_var]).copy()
        if len(clean) == 0:
            ax.text(0.5, 0.5, f'No data for {dv}', ha='center', va='center',
                    transform=ax.transAxes)
            continue

        # Default level order: descending n.  Purely a display choice -- both
        # Kruskal-Wallis and Mann-Whitney are invariant to group order -- chosen
        # so the largest, most interpretable cells sit on the left.
        if spec.get('order') is not None:
            groups = list(spec['order'])
        else:
            groups = list(clean[group_var].value_counts().index)
        n_groups = len(groups)

        # Dark -> electric-blue house gradient, one step per level.
        palette = build_linear_palette(
            ELECTRIC_BLUE_SPEC['vibrant_rgb'],
            dark_rgb=ELECTRIC_BLUE_SPEC['dark_rgb'],
            n_levels=max(n_groups, 2),
        )[:n_groups]

        # seaborn 0.13.2 needs a str-keyed palette dict and a str hue column
        # (see the note in multipanel_boxplot_grid: a list palette trips an
        # UnboundLocalError on 'boxprops' inside _configure_legend).
        hue_col = f'__hue_{group_var}'
        clean[hue_col] = clean[group_var].astype(str)
        groups = [str(g) for g in groups]
        pal_dict = dict(zip(groups, palette))

        # ── boxes ────────────────────────────────────────────────────────────
        sns.boxplot(
            data=clean, x=hue_col, y=dv, hue=hue_col,
            palette=pal_dict, order=groups, ax=ax,
            showfliers=show_outliers, showcaps=showcaps, width=box_width,
            whiskerprops={'linewidth': line_width},
            capprops={'linewidth': line_width},
            medianprops={'linewidth': line_width},
        )
        if ax.get_legend() is not None:
            ax.get_legend().remove()
        for patch in ax.patches:
            patch.set_alpha(box_alpha)
            patch.set_edgecolor('none')

        # ── points ───────────────────────────────────────────────────────────
        for i, group in enumerate(groups):
            gdf = clean[clean[hue_col] == group]
            sns.stripplot(
                data=gdf, x=hue_col, y=dv, order=groups,
                color=palette[i], alpha=strip_alpha, size=strip_size,
                jitter=0.2, ax=ax,
            )

        # ── mean markers + quartile rules ────────────────────────────────────
        half_width = box_width / 2
        for i, group in enumerate(groups):
            vals = clean[clean[hue_col] == group][dv]
            ax.plot(i, vals.mean(), marker=mean_marker, markersize=mean_marker_size,
                    color=palette[i], markeredgecolor=palette[i],
                    markeredgewidth=line_width / 2, zorder=3, alpha=1)
            ax.hlines([vals.quantile(0.25), vals.quantile(0.75)],
                      i - half_width, i + half_width,
                      colors=palette[i], linewidth=line_width, alpha=1, zorder=3)

        # Recolour whiskers/median lines to match their box.
        for line in ax.lines:
            xdata = line.get_xdata()
            if xdata is None or len(xdata) == 0:
                continue
            idx = int(np.clip(round(float(np.mean(xdata))), 0, n_groups - 1))
            line.set_color(palette[idx])
            line.set_alpha(1)

        # ── the test ─────────────────────────────────────────────────────────
        group_values = [clean[clean[hue_col] == g][dv].values for g in groups]
        group_ns     = [len(v) for v in group_values]

        # Post-hoc brackets and the stats annotation both want vertical space at
        # the top of the panel.  The test block records how many bracket levels
        # it needs; the y-limits are then set ONCE, afterwards, so the brackets
        # and the annotation cannot overlap each other (they did when each stage
        # grew the limits independently).
        sig_pairs = []
        pairwise_bracket_levels = 0

        if force_test is None:
            test = 'mann_whitney' if n_groups == 2 else 'kruskal_wallis'
        else:
            test = {'mannwhitney': 'mann_whitney', 'kruskal': 'kruskal_wallis'}[force_test]

        panel_key = f'{dv}_{group_var}'
        base = {
            'panel': panel_key, 'dv': dv, 'group_var': group_var,
            'test': test, 'n_groups': n_groups, 'sample_size': len(clean),
            'groups': '|'.join(groups),
            'group_ns': '|'.join(str(n) for n in group_ns),
        }
        results[panel_key] = dict(base)
        results[panel_key]['pairwise'] = {}

        if test == 'mann_whitney':
            if n_groups != 2:
                raise ValueError(
                    f'{panel_key}: Mann-Whitney requested but {n_groups} groups present'
                )
            v0, v1 = group_values
            stat, p_value = mannwhitneyu(v0, v1, alternative='two-sided')
            denom = len(v0) * len(v1)
            # Rank-biserial: +1 = every value in group 1 exceeds every value in
            # group 2; -1 = the reverse.  Same formula as multipanel_boxplot_grid.
            r_rb = ((2 * stat) / denom - 1) if denom > 0 else np.nan
            results[panel_key].update({
                'statistic': stat, 'p_value': p_value,
                'effect_size_name': 'rank_biserial_r',
                'effect_size': r_rb,
                'group_1_label': groups[0], 'group_1_n': len(v0),
                'group_1_median': float(np.median(v0)),
                'group_1_iqr_q1': float(np.quantile(v0, 0.25)),
                'group_1_iqr_q3': float(np.quantile(v0, 0.75)),
                'group_2_label': groups[1], 'group_2_n': len(v1),
                'group_2_median': float(np.median(v1)),
                'group_2_iqr_q1': float(np.quantile(v1, 0.25)),
                'group_2_iqr_q3': float(np.quantile(v1, 0.75)),
            })
            omnibus_p = p_value
            stat_label = f'U = {stat:,.0f}'
            es_label = f'r$_{{rb}}$ = {r_rb:.2f}'
        else:
            h_stat, p_value = kruskal(*group_values)
            n_total = len(clean)
            # epsilon^2 = H / (N - 1): proportion of rank variance explained,
            # bounded [0, 1].  Reported because a NULL result is the point here,
            # and a bare p-value cannot express "and the effect is tiny".
            eps2 = h_stat / (n_total - 1) if n_total > 1 else np.nan
            results[panel_key].update({
                'statistic': h_stat, 'p_value': p_value,
                'effect_size_name': 'epsilon_squared',
                'effect_size': eps2,
                'df': n_groups - 1,
            })
            omnibus_p = p_value
            stat_label = f'H({n_groups - 1}) = {h_stat:.2f}'
            es_label = f'$\\epsilon^2$ = {eps2:.3f}'

            # Post-hoc pairwise Mann-Whitney, every pair, always into the CSV.
            pairs = list(combinations(range(n_groups), 2))
            n_comparisons = len(pairs)
            for i, j in pairs:
                u_stat, p_pair = mannwhitneyu(
                    group_values[i], group_values[j], alternative='two-sided'
                )
                p_bonf = min(1.0, p_pair * n_comparisons)
                results[panel_key]['pairwise'][f'{groups[i]} vs {groups[j]}'] = {
                    'group_1': groups[i], 'group_2': groups[j],
                    'group_1_n': group_ns[i], 'group_2_n': group_ns[j],
                    'statistic': u_stat, 'p_value': p_pair,
                    'p_value_bonferroni': p_bonf,
                    'n_comparisons': n_comparisons,
                }
                if p_pair < 0.05:
                    sig_pairs.append((i, j, p_pair))

            # Brackets only when the omnibus test is significant AND the panel
            # is small enough to draw them legibly.
            if omnibus_p < 0.05 and sig_pairs and n_groups <= pairwise_line_max_groups:
                pairwise_bracket_levels = len(sig_pairs)

        # ── vertical layout: data, then brackets, then annotation ────────────
        # Budget the y-axis in data units so the three tiers never collide.
        # The bracket stack is anchored to the observed data maximum (NOT to
        # ax.get_ylim(), which seaborn has already padded), and the annotation
        # band is reserved above whatever the brackets end up needing.
        data_lo = float(clean[dv].min())
        data_hi = float(clean[dv].max())
        span = (data_hi - data_lo) or 1.0

        bracket_gap = pairwise_line_height_increment * span   # first bracket offset
        bracket_step = pairwise_line_height_increment * span  # spacing per level
        if pairwise_bracket_levels:
            # + a glyph allowance: the asterisks are drawn va='bottom' ON the
            # top bracket, so the stack is taller than its last line.
            bracket_top = (data_hi + bracket_gap
                           + bracket_step * (pairwise_bracket_levels - 1)
                           + 0.07 * span)
        else:
            bracket_top = data_hi
        # Three stacked bands above the data, bottom-to-top:
        #   brackets  ->  the per-group "n = ..." row  ->  the stats annotation.
        # Each is reserved explicitly so none can overlap another.
        n_band     = 0.11 * span if show_group_n else 0.0
        annot_band = 0.20 * span if annotate_stats else 0.03 * span

        n_row_y = bracket_top + 0.03 * span
        ax.set_ylim(data_lo - 0.05 * span, bracket_top + n_band + annot_band)

        for level, (i, j, p_pair) in enumerate(sig_pairs[:pairwise_bracket_levels]):
            y_line = data_hi + bracket_gap + bracket_step * level
            ax.plot([i, j], [y_line, y_line], color='black',
                    linewidth=line_width / 2, linestyle='--', zorder=10, alpha=0.6)
            ax.text((i + j) / 2, y_line, _stars(p_pair), ha='center', va='bottom',
                    fontsize=int(sig_marker_fontsize * 0.5), fontweight='bold')

        # Per-group N across the top of the panel, above each box.  Group
        # sizes are the first thing anyone asks of a free-text-derived
        # grouping, and they decide which levels are collapsible -- so they
        # belong on the figure, not only in the summary CSV.
        if show_group_n:
            n_fs = group_n_fontsize if group_n_fontsize is not None else \
                max(9, 16 - 0.4 * n_groups)
            for i, n_i in enumerate(group_ns):
                ax.text(i, n_row_y, f'n = {n_i}', ha='center', va='bottom',
                        fontsize=n_fs, color='#313a3d')

        # The omnibus statistic is printed in full rather than reduced to an
        # asterisk: these panels are expected to be null, and "H(14) = 12.31,
        # p = .58" is the reviewer-facing answer, whereas blank space is not.
        if annotate_stats:
            p_txt = 'p < .001' if omnibus_p < 0.001 else f'p = {omnibus_p:.3f}'.replace('0.', '.')
            stars = _stars(omnibus_p)
            ax.text(
                0.99, 0.99,
                f'{stat_label}, {p_txt}{"  " + stars if stars else ""}\n'
                f'{es_label},  N = {len(clean)}',
                transform=ax.transAxes, ha='right', va='top',
                fontsize=annot_fontsize, color='#313a3d',
            )

        # ── axes ─────────────────────────────────────────────────────────────
        # Both the y-label AND the y-axis stay on every panel: the DV changes
        # across panels, so the scales are not comparable and cannot be shared.
        # labelpad is in points, so it is scaled to the journal font sizes --
        # the old 10/12 were sized for 20 pt labels and leave a visible gap at
        # 9.2 pt.
        ax.set_ylabel(_wrap_label(ylabel, width=24),
                      fontsize=ylabel_fontsize, labelpad=4)

        # Whether this panel shows the x-axis at all.  With one shared grouping
        # variable only the bottom panel needs tick labels and an x-label;
        # repeating them on every row would say the same thing twice and cost
        # vertical space.
        show_x_axis = (not xaxis_bottom_row_only) or row == last_occupied_row

        # Long free-text-derived level names: wrap, then rotate.  The threshold
        # is 3 groups, not 5: level names here run to 27 characters
        # ("Over-Ear/On-Ear Wired Headphones"), and horizontal labels overrun
        # their slot well before the panel is full.
        wrapped = [_wrap_label(g, width=16) for g in groups]
        rotation = 45 if n_groups > 3 else 0
        ax.set_xticks(range(n_groups))
        if show_x_axis:
            ax.set_xticklabels(wrapped, rotation=rotation,
                               ha='right' if rotation else 'center')
            ax.set_xlabel(xlabel, fontsize=xlabel_fontsize, labelpad=6)
        else:
            ax.set_xticklabels([])
            ax.set_xlabel('')

        # tick_fontsize=None restores the old behaviour: shrink tick text as the
        # panel fills up, with a readability floor.  That existed for the wide
        # 11-17-level exploratory schemes; a fixed journal size is used instead
        # whenever the caller supplies one.
        if tick_fontsize is None:
            xtick_fs, ytick_fs = max(9, 17 - 0.45 * n_groups), 15
        else:
            xtick_fs = ytick_fs = tick_fontsize
        ax.tick_params(axis='x', labelsize=xtick_fs)
        ax.tick_params(axis='y', labelsize=ytick_fs)

        # Bold panel letter, upper-left, drawn by the same helper every
        # assembled manuscript figure uses (Arial bold, black, va='top').
        if panel_labels:
            _add_panel_label(ax, panel_labels[panel_idx], panel_label_fontsize,
                             'black', x_offset=0.01, y_offset=0.99)

        for side in ('top', 'right', 'bottom', 'left'):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis='both', which='both', length=0)

        print(f'  panel {panel_idx} [{dv} x {group_var}]  {test}  '
              f'{stat_label}, p = {omnibus_p:.4f}, N = {len(clean)}, k = {n_groups}')

    for panel_idx in range(n_panels, nrows * ncols):
        r, c = divmod(panel_idx, ncols)
        axes[r, c].axis('off')

    plt.tight_layout()

    summary_df = _results_to_long(results)

    if savepath:
        savepath = Path(savepath)
        savepath.parent.mkdir(parents=True, exist_ok=True)
        base_path = savepath.with_suffix('') if savepath.suffix else savepath
        # The scheme-1 figure is ~26 x 13 in; an uncompressed 300-dpi TIFF of it
        # is 120 MB.  LZW is lossless, universally readable, and accepted by
        # Scientific Reports, and brings that to a few MB.
        save_formats = (
            ('png',  {'dpi': dpi}),
            ('tiff', {'dpi': dpi, 'pil_kwargs': {'compression': 'tiff_lzw'}}),
            ('svg',  {}),
        )
        # bbox_inches='tight' crops to the artists, so the saved width is NOT
        # the figsize width.  This figure is built to a journal measure
        # (FIG_WIDTH_IN), so it is saved uncropped and tight_layout() above is
        # what keeps the labels inside the canvas.
        bbox = 'tight' if bbox_tight else None
        for ext, kw in save_formats:
            out = base_path.with_suffix(f'.{ext}')
            fig.savefig(out, format=ext, bbox_inches=bbox, facecolor='white', **kw)
            print(f'  saved -> {out}')

        summary_dir = base_path.parent / 'summary_results'
        summary_dir.mkdir(parents=True, exist_ok=True)
        summary_csv = summary_dir / f'{base_path.name}.csv'
        summary_df.to_csv(summary_csv, index=False)
        print(f'  saved -> {summary_csv}')

    plt.close(fig)
    return results, summary_df





def _wrap_label(label, width=16):
    """
    Wrap a hand-coded level name for use as an x-tick label.

    Breaks on spaces AND after '/', because these names come from free text and
    the informative boundary is often a slash: plain textwrap turns
    "No Headphones/Computer Audio" into "No Headphones/Co" + "mputer Audio" and
    "Other/Uncategorized" into "Other/Uncategori" + "zed".  Words are never
    split mid-token, so a line may exceed `width` when a single token does.
    """
    tokens = re.split(r'(?<=/)|(?<= )', str(label))
    lines, current = [], ''
    for tok in tokens:
        if current and len(current) + len(tok.rstrip()) > width:
            lines.append(current.rstrip())
            current = tok
        else:
            current += tok
    if current.rstrip():
        lines.append(current.rstrip())
    return '\n'.join(lines) if lines else str(label)


def _stars(p):
    """Conventional significance markers; '' when p >= .10 (matches the parent helper)."""
    if p < 0.001:
        return '***'
    if p < 0.01:
        return '**'
    if p < 0.05:
        return '*'
    if p < 0.1:
        return '~'
    return ''


def _results_to_long(results):
    """Flatten the per-panel results dict into one tidy row per test."""
    rows = []
    for panel_key, res in results.items():
        base = {k: res.get(k) for k in
                ('panel', 'dv', 'group_var', 'test', 'n_groups',
                 'sample_size', 'groups', 'group_ns')}
        overall = dict(base)
        overall.update({
            'row_type': 'overall',
            'comparison': '',
            'df': res.get('df'),
            'statistic': res.get('statistic'),
            'p_value': res.get('p_value'),
            'p_value_bonferroni': np.nan,   # omnibus tests are not corrected
            'effect_size_name': res.get('effect_size_name'),
            'effect_size': res.get('effect_size'),
            'group_1': res.get('group_1_label'),
            'group_2': res.get('group_2_label'),
            'group_1_n': res.get('group_1_n'),
            'group_2_n': res.get('group_2_n'),
            'group_1_median': res.get('group_1_median'),
            'group_1_iqr_q1': res.get('group_1_iqr_q1'),
            'group_1_iqr_q3': res.get('group_1_iqr_q3'),
            'group_2_median': res.get('group_2_median'),
            'group_2_iqr_q1': res.get('group_2_iqr_q1'),
            'group_2_iqr_q3': res.get('group_2_iqr_q3'),
        })
        rows.append(overall)

        for comparison, pw in res.get('pairwise', {}).items():
            row = dict(base)
            row.update({
                'row_type': 'pairwise_posthoc',
                'test': 'mann_whitney',
                'comparison': comparison,
                'df': np.nan,
                'statistic': pw['statistic'],
                'p_value': pw['p_value'],
                'p_value_bonferroni': pw['p_value_bonferroni'],
                'effect_size_name': '',
                'effect_size': np.nan,
                'group_1': pw['group_1'], 'group_2': pw['group_2'],
                'group_1_n': pw['group_1_n'], 'group_2_n': pw['group_2_n'],
            })
            rows.append(row)
    return pd.DataFrame(rows)


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    # The display class ships as a column of the analysis dataframe; it is not
    # derived here and not read from any side-car file.  See the DISPLAY CLASS
    # note above for where the coding is defined.
    df = pd.read_csv(most_recent_public_df(DATA_DIR), low_memory=False)
    print(f'  df shape: {df.shape}')

    # Integrity checks.  Any failure means the shipped dataframe is not the file
    # this script was written against -- stop rather than adapt silently.
    required = set(DV_COLS) | {GROUP_COL, RAW_FIELD}
    missing = sorted(required - set(df.columns))
    if missing:
        raise KeyError(
            f'The analysis dataframe is missing expected columns: {missing}. '
            f'{GROUP_COL!r} is a materialised column of df_public_*.csv; if it '
            'is absent the export needs regenerating, not this script patching.'
        )

    for dv in DV_COLS:
        if not pd.api.types.is_numeric_dtype(df[dv]):
            raise TypeError(f'{dv} is {df[dv].dtype}, expected numeric')
    if not df[DV_COLS[0]].isna().equals(df[DV_COLS[1]].isna()):
        raise ValueError(
            "d_prime_overall and vch_threshold no longer share the same missing "
            "rows; per-panel N would differ across columns. Investigate before "
            "regenerating."
        )
    n_dropped = int(df[DV_COLS[0]].isna().sum())
    print(f'  {n_dropped} rows have missing d\'/threshold (VCH-excluded) '
          f'-> N = {len(df) - n_dropped} before the non-responder exclusion')

    vc = df[GROUP_COL].value_counts(dropna=False)
    print(f'  {GROUP_COL}: {dict(vc)}')

    if DROP_NON_RESPONDERS:
        # Participants who left the hardware question blank are excluded: a
        # non-response is not a display class.  The shipped column already
        # carries NaN for them, so this is normally a no-op -- it is kept, and
        # asserted, so that a future export which labels non-responders instead
        # of nulling them cannot silently add a fourth group.
        df = df.copy()
        blank = df[RAW_FIELD].isna()
        n_blank_in_sample = int((blank & df[DV_COLS[0]].notna()).sum())
        already_null = int((blank & df[GROUP_COL].isna()).sum())
        if already_null != int(blank.sum()):
            raise ValueError(
                f'{GROUP_COL}: {int(blank.sum()) - already_null} rows have a '
                f'blank {RAW_FIELD} but a non-null display class. The coding and '
                'the raw field disagree; investigate the export before plotting.'
            )
        df.loc[blank, GROUP_COL] = np.nan
        print(f'  DROP_NON_RESPONDERS -> blank {RAW_FIELD}: {int(blank.sum())} '
              f'of {len(df)} rows dropped from the display-class panels '
              f'({n_blank_in_sample} of them had VCH data)')

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    group_col = GROUP_COL

    # One row per grouping variable, one column per DV.
    panel_specs = [{
        'dv': dv,
        'group_var': group_col,
        'ylabel': VARIABLE_REGISTRY[dv]['plot_label'],
        'xlabel': GROUP_LABELS[group_col],
        # Analyst decision 2026-08-28: a two-level column is tested with
        # Mann-Whitney rather than a two-group Kruskal-Wallis, so it
        # carries a U and a rank-biserial effect size.  Everything else is
        # Kruskal-Wallis.  `None` reproduces exactly that rule, but it is
        # spelled out per panel so the intent is not implicit.
        'force_test': 'mannwhitney' if df[group_col].nunique() == 2
                      else 'kruskal',
    } for dv in DV_COLS]

    print(f'\n=== {group_col} ===')
    _, summary = hardware_boxplot_grid(
        df,
        panel_specs,
        # One column, one DV per row: d' on top (a), threshold below (b).
        # DV_COLS order IS the row order.
        nrows=len(DV_COLS),
        ncols=1,
        fig_width=FIG_WIDTH_IN,
        panel_labels=[chr(ord('a') + i) for i in range(len(DV_COLS))],
        # Both rows share `group_col`, so only the bottom panel carries the
        # x-tick labels and the x-axis title.
        xaxis_bottom_row_only=True,
        savepath=OUT_DIR / 'supplementary_figure_s8',
    )

    # One grouping, one summary table. The per-panel statistics (omnibus plus
    # every pairwise post-hoc) are already written to
    # summary_results/supplementary_figure_s8.csv by hardware_boxplot_grid().
    omnibus = summary[summary.row_type == 'overall']
    print('\n--- omnibus tests ---')
    print(omnibus[['dv', 'group_var', 'test', 'n_groups',
                   'sample_size', 'statistic', 'p_value',
                   'effect_size_name', 'effect_size']].to_string(index=False))
    n_sig = int((omnibus.p_value < 0.05).sum())
    print(f'\n{n_sig} of {len(omnibus)} omnibus tests significant at p < .05 '
          f'(uncorrected across those {len(omnibus)}).')


if __name__ == '__main__':
    main()
