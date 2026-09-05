#!/usr/bin/env /usr/local/bin/python3.12
"""
beta_sigmoid_creator.py
-----------------------
Plots the unit-square sigmoid function used in the HGF response model for
two values of the action precision parameter β.

The function transforms the belief b_t ∈ (0, 1) into a probability P(u_t)
of making a subjective sensory observation u_t ∈ {0, 1} (signal absent/present):

    P(u_t) = b_t^β / (b_t^β + (1 − b_t)^β)

β controls the precision of this probability:
  β < 1  →  probabilities compressed toward 0.5 (more stochastic)
  β > 1  →  probabilities pushed toward 0 and 1 (more deterministic)

Two curves are plotted: β ∈ {0.25, 1.75}.

Staggered dotted reference lines are drawn at x = 0.3 and x = 0.7 for each
curve: vertical lines from y = 0 to the curve height, and horizontal lines
from x = 0 to the x reference position, in the respective curve color.

Reads:   nothing (pure math)
Outputs: results/supplement/hgf_figures/supplementary_figure_s1.png   (for manuscript drafting)
         results/supplement/hgf_figures/supplementary_figure_s1.tiff  (for journal submission)
         results/supplement/hgf_figures/supplementary_figure_s1.svg   (vector, for journal submission)

Output directory is gitignored — regenerate by running this script.

Run from any directory (paths resolve relative to this file):
    /usr/local/bin/python3.12 04_visualizations/supplement/beta_sigmoid_creator.py
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')

# ── Arial font (set before importing pyplot so rcParams take effect) ──────────
matplotlib.rcParams['font.family']     = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection
from pathlib import Path

# ── Output paths (resolved relative to this script's location) ───────────────
_THIS_DIR    = Path(__file__).resolve().parent          # .../04_visualizations/supplement/
_PROJECT_ROOT = _THIS_DIR.parent.parent                 # .../hppd_manuscript_public/
OUT_DIR      = _PROJECT_ROOT / 'results' / 'supplement' / 'hgf_figures'
OUT_PNG      = OUT_DIR / 'supplementary_figure_s1.png'
OUT_TIFF     = OUT_DIR / 'supplementary_figure_s1.tiff'
OUT_SVG      = OUT_DIR / 'supplementary_figure_s1.svg'

# ── Config ───────────────────────────────────────────────────────────────────

# β values to plot — β=1 (identity / linear) deliberately excluded.
BETAS  = [0.25, 1.75]

DPI      = 600

# Display toggles
SHOW_LEGEND = False   # set True to show per-β legend
SHOW_TICKS  = False   # set True to show x/y tick marks and labels

# Colors: electric blue (lowest β) → dark (highest β).
# Index order matches BETAS order.
COLORS  = ['#00BFFF', '#141414']
LWIDTHS = [2.0, 2.0]

# X positions at which to draw reference lines (unit-square domain).
# Both curves share the exact same x coordinate; they are visually distinguished
# by dot size: light blue uses thick dotted line (large dots) and dark gray uses
# thin dotted line (small dots). Both reference lines are drawn with alpha=0.5.
X_REFS       = [0.3, 0.7]
REF_LW_THIN  = 1.2   # linewidth for dark-gray reference lines  → small dots
REF_LW_THICK = 3.5   # linewidth for light-blue reference lines → large dots
REF_ALPHA    = 0.5   # transparency for all reference lines

# ── Math ─────────────────────────────────────────────────────────────────────

def unit_square_sigmoid(b, beta):
    """P(u_t) = b^β / (b^β + (1−b)^β). Safe for b ∈ (0, 1)."""
    b_b   = np.power(b, beta)
    inv_b = np.power(1.0 - b, beta)
    return b_b / (b_b + inv_b)

b_vals = np.linspace(0.001, 0.999, 1000)   # avoid exact 0 and 1

# ── Plot ─────────────────────────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(4, 4))

# Main sigmoid curves
for beta, color, lw in zip(BETAS, COLORS, LWIDTHS):
    p_vals = unit_square_sigmoid(b_vals, beta)
    ax.plot(b_vals, p_vals, color=color, linewidth=lw,
            label=f'β = {beta}', zorder=3)

# ── Reference lines at X_REFS ────────────────────────────────────────────────
# Both curves share the same x coordinate at each reference position.
# Visual distinction is by dot size:
#   • Light blue (COLORS[0], low β)  → thick dotted line (large dots, REF_LW_THICK)
#   • Dark gray  (COLORS[1], high β) → thin dotted line  (small dots, REF_LW_THIN)
# Vertical lines go from y = 0 to the curve's y-value at x_ref.
# Horizontal lines go from x = 0 to x_ref at the same y-value.
# Both are drawn at alpha = REF_ALPHA (0.5).
for x_ref in X_REFS:
    y_low  = unit_square_sigmoid(x_ref, BETAS[0])   # y on the light-blue curve
    y_high = unit_square_sigmoid(x_ref, BETAS[1])   # y on the dark-gray curve

    # Vertical dotted lines — each rises to its own curve's height at x_ref
    ax.plot([x_ref, x_ref], [0, y_low],
            color=COLORS[0], linewidth=REF_LW_THICK, linestyle=':', alpha=REF_ALPHA, zorder=2)
    ax.plot([x_ref, x_ref], [0, y_high],
            color=COLORS[1], linewidth=REF_LW_THIN,  linestyle=':', alpha=REF_ALPHA, zorder=2)

    # Horizontal dotted lines — run from the y-axis to x_ref at each curve's y-value
    ax.plot([0, x_ref], [y_low,  y_low],
            color=COLORS[0], linewidth=REF_LW_THICK, linestyle=':', alpha=REF_ALPHA, zorder=2)
    ax.plot([0, x_ref], [y_high, y_high],
            color=COLORS[1], linewidth=REF_LW_THIN,  linestyle=':', alpha=REF_ALPHA, zorder=2)

# ── Axes limits and ticks ─────────────────────────────────────────────────────

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)

if SHOW_TICKS:
    ax.set_xticks([0, 0.5, 1])
    ax.set_yticks([0, 0.5, 1])
    ax.tick_params(labelsize=11)
else:
    ax.set_xticks([])
    ax.set_yticks([])

# ── Axis labels (always shown, Arial via rcParams) ────────────────────────────
ax.set_xlabel("P(Target | Cue)",      fontsize=11)
ax.set_ylabel("P(Reporting 'Yes')", fontsize=11)

if SHOW_LEGEND:
    ax.legend(fontsize=11, framealpha=0.9, loc='upper left')

# ── Gradient arrow showing direction of increasing β (blue → dark gray) ──────
# Drawn in axes-fraction coordinates via a LineCollection so it floats above
# the plot area independently of the data scale.
ARROW_X0, ARROW_X1 = 0.57, 0.83   # left/right extent in axes fraction
ARROW_Y             = 0.91         # vertical position in axes fraction
N_SEG               = 200          # gradient resolution

_xs = np.linspace(ARROW_X0, ARROW_X1, N_SEG + 1)
_ys = np.full(N_SEG + 1, ARROW_Y)
_pts  = np.column_stack([_xs, _ys]).reshape(-1, 1, 2)
_segs = np.concatenate([_pts[:-1], _pts[1:]], axis=1)

_cmap = mcolors.LinearSegmentedColormap.from_list(
    'arrow_grad', [COLORS[0], COLORS[-1]]   # blue → dark
)
_lc = LineCollection(_segs, cmap=_cmap, transform=ax.transAxes,
                     linewidth=1.6, zorder=5, clip_on=False)
_lc.set_array(np.linspace(0, 1, N_SEG))
ax.add_collection(_lc)

# Dark arrowhead at the right end of the gradient arrow
ax.annotate(
    '',
    xy=(ARROW_X1, ARROW_Y), xycoords='axes fraction',
    xytext=(ARROW_X1 - 0.001, ARROW_Y), textcoords='axes fraction',
    arrowprops=dict(arrowstyle='->', color=COLORS[-1], lw=1.6),
    zorder=6,
)
ax.text(
    (ARROW_X0 + ARROW_X1) / 2, ARROW_Y + 0.045,
    'increasing β',
    transform=ax.transAxes,
    ha='center', va='bottom', fontsize=9.5, color='#444444',
)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(False)

fig.tight_layout()
OUT_DIR.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_PNG,  dpi=DPI, bbox_inches='tight')
fig.savefig(OUT_TIFF, dpi=DPI, bbox_inches='tight')
fig.savefig(OUT_SVG,  dpi=DPI, bbox_inches='tight')
plt.close(fig)
print(f'Saved: {OUT_PNG.resolve()}')
print(f'Saved: {OUT_TIFF.resolve()}')
