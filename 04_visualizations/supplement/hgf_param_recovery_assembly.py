#!/usr/bin/env /usr/local/bin/python3.12
"""
hgf_param_recovery_assembly.py

Assembles four prior-based HGF parameter recovery scatter plots into a
single stacked figure arranged in two pairs:

    ─── 2-Level ───  (centered horizontally — 3 scatter subplots, narrower)
    [2-level empiric]    ← "Empiric\nLikelihood" on left
    [2-level nominal]    ← "Nominal\nLikelihood" on left

    ─── 3-Level ───  (full content width — 4 scatter subplots, wider)
    [3-level empiric]    ← "Empiric\nLikelihood" on left
    [3-level nominal]    ← "Nominal\nLikelihood" on left

Reads
-----
    the Julia HGF output tree (02_hgf_modeling/julia_outputs/):
        param_recovery/prior_based_mcmc/figures/
        prior_recovery_scatter_{x}level_{type}.png
    for x ∈ {2, 3}, type ∈ {'empiric', 'nominal'}

Writes
------
    {repo_root}/results/supplement/hgf_figures/
        param_recovery_scatterplots_all_models.png

Layout details
--------------
    All four panels are rendered at the same height, determined by the
    3-level image width (which defines the full content width).
    The narrower 2-level panels are centered horizontally within the
    3-level content width.

    Pair titles "2-Level" / "3-Level" appear above each pair.
    Row labels "Empiric\\nLikelihood" / "Nominal\\nLikelihood" appear to the
    left of each panel, aligned to the 3-level left edge for visual
    consistency (not relative to each panel's own left edge).

    Shared axis labels (bold, fontsize=20):
      - "Generative parameter" at the figure bottom (centered)
      - "MAP recovered" at the figure left (rotated 90°, spanning all panels)

Notes
-----
    Per-image titles, X labels, and Y labels are stripped by a
    content-cluster crop algorithm.  See crop_image() / _find_clusters_1d().
    Positioning uses fig.add_axes() with inch-based geometry so that the
    2-level centering offset and equal-height constraint are exact.
"""

from pathlib import Path
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from PIL import Image

# Use Arial throughout; fall back to Helvetica then DejaVu Sans if not found.
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

# ── Paths ────────────────────────────────────────────────────────────────────
SCRIPT_DIR  = Path(__file__).resolve().parent           # 04_visualizations/supplement/
REPO_ROOT   = SCRIPT_DIR.parent.parent                  # hppd_manuscript_public/
# Vendored HGF pipeline outputs (see "Script reference" in 02_hgf_modeling/README.md).
FIGURES_DIR = (REPO_ROOT / '02_hgf_modeling' / 'julia_outputs'
               / 'param_recovery' / 'prior_based_mcmc' / 'figures')
OUT_DIR = REPO_ROOT / 'results' / 'supplement' / 'hgf_figures'
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Panel order (top to bottom in the final figure) ──────────────────────────
# Pair 1 (top):    2-level empiric, 2-level nominal   (centered, narrower)
# Pair 2 (bottom): 3-level empiric, 3-level nominal   (full content width)
PANELS = [
    ('2', 'empiric'),   # index 0
    ('2', 'nominal'),   # index 1
    ('3', 'empiric'),   # index 2
    ('3', 'nominal'),   # index 3
]
ROW_LABELS = [
    'Empiric\nLikelihood',
    'Nominal\nLikelihood',
    'Empiric\nLikelihood',
    'Nominal\nLikelihood',
]
PAIR_TITLES    = ['2-Level', '3-Level']   # placed above panel 0 and panel 2
PAIR_TOP_IDX   = [0, 2]
SHARED_X_LABEL = 'Generative parameter'
SHARED_Y_LABEL = 'MAP recovered'


# ── Crop helpers ─────────────────────────────────────────────────────────────

def _find_clusters_1d(mask: np.ndarray, min_gap: int = 8) -> list[tuple[int, int]]:
    """
    Find contiguous runs of True in a 1-D bool array, merging runs whose
    intervening gap is shorter than min_gap.

    Returns a list of (start, end) tuples where end is exclusive.
    """
    padded = np.concatenate([[False], mask.astype(bool), [False]])
    diff   = np.diff(padded.astype(np.int8))
    starts = np.where(diff ==  1)[0]
    ends   = np.where(diff == -1)[0]
    merged: list[tuple[int, int]] = []
    for s, e in zip(starts, ends):
        if merged and (s - merged[-1][1]) < min_gap:
            merged[-1] = (merged[-1][0], e)
        else:
            merged.append((s, e))
    return merged


def crop_image(img_arr: np.ndarray, white_thresh: int = 245) -> np.ndarray:
    """
    Remove the outer title, x-axis label, y-axis label, and right whitespace
    from a matplotlib-rendered PNG.

    Strategy: detect "content clusters" — contiguous bands of non-white rows
    or columns.  matplotlib figures have the structure:

        rows:  [title text] … [axes + scatter + tick labels] … [x-axis label]
        cols:  [y-axis label (rotated)] … [tick labels + axes content] … [whitespace]

    For rows: skip the first cluster (title) and the last cluster (x-label).
    For cols: skip the first cluster (y-label) and keep through the last cluster.

    Parameters
    ----------
    img_arr : np.ndarray (H, W, 3)  uint8
    white_thresh : int
        Pixels with ANY channel below this are considered non-white (content).

    Returns
    -------
    np.ndarray  Cropped image.
    """
    rgb = img_arr[:, :, :3]
    h, w = rgb.shape[:2]
    dark         = np.any(rgb < white_thresh, axis=2)
    row_has_dark = dark.any(axis=1)
    col_has_dark = dark.any(axis=0)
    row_clusters = _find_clusters_1d(row_has_dark, min_gap=8)
    col_clusters = _find_clusters_1d(col_has_dark, min_gap=8)
    print(f'  row clusters ({len(row_clusters)}): {row_clusters}')
    print(f'  col clusters ({len(col_clusters)}): {col_clusters}')

    if len(row_clusters) >= 3:
        top    = row_clusters[1][0]
        bottom = row_clusters[-2][1]
    elif len(row_clusters) == 2:
        top    = row_clusters[1][0]
        bottom = row_clusters[1][1]
    else:
        top    = row_clusters[0][0] if row_clusters else 0
        bottom = row_clusters[0][1] if row_clusters else h

    if len(col_clusters) >= 2:
        left  = col_clusters[1][0]
        right = col_clusters[-1][1]
    else:
        left  = col_clusters[0][0] if col_clusters else 0
        right = col_clusters[0][1] if col_clusters else w

    pad = 4
    top    = max(0, top    - pad)
    bottom = min(h, bottom + pad)
    left   = max(0, left   - pad)
    right  = min(w, right  + pad)
    print(f'  → crop: top={top} bottom={bottom} left={left} right={right}  '
          f'(original {h}×{w} → cropped {bottom-top}×{right-left})')
    return img_arr[top:bottom, left:right]


# ── Load and crop all panels ─────────────────────────────────────────────────
imgs: list[np.ndarray] = []
for x, kind in PANELS:
    fname = f'prior_recovery_scatter_{x}level_{kind}.png'
    fpath = FIGURES_DIR / fname
    print(f'Loading: {fname}')
    img = np.array(Image.open(fpath).convert('RGB'))
    imgs.append(crop_image(img))

# ── Image dimensions ─────────────────────────────────────────────────────────
w_px  = [img.shape[1] for img in imgs]    # pixel widths: [1805, 1805, 2440, 2440]
h_px  = [img.shape[0] for img in imgs]    # pixel heights: all ≈ 483
max_w = max(w_px)                          # 3-level width defines content width
img_h = max(h_px)                          # canonical image height (px), same for all

# ── Figure geometry (all measurements in inches) ──────────────────────────────
FIGW = 14.0        # total figure width

# Outer margins
LEFT_IN   = 2.0    # left margin: shared y-label + row labels
RIGHT_IN  = 0.15   # right margin
TOP_IN    = 0.15   # top margin
BOTTOM_IN = 0.75   # bottom margin for shared x-label

content_w_in = FIGW - LEFT_IN - RIGHT_IN   # width available for image content

# Scale: widest images (3-level, max_w px) fill the full content width.
# This guarantees all four panels render at the same physical height.
scale    = content_w_in / max_w            # inches per pixel
img_h_in = img_h * scale                   # image height in inches (equal for all panels)
w3_in    = max_w   * scale                 # 3-level panel width  (= content_w_in)
w2_in    = w_px[0] * scale                 # 2-level panel width  (narrower)

# Horizontal positions: 3-level fills LEFT_IN→RIGHT edge; 2-level is centered within that
left3_in = LEFT_IN
left2_in = LEFT_IN + (w3_in - w2_in) / 2   # centers 2-level over 3-level

# Vertical spacing
PAIR_GAP_IN  = 0.18   # gap between empiric and nominal within a pair
TITLE_H_IN   = 0.28   # height reserved for each pair title text
TITLE_PAD_IN = 0.10   # gap between an image top and the pair title above it
BIG_GAP_IN   = 0.30   # extra gap between "3-Level" title and 2-level nominal bottom

# Bottom-up layout (y = distance from figure bottom in inches):
#
#   BOTTOM_IN   ←  shared x-label lives here
#   y_nom3      →  3-level nominal bottom
#   + img_h_in  →  3-level nominal top
#   + PAIR_GAP  →  y_emp3
#   + img_h_in  →  3-level empiric top
#   + PAD       →  y_t3 ("3-Level" title bottom)
#   + TITLE_H   →  "3-Level" title top
#   + BIG_GAP   →  y_nom2 (2-level nominal bottom)
#   + img_h_in  →  2-level nominal top
#   + PAIR_GAP  →  y_emp2
#   + img_h_in  →  2-level empiric top
#   + PAD       →  y_t2 ("2-Level" title bottom)
#   + TITLE_H   →  "2-Level" title top
#   + TOP_IN    →  figure top
#
y_nom3 = BOTTOM_IN
y_emp3 = y_nom3 + img_h_in + PAIR_GAP_IN
y_t3   = y_emp3 + img_h_in + TITLE_PAD_IN           # "3-Level" title bottom edge
y_nom2 = y_t3   + TITLE_H_IN + BIG_GAP_IN
y_emp2 = y_nom2 + img_h_in + PAIR_GAP_IN
y_t2   = y_emp2 + img_h_in + TITLE_PAD_IN           # "2-Level" title bottom edge
FIGH   = y_t2   + TITLE_H_IN + TOP_IN               # total figure height

print(f'\nFigure: {FIGW:.1f} × {FIGH:.2f} in')
print(f'Image height: {img_h_in:.3f} in ({img_h} px)')
print(f'2-level width: {w2_in:.3f} in  |  3-level width: {w3_in:.3f} in')
print(f'2-level centering offset: {(left2_in - LEFT_IN):.3f} in from 3-level left edge')

# ── Build figure ──────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(FIGW, FIGH), facecolor='white')


def _frac(x_in: float, y_in: float, w_in: float, h_in: float) -> list[float]:
    """Convert inch coordinates to figure [left, bottom, width, height] fractions."""
    return [x_in / FIGW, y_in / FIGH, w_in / FIGW, h_in / FIGH]


# (bottom_in, left_in, width_in) for each panel in PANELS order (top → bottom)
panel_geom = [
    (y_emp2, left2_in, w2_in),   # 0: 2-level empiric
    (y_nom2, left2_in, w2_in),   # 1: 2-level nominal
    (y_emp3, left3_in, w3_in),   # 2: 3-level empiric
    (y_nom3, left3_in, w3_in),   # 3: 3-level nominal
]

axes: list[plt.Axes] = []
for y_b, x_l, w in panel_geom:
    ax = fig.add_axes(_frac(x_l, y_b, w, img_h_in))
    axes.append(ax)

# ── Display images ────────────────────────────────────────────────────────────
for ax, img in zip(axes, imgs):
    ax.imshow(img, aspect='auto', interpolation='lanczos')
    ax.axis('off')

# ── Row labels: "Empiric\nLikelihood" / "Nominal\nLikelihood" ─────────────────
# All four labels are placed at the same x position — aligned to the 3-level
# left edge offset — so they form a consistent vertical column regardless of
# whether the panel to their right is 2-level (centered) or 3-level (full width).
row_label_x_in = left3_in - 0.04 * w3_in    # same fractional offset as the widest panels
panel_y_centers = [y_b + img_h_in / 2 for y_b, _, _ in panel_geom]
for y_c, label in zip(panel_y_centers, ROW_LABELS):
    fig.text(
        row_label_x_in / FIGW, y_c / FIGH, label,
        ha='right', va='center',
        fontsize=11, multialignment='center', linespacing=1.4,
    )

# ── Pair titles: "2-Level" and "3-Level" ──────────────────────────────────────
# Both are centered at the same horizontal midpoint (LEFT_IN + w3_in/2), because
# the 2-level panels are centered within the 3-level width.
title_cx_in = LEFT_IN + w3_in / 2
title_y_centers = [y_t2 + TITLE_H_IN / 2, y_t3 + TITLE_H_IN / 2]   # "2-Level", "3-Level"
for title, y_c in zip(PAIR_TITLES, title_y_centers):
    fig.text(
        title_cx_in / FIGW, y_c / FIGH, title,
        ha='center', va='center',
        fontsize=12, fontweight='bold',
    )

# ── Shared X label ────────────────────────────────────────────────────────────
fig.text(
    (LEFT_IN + w3_in / 2) / FIGW,
    (BOTTOM_IN * 0.45) / FIGH,
    SHARED_X_LABEL,
    ha='center', va='center',
    fontsize=20, fontweight='bold',
)

# ── Shared Y label ────────────────────────────────────────────────────────────
# Centered vertically across all four panels
all_panels_cy = (y_nom3 + y_emp2 + img_h_in) / 2
fig.text(
    0.40 / FIGW, all_panels_cy / FIGH,
    SHARED_Y_LABEL,
    ha='center', va='center',
    fontsize=20, fontweight='bold',
    rotation=90,
)

# ── Save ─────────────────────────────────────────────────────────────────────
out_path = OUT_DIR / 'param_recovery_scatterplots_all_models.png'
fig.savefig(out_path, dpi=200, bbox_inches='tight', facecolor='white')
tiff_path = OUT_DIR / 'param_recovery_scatterplots_all_models.tiff'
fig.savefig(tiff_path, dpi=200, bbox_inches='tight', facecolor='white')
svg_path = OUT_DIR / 'param_recovery_scatterplots_all_models.svg'
fig.savefig(svg_path, dpi=200, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'\nSaved → {out_path}')
print(f'Saved → {tiff_path}')
