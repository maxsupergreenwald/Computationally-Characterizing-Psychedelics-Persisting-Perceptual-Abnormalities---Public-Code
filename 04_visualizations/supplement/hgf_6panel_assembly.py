#!/usr/bin/env /usr/local/bin/python3.12
"""
hgf_6panel_assembly.py

Assembles six HGF validation panels into a single 2-row × 3-column
publication-ready figure.

Layout
------
Row 1 (top):
    a  param_recovery_scatterplots_all_models.png  — prior-based scatter (4 models)
    b  model_identifiability_confusion.png          — 4-way confusion matrix
    c  bms_ef_pxp_modified.png                     — RFX-BMS Ef + PXP bars

Row 2 (bottom):
    d  corr_gen_vs_rec_2level.png                  — posterior recovery scatter (ν,β,ω)
    e  pair_beta_nu_2level.png                     — pair plot (generative × recovered)
    f  ppc_2level_stacked.png                      — PPC conditions + blocks

Each row is set at a uniform height determined so that all three panels
fill the full content width at their native pixel aspect ratios (no
stretching or cropping).  The two row heights differ slightly because
the sums of the two rows' aspect ratios differ.

Reads
-----
    {repo_root}/results/supplement/hgf_figures/
        param_recovery_scatterplots_all_models.png
        bms_ef_pxp_modified.png
        ppc_2level_stacked.png
    julia_hgf_ch/param_recovery/prior_based_mcmc/figures/aic_bic/
        model_identifiability_confusion_bic.png
    julia_hgf_ch/param_recovery/figures/
        corr_gen_vs_rec_2level.png
        pair_beta_nu_2level.png

Writes
------
    {repo_root}/results/supplement/hgf_figures/supplementary_figure_s3.png

Common things to change
-----------------------
    FIGW          — total figure width in inches (default 18).
    COL_GAP_IN    — horizontal gap between adjacent panels (default 0.15 in).
    ROW_GAP_IN    — vertical gap between the two rows (default 0.25 in).
    LEFT/RIGHT_IN — outer horizontal margins (default 0.30 in each).
    TOP/BOT_IN    — outer vertical margins (default 0.20 in each).
    LABEL_SIZE    — font size for panel labels A–F (default 14).
"""

from pathlib import Path
import subprocess
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from PIL import Image

matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']


# ── Per-panel preprocessing ────────────────────────────────────────────────────
# Some source images carry a figure-level title that looks inconsistent in an
# assembled multi-panel figure (other panels have no title).  We crop those
# titles off before placing the panel.
#
# Strategy: detect content clusters along rows (contiguous bands of non-white
# pixels).  If the first cluster is a narrow title band, skip it and keep
# everything from the second cluster downward.

def _find_clusters_1d(mask: np.ndarray, min_gap: int = 8) -> list:
    """Contiguous runs of True in a 1-D bool array, merging gaps < min_gap."""
    padded = np.concatenate([[False], mask.astype(bool), [False]])
    diff   = np.diff(padded.astype(np.int8))
    starts = np.where(diff ==  1)[0]
    ends   = np.where(diff == -1)[0]
    merged = []
    for s, e in zip(starts, ends):
        if merged and (s - merged[-1][1]) < min_gap:
            merged[-1] = (merged[-1][0], e)
        else:
            merged.append((s, e))
    return merged


def crop_top_title(img_arr: np.ndarray, white_thresh: int = 245) -> np.ndarray:
    """
    Remove a title band from the top of an image.

    Detects row-wise content clusters.  If the topmost cluster spans fewer than
    15 % of the total image height (consistent with a title line vs. the main
    axes content), it is treated as the title and dropped.  Everything from the
    start of the second cluster downward is kept.

    Returns the original array unchanged if there is only one content cluster
    (no detectable title separation).
    """
    rgb = img_arr[:, :, :3]
    dark = np.any(rgb < white_thresh, axis=2)
    clusters = _find_clusters_1d(dark.any(axis=1), min_gap=8)
    if len(clusters) >= 2:
        first_h = clusters[0][1] - clusters[0][0]
        if first_h < 0.15 * img_arr.shape[0]:   # title is < 15 % of image height
            print(f'    crop_top_title: removed rows 0–{clusters[1][0]} '
                  f'(title cluster was {first_h} px tall)')
            return img_arr[clusters[1][0]:]
    print('    crop_top_title: no separable title found; keeping full image')
    return img_arr

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT  = SCRIPT_DIR.parent.parent
SUPP_FIG   = REPO_ROOT / 'results' / 'supplement' / 'hgf_figures'
# Vendored HGF pipeline outputs (see "Script reference" in 02_hgf_modeling/README.md).
JUL_PAR    = REPO_ROOT / '02_hgf_modeling' / 'julia_outputs' / 'param_recovery'
OUT_DIR    = SUPP_FIG
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Run prerequisite scripts ───────────────────────────────────────────────────
# hgf_param_recovery_assembly.py, hgf_bms_modified.py, and hgf_ppc_assembly.py
# must all run before panel images are loaded below.
_PREREQS = [
    SCRIPT_DIR / 'hgf_param_recovery_assembly.py',
    SCRIPT_DIR / 'hgf_bms_modified.py',
    SCRIPT_DIR / 'hgf_ppc_assembly.py',
]
for _prereq in _PREREQS:
    print(f'\n[6panel] Running prerequisite: {_prereq.name}')
    try:
        subprocess.run([sys.executable, str(_prereq)], check=True)
    except subprocess.CalledProcessError:
        # Allow failure if the output already exists from a prior run
        print(f'  WARNING: {_prereq.name} failed — using cached output if available')

PANEL_PATHS = [
    # Row 1
    SUPP_FIG / 'param_recovery_scatterplots_all_models.png',     # A
    JUL_PAR  / 'prior_based_mcmc/figures/aic_bic/model_identifiability_confusion_bic.png',  # B
    SUPP_FIG / 'bms_ef_pxp_modified.png',                        # C
    # Row 2
    JUL_PAR  / 'figures/corr_gen_vs_rec_2level_empiric.png',      # D
    JUL_PAR  / 'figures/pair_beta_nu_2level_empiric.png',         # E
    SUPP_FIG / 'ppc_2level_stacked.png',                          # F
]
PANEL_LABELS = ['a', 'b', 'c', 'd', 'e', 'f']

# ── Figure geometry ────────────────────────────────────────────────────────────
FIGW       = 18.0   # total figure width (inches)
LEFT_IN    = 0.30   # outer left margin
RIGHT_IN   = 0.30   # outer right margin
TOP_IN     = 0.20   # outer top margin
BOT_IN     = 0.20   # outer bottom margin
COL_GAP_IN = 0.15   # gap between panels in the same row
ROW_GAP_IN = 0.25   # gap between row 1 and row 2
LABEL_SIZE = 16     # panel label font size

# Panels whose figure-level title should be cropped before assembly.
# Key = index into PANEL_PATHS (0-based).  Value = crop function to apply.
PANEL_PREPROCESS = {
    1: crop_top_title,   # B — model_identifiability_confusion.png has a title
}

# ── Load all images ────────────────────────────────────────────────────────────
imgs = []
for idx, p in enumerate(PANEL_PATHS):
    arr = np.array(Image.open(p).convert('RGB'))
    if idx in PANEL_PREPROCESS:
        print(f'  Preprocessing panel index {idx} ({p.name}):')
        arr = PANEL_PREPROCESS[idx](arr)
    imgs.append(arr)
    print(f'  {p.name}: {arr.shape[1]}×{arr.shape[0]} px  '
          f'aspect={arr.shape[1]/arr.shape[0]:.4f}')

aspects = [img.shape[1] / img.shape[0] for img in imgs]  # W/H

# ── Compute row heights ────────────────────────────────────────────────────────
# Content width available for panels (excluding outer margins and column gaps)
content_w = FIGW - LEFT_IN - RIGHT_IN - 2 * COL_GAP_IN

# Row 1: panels A, B, C (indices 0-2)
r1_sum = sum(aspects[:3])   # sum of aspect ratios
H1 = content_w / r1_sum     # height that makes all three fill content_w exactly

# Row 2: panels D, E, F (indices 3-5)
r2_sum = sum(aspects[3:])
H2 = content_w / r2_sum

FIGH = BOT_IN + H2 + ROW_GAP_IN + H1 + TOP_IN
print(f'\nFIGW={FIGW:.1f}  FIGH={FIGH:.2f}  H1={H1:.3f}  H2={H2:.3f} in')

# ── Panel positions (inches from figure bottom-left) ──────────────────────────
def panel_positions(row_aspects, H, y_bot):
    """Return list of (x_left, y_bot, width, height) in inches for each panel."""
    positions = []
    x = LEFT_IN
    for asp in row_aspects:
        w = asp * H
        positions.append((x, y_bot, w, H))
        x += w + COL_GAP_IN
    return positions

row1_pos = panel_positions(aspects[:3], H1, y_bot=BOT_IN + H2 + ROW_GAP_IN)
row2_pos = panel_positions(aspects[3:], H2, y_bot=BOT_IN)
all_pos  = row1_pos + row2_pos

# ── Build figure ───────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(FIGW, FIGH), facecolor='white')

def _frac(x_in, y_in, w_in, h_in):
    return [x_in/FIGW, y_in/FIGH, w_in/FIGW, h_in/FIGH]

for i, (img, pos, label) in enumerate(zip(imgs, all_pos, PANEL_LABELS)):
    x, y, w, h = pos
    ax = fig.add_axes(_frac(x, y, w, h))
    ax.imshow(img, aspect='auto', interpolation='lanczos')
    ax.axis('off')
    # Panel label at top-left inside the panel.
    # Panel E (index 4) is nudged right to clear the pair-plot histogram.
    # All labels shifted 2% leftward relative to original positions.
    label_x = -0.055 if i == 4 else -0.01
    ax.text(label_x, 0.99, label,
            transform=ax.transAxes,
            fontsize=LABEL_SIZE, fontweight='bold',
            va='top', ha='left',
            color='black')

# ── Save ──────────────────────────────────────────────────────────────────────
out_path = OUT_DIR / 'supplementary_figure_s3.png'
fig.savefig(out_path, dpi=200, bbox_inches='tight', facecolor='white')
tiff_path = OUT_DIR / 'supplementary_figure_s3.tiff'
fig.savefig(tiff_path, dpi=200, bbox_inches='tight', facecolor='white')
svg_path = OUT_DIR / 'supplementary_figure_s3.svg'
fig.savefig(svg_path, dpi=200, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'\nSaved → {out_path}')
print(f'Saved → {tiff_path}')
