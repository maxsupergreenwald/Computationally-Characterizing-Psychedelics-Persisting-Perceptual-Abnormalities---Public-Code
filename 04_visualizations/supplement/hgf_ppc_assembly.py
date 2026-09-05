#!/usr/bin/env /usr/local/bin/python3.12
"""
hgf_ppc_assembly.py

Assembles two PPC (posterior predictive check) figures into a single
stacked multipanel PNG:

    [ppc_classic_2level_empiric_conditions.png]   ← top panel
    [ppc_classic_2level_empiric_blocks.png]        ← bottom panel

Edits applied before assembly
------------------------------
    Bottom panel — the legend box (top-right corner of the plot) is removed
                   by filling it with white (top panel already has legend).

No other labels or annotations are added (unlabeled multipanel figure).
The source images already have correct axis labels
("QUEST-Derived Stimulus Intensity (% Detection Probability)" / "Empiric Detection Probability").

Reads
-----
    the Julia HGF output tree (02_hgf_modeling/julia_outputs/):
        param_recovery/ppc_classic/figures/2level_empiric/
        ppc_classic_2level_empiric_conditions.png   (750 × 1050 px)
        ppc_classic_2level_empiric_blocks.png       (750 × 1500 px)

Writes
------
    {repo_root}/results/supplement/hgf_figures/ppc_2level_stacked.png

Layout
------
    Both panels are 750 px tall.  The narrower top panel (1050 px wide) is
    centered horizontally over the wider bottom panel (1500 px wide), so
    the x-axes align by visual center.  A small gap separates the panels.

Notes on crop coordinates (pixel space)
-----------------------------------------
    Bottom panel legend box: rows 32–135, cols 995–1467.
      Detected via pixel scanning of the top-right of the plot area.
      Padding of ~4 px is included.

Common things to change
-----------------------
    GAP_PX: pixel gap between the two panels (default 30 px).
    LEGEND_BOX: bounding box of the legend in the bottom panel, as
        (row_start, row_end, col_start, col_end).  Update if the source
        images are re-rendered at a different resolution or style.
"""

import sys
from pathlib import Path
import numpy as np
from PIL import Image

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent           # 04_visualizations/supplement/
REPO_ROOT  = SCRIPT_DIR.parent.parent                  # hppd_manuscript_public/
# Vendored HGF pipeline outputs (see "Script reference" in 02_hgf_modeling/README.md).
SRC_DIR    = (REPO_ROOT / '02_hgf_modeling' / 'julia_outputs'
              / 'param_recovery' / 'ppc_classic' / 'figures' / '2level_empiric')
OUT_DIR = REPO_ROOT / 'results' / 'supplement' / 'hgf_figures'
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Edit coordinates (pixel space, 0-indexed) ─────────────────────────────────
# Bottom panel: bounding box of legend to remove
LEGEND_BOX = (32, 135, 995, 1467)   # (row_start, row_end, col_start, col_end)

# ── Gap between panels ─────────────────────────────────────────────────────────
GAP_PX = 30   # pixels of white space between the two panels


# ── Load images ───────────────────────────────────────────────────────────────
top_path = SRC_DIR / 'ppc_classic_2level_empiric_conditions.png'
bot_path = SRC_DIR / 'ppc_classic_2level_empiric_blocks.png'

top_src = np.array(Image.open(top_path).convert('RGB'))
bot_src = np.array(Image.open(bot_path).convert('RGB'))

print(f'Top panel:    {top_src.shape[0]} × {top_src.shape[1]} px')
print(f'Bottom panel: {bot_src.shape[0]} × {bot_src.shape[1]} px')


# ── Edit bottom panel: remove legend ──────────────────────────────────────────
# Top panel kept as-is (correct axis labels from ppc_classic, legend visible).
top_edited = top_src

bot_edit = bot_src.copy()
lr0, lr1, lc0, lc1 = LEGEND_BOX
bot_edit[lr0:lr1, lc0:lc1, :] = 255  # white out legend rectangle
bot_edited = bot_edit
print(f'Bottom panel: removed legend box [{lr0}:{lr1}, {lc0}:{lc1}]')


# ── Assemble: top (1050 px) centered over bottom (1500 px) ────────────────────
top_h, top_w = top_edited.shape[:2]   # 750, 1050
bot_h, bot_w = bot_edited.shape[:2]   # 750, 1500
canvas_w     = bot_w                  # full width set by the wider panel
canvas_h     = top_h + GAP_PX + bot_h

canvas = Image.new('RGB', (canvas_w, canvas_h), (255, 255, 255))

top_left_px = (canvas_w - top_w) // 2   # center the narrower top panel
canvas.paste(Image.fromarray(top_edited), (top_left_px, 0))
canvas.paste(Image.fromarray(bot_edited), (0, top_h + GAP_PX))

print(f'\nCanvas: {canvas_w} × {canvas_h} px  '
      f'(top offset: {top_left_px} px from left)')

out_path = OUT_DIR / 'ppc_2level_stacked.png'
canvas.save(out_path)
tiff_path = OUT_DIR / 'ppc_2level_stacked.tiff'
canvas.save(tiff_path)
# SVG: the canvas is a raster composite, so the SVG embeds it losslessly rather
# than vectorising it (figure_assembly.raster_to_svg).
svg_path = OUT_DIR / 'ppc_2level_stacked.svg'
sys.path.insert(0, str(REPO_ROOT / 'modules'))
from figure_assembly import raster_to_svg
raster_to_svg(out_path, svg_path)
print(f'Saved → {out_path}')
print(f'Saved → {tiff_path}')
print(f'Saved → {svg_path}')
