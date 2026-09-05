#!/usr/bin/env /usr/local/bin/python3.12
"""
hgf_bms_modified.py

Re-renders the RFX-BMS bar chart (Ef + PXP) from the pre-computed
bms_summary.csv with two cosmetic changes:
  1. Title removed (the original suptitle "RFX-BMS | n = … | BOR = …").
  2. X-axis tick labels capitalised: "level" → "Level",
     "empirical" → "Empirical", "nominal" → "Nominal".

Reads
-----
    the Julia HGF output tree (02_hgf_modeling/julia_outputs/):
        model_comparison/bms/bms_summary.csv

Writes
------
    {repo_root}/results/supplement/hgf_figures/bms_ef_pxp_modified.png

Derivation
----------
Plotting code surgically extracted from
julia_hgf_ch/bms_vch.py (DO_PLOT section, lines 239–295).
Only MODEL_LABELS and the suptitle call are changed.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR  = Path(__file__).resolve().parent
REPO_ROOT   = SCRIPT_DIR.parent.parent
# Written by 02_hgf_modeling/bms_vch.py (DO_BMS stage), which consumes the
# per-subject *_lme.csv files produced by bms_vch.jl on the cluster.
BMS_CSV     = (REPO_ROOT / '02_hgf_modeling' / 'model_comparison' / 'bms'
               / 'bms_summary.csv')
OUT_DIR = REPO_ROOT / 'results' / 'supplement' / 'hgf_figures'
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Config (kept in sync with bms_vch.py) ─────────────────────────────────────
MODEL_KEYS = ["2level_empiric", "3level_empiric", "2level_nominal", "3level_nominal"]

# Capitalised versions of the original labels (only capitalisation changed)
MODEL_LABELS = {
    "2level_empiric"  : "2-Level\n(Empirical)",
    "3level_empiric"  : "3-Level\n(Empirical)",
    "2level_nominal"  : "2-Level\n(Nominal)",
    "3level_nominal"  : "3-Level\n(Nominal)",
}

_PALETTE = ["#4472C4", "#ED7D31", "#70AD47", "#FFC000",
            "#5B9BD5", "#C55A11", "#375623"]

# ── Plot (DO_PLOT section from bms_vch.py, title removed) ─────────────────────
matplotlib.rcParams.update({
    "font.family":     "Arial",
    "axes.labelsize":  20,
    "xtick.labelsize": 13,
    "ytick.labelsize": 14,
})

summary  = pd.read_csv(BMS_CSV)
n_models = len(summary)
bor      = summary["BOR"].iloc[0]
n_sub    = summary["n_subjects"].iloc[0]
x        = np.arange(n_models)
colors   = [_PALETTE[i % len(_PALETTE)] for i in range(n_models)]
labels   = [MODEL_LABELS.get(k, k) for k in summary["model"]]
bar_w    = max(0.3, min(0.6, 2.4 / n_models))
fig_w    = max(9, 2.5 * n_models)

fig, axes = plt.subplots(1, 2, figsize=(fig_w, 4))

# ── Ef ────────────────────────────────────────────────────────────────────────
axes[0].bar(x, summary["Ef"], color=colors, width=bar_w, edgecolor="none")
axes[0].axhline(1 / n_models, color="gray", lw=1.2, linestyle="--", alpha=0.7)
axes[0].set_xticks(x)
axes[0].set_xticklabels(labels)
axes[0].set_ylabel("Ef")
axes[0].set_ylim(0, 1)
for xi, ef in zip(x, summary["Ef"]):
    axes[0].text(xi, ef + 0.02, f"{ef:.3f}", ha="center", va="bottom",
                 fontsize=11, fontweight="bold")
sns.despine(ax=axes[0])

# ── PXP ───────────────────────────────────────────────────────────────────────
axes[1].bar(x, summary["PXP"], color=colors, width=bar_w, edgecolor="none")
axes[1].axhline(1 / n_models, color="gray", lw=1.2, linestyle="--", alpha=0.7)
axes[1].set_xticks(x)
axes[1].set_xticklabels(labels)
axes[1].set_ylabel("PXP")
axes[1].set_ylim(0, 1)
for xi, pxp_v in zip(x, summary["PXP"]):
    axes[1].text(xi, pxp_v + 0.02, f"{pxp_v:.3f}", ha="center",
                 va="bottom", fontsize=11, fontweight="bold")
sns.despine(ax=axes[1])

# Title intentionally omitted (was: fig.suptitle(...))
fig.tight_layout()

out_path = OUT_DIR / "bms_ef_pxp_modified.png"
fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="white")
tiff_path = OUT_DIR / "bms_ef_pxp_modified.tiff"
fig.savefig(tiff_path, dpi=200, bbox_inches="tight", facecolor="white")
svg_path = OUT_DIR / "bms_ef_pxp_modified.svg"
fig.savefig(svg_path, dpi=200, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"Saved → {out_path}")
print(f"Saved → {tiff_path}")
