"""
prior_recovery_aic_bic.py  ─  AIC/BIC confusion matrices for 4-way prior-based recovery
=========================================================================================

Reads the compiled prior-based MCMC recovery results (prior_recovery_mcmc_summary.csv,
which includes aic_winner and bic_winner columns added by prior_recovery_aic_bic.jl),
then produces two confusion matrix figures — one for BIC and one for AIC — in the
same visual style as model_identifiability_confusion.png (bridge sampling LME version).

Figures are saved to:
  param_recovery/prior_based_mcmc/figures/aic_bic/
    model_identifiability_confusion_bic.png
    model_identifiability_confusion_aic.png

The summary CSV is recompiled from the per-simulation CSVs on each run to pick up
any newly added columns (e.g. aic_winner / bic_winner from the patch script).
"""

################################################################################
# ─── ROLE IN THE MANUSCRIPT ───────────────────────────────────────────────────
#
# PRIOR-BASED RECOVERY, plotting step.  Produces Supplementary Fig. S3b —
# the published model-identifiability confusion matrix.
#
#     prior_recovery_aic_bic.jl ──▶ THIS SCRIPT
#
# How to read the figure it makes
#     Rows are the model that generated the data; columns are the model that won
#     after inversion.  Cells are row-normalised percentages, so a perfectly
#     identifiable model set would be a bright leading diagonal.
#
#     What was actually observed is that the 2-level models won regardless of
#     which model generated the data and regardless of stimulus convention —
#     consistent with the third level carrying no recoverable information, and
#     with omega3 being unrecoverable (Fig. S3a).
#
#     The AIC version is emitted alongside the BIC one for completeness; BIC is
#     what the manuscript reports.
################################################################################

import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path

################################################################################
# ─── PATHS ────────────────────────────────────────────────────────────────────
################################################################################

# Resolved from this file's own location — see hgf_pipeline.py LOCAL PATHS.
LOCAL_JULIA_CH_DIR = os.path.dirname(os.path.abspath(__file__))

LOCAL_PRIOR_MCMC_DIR = os.path.join(LOCAL_JULIA_CH_DIR, "param_recovery", "prior_based_mcmc")
RESULTS_DIR  = os.path.join(LOCAL_PRIOR_MCMC_DIR, "results")
FIGURES_DIR  = os.path.join(LOCAL_PRIOR_MCMC_DIR, "figures", "aic_bic")

os.makedirs(FIGURES_DIR, exist_ok=True)

################################################################################
# ─── MODEL METADATA (identical to prior_recovery_vch_mcmc.py) ─────────────────
################################################################################

N_SIM       = 500
MODEL_TYPES = ["2level_empiric", "3level_empiric",
               "2level_nominal",  "3level_nominal"]

MODEL_LABELS = {
    "2level_empiric" : "2-level\n(empiric)",
    "3level_empiric" : "3-level\n(empiric)",
    "2level_nominal" : "2-level\n(nominal)",
    "3level_nominal" : "3-level\n(nominal)",
}

################################################################################
# ─── COMPILE SUMMARY CSV ──────────────────────────────────────────────────────
# Recompile from per-sim CSVs so we always pick up the latest columns.
################################################################################

dfs     = []
n_found = 0
n_miss  = 0

for true_model in MODEL_TYPES:
    for sim_index in range(1, N_SIM + 1):
        fpath = os.path.join(RESULTS_DIR, f"sim{sim_index}_{true_model}.csv")
        if not os.path.exists(fpath):
            n_miss += 1
            continue
        try:
            dfs.append(pd.read_csv(fpath))
            n_found += 1
        except Exception as exc:
            print(f"  Warning: could not read {fpath}: {exc}")
            n_miss += 1

if not dfs:
    raise FileNotFoundError(f"No result CSVs found in {RESULTS_DIR}")

summary = (pd.concat(dfs, ignore_index=True)
             .sort_values(["true_model", "sim_index"])
             .reset_index(drop=True))

expected = N_SIM * len(MODEL_TYPES)
print(f"Compiled {n_found}/{expected} result files ({n_miss} missing)")

# Check AIC/BIC columns are present
for col in ("aic_winner", "bic_winner"):
    missing_col = summary[col].isna().sum()
    if missing_col > 0:
        print(f"  WARNING: {missing_col} rows have missing {col}")
    else:
        print(f"  {col}: all {len(summary)} rows present")

# Save updated summary
summary_path = os.path.join(RESULTS_DIR, "prior_recovery_mcmc_summary.csv")
summary.to_csv(summary_path, index=False)
print(f"Saved: {summary_path}\n")

################################################################################
# ─── CONFUSION MATRIX FUNCTION ────────────────────────────────────────────────
# Reproduces the exact figure style from prior_recovery_vch_mcmc.py:
#   Figure 2 — 4×4 confusion matrix heatmap
################################################################################

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

matplotlib.rcParams.update({
    "font.family"    : "Arial",
    "axes.labelsize" : 20,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 13,
})


def plot_confusion_matrix(summary_df, winner_col, criterion_label, out_path):
    """
    Reproduce the bridge-sampling confusion matrix figure exactly, substituting
    `winner_col` (e.g. 'aic_winner' or 'bic_winner') for 'winner', and updating
    the axis label and title to reflect the criterion.

    Parameters
    ----------
    summary_df    : compiled summary DataFrame (one row per simulation)
    winner_col    : column name containing the selected-model string
    criterion_label : short label for axis/title (e.g. 'AIC', 'BIC')
    out_path      : full path for the saved PNG
    """
    n_models = len(MODEL_TYPES)

    # ── Build count matrix ────────────────────────────────────────────────────
    conf = np.zeros((n_models, n_models))
    for i, true_m in enumerate(MODEL_TYPES):
        sub = summary_df[summary_df["true_model"] == true_m]
        for j, win_m in enumerate(MODEL_TYPES):
            conf[i, j] = (sub[winner_col] == win_m).sum()

    row_totals = conf.sum(axis=1, keepdims=True)
    conf_pct   = np.where(row_totals > 0, 100 * conf / row_totals, 0)

    # ── Print text summary ────────────────────────────────────────────────────
    col_w  = 18
    header = f"{'':25s}" + "".join(f"{'winner='+m:>{col_w}s}" for m in MODEL_TYPES)
    print(f"Model identifiability confusion matrix ({criterion_label}):")
    print("-" * (25 + col_w * len(MODEL_TYPES)))
    print(header)
    for i, true_m in enumerate(MODEL_TYPES):
        n   = int(row_totals[i, 0])
        row = f"  true={true_m:<20s}"
        for j in range(n_models):
            cnt = int(conf[i, j])
            pct = conf_pct[i, j]
            row += f"  {cnt:3d}/{n} ({pct:5.1f}%){'':<1s}"
        print(row)

    valid   = summary_df[winner_col].isin(MODEL_TYPES)
    correct = (summary_df.loc[valid, "true_model"] ==
               summary_df.loc[valid, winner_col]).sum()
    total   = valid.sum()
    print(f"\n  Overall correct: {correct}/{total} ({100*correct/total:.1f}%)\n")

    # ── Figure (identical geometry to bridge sampling version) ────────────────
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(conf_pct, cmap="Blues", vmin=0, vmax=100)

    x_labels = [MODEL_LABELS[m] for m in MODEL_TYPES]
    y_labels  = [MODEL_LABELS[m] for m in MODEL_TYPES]

    ax.set_xticks(range(n_models))
    ax.set_yticks(range(n_models))
    ax.set_xticklabels(x_labels, fontsize=12)
    ax.set_yticklabels(y_labels, fontsize=12)
    ax.set_xlabel(f"Model selected by {criterion_label}", fontsize=16)
    ax.set_ylabel("Data-generating model", fontsize=16)

    for i in range(n_models):
        for j in range(n_models):
            txt   = f"{conf_pct[i,j]:.1f}%\n(n={int(conf[i,j])})"
            color = "white" if conf_pct[i, j] > 60 else "black"
            ax.text(j, i, txt, ha="center", va="center",
                    fontsize=11, color=color,
                    fontweight="bold" if i == j else "normal")

    fig.colorbar(im, ax=ax, label="% simulations", shrink=0.85)
    fig.suptitle(f"Model identifiability — {criterion_label} (MCMC)", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {os.path.basename(out_path)}")


################################################################################
# ─── GENERATE FIGURES ─────────────────────────────────────────────────────────
# BIC first (preferred criterion for model selection; penalises complexity more),
# then AIC.
################################################################################

plot_confusion_matrix(
    summary,
    winner_col      = "bic_winner",
    criterion_label = "BIC",
    out_path        = os.path.join(FIGURES_DIR, "model_identifiability_confusion_bic.png"),
)

plot_confusion_matrix(
    summary,
    winner_col      = "aic_winner",
    criterion_label = "AIC",
    out_path        = os.path.join(FIGURES_DIR, "model_identifiability_confusion_aic.png"),
)

print(f"\nAll figures saved to: {FIGURES_DIR}")
