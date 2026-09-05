#!/usr/bin/env python3
"""
hgf_ppc_oos_assembly.py
========================
Three-panel supplementary figure (a, b, c) for the HPPD manuscript:

  a. Out-of-set empirical condition-by-condition detection rates (SP users).
     Source image from julia_hgf_ch/out_of_set_data/figures/.

  b. PPC conditions: stacked pair — empiric likelihood (top) over nominal
     likelihood (bottom).  Pre-rendered images from ppc_classic, with the
     x-axis label removed from the top subplot and a single shared y-axis
     label centered vertically.

  c. Out-of-set non-hallucinator group (n=29): per-block empirical detection
     rates at 75% contrast intensity, with 94% CI. Block is derived from trial
     number (block = ceil(trial / 30), 12 blocks of 30 trials). No legend.
     Note: condition-75 trials are sparse in later blocks due to QUEST adaptive
     staircasing (~2 trials/block in blocks 7-12), so the group mean is noisy.

Layout: 3 columns (a | b | c), assembled via figure_assembly.py with
standard 10pt bold Arial panel labels.

Reads:
    julia_hgf_ch/out_of_set_data/figures/empirical_condition_by_condition_spusers.png
    julia_hgf_ch/param_recovery/ppc_classic/figures/2level_empiric/ppc_classic_2level_empiric_conditions.png
    julia_hgf_ch/param_recovery/ppc_classic/figures/2level_nominal/ppc_classic_2level_nominal_conditions.png
    julia_hgf_ch/out_of_set_data/behavioral_data_OUT_OF_SET_with_metadata.csv

Writes:
    results/supplement/hgf_figures/supplementary_figure_s2.png
    results/supplement/hgf_figures/supplementary_figure_s2.tiff

Usage:
    cd <repo root>
    /usr/local/bin/python3.12 04_visualizations/supplement/hgf_ppc_oos_assembly.py

Author: Max Greenwald
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats as sp_stats
import seaborn as sns

# ─── PATHS ───────────────────────────────────────────────────────────────────
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
# Vendored HGF pipeline outputs (see "Script reference" in 02_hgf_modeling/README.md).
JULIA_BASE   = os.path.join(PROJECT_ROOT, '02_hgf_modeling', 'julia_outputs')

# Import figure_assembly from modules/
sys.path.insert(0, os.path.join(PROJECT_ROOT, 'modules'))
from figure_assembly import assemble_manuscript_figure

PANEL_A_IMG  = os.path.join(JULIA_BASE, 'out_of_set_data', 'figures',
                            'empirical_condition_by_condition_spusers.png')
# Panel B: generated from compiled PPC CSVs (empiric + nominal)
PANEL_B_EMP_CSV = os.path.join(JULIA_BASE, 'param_recovery', 'ppc_classic', 'results',
                               '2level_empiric', 'ppc_classic_2level_empiric_all.csv')
PANEL_B_NOM_CSV = os.path.join(JULIA_BASE, 'param_recovery', 'ppc_classic', 'results',
                               '2level_nominal', 'ppc_classic_2level_nominal_all.csv')
OOS_CSV      = os.path.join(JULIA_BASE, 'out_of_set_data',
                            'behavioral_data_OUT_OF_SET_with_metadata.csv')

OUT_DIR = os.path.join(PROJECT_ROOT, 'results', 'supplement', 'hgf_figures')
os.makedirs(OUT_DIR, exist_ok=True)
OUT_PATH = os.path.join(OUT_DIR, 'supplementary_figure_s2.png')

# Temporary panel PNGs (cleaned up after assembly)
_TMP_B = os.path.join(OUT_DIR, '_tmp_panel_b.png')
_TMP_C = os.path.join(OUT_DIR, '_tmp_panel_c.png')

# ─── CONFIG ──────────────────────────────────────────────────────────────────
PANEL_DPI      = 150       # DPI for generated panel C
CI_LEVEL       = 0.94      # 94% confidence interval
Z_94           = sp_stats.norm.ppf(1 - (1 - CI_LEVEL) / 2)  # ~1.88
N_BLOCKS       = 12
TRIALS_PER_BLK = 30   # 360 trials / 12 blocks = 30 trials per block

BLOCK_LABELS = [str(b) for b in range(1, N_BLOCKS + 1)]

# Panel B config — generated from scratch via matplotlib
PANEL_B_YLABEL    = 'Empiric Detection Probability'
PANEL_B_XLABEL    = 'QUEST-Derived Stimulus Intensity\n(% Detection Probability)'
N_BOOTSTRAP       = 10000   # bootstrap resamples for CI (matches ppc_classic_vch.py)
COND_COLS         = ["det_rate_0.0", "det_rate_0.25", "det_rate_0.5", "det_rate_0.75"]
COND_LABELS       = ["0%\nDetection\nProbability",
                     "25%\nDetection\nProbability",
                     "50%\nDetection\nProbability",
                     "75%\nDetection\nProbability"]

# ─── FONT SETUP ──────────────────────────────────────────────────────────────
matplotlib.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.labelsize":  18,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
})


# ─── HELPERS ─────────────────────────────────────────────────────────────────

def _bootstrap_group_ci(data_2d, rng, n_bootstrap=N_BOOTSTRAP):
    """
    Bootstrap CI for the group mean across participants.
    Mirrors ppc_classic_vch.py bootstrap_group_ci() exactly.
    """
    alpha_lo = (100 - CI_LEVEL * 100) / 2 / 100
    alpha_hi = 1.0 - alpha_lo
    n_pts = data_2d.shape[0]
    group_mean = np.nanmean(data_2d, axis=0)
    boot_means = np.full((n_bootstrap, data_2d.shape[1]), np.nan)
    for i in range(n_bootstrap):
        idx = rng.integers(0, n_pts, size=n_pts)
        boot_means[i] = np.nanmean(data_2d[idx], axis=0)
    lo = np.nanquantile(boot_means, alpha_lo, axis=0)
    hi = np.nanquantile(boot_means, alpha_hi, axis=0)
    return group_mean, lo, hi


def _plot_ppc_conditions(ax, df_all, rng, show_legend=True, show_xaxis=True):
    """
    Plot PPC conditions on a given axes.  Mirrors ppc_classic_vch.py _plot_ppc()
    + _save_ppc_fig() logic exactly: spaghetti lines, group mean ± CI for both
    empirical and simulated series.
    """
    emp = df_all[df_all["source"] == "empirical"].copy()
    sim = df_all[df_all["source"] == "sim"].copy()
    xs = list(range(len(COND_COLS)))
    ci_pct = int(CI_LEVEL * 100)

    # Spaghetti: individual simulated participant traces
    jitter_rng = np.random.default_rng(0)
    for _, row in sim[COND_COLS].iterrows():
        jitter = jitter_rng.uniform(-0.15, 0.15)
        ax.plot([x + jitter for x in xs], row.values.astype(float),
                color="steelblue", alpha=0.03, linewidth=0.6)

    # Empirical group mean + CI
    emp_arr = emp[COND_COLS].values.astype(float)
    emp_mean, emp_lo, emp_hi = _bootstrap_group_ci(emp_arr, rng)
    ax.fill_between(xs, emp_lo, emp_hi, alpha=0.30, color="black")
    ax.plot(xs, emp_mean, color="black", linewidth=2,
            marker="o", markersize=6,
            label=f"Empirical (mean ± {ci_pct}% CI)")

    # Simulated group mean + CI
    sim_arr = sim[COND_COLS].values.astype(float)
    sim_mean, sim_lo, sim_hi = _bootstrap_group_ci(sim_arr, rng)
    ax.fill_between(xs, sim_lo, sim_hi, alpha=0.30, color="steelblue")
    ax.plot(xs, sim_mean, color="steelblue", linewidth=2, linestyle="--",
            marker="o", markersize=6,
            label=f"Simulated (mean ± {ci_pct}% CI)")

    ax.set_ylim(-0.05, 1.05)
    sns.despine(ax=ax)

    if show_xaxis:
        ax.set_xticks(xs)
        ax.set_xticklabels(COND_LABELS, fontsize=8)
        ax.set_xlabel(PANEL_B_XLABEL, fontsize=18)
    else:
        # Hide x-axis entirely: ticks, tick labels, spine, label
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.set_xlabel("")
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(axis="x", length=0)

    if show_legend:
        ax.legend(loc="upper left", fontsize=10)


# ═══════════════════════════════════════════════════════════════════════════════
# PANEL B — Stacked PPC conditions: generated from scratch
# ═══════════════════════════════════════════════════════════════════════════════
print("[Panel B] Generating empiric + nominal PPC conditions from CSV data …")

df_emp = pd.read_csv(PANEL_B_EMP_CSV)
df_nom = pd.read_csv(PANEL_B_NOM_CSV)
print(f"  Empiric CSV: {df_emp['source'].value_counts().to_dict()}")
print(f"  Nominal CSV: {df_nom['source'].value_counts().to_dict()}")

rng_b = np.random.default_rng(42)   # same seed as ppc_classic_vch.py

fig_b, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(10, 5.5), sharex=False)
fig_b.subplots_adjust(hspace=0.12)

# Top subplot: empiric likelihood — legend, NO x-axis
_plot_ppc_conditions(ax_top, df_emp, rng_b, show_legend=True, show_xaxis=False)

# Bottom subplot: nominal likelihood — x-axis, NO legend
_plot_ppc_conditions(ax_bot, df_nom, rng_b, show_legend=False, show_xaxis=True)

# Shared y-axis label
fig_b.text(0.01, 0.5, PANEL_B_YLABEL, va='center', ha='left',
           rotation='vertical', fontsize=18)
# Remove per-subplot y-labels (shared label covers both)
ax_top.set_ylabel("")
ax_bot.set_ylabel("")

fig_b.savefig(_TMP_B, dpi=PANEL_DPI, bbox_inches='tight', facecolor='white')
plt.close(fig_b)
print(f"  Saved temp panel B: {_TMP_B}")

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL C — Out-of-set non-hallucinator group: per-block detection rate at 75%
# ═══════════════════════════════════════════════════════════════════════════════
print("\n[Panel C] Loading out-of-set trial-level data …")
df_oos = pd.read_csv(OOS_CSV)

# VCH modality, QC-filtered (same filter as empirical_condition_analysis.py)
df_vch = df_oos[(df_oos['modality'] == 'v')
                & df_oos['six_flag_qc_pass']
                & (df_oos['responseTime'] > 0)].copy()

# Non-hallucinator subset: max_vh_freq == 0 AND max_ah_freq == 0, NaN → 0
# (see README_empirical_condition_redetermination.md for rationale)
df_vch['_vh'] = df_vch['max_vh_freq'].fillna(0)
df_vch['_ah'] = df_vch['max_ah_freq'].fillna(0)
df_nonhall = df_vch[(df_vch['_vh'] == 0) & (df_vch['_ah'] == 0)].copy()
n_nonhall = df_nonhall['sudo_rec'].nunique()
print(f"  Non-hallucinator VCH trials: {len(df_nonhall)} from {n_nonhall} subjects")

# Condition = 75 only
df_panelc = df_nonhall[df_nonhall['condition'] == 75].copy()
print(f"  Condition-75 trials: {len(df_panelc)}")

# Derive block from trial number: block = ceil(trial / 30)
df_panelc['block'] = ((df_panelc['trial'] - 1) // TRIALS_PER_BLK) + 1

# Per-subject, per-block detection rate (mean of binary response)
subj_block = df_panelc.groupby(['sudo_rec', 'block'])['response'].mean().reset_index()
subj_block.rename(columns={'response': 'det_rate'}, inplace=True)

# Print trial counts per block to document sparsity
trial_counts = df_panelc.groupby('block')['trial'].count()
print(f"  Condition-75 trials per block (total across {n_nonhall} subjects):")
for b in range(1, N_BLOCKS + 1):
    cnt = trial_counts.get(b, 0)
    print(f"    Block {b:2d}: {cnt:4d} trials")

# Group mean and SEM across subjects, per block
block_means = []
block_ses   = []
for b in range(1, N_BLOCKS + 1):
    rates = subj_block[subj_block['block'] == b]['det_rate'].values
    block_means.append(np.mean(rates) if len(rates) > 0 else np.nan)
    block_ses.append(sp_stats.sem(rates) if len(rates) > 1 else 0.0)
block_means = np.array(block_means)
block_ses   = np.array(block_ses)

fig_c, ax_c = plt.subplots(figsize=(7, 5.5))

ax_c.errorbar(range(N_BLOCKS), block_means, yerr=block_ses * Z_94,
              color="black", linewidth=2, marker="o", markersize=5,
              capsize=3)

ax_c.set_xticks(list(range(N_BLOCKS)))
ax_c.set_xticklabels(BLOCK_LABELS)
ax_c.set_xlabel("Block")
ax_c.set_ylabel("Empiric Detection Probability\n(75% Intensity Condition)")
ax_c.set_ylim(-0.05, 1.05)
sns.despine(ax=ax_c)
fig_c.tight_layout()
fig_c.savefig(_TMP_C, dpi=PANEL_DPI, bbox_inches='tight', facecolor='white')
plt.close(fig_c)
print(f"  Saved temp panel C: {_TMP_C}")

# ═══════════════════════════════════════════════════════════════════════════════
# ASSEMBLE — 3-column layout via figure_assembly.py
# ═══════════════════════════════════════════════════════════════════════════════
print("\n[Assembly] Composing 3-column figure …")

config = {
    "output_path": OUT_PATH,
    "fig_width_inches": 14.0,
    "dpi": 300,
    "gap_inches": 0.08,
    "label_fontsize": 10,
    "label_color": "black",
    "rows": [
        {
            "panels": [
                {"path": PANEL_A_IMG, "label": "a"},
                {"path": _TMP_B,      "label": "b"},
                {"path": _TMP_C,      "label": "c"},
            ]
        }
    ]
}

result_path = assemble_manuscript_figure(config, base_dir=PROJECT_ROOT)
print(f"  Saved: {result_path}")

# Also save as TIFF using the same assembly function
tiff_config = dict(config)
tiff_config["output_path"] = OUT_PATH.replace('.png', '.tiff')
tiff_result_path = assemble_manuscript_figure(tiff_config, base_dir=PROJECT_ROOT)
print(f"  Saved: {tiff_result_path}")

# SVG: the assembled figure is a raster composite, so the SVG embeds it
# losslessly rather than vectorising it.
from figure_assembly import raster_to_svg
svg_out = OUT_PATH.replace('.png', '.svg')
raster_to_svg(result_path, svg_out)
print(f"  Saved: {svg_out}")

# Clean up temp files
for tmp in (_TMP_B, _TMP_C):
    if os.path.exists(tmp):
        os.remove(tmp)

print("Done.")
