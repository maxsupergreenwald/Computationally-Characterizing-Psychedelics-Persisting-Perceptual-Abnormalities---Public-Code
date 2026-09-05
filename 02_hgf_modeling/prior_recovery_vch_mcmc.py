"""
prior_recovery_vch_mcmc.py  ─  VCH HGF 4-Way Prior-Based Recovery (MCMC + Bridge Sampling)
============================================================================================

Master script for MCMC-based prior parameter and model recovery, extending the
MAP-based analysis in prior_recovery_vch.py with full MCMC fitting (NUTS, MAP
init, 4 chains × 1000 iterations) and bridge sampling LME (Meng & Wong, 1996;
Gronau et al., 2017) for model comparison.

For each of N_SIM=500 simulation indices × 4 true models, prior_recovery_vch_mcmc.jl:
  · samples generative parameters from the prior of the true_model
  · generates a 360-trial VCH sequence (fixed block multisets, within-block shuffled)
  · maps condition labels to intensities using the TRUE model's convention
  · simulates binary responses via the manual forward pass
  · fits ALL FOUR candidate models via MCMC (NUTS, MAP init, 4 chains × 1000),
    each with ITS OWN condition mapping
  · computes bridge sampling LME (MVN proposal in unconstrained space) per model
  · saves one result CSV per (sim_index, true_model)

Sections controlled by toggle flags:

  DO_GENERATE  Build dSQ job array + print transfer/submission commands
  DO_PULL      Print scp commands to pull results back from HPC
  DO_COMPILE   Stack sim*.csv → prior_recovery_mcmc_summary.csv; print confusion matrix
  DO_PLOT      Parameter recovery scatter + confusion heatmap + r comparison chart
  DO_ALL       Shortcut to enable all four sections
"""

################################################################################
# ─── ROLE IN THE MANUSCRIPT ───────────────────────────────────────────────────
#
# PRIOR-BASED RECOVERY, part 2 of 3.  Produces Supplementary Fig. S3a.
#
#     prior_recovery_vch_mcmc.jl (HPC) ──▶ THIS SCRIPT (compile + plot)
#
# Responsibilities
#     Builds the 2000-task dSQ array, prints the transfer/submit/pull commands,
#     stacks the per-simulation CSVs into one summary, and draws the
#     generative-vs-recovered scatter panels.
#
# Which figure was published
#     The scatter panels here are the source of Fig. S3a (assembled into a
#     2×2 layout by hppd_manuscript_master/03_visualization/supplement/
#     hgf_param_recovery_assembly.py).
#
#     The confusion matrix produced by THIS script uses bridge-sampling evidence.
#     The manuscript reports the BIC version instead — see
#     prior_recovery_aic_bic.py.  Both are computed from identical chains; the BIC
#     variant is the one that appears as Fig. S3b.
################################################################################

import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path

################################################################################
# ─── ACTIVE CONFIG ────────────────────────────────────────────────────────────
################################################################################

N_SIM       = 500                                          # simulations per true model
MODEL_TYPES = ["2level_empiric", "3level_empiric",
               "2level_nominal",  "3level_nominal"]        # all four true models

DO_ALL      = False
DO_GENERATE = True     # first run: build the 2000-task dSQ array
DO_PULL     = False
DO_COMPILE  = False    # enable after pulling results
DO_PLOT     = False    # enable after pulling results

if DO_ALL:
    DO_GENERATE = DO_PULL = DO_COMPILE = DO_PLOT = True

################################################################################
# ─── PATHS ────────────────────────────────────────────────────────────────────
################################################################################

# Resolved from this file's own location — see hgf_pipeline.py LOCAL PATHS.
LOCAL_JULIA_CH_DIR = os.path.dirname(os.path.abspath(__file__))
HPPD_PROJ_DIR      = os.path.dirname(LOCAL_JULIA_CH_DIR)

sys.path.insert(0, LOCAL_JULIA_CH_DIR)
sys.path.insert(0, os.path.join(HPPD_PROJ_DIR, "modules"))

from hgf_pipeline import (
    HPC_JULIA_CH_DIR,
    HPC_TRANSFER_HOST,
    HPC_LOGIN_HOST,
    HPC_USER,
    JULIA_MODULE_VCH,
)

LOCAL_PRIOR_MCMC_DIR = os.path.join(LOCAL_JULIA_CH_DIR, "param_recovery", "prior_based_mcmc")
COMMANDS_DIR = os.path.join(LOCAL_PRIOR_MCMC_DIR, "commands")
RESULTS_DIR  = os.path.join(LOCAL_PRIOR_MCMC_DIR, "results")
FIGURES_DIR  = os.path.join(LOCAL_PRIOR_MCMC_DIR, "figures")

HPC_PRIOR_MCMC_DIR = f"{HPC_JULIA_CH_DIR}/param_recovery/prior_based_mcmc"

# Walltime per dSQ task — each task runs 1 Julia call that fits 4 models
# via full MCMC (4 chains × 1000 iterations per model).
# Estimated ~60 min per call; 240:00 for safety margin.
WALLTIME = "240:00"

################################################################################
# ─── PARAMETER / LABEL METADATA ───────────────────────────────────────────────
################################################################################

MODEL_LABELS = {
    "2level_empiric" : "2-level\n(empiric)",
    "3level_empiric" : "3-level\n(empiric)",
    "2level_nominal" : "2-level\n(nominal)",
    "3level_nominal" : "3-level\n(nominal)",
}

MODEL_PARAMS = {
    "2level_empiric" : ["nu", "beta", "omega"],
    "3level_empiric" : ["nu", "beta", "omega", "omega3"],
    "2level_nominal" : ["nu", "beta", "omega"],
    "3level_nominal" : ["nu", "beta", "omega", "omega3"],
}
SHARED_PARAMS = ["nu", "beta", "omega"]

PARAM_SYMBOLS = {
    "nu"    : r"$\nu$",
    "beta"  : r"$\beta$",
    "omega" : r"$\omega$",
    "omega3": r"$\omega_3$",
}
PARAM_LABELS = {
    "nu"    : r"$\nu$ (prior weight)",
    "beta"  : r"$\beta$ (action precision)",
    "omega" : r"$\omega$ (volatility)",
    "omega3": r"$\omega_3$ (meta-volatility)",
}

SAME_LEVEL_PAIR = {
    "2level_empiric" : ("2level_empiric", "2level_nominal"),
    "3level_empiric" : ("3level_empiric", "3level_nominal"),
    "2level_nominal" : ("2level_nominal", "2level_empiric"),
    "3level_nominal" : ("3level_nominal", "3level_empiric"),
}

################################################################################
# ─── REFERENCE COMMANDS (always printed) ──────────────────────────────────────
################################################################################

_W = 70
print("=" * _W)
print("PRIOR-BASED RECOVERY — MCMC + BRIDGE SAMPLING (4-way) — KEY COMMANDS")
print("=" * _W)
print(f"""
  PUSH script + job array to HPC:
    scp {LOCAL_JULIA_CH_DIR}/prior_recovery_vch_mcmc.jl \\
        {HPC_USER}@{HPC_TRANSFER_HOST}:{HPC_JULIA_CH_DIR}/
    scp {COMMANDS_DIR}/prior_recovery_vch_mcmc.txt \\
        {HPC_USER}@{HPC_TRANSFER_HOST}:{HPC_PRIOR_MCMC_DIR}/commands/

  SUBMIT (run on bouchet):
    mkdir -p {HPC_PRIOR_MCMC_DIR}/{{commands,results,logs}}
    cd {HPC_PRIOR_MCMC_DIR}/commands
    module load dSQ
    dsq --job-file prior_recovery_vch_mcmc.txt --mem-per-cpu 4g -t {WALLTIME} --mail-type ALL
    sbatch dsq-prior_recovery_vch_mcmc-<date>.sh

  PULL results to local:
    scp -r {HPC_USER}@{HPC_TRANSFER_HOST}:{HPC_PRIOR_MCMC_DIR}/results/ \\
        {RESULTS_DIR}/
""")
print("=" * _W)
print()

################################################################################
# ─── DO_GENERATE: Build job array + print HPC commands ────────────────────────
################################################################################

if DO_GENERATE:
    os.makedirs(COMMANDS_DIR, exist_ok=True)

    # One Julia command per (sim_index × true_model) — 2000 total
    # NO bundling: each MCMC call takes ~60 min (4 full MCMC fits)
    job_lines = []
    for true_model in MODEL_TYPES:
        for sim_index in range(1, N_SIM + 1):
            cmd = (
                f"module load {JULIA_MODULE_VCH}; "
                f"cd {HPC_JULIA_CH_DIR}; "
                f"export GKSwstype=nul; "
                f"export JULIA_CPU_TARGET=generic; "
                f"julia --project=. --threads 4 prior_recovery_vch_mcmc.jl "
                f"{sim_index} {true_model}"
            )
            job_lines.append(cmd)

    fname = "prior_recovery_vch_mcmc.txt"
    fpath = os.path.join(COMMANDS_DIR, fname)
    dsq_cmd = (
        f"module load dSQ; dsq --job-file {fname} "
        f"--mem-per-cpu 4g -t {WALLTIME} --mail-type ALL"
    )
    with open(fpath, "w") as f:
        f.write(f"# to run: {dsq_cmd}\n")
        for line in job_lines:
            f.write(line + "\n")

    n_tasks = len(job_lines)
    print(f"Wrote: {fpath}")
    print(f"  ({N_SIM} sims × {len(MODEL_TYPES)} true models = {n_tasks} Julia calls "
          f"→ {n_tasks} SLURM tasks, walltime {WALLTIME} each)")

    print("\n" + "-" * 70)
    print("TRANSFER COMMANDS — copy-paste in order:")
    print("-" * 70)
    print(f"\n# ── Push updated prior_recovery_vch_mcmc.jl to HPC ─────────────────")
    print(f"scp {LOCAL_JULIA_CH_DIR}/prior_recovery_vch_mcmc.jl "
          f"{HPC_USER}@{HPC_TRANSFER_HOST}:{HPC_JULIA_CH_DIR}/")
    print(f"\n# ── Push job array ───────────────────────────────────────────────────")
    print(f"scp {COMMANDS_DIR}/{fname} "
          f"{HPC_USER}@{HPC_TRANSFER_HOST}:{HPC_PRIOR_MCMC_DIR}/commands/")
    print("\n" + "-" * 70)
    print("SUBMISSION COMMANDS — SSH in and run:")
    print("-" * 70)
    print(f"\nssh {HPC_USER}@{HPC_LOGIN_HOST}")
    print(f"mkdir -p {HPC_PRIOR_MCMC_DIR}/commands {HPC_PRIOR_MCMC_DIR}/results "
          f"{HPC_PRIOR_MCMC_DIR}/logs")
    print(f"cd {HPC_PRIOR_MCMC_DIR}/commands")
    print("module load dSQ")
    print(f"dsq --job-file {fname} --mem-per-cpu 4g -t {WALLTIME} --mail-type ALL")
    print(f"sbatch dsq-{fname.replace('.txt','')}-<date>.sh")

################################################################################
# ─── DO_PULL: Print scp commands to retrieve results ──────────────────────────
################################################################################

if DO_PULL:
    print("\n" + "=" * 70)
    print("PULL COMMANDS — run after HPC jobs complete:")
    print("=" * 70)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    hpc_src = f"{HPC_PRIOR_MCMC_DIR}/results/"
    print(f"\nscp -r {HPC_USER}@{HPC_TRANSFER_HOST}:{hpc_src} {RESULTS_DIR}/")

################################################################################
# ─── DO_COMPILE: Stack results + print 4×4 confusion matrix ──────────────────
################################################################################

if DO_COMPILE:
    if not os.path.isdir(RESULTS_DIR):
        print(f"\nResults directory not found: {RESULTS_DIR}")
        print("Run DO_PULL first to retrieve results from HPC.")
    else:
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
            print("\nNo result files found — have you pulled results from HPC?")
        else:
            summary = (pd.concat(dfs, ignore_index=True)
                         .sort_values(["true_model", "sim_index"])
                         .reset_index(drop=True))

            expected = N_SIM * len(MODEL_TYPES)
            print(f"\nCompiled {n_found}/{expected} expected result files "
                  f"({n_miss} missing)")

            out_path = os.path.join(RESULTS_DIR, "prior_recovery_mcmc_summary.csv")
            summary.to_csv(out_path, index=False)
            print(f"Saved: {out_path}")

            # ── 4×4 confusion matrix (counts + %) ─────────────────────────
            col_w  = 18
            header = f"{'':25s}" + "".join(f"{'winner='+m:>{col_w}s}" for m in MODEL_TYPES)
            print("\nModel identifiability confusion matrix (bridge sampling LME):")
            print("-" * (25 + col_w * len(MODEL_TYPES)))
            print(header)
            for true_model in MODEL_TYPES:
                sub = summary[summary["true_model"] == true_model]
                n   = len(sub)
                if n == 0:
                    continue
                row_str = f"  true={true_model:<20s}"
                for win_m in MODEL_TYPES:
                    cnt = (sub["winner"] == win_m).sum()
                    pct = 100 * cnt / n
                    row_str += f"  {cnt:3d}/{n} ({pct:5.1f}%){'':<1s}"
                print(row_str)

            valid   = summary["winner"].isin(MODEL_TYPES)
            correct = (summary.loc[valid, "true_model"] ==
                       summary.loc[valid, "winner"]).sum()
            total   = valid.sum()
            print(f"\n  Overall correct: {correct}/{total} ({100*correct/total:.1f}%)")

            # ── Convergence summary ──────────────────────────────────────
            print("\nMCMC convergence summary:")
            for fitted_model in MODEL_TYPES:
                col = f"converged_{fitted_model}"
                if col in summary.columns:
                    n_conv = summary[col].sum()
                    n_tot  = summary[col].notna().sum()
                    print(f"  {fitted_model:<20s}: {n_conv}/{n_tot} converged "
                          f"({100*n_conv/n_tot:.1f}%)")

################################################################################
# ─── DO_PLOT ──────────────────────────────────────────────────────────────────
################################################################################

if DO_PLOT:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.stats import pearsonr

    matplotlib.rcParams.update({
        "font.family"    : "Arial",
        "axes.labelsize" : 20,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "legend.fontsize": 13,
    })

    _PALETTE = ["#4472C4", "#ED7D31", "#70AD47", "#FFC000"]
    _COLOR_EMPIRIC = "#4472C4"
    _COLOR_NOMINAL = "#ED7D31"

    summary_path = os.path.join(RESULTS_DIR, "prior_recovery_mcmc_summary.csv")
    if not os.path.exists(summary_path):
        print("\nSummary CSV not found. Run DO_COMPILE first.")
    else:
        os.makedirs(FIGURES_DIR, exist_ok=True)
        summary = pd.read_csv(summary_path)
        print(f"\nPlotting {len(summary)} simulation results...")

        # ── Figure 1: Parameter recovery scatter per true_model ──────────
        for true_model in MODEL_TYPES:
            sub     = summary[summary["true_model"] == true_model].copy()
            params  = MODEL_PARAMS[true_model]
            mcmc_pfx = f"mcmc_{true_model}"

            fig, axes = plt.subplots(1, len(params), figsize=(4.5 * len(params), 4.5))
            if len(params) == 1:
                axes = [axes]

            for ax, p in zip(axes, params):
                gen_col = f"gen_{p}"
                rec_col = f"{mcmc_pfx}_{p}"
                if gen_col not in sub.columns or rec_col not in sub.columns:
                    ax.set_visible(False)
                    continue

                data = sub[[gen_col, rec_col]].dropna()
                if len(data) < 2:
                    ax.set_visible(False)
                    continue

                ax.scatter(data[gen_col], data[rec_col],
                           alpha=0.4, s=18, color="#1976D2", linewidths=0)

                vmin = min(data[gen_col].min(), data[rec_col].min())
                vmax = max(data[gen_col].max(), data[rec_col].max())
                ax.plot([vmin, vmax], [vmin, vmax], color="gray", lw=1.5, alpha=0.5)

                r, _ = pearsonr(data[gen_col], data[rec_col])
                ax.text(0.95, 0.05, f"$r = {r:.3f}$",
                        transform=ax.transAxes, ha="right", va="bottom", fontsize=14)
                ax.text(0.05, 0.95, PARAM_SYMBOLS.get(p, p),
                        transform=ax.transAxes, ha="left", va="top", fontsize=24)
                sns.despine(ax=ax)

            fig.supxlabel("Generative parameter", fontsize=20, y=0.01)
            fig.supylabel("MCMC posterior median", fontsize=20, x=0.01)
            title_label = MODEL_LABELS[true_model].replace("\n", " ")
            fig.suptitle(f"Parameter recovery (MCMC) — true model: {title_label}",
                         fontsize=14, y=1.02)
            fig.tight_layout(rect=[0.06, 0.06, 1, 1])
            out = os.path.join(FIGURES_DIR, f"prior_recovery_scatter_{true_model}.png")
            fig.savefig(out, dpi=150, bbox_inches="tight")
            plt.close(fig)
            print(f"  Saved: {os.path.basename(out)}")

        # ── Figure 2: 4×4 Confusion matrix heatmap ──────────────────────
        n_models = len(MODEL_TYPES)
        conf     = np.zeros((n_models, n_models))
        for i, true_m in enumerate(MODEL_TYPES):
            sub = summary[summary["true_model"] == true_m]
            for j, win_m in enumerate(MODEL_TYPES):
                conf[i, j] = (sub["winner"] == win_m).sum()

        row_totals = conf.sum(axis=1, keepdims=True)
        conf_pct   = np.where(row_totals > 0, 100 * conf / row_totals, 0)

        fig, ax = plt.subplots(figsize=(7, 6))
        im = ax.imshow(conf_pct, cmap="Blues", vmin=0, vmax=100)

        x_labels = [MODEL_LABELS[m] for m in MODEL_TYPES]
        y_labels  = [MODEL_LABELS[m] for m in MODEL_TYPES]

        ax.set_xticks(range(n_models))
        ax.set_yticks(range(n_models))
        ax.set_xticklabels(x_labels, fontsize=12)
        ax.set_yticklabels(y_labels, fontsize=12)
        ax.set_xlabel("Model selected by bridge sampling LME", fontsize=16)
        ax.set_ylabel("Data-generating model", fontsize=16)

        for i in range(n_models):
            for j in range(n_models):
                txt   = f"{conf_pct[i,j]:.1f}%\n(n={int(conf[i,j])})"
                color = "white" if conf_pct[i, j] > 60 else "black"
                ax.text(j, i, txt, ha="center", va="center",
                        fontsize=11, color=color,
                        fontweight="bold" if i == j else "normal")

        fig.colorbar(im, ax=ax, label="% simulations", shrink=0.85)
        fig.suptitle("Model identifiability — bridge sampling LME (MCMC)",
                     fontsize=13)
        fig.tight_layout()
        out = os.path.join(FIGURES_DIR, "model_identifiability_confusion.png")
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {os.path.basename(out)}")

        # ── Figure 3: Cross-condition recoverability comparison ──────────
        fig, axes = plt.subplots(2, 2, figsize=(11, 9), sharey=True)
        axes_flat = axes.flat

        for ax, true_model in zip(axes_flat, MODEL_TYPES):
            sub = summary[summary["true_model"] == true_model].copy()
            emp_key, nom_key = SAME_LEVEL_PAIR[true_model]
            params = MODEL_PARAMS[true_model]

            x      = np.arange(len(params))
            bar_w  = 0.35

            rs_emp = []
            rs_nom = []
            for p in params:
                gen_col  = f"gen_{p}"
                emp_col  = f"mcmc_{emp_key}_{p}"
                nom_col  = f"mcmc_{nom_key}_{p}"

                if gen_col in sub.columns and emp_col in sub.columns:
                    d = sub[[gen_col, emp_col]].dropna()
                    rs_emp.append(pearsonr(d[gen_col], d[emp_col])[0] if len(d) > 1 else np.nan)
                else:
                    rs_emp.append(np.nan)

                if gen_col in sub.columns and nom_col in sub.columns:
                    d = sub[[gen_col, nom_col]].dropna()
                    rs_nom.append(pearsonr(d[gen_col], d[nom_col])[0] if len(d) > 1 else np.nan)
                else:
                    rs_nom.append(np.nan)

            emp_label = MODEL_LABELS[emp_key].replace("\n", " ")
            nom_label = MODEL_LABELS[nom_key].replace("\n", " ")

            ax.bar(x - bar_w / 2, rs_emp, bar_w,
                   color=_COLOR_EMPIRIC, label=emp_label, edgecolor="none")
            ax.bar(x + bar_w / 2, rs_nom, bar_w,
                   color=_COLOR_NOMINAL, label=nom_label, edgecolor="none")

            ax.axhline(0, color="black", lw=0.8)
            ax.set_xticks(x)
            ax.set_xticklabels([PARAM_SYMBOLS.get(p, p) for p in params], fontsize=14)
            ax.set_ylim(-0.1, 1.05)
            ax.set_title(f"True: {MODEL_LABELS[true_model].replace(chr(10), ' ')}",
                         fontsize=13)
            ax.legend(fontsize=11, loc="lower right")
            sns.despine(ax=ax)

        fig.supylabel("Pearson r (gen vs MCMC median)", fontsize=16, x=0.01)
        fig.suptitle("Parameter recoverability (MCMC): empiric vs. nominal mapping",
                     fontsize=14)
        fig.tight_layout(rect=[0.04, 0, 1, 1])
        out = os.path.join(FIGURES_DIR, "param_recovery_by_mapping.png")
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {os.path.basename(out)}")

        # ── Copy key figures to model_comparison/figures/ ────────────────
        top_figures_dir = os.path.join(LOCAL_JULIA_CH_DIR, "model_comparison", "figures")
        os.makedirs(top_figures_dir, exist_ok=True)
        for fname_fig in [
            "prior_recovery_scatter_2level_empiric.png",
            "prior_recovery_scatter_3level_empiric.png",
            "prior_recovery_scatter_2level_nominal.png",
            "prior_recovery_scatter_3level_nominal.png",
            "model_identifiability_confusion.png",
            "param_recovery_by_mapping.png",
        ]:
            # Prefix with "mcmc_" to avoid overwriting MAP-based figures
            src = os.path.join(FIGURES_DIR, fname_fig)
            dst = os.path.join(top_figures_dir, f"mcmc_{fname_fig}")
            if os.path.exists(src):
                import shutil
                shutil.copy2(src, dst)
                print(f"  Copied: {fname_fig} → model_comparison/figures/mcmc_{fname_fig}")

        print(f"\nAll figures saved to: {FIGURES_DIR}")
