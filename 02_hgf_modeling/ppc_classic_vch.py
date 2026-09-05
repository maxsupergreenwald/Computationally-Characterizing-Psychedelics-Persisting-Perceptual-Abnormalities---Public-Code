### ppc_classic_vch.py  ─  VCH Classic PPC master script (local)
###
### Simulates ONE response sequence per participant at their posterior median
### parameters (no MAP, no Hessian, no covariance matrix), then plots group-mean
### detection rates with bootstrapped confidence intervals for both the real and
### simulated datasets.
###
### Flags:
###   DO_RUN     — load wide df, write medians CSVs, call ppc_classic_vch.jl
###   DO_COMPILE — stack per-participant CSVs into one summary file per model type
###   DO_PLOT    — generate figures (6 per model type)
###
### Julia skips participants whose output file already exists, so DO_RUN is safe
### to re-run after an interruption.

################################################################################
# ─── ROLE IN THE MANUSCRIPT ───────────────────────────────────────────────────
#
# POSTERIOR PREDICTIVE CHECKS, part 2 of 2.  Produces Supplementary Fig. S3f and
# the two PPC panels of Supplementary Fig. S2b.
#
#     ppc_classic_vch.jl ──▶ THIS SCRIPT (run / compile / plot)
#
# Responsibilities
#     Pulls each participant's posterior medians from the wide df, writes the
#     medians CSV that the Julia script consumes, invokes it once per model type,
#     then plots observed against simulated detection rates with bootstrapped
#     94% intervals plus per-participant spaghetti lines.
#
#     Six figures are produced per model type: per condition, per block pooled,
#     and per block within each of the four conditions.
#
# The trap when running the nominal variant
#     For "2level_nominal", PARAM_COLS points at the SAME wide-df columns as
#     "2level_empiric" (vch_beta / vch_nu / vch_omega), because
#     import_hgf_results_unified.py writes nominal results under the empiric
#     column names.  Those columns therefore hold whichever variant was imported
#     most recently.  Run one convention at a time — import, then PPC, then
#     import the other — or the nominal figures will be drawn from empiric
#     parameter estimates without any error being raised.
#
#     Published panels: 2level_empiric conditions + blocks are stacked into
#     Fig. S3f; 2level_empiric and 2level_nominal conditions are paired in
#     Fig. S2b.  Both assemblies live in hppd_manuscript_master/03_visualization/
#     supplement/.
################################################################################

import os
import sys
import subprocess
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

################################################################################
# ─── CONFIG ───────────────────────────────────────────────────────────────────
################################################################################

CI_PCT      = 94     # bootstrap confidence interval width (percent); easy to change
N_BOOTSTRAP = 10000  # bootstrap resamples for CI

# Model types to run.  For "2level_nominal", the wide-df columns (vch_beta,
# vch_nu, vch_omega) are shared with "2level_empiric" — they reflect whichever
# variant was last imported via import_hgf_results_unified.py.  Run each type
# separately if both nominal and empiric results have been imported.
MODEL_TYPES = ["2level_empiric", "2level_nominal", "3level_empiric"]

DO_RUN     = True
DO_COMPILE = True
DO_PLOT    = True

################################################################################
# ─── PATHS ────────────────────────────────────────────────────────────────────
################################################################################

SCRIPT_DIR    = os.path.dirname(os.path.abspath(__file__))
HPPD_PROJ_DIR = os.path.dirname(SCRIPT_DIR)   # public repo root

sys.path.insert(0, SCRIPT_DIR)

STIM_DIR     = os.path.join(SCRIPT_DIR, "data_n_cmnds", "hppd_manuscript", "hppd", "vch_data")
RESULTS_BASE = os.path.join(SCRIPT_DIR, "param_recovery", "ppc_classic", "results")
FIGURES_DIR  = os.path.join(SCRIPT_DIR, "param_recovery", "ppc_classic", "figures")
MEDIANS_DIR  = os.path.join(SCRIPT_DIR, "param_recovery", "ppc_classic", "medians")

os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(MEDIANS_DIR, exist_ok=True)

################################################################################
# ─── PARAMETER MAPPINGS ───────────────────────────────────────────────────────
# Julia parameter name → wide-df column name.
# Confirmed against import_hgf_results_unified.py and param_recovery_vch.py.
#
# Note on 2level_nominal: import_hgf_results_unified.py writes nominal results
# to the same wide-df column names as 2level_empiric (vch_beta/nu/omega) — only
# the output filename differs.  Chain imports carefully to avoid overwriting.
################################################################################

PARAM_COLS = {
    "2level_empiric": {
        "action_precision":       "vch_beta",
        "prior_posterior_weight": "vch_nu",
        "xprob_volatility":       "vch_omega",
    },
    "2level_nominal": {
        "action_precision":       "vch_beta",
        "prior_posterior_weight": "vch_nu",
        "xprob_volatility":       "vch_omega",
    },
    "3level_empiric": {
        "action_precision":       "vch_beta_3lev",
        "prior_posterior_weight": "vch_nu_3lev",
        "xprob_volatility":       "vch_omega_3lev",
        "xvol_volatility":        "vch_omega3_3lev",
    },
}

################################################################################
# ─── COLUMN SETS ──────────────────────────────────────────────────────────────
################################################################################

COND_COLS       = ["det_rate_0.0", "det_rate_0.25", "det_rate_0.5", "det_rate_0.75"]
BLOCK_COLS      = [f"block_{b}" for b in range(1, 13)]
BLOCK_COND0_COLS  = [f"block_{b}_cond0"  for b in range(1, 13)]
BLOCK_COND25_COLS = [f"block_{b}_cond25" for b in range(1, 13)]
BLOCK_COND50_COLS = [f"block_{b}_cond50" for b in range(1, 13)]
BLOCK_COND75_COLS = [f"block_{b}_cond75" for b in range(1, 13)]

################################################################################
# ─── DO_RUN ───────────────────────────────────────────────────────────────────
################################################################################

if DO_RUN:
    # Public, PII-free wide df read directly.  load_and_prepare_data() is NOT
    # called — the public file is already that function's output, and the fitted
    # HGF medians simulated from here are among the columns it contains.
    from hgf_pipeline import load_public_wide_df

    df, _wide_path = load_public_wide_df(project="hppd_manuscript")
    print(f"Loaded wide df: {len(df)} participants "
          f"({os.path.basename(_wide_path)})")

    # Availability check, PER MODEL TYPE.  A variant whose parameter columns are
    # absent from this wide df is reported by name and skipped; the remaining
    # variants still run.  Checking the union up front would let one unavailable
    # variant block every other one.
    _unavailable = {
        mt: [c for c in PARAM_COLS[mt].values() if c not in df.columns]
        for mt in MODEL_TYPES
    }
    _unavailable = {mt: cols for mt, cols in _unavailable.items() if cols}
    for mt, cols in _unavailable.items():
        print(f"  SKIPPING model type '{mt}': wide df has no {cols}.\n"
              f"    Run import_hgf_results_unified.py with MODEL_TYPE='{mt}' "
              f"to add them.")
    MODEL_TYPES = [mt for mt in MODEL_TYPES if mt not in _unavailable]
    if not MODEL_TYPES:
        raise KeyError(
            "No model type in MODEL_TYPES has its parameter columns present in "
            f"{os.path.basename(_wide_path)}; nothing to do."
        )

    staged_ids = {
        int(os.path.splitext(os.path.basename(f))[0])
        for f in glob.glob(os.path.join(STIM_DIR, "*.csv"))
    }
    print(f"Found {len(staged_ids)} staged stimulus CSVs in: {STIM_DIR}")

    for mt in MODEL_TYPES:
        param_cols    = PARAM_COLS[mt]
        required_cols = list(param_cols.values())

        # Filter to participants with complete medians AND a staged stimulus CSV
        has_params  = df[required_cols].notna().all(axis=1)
        has_stim    = df["record_id"].isin(staged_ids)
        eligible_df = df[has_params & has_stim].copy()

        n_no_params = (~has_params).sum()
        n_no_stim   = (has_params & ~has_stim).sum()
        print(f"\n[{mt}] {len(eligible_df)} eligible participants "
              f"({n_no_params} missing params, {n_no_stim} missing stimulus CSV)")

        if eligible_df.empty:
            print(f"  [WARNING] No eligible participants — skipping {mt}.")
            continue

        # Write medians CSV for Julia (columns: record_id + Julia param names)
        medians_rows = []
        for _, row in eligible_df.iterrows():
            mrow = {"record_id": int(row["record_id"])}
            for julia_name, df_col in param_cols.items():
                mrow[julia_name] = float(row[df_col])
            medians_rows.append(mrow)

        medians_df  = pd.DataFrame(medians_rows)
        medians_csv = os.path.join(MEDIANS_DIR, f"medians_{mt}.csv")
        medians_df.to_csv(medians_csv, index=False)
        print(f"  Wrote medians CSV: {medians_csv} ({len(medians_df)} rows)")

        # Create results directory
        os.makedirs(os.path.join(RESULTS_BASE, mt), exist_ok=True)

        # Call Julia
        cmd = [
            "julia", "--project=.", "--threads", "1",
            "ppc_classic_vch.jl", mt, medians_csv,
        ]
        print(f"\n{'='*60}")
        print(f"Running PPC classic — model={mt}")
        print(f"  {' '.join(cmd)}")
        print("=" * 60)
        result = subprocess.run(cmd, cwd=SCRIPT_DIR)
        if result.returncode != 0:
            print(f"[ERROR] ppc_classic_vch.jl {mt} exited with code {result.returncode}",
                  file=sys.stderr)

################################################################################
# ─── DO_COMPILE ───────────────────────────────────────────────────────────────
################################################################################

if DO_COMPILE:
    for mt in MODEL_TYPES:
        results_dir = os.path.join(RESULTS_BASE, mt)
        csv_files   = sorted(glob.glob(os.path.join(results_dir, "*_ppc_classic.csv")))
        if not csv_files:
            print(f"[COMPILE] No result files found in {results_dir}")
            continue
        combined = pd.concat([pd.read_csv(f) for f in csv_files], ignore_index=True)
        out_path = os.path.join(results_dir, f"ppc_classic_{mt}_all.csv")
        combined.to_csv(out_path, index=False)
        n_sub = combined["record_id"].nunique()
        n_emp = (combined["source"] == "empirical").sum()
        n_sim = (combined["source"] == "sim").sum()
        print(f"[COMPILE] {mt}: {n_sub} participants, "
              f"{n_emp} empirical + {n_sim} sim rows → {out_path}")

################################################################################
# ─── DO_PLOT ──────────────────────────────────────────────────────────────────
# For both real and simulated data:
#   group mean  = mean across all participants (NaN-ignoring)
#   CI          = bootstrapped by resampling participants with replacement
#
# Simulated data: one row per participant (source == "sim"), no sample_idx.
# Empirical data: one row per participant (source == "empirical").
################################################################################

if DO_PLOT:
    matplotlib.rcParams.update({
        "font.family":     "Arial",
        "axes.labelsize":  20,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "legend.fontsize": 14,
    })

    alpha_lo = (100 - CI_PCT) / 2 / 100   # e.g. 0.03 for 94% CI
    alpha_hi = 1.0 - alpha_lo             # e.g. 0.97

    _rng = np.random.default_rng(42)

    def bootstrap_group_ci(data_2d):
        """
        Bootstrap CI for the group mean across participants.

        Parameters
        ----------
        data_2d : (n_participants, n_cols) float array; NaN-safe.

        Returns
        -------
        group_mean, lo, hi — each shape (n_cols,)
        """
        n_pts = data_2d.shape[0]
        group_mean = np.nanmean(data_2d, axis=0)
        boot_means = np.full((N_BOOTSTRAP, data_2d.shape[1]), np.nan)
        for i in range(N_BOOTSTRAP):
            idx = _rng.integers(0, n_pts, size=n_pts)
            boot_means[i] = np.nanmean(data_2d[idx], axis=0)
        lo = np.nanquantile(boot_means, alpha_lo, axis=0)
        hi = np.nanquantile(boot_means, alpha_hi, axis=0)
        return group_mean, lo, hi

    def _plot_ppc(ax, cols, emp, sim):
        """Add simulated spaghetti + group mean + CI bands for both series."""
        xs = list(range(len(cols)))

        # ── Simulated spaghetti (behind everything) ────────────────────────────
        jitter_rng = np.random.default_rng(0)
        for _, row in sim[cols].iterrows():
            jitter = jitter_rng.uniform(-0.15, 0.15)
            ax.plot([x + jitter for x in xs], row.values.astype(float),
                    color="steelblue", alpha=0.03, linewidth=0.6)

        # ── Empirical group mean + CI ──────────────────────────────────────────
        emp_arr = emp[cols].values.astype(float)
        emp_mean, emp_lo, emp_hi = bootstrap_group_ci(emp_arr)
        ax.fill_between(xs, emp_lo, emp_hi, alpha=0.30, color="black")
        ax.plot(xs, emp_mean, color="black", linewidth=2,
                marker="o", markersize=6, label=f"Empirical (mean ± {CI_PCT}% CI)")

        # ── Simulated group mean + CI ──────────────────────────────────────────
        sim_arr = sim[cols].values.astype(float)
        sim_mean, sim_lo, sim_hi = bootstrap_group_ci(sim_arr)
        ax.fill_between(xs, sim_lo, sim_hi, alpha=0.30, color="steelblue")
        ax.plot(xs, sim_mean, color="steelblue", linewidth=2,
                linestyle="--", marker="o", markersize=6,
                label=f"Simulated (mean ± {CI_PCT}% CI)")

    def _save_ppc_fig(fig, ax, x_range, tick_labels, xlabel, ylabel,
                      legend_loc, fig_path):
        ax.set_xticks(list(x_range))
        ax.set_xticklabels(tick_labels, fontsize=8)
        ax.set_xlabel(xlabel, fontsize=10)
        ax.set_ylabel(ylabel)
        ax.set_ylim(-0.05, 1.05)
        ax.legend(loc=legend_loc)
        sns.despine(ax=ax)
        fig.tight_layout()
        fig.savefig(fig_path, dpi=150)
        plt.close(fig)
        print(f"[PLOT] Saved {fig_path}")

    COND_LABELS  = ["0%\nDetection\nProbability",
                    "25%\nDetection\nProbability",
                    "50%\nDetection\nProbability",
                    "75%\nDetection\nProbability"]
    BLOCK_LABELS = [str(b) for b in range(1, 13)]

    for mt in MODEL_TYPES:
        results_dir   = os.path.join(RESULTS_BASE, mt)
        combined_path = os.path.join(results_dir, f"ppc_classic_{mt}_all.csv")
        if not os.path.exists(combined_path):
            print(f"[PLOT] Compiled file not found: {combined_path} — run DO_COMPILE first")
            continue

        df_all  = pd.read_csv(combined_path)
        emp     = df_all[df_all["source"] == "empirical"].copy()
        sim     = df_all[df_all["source"] == "sim"].copy()
        fig_dir = os.path.join(FIGURES_DIR, mt)
        os.makedirs(fig_dir, exist_ok=True)
        print(f"\n[{mt}] Plotting {emp['record_id'].nunique()} participants → {fig_dir}")

        # ── Figure 1: per-condition detection rates ───────────────────────────
        fig, ax = plt.subplots(figsize=(7, 5))
        _plot_ppc(ax, COND_COLS, emp, sim)
        _save_ppc_fig(fig, ax, range(4), COND_LABELS,
                      "QUEST-Derived Stimulus Intensity (% Detection Probability)",
                      "Empiric Detection Probability",
                      "upper left",
                      os.path.join(fig_dir, f"ppc_classic_{mt}_conditions.png"))

        # ── Figure 2: per-block detection rates (all conditions pooled) ───────
        fig, ax = plt.subplots(figsize=(10, 5))
        _plot_ppc(ax, BLOCK_COLS, emp, sim)
        _save_ppc_fig(fig, ax, range(12), BLOCK_LABELS,
                      "Block", "Empiric Detection Probability",
                      "upper right",
                      os.path.join(fig_dir, f"ppc_classic_{mt}_blocks.png"))

        # ── Figures 3–6: per-block, condition-specific ────────────────────────
        for cond_label, cond_cols, stem in [
            ("0%",  BLOCK_COND0_COLS,  "blocks_cond0"),
            ("25%", BLOCK_COND25_COLS, "blocks_cond25"),
            ("50%", BLOCK_COND50_COLS, "blocks_cond50"),
            ("75%", BLOCK_COND75_COLS, "blocks_cond75"),
        ]:
            if not all(c in df_all.columns for c in cond_cols):
                print(f"[PLOT] {cond_label} columns missing in {combined_path} — skipping")
                continue
            fig, ax = plt.subplots(figsize=(10, 5))
            _plot_ppc(ax, cond_cols, emp, sim)
            _save_ppc_fig(fig, ax, range(12), BLOCK_LABELS,
                          "Block",
                          f"Empiric Detection Probability ({cond_label})",
                          "upper right",
                          os.path.join(fig_dir, f"ppc_classic_{mt}_{stem}.png"))

print("\nDone.")
