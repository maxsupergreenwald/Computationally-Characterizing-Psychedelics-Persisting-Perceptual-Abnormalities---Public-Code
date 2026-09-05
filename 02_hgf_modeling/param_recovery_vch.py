"""
param_recovery_vch.py  ─  VCH HGF Parameter Recovery
======================================================

Master script for VCH HGF parameter recovery for the HPPD manuscript.

For each participant with non-missing HGF estimates (vch_nu, vch_beta, vch_omega
for 2-level; vch_nu_3lev, vch_beta_3lev, vch_omega_3lev, vch_omega3_3lev for
3-level), this script uses those fitted parameter values as "generative" ground
truth, simulates synthetic VCH responses 10 times per participant, and refits the
HGF to each simulated dataset.

Five sections controlled by toggle flags:

  DO_GENERATE  Build dSQ job array .txt files + print transfer/submission commands
  DO_PULL      Print scp commands to pull results back from HPC
  DO_COMPILE   Read *_medians.csv + *_rhats.csv → summary CSV; compute per-run
               beta/nu posterior correlation from chains CSV; flag bad-rhat runs
  DO_PLOT      Generate figures (only converged runs are plotted):
                 · Scatter: generative vs recovered per parameter (all iterations)
                 · Box plot: pooled beta/nu posterior correlation coefficients
                 · 3D scatter: 3 plots of generative vs recovered nu/beta
  DO_ALL       Shortcut to enable all four sections

Parameter column names (confirmed against actual HPC output files):

  _medians.csv  : action_precision, prior_posterior_weight,
                  xprob_volatility [, xvol_volatility]
  _rhats.csv    : action_precision.session[1], prior_posterior_weight.session[1],
                  xprob_volatility.session[1] [, xvol_volatility.session[1]]
  chains .csv   : iteration, chain,
                  action_precision.session[1], prior_posterior_weight.session[1],
                  xprob_volatility.session[1] [, xvol_volatility.session[1]], lp, ...

Wide-df → Julia parameter mapping:
  2-level: vch_nu → prior_posterior_weight,  vch_beta → action_precision,
           vch_omega → xprob_volatility
  3-level: vch_nu_3lev, vch_beta_3lev, vch_omega_3lev, vch_omega3_3lev

HOW TO USE
----------
1. Run DO_GENERATE (default) to build job arrays and print HPC commands.
2. Transfer Julia script + job arrays to HPC, submit, wait for completion.
3. Pull results with DO_PULL commands.
4. Set DO_COMPILE=True to build summary CSV (includes rhat QC).
5. Set DO_PLOT=True to generate figures (bad-rhat runs automatically excluded).
"""

################################################################################
# ─── ROLE IN THE MANUSCRIPT ───────────────────────────────────────────────────
#
# POSTERIOR-BASED RECOVERY, part 2 of 2.  Produces Supplementary Fig. S3d and S3e.
#
#     param_recovery_vch.jl (HPC) ──▶ THIS SCRIPT (generate / compile / plot)
#
# Responsibilities
#     Reads each participant's fitted medians from the wide df, emits the dSQ job
#     array that simulates and refits from them, then compiles the recovered
#     medians and draws the figures.
#
#     Panel d of the assembled supplementary figure is corr_gen_vs_rec_2level_
#     empiric.png; panel e is pair_beta_nu_2level_empiric.png.
#
# Convergence gate — note the asymmetry
#     A recovery run with ANY parameter R-hat outside 0.9–1.1 is written to the
#     summary CSV flagged converged = False, and then excluded from every figure.
#     The flag is preserved in the CSV so the exclusion rate can be inspected; the
#     figures alone will not reveal it.
################################################################################

import os
import sys
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

################################################################################
# ─── ACTIVE CONFIG — edit these on each run ───────────────────────────────────
################################################################################

MODEL_TYPES  = ["2level_empiric", "3level_empiric"]   # Julia HGF model variants to recover
N_ITERATIONS = 10                     # recovery iterations per participant

# Rhat convergence threshold (consistent with import_hgf_results_unified.py)
RHAT_LOWER = 0.9
RHAT_UPPER = 1.1

DO_ALL      = False   # shortcut: enables all sections below
DO_GENERATE = True    # generate job arrays + print HPC commands
DO_PULL     = False   # print scp pull commands
DO_COMPILE  = False   # compile results → summary CSV with rhat QC (after pull)
DO_PLOT     = False   # generate figures (only converged runs plotted; after pull)

if DO_ALL:
    DO_GENERATE = DO_PULL = DO_COMPILE = DO_PLOT = True

################################################################################
# ─── PATHS ────────────────────────────────────────────────────────────────────
################################################################################

# Resolved from this file's own location — see hgf_pipeline.py LOCAL PATHS.
LOCAL_JULIA_CH_DIR = os.path.dirname(os.path.abspath(__file__))
HPPD_PROJ_DIR      = os.path.dirname(LOCAL_JULIA_CH_DIR)

sys.path.insert(0, LOCAL_JULIA_CH_DIR)

from hgf_pipeline import (
    PROJECT_CONFIGS,
    HPC_JULIA_CH_DIR,
    HPC_TRANSFER_HOST,
    HPC_LOGIN_HOST,
    HPC_USER,
    JULIA_MODULE_VCH,
    load_public_wide_df,
)

cfg       = PROJECT_CONFIGS["hppd_manuscript"]
PROJ_KEY  = cfg["project_hpc_key"]   # "hppd_manuscript"
TIMEPOINT = cfg["timepoints_vch"][0]  # "hppd"

LOCAL_PARAM_RECOVERY_DIR = os.path.join(LOCAL_JULIA_CH_DIR, "param_recovery")
COMMANDS_DIR = os.path.join(LOCAL_PARAM_RECOVERY_DIR, "commands")
RESULTS_DIR  = os.path.join(LOCAL_PARAM_RECOVERY_DIR, "results")
FIGURES_DIR  = os.path.join(LOCAL_PARAM_RECOVERY_DIR, "figures")

HPC_PARAM_RECOVERY_DIR = f"{HPC_JULIA_CH_DIR}/param_recovery"

# Walltime per dSQ task (2 bundled Julia calls) — mirrors main-analysis limits
VCH_TIMELIMITS = {
    "2level_empiric": "60:00",
    "3level_empiric": "180:00",
}

################################################################################
# ─── REFERENCE COMMANDS (always printed regardless of flags) ──────────────────
################################################################################

_W = 70
print("=" * _W)
print("POSTERIOR-BASED RECOVERY — KEY COMMANDS")
print("=" * _W)
print(f"""
  PUSH script + job arrays to HPC:
    scp {LOCAL_JULIA_CH_DIR}/param_recovery_vch.jl \\
        {HPC_USER}@{HPC_TRANSFER_HOST}:{HPC_JULIA_CH_DIR}/
    scp {COMMANDS_DIR}/param_recovery_2level.txt \\
        {HPC_USER}@{HPC_TRANSFER_HOST}:{HPC_PARAM_RECOVERY_DIR}/commands/
    scp {COMMANDS_DIR}/param_recovery_3level.txt \\
        {HPC_USER}@{HPC_TRANSFER_HOST}:{HPC_PARAM_RECOVERY_DIR}/commands/

  SUBMIT (run on bouchet):
    mkdir -p {HPC_PARAM_RECOVERY_DIR}/{{commands,results,logs}}
    cd {HPC_PARAM_RECOVERY_DIR}/commands
    module load dSQ
    dsq --job-file param_recovery_2level.txt --mem-per-cpu 4g -t 60:00 --mail-type ALL
    sbatch dsq-param_recovery_2level-<date>.sh
    dsq --job-file param_recovery_3level.txt --mem-per-cpu 4g -t 180:00 --mail-type ALL
    sbatch dsq-param_recovery_3level-<date>.sh

  PULL results to local:
    scp -r {HPC_USER}@{HPC_TRANSFER_HOST}:{HPC_PARAM_RECOVERY_DIR}/results/2level/ \\
        {RESULTS_DIR}/2level/
    scp -r {HPC_USER}@{HPC_TRANSFER_HOST}:{HPC_PARAM_RECOVERY_DIR}/results/3level/ \\
        {RESULTS_DIR}/3level/
""")
print("=" * _W)
print()

################################################################################
# ─── PARAMETER MAPPINGS ───────────────────────────────────────────────────────
# Confirmed against:
#   · results/vch/{2level,3level}/hppd/vch/*_medians.csv
#   · import_hgf_results_unified.py  PARAM_NAMES / RHAT_PARAM_NAMES
#   · variable_labels.py  iv_type_dict["vch_computations"] / ["vch_comps_3lev"]
################################################################################

# Julia parameter name → wide-df column name (source of generative values)
PARAM_COLS = {
    "2level_empiric": {
        "prior_posterior_weight": "vch_nu",
        "action_precision":       "vch_beta",
        "xprob_volatility":       "vch_omega",
    },
    "3level_empiric": {
        "prior_posterior_weight": "vch_nu_3lev",
        "action_precision":       "vch_beta_3lev",
        "xprob_volatility":       "vch_omega_3lev",
        "xvol_volatility":        "vch_omega3_3lev",
    },
}

# Julia parameter name → short label for summary CSV columns and plot axes
PARAM_SHORT = {
    "prior_posterior_weight": "nu",
    "action_precision":       "beta",
    "xprob_volatility":       "omega",
    "xvol_volatility":        "omega3",
}

# Chains CSV column → short label (for posterior correlation computation)
CHAINS_COL_MAP = {
    "action_precision.session[1]":       "beta",
    "prior_posterior_weight.session[1]": "nu",
    "xprob_volatility.session[1]":       "omega",
    "xvol_volatility.session[1]":        "omega3",
}

# Rhat CSV column → short label
RHAT_COL_MAP = {
    "action_precision.session[1]":       "beta",
    "prior_posterior_weight.session[1]": "nu",
    "xprob_volatility.session[1]":       "omega",
    "xvol_volatility.session[1]":        "omega3",
}

################################################################################
# ─── LOAD WIDE DF ─────────────────────────────────────────────────────────────
# Required for DO_GENERATE (participant selection) and DO_COMPILE (ground truth).
#
# Reads the public, PII-free wide df directly.  load_and_prepare_data() is NOT
# called: the public file is already that function's output, so the derived
# columns exist and re-deriving them is neither possible nor correct.  The fitted
# HGF medians used here as generative ground truth are among those columns.
#
# `require` lists the parameter columns every selected model type needs, so a
# missing column fails immediately, naming itself, rather than surfacing later
# as an opaque KeyError during participant selection.
################################################################################

need_df = DO_GENERATE or DO_COMPILE
if need_df:
    df, _wide_path = load_public_wide_df(project="hppd_manuscript")
    print(f"Loaded wide df: {len(df)} participants "
          f"({os.path.basename(_wide_path)})")

    # Availability check, PER MODEL TYPE — see the equivalent block in
    # ppc_classic_vch.py.  A variant with absent parameter columns is named and
    # skipped rather than aborting the whole run.
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

################################################################################
# ─── DO_GENERATE: Build job arrays and print HPC commands ─────────────────────
################################################################################

if DO_GENERATE:
    os.makedirs(COMMANDS_DIR, exist_ok=True)

    # Staged stimulus CSVs — written by hgf_pipeline.py step 2, already on HPC.
    # We check locally that the staging dir exists so we can list eligible IDs.
    stim_dir = os.path.join(
        LOCAL_JULIA_CH_DIR, "data_n_cmnds", PROJ_KEY, TIMEPOINT, "vch_data"
    )
    if not os.path.isdir(stim_dir):
        raise FileNotFoundError(
            f"VCH staging directory not found: {stim_dir}\n"
            "Run hgf_pipeline.py (step 2) first to stage the stimulus CSVs."
        )

    staged_ids = {int(Path(f).stem) for f in Path(stim_dir).glob("*.csv")}
    print(f"\nFound {len(staged_ids)} staged VCH stimulus CSVs in:\n  {stim_dir}")

    for model_type in MODEL_TYPES:
        param_cols    = PARAM_COLS[model_type]       # julia_name → df_col_name
        required_cols = list(param_cols.values())

        # Participants with complete (non-NaN) parameter estimates
        has_params   = df[required_cols].notna().all(axis=1)
        eligible_df  = df[has_params].copy()
        print(f"\n[{model_type}]  {has_params.sum()}/{len(df)} participants "
              f"have complete parameters: {required_cols}")

        # Further filter: must have a staged stimulus CSV on disk
        has_stim    = eligible_df["record_id"].isin(staged_ids)
        n_no_stim   = (~has_stim).sum()
        eligible_df = eligible_df[has_stim].copy()
        if n_no_stim > 0:
            print(f"  ({n_no_stim} excluded — no staged stimulus CSV found)")
        print(f"  → {len(eligible_df)} participants eligible for recovery")

        if eligible_df.empty:
            print(f"  [WARNING] No eligible participants — skipping {model_type}.")
            continue

        # Build one Julia command per (participant × iteration)
        job_lines = []
        for _, row in eligible_df.iterrows():
            rid = int(row["record_id"])
            nu    = row[param_cols["prior_posterior_weight"]]
            beta  = row[param_cols["action_precision"]]
            omega = row[param_cols["xprob_volatility"]]

            for iteration in range(1, N_ITERATIONS + 1):
                if model_type == "3level_empiric":
                    omega3 = row[param_cols["xvol_volatility"]]
                    args = (f"{rid} {model_type} {iteration} "
                            f"{nu:.10f} {beta:.10f} {omega:.10f} {omega3:.10f}")
                else:
                    args = (f"{rid} {model_type} {iteration} "
                            f"{nu:.10f} {beta:.10f} {omega:.10f}")

                cmd = (
                    f"module load {JULIA_MODULE_VCH}; "
                    f"cd {HPC_JULIA_CH_DIR}; "
                    f"export GKSwstype=nul; "
                    f"julia --project=. --threads 4 param_recovery_vch.jl {args}"
                )
                job_lines.append(cmd)

        # Write job array: 2 Julia calls bundled per SLURM task (same pattern
        # as hgf_pipeline.py step 3)
        fname     = f"param_recovery_{model_type}.txt"
        fpath     = os.path.join(COMMANDS_DIR, fname)
        timelimit = VCH_TIMELIMITS.get(model_type, "120:00")
        dsq_cmd   = (
            f"module load dSQ; dsq --job-file {fname} "
            f"--mem-per-cpu 4g -t {timelimit} --mail-type ALL"
        )

        with open(fpath, "w") as f:
            f.write(f"# to run: {dsq_cmd}\n")
            for i in range(0, len(job_lines), 2):
                bundle = job_lines[i : i + 2]
                f.write("; ".join(bundle) + "\n")

        n_jobs  = len(job_lines)
        n_tasks = (n_jobs + 1) // 2
        print(f"  Wrote: {fpath}")
        print(f"  ({len(eligible_df)} participants × {N_ITERATIONS} iters "
              f"= {n_jobs} Julia calls → {n_tasks} SLURM tasks)")

    # ── Transfer commands ──────────────────────────────────────────────────────
    print("\n" + "-" * 70)
    print("TRANSFER COMMANDS — copy-paste in order:")
    print("-" * 70)

    print(f"\n# ── Transfer param_recovery_vch.jl to HPC ──────────────────────────")
    print(f"scp {LOCAL_JULIA_CH_DIR}/param_recovery_vch.jl "
          f"{HPC_USER}@{HPC_TRANSFER_HOST}:{HPC_JULIA_CH_DIR}/")

    print(f"\n# ── Transfer job array files ─────────────────────────────────────────")
    for model_type in MODEL_TYPES:
        fname = f"param_recovery_{model_type}.txt"
        print(f"scp {COMMANDS_DIR}/{fname} "
              f"{HPC_USER}@{HPC_TRANSFER_HOST}:{HPC_PARAM_RECOVERY_DIR}/commands/")

    # ── Submission commands ────────────────────────────────────────────────────
    print("\n" + "-" * 70)
    print("SUBMISSION COMMANDS — SSH in and run:")
    print("-" * 70)
    print(f"\nssh {HPC_USER}@{HPC_LOGIN_HOST}")
    print(f"mkdir -p {HPC_PARAM_RECOVERY_DIR}/commands {HPC_PARAM_RECOVERY_DIR}/results")
    print(f"cd {HPC_PARAM_RECOVERY_DIR}/commands")
    for model_type in MODEL_TYPES:
        fname     = f"param_recovery_{model_type}.txt"
        timelimit = VCH_TIMELIMITS.get(model_type, "120:00")
        print(f"\n# ── {model_type} ───────────────────────────────────────────────────")
        print(f"module load dSQ")
        print(f"dsq --job-file {fname} --mem-per-cpu 4g -t {timelimit} --mail-type ALL")

################################################################################
# ─── DO_PULL: Print scp commands to retrieve results ──────────────────────────
################################################################################

if DO_PULL:
    print("\n" + "=" * 70)
    print("PULL COMMANDS — run these AFTER HPC jobs complete:")
    print("=" * 70)
    for model_type in MODEL_TYPES:
        local_dst = os.path.join(RESULTS_DIR, model_type)
        os.makedirs(local_dst, exist_ok=True)
        hpc_src = f"{HPC_PARAM_RECOVERY_DIR}/results/{model_type}/"
        print(f"\n# ── {model_type} ───────────────────────────────────────────────────")
        print(f"scp -r {HPC_USER}@{HPC_TRANSFER_HOST}:{hpc_src} {local_dst}/")

################################################################################
# ─── DO_COMPILE: Build summary CSV with rhat QC ───────────────────────────────
#
# For each recovered run (record_id × iteration):
#   · Reads _medians.csv → recovered parameter medians
#   · Reads _rhats.csv   → checks all params pass RHAT_LOWER < r < RHAT_UPPER
#   · Reads chains .csv  → computes Pearson r(beta, nu) across all MCMC samples
#   · Joins ground-truth values from wide df
#   · Flags run as converged = True/False
#
# Non-converged runs are included in the summary CSV (with converged=False)
# so you can inspect them, but DO_PLOT silently excludes them.
#
# Summary CSV columns:
#   record_id, model_type, iteration,
#   gen_nu, gen_beta, gen_omega [, gen_omega3],
#   rec_nu, rec_beta, rec_omega [, rec_omega3],
#   rhat_nu, rhat_beta, rhat_omega [, rhat_omega3],
#   beta_nu_chain_r,
#   converged                   ← False if ANY param has rhat outside threshold
################################################################################

if DO_COMPILE:
    from scipy.stats import pearsonr

    for model_type in MODEL_TYPES:
        results_subdir = os.path.join(RESULTS_DIR, model_type)
        if not os.path.isdir(results_subdir):
            print(f"\n[{model_type}] Results directory not found: {results_subdir}")
            print(f"  Have you pulled results from HPC? Skipping.")
            continue

        param_cols  = PARAM_COLS[model_type]
        julia_names = list(param_cols.keys())

        # All _medians.csv files in this model type's results directory
        medians_files = sorted(Path(results_subdir).glob("*_iter*_medians.csv"))
        if not medians_files:
            print(f"\n[{model_type}] No *_iter*_medians.csv files found in "
                  f"{results_subdir}")
            print(f"  Have you pulled results from HPC? Skipping.")
            continue

        print(f"\n[{model_type}] Compiling {len(medians_files)} medians files...")

        rows = []
        n_bad_rhat = 0

        for mfile in medians_files:
            try:
                mdf = pd.read_csv(mfile)
            except Exception as exc:
                print(f"  Warning: could not read {mfile.name}: {exc}")
                continue

            if mdf.empty:
                continue

            record_id = int(mdf["record_id"].iloc[0])
            iteration = int(mdf["iteration"].iloc[0])
            row       = {
                "record_id":  record_id,
                "model_type": model_type,
                "iteration":  iteration,
            }

            # ── Recovered medians ──────────────────────────────────────────────
            # Medians CSV columns: action_precision, prior_posterior_weight,
            # xprob_volatility [, xvol_volatility]
            for julia_name in julia_names:
                short = PARAM_SHORT[julia_name]
                row[f"rec_{short}"] = (
                    float(mdf[julia_name].iloc[0])
                    if julia_name in mdf.columns else np.nan
                )

            # ── Ground-truth from wide df ──────────────────────────────────────
            wide_row = df[df["record_id"] == record_id]
            if wide_row.empty:
                print(f"  Warning: record_id {record_id} not found in wide df — "
                      f"skipping.")
                continue
            for julia_name, df_col in param_cols.items():
                short = PARAM_SHORT[julia_name]
                row[f"gen_{short}"] = (
                    float(wide_row[df_col].iloc[0])
                    if df_col in wide_row.columns else np.nan
                )

            # ── Rhat QC ────────────────────────────────────────────────────────
            # Rhat CSV columns: action_precision.session[1],
            # prior_posterior_weight.session[1], xprob_volatility.session[1]
            # [, xvol_volatility.session[1]], record_id, iteration
            rhat_file = mfile.parent / mfile.name.replace("_medians.csv", "_rhats.csv")
            converged = True   # assume good until proven otherwise

            if rhat_file.exists():
                try:
                    rdf = pd.read_csv(rhat_file)
                    for rhat_col, short in RHAT_COL_MAP.items():
                        if rhat_col not in rdf.columns:
                            continue
                        rhat_val = float(rdf[rhat_col].iloc[0])
                        row[f"rhat_{short}"] = rhat_val
                        if not (RHAT_LOWER <= rhat_val <= RHAT_UPPER):
                            converged = False
                except Exception as exc:
                    print(f"  Warning: could not read rhats for "
                          f"record_id={record_id} iter={iteration}: {exc}")
                    converged = False
            else:
                print(f"  Warning: no rhat file for "
                      f"record_id={record_id} iter={iteration} — marking as "
                      f"non-converged.")
                converged = False

            if not converged:
                n_bad_rhat += 1

            row["converged"] = converged

            # ── Posterior beta / nu correlation ────────────────────────────────
            # Chains CSV columns: action_precision.session[1],
            # prior_posterior_weight.session[1] (+ others we don't need here)
            chains_csv = mfile.parent / mfile.name.replace("_medians.csv", ".csv")
            row["beta_nu_chain_r"] = np.nan

            if chains_csv.exists():
                try:
                    cdf = pd.read_csv(
                        chains_csv,
                        usecols=[
                            "action_precision.session[1]",
                            "prior_posterior_weight.session[1]",
                        ],
                    )
                    beta_s = cdf["action_precision.session[1]"].dropna()
                    nu_s   = cdf["prior_posterior_weight.session[1]"].dropna()
                    if len(beta_s) > 2 and len(nu_s) > 2:
                        r, _ = pearsonr(beta_s, nu_s)
                        row["beta_nu_chain_r"] = float(r)
                except Exception as exc:
                    print(f"  Warning: chains correlation failed for "
                          f"record_id={record_id} iter={iteration}: {exc}")
            else:
                print(f"  Warning: chains CSV not found for "
                      f"record_id={record_id} iter={iteration} — "
                      f"beta_nu_chain_r will be NaN.")

            rows.append(row)

        if not rows:
            print(f"[{model_type}] No valid rows compiled.")
            continue

        summary_df = (
            pd.DataFrame(rows)
            .sort_values(["record_id", "iteration"])
            .reset_index(drop=True)
        )

        n_total     = len(summary_df)
        n_converged = summary_df["converged"].sum()
        print(f"  Compiled {n_total} runs: "
              f"{n_converged} converged, {n_bad_rhat} bad-rhat (excluded from plots)")

        out_path = os.path.join(
            RESULTS_DIR, model_type,
            f"param_recovery_summary_{model_type}.csv"
        )
        summary_df.to_csv(out_path, index=False)
        print(f"  Summary saved: {out_path}")

################################################################################
# ─── DO_PLOT: Generate figures (converged runs only) ──────────────────────────
#
# Figures per model_type:
#   1. corr_gen_vs_rec_{model_type}.png
#        Scatter: generative vs recovered for each parameter.
#        All 10 iterations as separate points, colored by participant.
#        Identity line + regression line + Pearson r in title.
#   2. posterior_corr_boxplot_beta_nu_{model_type}.png
#        Box plot of beta/nu posterior correlation coefficients (pooled across
#        all converged participant × iteration runs).
#   3. 3d_gennu_recnu_recbeta_{model_type}.png     [x=gen_nu, y=rec_nu, z=rec_beta]
#   4. 3d_genbeta_recbeta_recnu_{model_type}.png   [x=gen_beta, y=rec_beta, z=rec_nu]
#   5. 3d_gennu_recnu_genbeta_{model_type}.png     [x=gen_nu, y=rec_nu, z=gen_beta]
#
# Non-converged runs (converged=False) are silently excluded from all figures.
################################################################################

if DO_PLOT:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    import plotly.graph_objects as go
    from scipy.stats import pearsonr

    matplotlib.rcParams.update({
        "font.family":     "Arial",
        "axes.labelsize":  20,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "legend.fontsize": 14,
    })

    # Italic parameter symbols for subplot corner annotations
    PARAM_SYMBOLS = {
        "nu":     r"$\nu$",
        "beta":   r"$\beta$",
        "omega":  r"$\omega$",
        "omega3": r"$\omega_3$",
    }

    os.makedirs(FIGURES_DIR, exist_ok=True)

    for model_type in MODEL_TYPES:
        summary_path = os.path.join(
            RESULTS_DIR, model_type,
            f"param_recovery_summary_{model_type}.csv"
        )
        if not os.path.exists(summary_path):
            print(f"\n[{model_type}] Summary CSV not found — run DO_COMPILE first.")
            continue

        full_summary = pd.read_csv(summary_path)
        n_total      = len(full_summary)
        summary      = full_summary[full_summary["converged"] == True].copy()
        n_excluded   = n_total - len(summary)

        print(f"\n[{model_type}] Plotting {len(summary)}/{n_total} converged runs "
              f"({n_excluded} excluded for bad rhat)")

        if summary.empty:
            print(f"  [WARNING] No converged runs to plot — skipping {model_type}.")
            continue

        param_julia_names = list(PARAM_COLS[model_type].keys())
        short_names       = [PARAM_SHORT[p] for p in param_julia_names]
        n_params          = len(short_names)

        # Unique participant colour map
        part_ids  = sorted(summary["record_id"].unique())
        cmap      = plt.get_cmap("viridis", len(part_ids))
        id_colour = {pid: cmap(i) for i, pid in enumerate(part_ids)}

        # ── Figure 1: Generative vs Recovered scatter ────────────────────────
        fig, axes = plt.subplots(1, n_params, figsize=(5 * n_params, 5))
        if n_params == 1:
            axes = [axes]

        for ax, short in zip(axes, short_names):
            gen_col = f"gen_{short}"
            rec_col = f"rec_{short}"
            if gen_col not in summary.columns or rec_col not in summary.columns:
                ax.set_visible(False)
                continue

            sub = summary[[gen_col, rec_col, "record_id"]].dropna()
            if sub.empty:
                ax.set_visible(False)
                continue

            colours = [id_colour[pid] for pid in sub["record_id"]]
            ax.scatter(sub[gen_col], sub[rec_col],
                       c=colours, alpha=0.1, s=3, linewidths=0)

            vmin = min(sub[gen_col].min(), sub[rec_col].min())
            vmax = max(sub[gen_col].max(), sub[rec_col].max())
            ax.plot([vmin, vmax], [vmin, vmax], "k--", lw=1, label="Identity")

            coeffs = np.polyfit(sub[gen_col], sub[rec_col], 1)
            xs = np.linspace(vmin, vmax, 200)
            ax.plot(xs, np.polyval(coeffs, xs), "r-", lw=1.5, label="OLS")

            r, _ = pearsonr(sub[gen_col], sub[rec_col])
            ax.text(0.95, 0.05, f"$r = {r:.3f}$",
                    transform=ax.transAxes, ha="right", va="bottom", fontsize=14)
            ax.text(0.05, 0.95, PARAM_SYMBOLS.get(short, short),
                    transform=ax.transAxes, ha="left", va="top", fontsize=24)

            ax.set_xlabel(f"Generative {short}", fontsize=16)
            ax.set_ylabel(f"Recovered {short}", fontsize=16)
            ax.legend(fontsize=10, framealpha=0.5)
            sns.despine(ax=ax)

        fig.tight_layout()
        out = os.path.join(FIGURES_DIR, f"corr_gen_vs_rec_{model_type}.png")
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {os.path.basename(out)}")

        # ── Figure 2: pair plot (all parameters: nu, beta, omega[, omega3]) ──────
        # Build columns from short_names so omega (and omega3 for 3-level) are
        # included automatically. Ordering: all gen_ columns, then all rec_ columns.
        PAIR_LABEL_MAP = {
            "gen_nu":     r"Gen. $\nu$",
            "gen_beta":   r"Gen. $\beta$",
            "gen_omega":  r"Gen. $\omega$",
            "gen_omega3": r"Gen. $\omega_3$",
            "rec_nu":     r"Rec. $\nu$",
            "rec_beta":   r"Rec. $\beta$",
            "rec_omega":  r"Rec. $\omega$",
            "rec_omega3": r"Rec. $\omega_3$",
        }
        pair_cols = (
            [f"gen_{s}" for s in short_names] +
            [f"rec_{s}" for s in short_names]
        )
        # Keep only columns present in summary (guard against missing omega3 etc.)
        pair_cols = [c for c in pair_cols if c in summary.columns]
        pair_data = summary[pair_cols].dropna().copy()
        pair_data.columns = [PAIR_LABEL_MAP.get(c, c) for c in pair_cols]

        n_vars    = len(pair_cols)
        cell_size = 2.5
        g = sns.PairGrid(pair_data, diag_sharey=False, height=cell_size)
        g.map_diag(sns.histplot, bins=20, color="#4472C4", alpha=0.7)
        g.map_lower(plt.scatter, alpha=0.1, s=3, color="#4472C4", linewidths=0)
        g.map_upper(plt.scatter, alpha=0.1, s=3, color="#4472C4", linewidths=0)

        # Pearson r in each off-diagonal cell
        cols = pair_data.columns.tolist()
        for i in range(len(cols)):
            for j in range(len(cols)):
                if i == j:
                    continue
                ax = g.axes[i, j]
                r, _ = pearsonr(pair_data[cols[j]], pair_data[cols[i]])
                ax.text(0.05, 0.95, f"$r = {r:.2f}$",
                        transform=ax.transAxes, ha="left", va="top",
                        fontsize=10, color="black")

        for ax in g.axes.flat:
            ax.tick_params(labelsize=11)
            sns.despine(ax=ax)
        g.figure.set_size_inches(cell_size * n_vars, cell_size * n_vars)
        g.figure.tight_layout()
        out = os.path.join(FIGURES_DIR, f"pair_beta_nu_{model_type}.png")
        g.figure.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(g.figure)
        print(f"  Saved: {os.path.basename(out)}")

        # ── 3D interactive HTML plots (Plotly) ────────────────────────────────
        required_3d = ["gen_nu", "rec_nu", "gen_beta", "rec_beta"]
        missing_3d  = [c for c in required_3d if c not in summary.columns]
        if missing_3d:
            print(f"  [WARNING] Missing columns for 3D plots: {missing_3d} — skipping.")
        else:
            sub3d = summary[required_3d + ["record_id"]].dropna()

            def _make_3d_html(xs, ys, zs, xlabel, ylabel, zlabel, outname):
                fig_3d = go.Figure(data=[go.Scatter3d(
                    x=xs, y=ys, z=zs,
                    mode="markers",
                    marker=dict(
                        size=2,
                        color=ys,
                        colorscale="Viridis",
                        opacity=0.5,
                        colorbar=dict(title=ylabel, thickness=15),
                    ),
                )])
                fig_3d.update_layout(
                    scene=dict(
                        xaxis_title=xlabel,
                        yaxis_title=ylabel,
                        zaxis_title=zlabel,
                    ),
                    font=dict(family="Arial", size=13),
                    margin=dict(l=0, r=0, b=0, t=30),
                )
                path = os.path.join(FIGURES_DIR, outname)
                fig_3d.write_html(path)
                print(f"  Saved: {outname}")

            _make_3d_html(
                sub3d["gen_nu"], sub3d["rec_nu"], sub3d["rec_beta"],
                "Generative ν", "Recovered ν", "Recovered β",
                f"3d_gennu_recnu_recbeta_{model_type}.html",
            )
            _make_3d_html(
                sub3d["gen_beta"], sub3d["rec_beta"], sub3d["rec_nu"],
                "Generative β", "Recovered β", "Recovered ν",
                f"3d_genbeta_recbeta_recnu_{model_type}.html",
            )
            _make_3d_html(
                sub3d["gen_nu"], sub3d["rec_nu"], sub3d["gen_beta"],
                "Generative ν", "Recovered ν", "Generative β",
                f"3d_gennu_recnu_genbeta_{model_type}.html",
            )

    print(f"\nAll figures saved to: {FIGURES_DIR}")
