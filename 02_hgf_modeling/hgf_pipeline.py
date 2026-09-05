"""
hgf_pipeline.py  ─  HGF Analysis Pipeline (VCH task, hppd_manuscript project)
==============================================================================

PURPOSE
-------
Five-step pipeline for running Hierarchical Gaussian Filter (HGF) analyses on
this manuscript's VCH (visual contrast-detection) task data.

  Step 1  Parse compressed JSON task data from the project-level wide df
          → per-timepoint trial-level long DataFrame
  Step 2  Export per-participant CSVs + master long DataFrame to local
          staging directories ready for transfer to HPC
  Step 3  Generate SLURM dSQ job-array .txt files; print exact terminal
          commands for scp transfer, SSH, and dSQ submission
  Step 4  Print exact scp commands to pull HPC results back locally
  Step 5  Parse HPC outputs → summary metrics → merge into the wide df.
          Handled by import_hgf_results_unified.py, NOT by this script; run it
          after pulling results (see the printed Step 4 commands).

HOW TO USE (quick start)
------------------------
1. Set the variables in ACTIVE CONFIG (below) for your current run.
2. Run:
       python3 hgf_pipeline.py
   from any directory (the script uses absolute paths).
3. Follow the printed terminal commands to transfer, submit, and pull.

NEW-USER BACKGROUND
--------------------
Input data
  The pipeline reads a wide-format CSV ("wide df") where every row is one
  participant and the VCH task is stored as a single BASE64-encoded, gzip-
  compressed JSON blob in one column.  The JSON contains trial-by-trial
  responses and stimulus values.

VCH (visual contrast-detection)
  JSON key "component_1" … "component_4" → four blocks per participant.
  Stimulus intensity is ranked to derive a nominal condition (0/25/50/75 %)
  which is then mapped to an empirically-derived hit-rate (empirical_condition)
  from a COPE normative sample.  The HGF uses this empirical hit-rate as the
  trial-by-trial stimulus strength input.
  Julia script: model_fitting_singleagent_forarray.jl  (binary 2-level or
  3-level HGF; 4 MCMC chains; MCMCSerial; MAP init)

  MODEL_TYPE selects the HGF variant.  Empiric models use empirically-derived
  hit-rates as stimulus strength; nominal models use the raw condition
  proportions {0.0, 0.25, 0.5, 0.75}.  Results land in results/{MODEL_TYPE}/.
  "2level_nominal" / "3level_nominal" are sensitivity checks (Supplementary
  Table S5).

HPC cluster
  Yale Bouchet cluster (SLURM).  Transfer via scp to the transfer node;
  submit via SSH to the login node.  dSQ bundles lines from a .txt file
  into a SLURM array.

Required local repo structure
  02_hgf_modeling/       ← LOCAL_JULIA_CH_DIR
    model_fitting_singleagent_forarray.jl
    Project.toml, Manifest.toml
    data/                ← master long CSVs (gitignored)
    data_n_cmnds/        ← per-subject CSVs + job arrays (gitignored)
      array_scripts/
      {project_hpc_key}/{timepoint}/vch_data/{record_id}.csv
"""

################################################################################
# ─── ROLE IN THE MANUSCRIPT ───────────────────────────────────────────────────
#
# STAGE 1 OF 3 in the HGF pipeline.  Nothing in this directory produces a
# parameter estimate until this script has run.
#
#     THIS SCRIPT ──▶ model_fitting_singleagent_forarray.jl ──▶
#                     import_hgf_results_unified.py
#
# What it is responsible for
#     Turning the compressed JSON blob stored in each participant's REDCap row
#     into one tidy trial-level CSV per participant, laid out exactly where the
#     Julia fitting script expects to find it, plus the dSQ job file that fans
#     those participants out across a SLURM array.
#
#     It also fixes the stimulus intensity used by the HGF, via
#     VCH_STAGING_MAPPINGS.  Those values are duplicated (deliberately) in the
#     Julia fitting script, which recomputes them at runtime as a guard — see §4
#     of model_fitting_singleagent_forarray.jl.  If you change the mapping in one
#     place you MUST change it in both, and re-stage.
#
# Inputs and outputs
#     Reads   the project wide df, data/final/df_bl_<date>.csv
#     Writes  data_n_cmnds/<project>/<timepoint>/vch_data/<record_id>.csv
#             data_n_cmnds/array_scripts/<jobfile>.txt
#             data/vch_master_<timepoint>_<date>.csv
#     Prints  the exact scp / dsq / sbatch commands to run next.  Use the printed
#             commands rather than reconstructing paths by hand.
################################################################################

import os
import json
import gzip
import base64
import pandas as pd
import numpy as np
from datetime import datetime
from pathlib import Path


################################################################################
# ─── ACTIVE CONFIG — edit these on each run ───────────────────────────────────
################################################################################

# Which project.  Only "hppd_manuscript" is configured in this repository.
ACTIVE_PROJECT = "hppd_manuscript"

# Which tasks to process.  Only "vch" is configured in this repository.
ACTIVE_TASKS = ["vch"]

# "skip"  → skip subjects whose output file already exists on HPC
# anything else → overwrite / rerun all subjects
SKIP_EXISTING = "skip"

# True  → only process the first 3 record_id values (dry-run / sanity check)
# False → process all subjects
TEST_MODE = False

# HGF model variant.  One of:
#   "2level_empiric"             — 2-level HGF, empirical stimulus strengths (primary analysis)
#   "2level_nominal"             — 2-level HGF, raw condition proportions (sensitivity check)
#   "3level_empiric"             — 3-level HGF, empirical stimulus strengths
#   "3level_nominal"             — 3-level HGF, raw condition proportions (sensitivity check)
MODEL_TYPE = "2level_empiric"

# Optional: supply an explicit wide-df path, e.g. to target a specific dated file.
# Set to None to auto-select the most-recently-modified file with the correct prefix.
WIDE_DF_PATH = None


################################################################################
# ─── LOCAL PATHS ──────────────────────────────────────────────────────────────
################################################################################

# Resolved from this file's own location so the repository can be cloned or
# moved anywhere without editing paths.
#   LOCAL_JULIA_CH_DIR  — this directory (02_hgf_modeling).  All Julia scripts,
#                         staging directories, and generated results live here.
#   PUBLIC_REPO_DIR     — the repository root, one level up.  Holds data/final/.
LOCAL_JULIA_CH_DIR  = os.path.dirname(os.path.abspath(__file__))
PUBLIC_REPO_DIR     = os.path.dirname(LOCAL_JULIA_CH_DIR)


################################################################################
# ─── HPC PATHS — change only if the cluster layout changes ────────────────────
################################################################################

HPC_BASE          = "/nfs/roberts/scratch/pi_arp29/msg74"
HPC_JULIA_CH_DIR  = f"{HPC_BASE}/julia_hgf_ch"
HPC_TRANSFER_HOST = "transfer-bouchet.ycrc.yale.edu"
HPC_LOGIN_HOST    = "bouchet.ycrc.yale.edu"
HPC_USER          = "msg74"
JULIA_MODULE_VCH  = "Julia/1.11.4-linux-x86_64"   # for model_fitting_singleagent_forarray.jl
# ActionModels 0.7 requires Julia 1.11+.

# SLURM time limits per model type (HH:MM format for dSQ -t flag)
VCH_TIMELIMITS = {
    "2level_empiric"            : "15:00",
    "2level_nominal"            : "15:00",
    "3level_empiric"            : "35:00",
    "3level_nominal"            : "35:00",
}


################################################################################
# ─── PROJECT CONFIG ────────────────────────────────────────────────────────────
# Key fields:
#
#   proj_dir          Local project root; wide df lives at proj_dir/data/final/
#   proj_suffix       Wide-df filename prefix component: df_{proj_suffix}_{date}.csv
#   project_hpc_key   String used in HPC staging + results paths for project isolation
#                       e.g.  results/{model_type}/{project_hpc_key}/{timepoint}/vch/
#
#   timepoints_vch    Ordered list of timepoints for the VCH HGF.
#
#   task_cols         Maps "vch" → list of wide-df column names, one entry per
#                     timepoint (same order as timepoints_vch).
#
#   long_df_dir       Where to save combined long DataFrames locally.
#   julia_data_dir    Where to save dated master CSVs.
#   array_scripts_dir Where to save dSQ job-array .txt files.
#
#   results_local     Maps task key → local directory for pulled HPC results.
#
# ─── HPC staging layout ────────────────────────────────────────────────────────
#
#   Per-subject CSVs  → 02_hgf_modeling/data_n_cmnds/{proj_key}/{tp}/vch_data/{rid}.csv
#   Job arrays        → 02_hgf_modeling/data_n_cmnds/array_scripts/{proj_key}_{tp}_vch_hgf_{model_type}.txt
#
# ─── HPC results layout ────────────────────────────────────────────────────────
#
#   → 02_hgf_modeling/results/vch/{model_type}/{tp}/vch/{rid}.*  (after pull)
################################################################################

PROJECT_CONFIGS = {

    # ── Cross-sectional: HPPD Manuscript ──────────────────────────────────────
    # Single imaging/drug session; the timepoint label "hppd" is used throughout
    # so that Julia results land in a consistently named directory.
    "hppd_manuscript": {
        "proj_dir"        : PUBLIC_REPO_DIR,

        # Wide-df filename stem: the pipeline auto-selects the most recent
        # data/final/df_{proj_suffix}_*.csv.  The public, PII-free wide df is
        # df_public_<date>.csv, so the suffix is "public" (it was "bl" when this
        # pipeline read the private df_bl_<date>.csv).
        "proj_suffix"     : "public",
        "project_hpc_key" : "hppd_manuscript",

        "timepoints_vch"  : ["hppd"],

        "task_cols": {
            "vch": ["task_data_vch_short_psychedelic_bl"],
        },

        "long_df_dir"      : os.path.join(PUBLIC_REPO_DIR, "data", "final"),
        "julia_data_dir"   : os.path.join(LOCAL_JULIA_CH_DIR, "data"),
        "array_scripts_dir": os.path.join(LOCAL_JULIA_CH_DIR, "data_n_cmnds", "array_scripts"),

        # Pulled HPC results land inside this module, so the HGF stage is
        # self-contained.  This directory is gitignored.
        "results_local": {
            "vch": os.path.join(LOCAL_JULIA_CH_DIR, "results", "vch"),
        },
    },
}


################################################################################
# ─── EMPIRICAL CONDITION MAPPINGS ─────────────────────────────────────────────
# Maps nominal condition proportion (0.0 – 0.75) to the empirically-derived
# hit rate used as the HGF stimulus-strength input.
# Keys are proportions (percent / 100).
#
# ── 2026-06-17 update ─────────────────────────────────────────────────────────
# Values for 25 %, 50 %, and 75 % updated to reflect detection rates from the
# NON-HALLUCINATOR subset of the out-of-set (COPE) normative cohort: subjects
# with max_vh_freq == 0 AND max_ah_freq == 0 (NaN treated as 0 per
# out_of_set_data/README), n=29 VCH-passing. Rationale: the HGF likelihood
# should be grounded in perceptual behavior of people without hallucination
# histories.
# Analysis and justification:
#   out_of_set_data/README_empirical_condition_redetermination.md
#   out_of_set_data/figures/empirical_condition_by_condition.png
#   out_of_set_data/empirical_condition_full_vs_nonhall.csv
#
# IMPORTANT — 0 % condition is ALWAYS 0.0, never the empirical false-alarm rate.
# At 0 % contrast there is literally no stimulus; the ground-truth detection
# probability is exactly 0.0.  Using the empirical FA rate here would
# misrepresent the task structure to the HGF.
#
# Previous values (QC-filtered full sample, n=114; used through 2026-06-16):
#   0.00: 0.0   (ground-truth: no stimulus → 0.0)
#   0.25: 0.44184090362893536
#   0.50: 0.7010130961205832
#   0.75: 0.8690302144249513
################################################################################

VCH_CONDITION_TO_EMPIRICAL = {
    0.75: 0.8994252873563219,
    0.50: 0.7115104419621175,
    0.25: 0.4180444024563061,
    0.00: 0.0,  # ground-truth: no stimulus at 0 % contrast → 0.0 (NOT the empirical FA rate)
}

# ── Staging mapping lookup ─────────────────────────────────────────────────────
# Maps MODEL_TYPE → the condition→empirical_condition dict written into staging
# CSVs by parse_vch_data().  Scripts that read empirical_condition directly from
# staged CSVs (bms_vch.jl, param_recovery_vch.jl, ppc_classic_vch.jl) rely on this.
# model_fitting_singleagent_forarray.jl recomputes empirical_condition at runtime
# from condition using its own VCH_STIM_MAPPINGS dict (cross-contamination guard),
# so the staging value is overridden there.
# Nominal model types default to the corrected empirical mapping in staging
# CSVs; Julia ignores/overrides that column for those types anyway.
VCH_STAGING_MAPPINGS = {
    "2level_empiric":             VCH_CONDITION_TO_EMPIRICAL,
    "3level_empiric":             VCH_CONDITION_TO_EMPIRICAL,
    # Nominal: Julia reads condition directly; empiric fallback here.
    "2level_nominal":             VCH_CONDITION_TO_EMPIRICAL,
    "3level_nominal":             VCH_CONDITION_TO_EMPIRICAL,
}


################################################################################
# ─── HELPER FUNCTIONS ─────────────────────────────────────────────────────────
################################################################################

def load_most_recent_csv(directory, prefix):
    """Return the path of the most-recently-modified CSV in *directory* whose
    filename starts with *prefix* (ignores files containing 'LABELS')."""
    files = [
        os.path.join(directory, f)
        for f in os.listdir(directory)
        if f.endswith(".csv") and f.startswith(prefix) and "LABELS" not in f
    ]
    if not files:
        raise FileNotFoundError(
            f"No CSV with prefix '{prefix}' found in {directory}"
        )
    files.sort(key=os.path.getmtime, reverse=True)
    print(f"  Auto-selected: {files[0]}")
    return files[0]


def load_public_wide_df(project=None, df_path=None, require=()):
    """Load the public, PII-free wide dataframe for *project*.

    This is the ONLY supported way for scripts in this directory to obtain the
    participant-level dataframe.

    Why this exists rather than load_and_prepare_data()
        In the private analysis repository the wide df is built by calling
        load_and_prepare_data(), which derives dozens of columns that the raw
        REDCap export does not contain.  The public file distributed here is the
        OUTPUT of that function, with identifiable fields removed — so those
        derived columns are already present and load_and_prepare_data() must NOT
        be run on it again.  Doing so would attempt to re-derive columns from
        source fields that are no longer in the file.

    Parameters
    ----------
    project : str or None
        Key into PROJECT_CONFIGS.  Defaults to ACTIVE_PROJECT.
    df_path : str or None
        Explicit path to a wide df, bypassing auto-selection.  Use this to pin a
        specific dated file when file identity matters.
    require : iterable of str
        Column names that must be present.  Any that are missing raise KeyError
        naming them explicitly, rather than failing later with an opaque error.

    Returns
    -------
    (DataFrame, str)
        The dataframe, and the path it was loaded from.
    """
    cfg_ = PROJECT_CONFIGS[project or ACTIVE_PROJECT]
    if df_path is None:
        df_path = load_most_recent_csv(
            cfg_["long_df_dir"], f"df_{cfg_['proj_suffix']}_"
        )
    df_ = pd.read_csv(df_path, low_memory=False)
    df_["record_id"] = df_["record_id"].astype(int)

    missing = [c for c in require if c not in df_.columns]
    if missing:
        raise KeyError(
            f"Wide df {os.path.basename(df_path)} is missing required "
            f"column(s): {missing}.\n"
            "These are HGF parameter columns produced by "
            "import_hgf_results_unified.py.  Run that import for the "
            "corresponding model type first, or point df_path at a wide df "
            "that already contains them."
        )
    return df_, df_path


def _decompress_json(cell):
    """Decompress a base64 + gzip + JSON cell and return the parsed Python object."""
    decoded      = base64.b64decode(cell)
    decompressed = gzip.decompress(decoded).decode("utf-8")
    return json.loads(decompressed)


def parse_vch_data(df, task_col, record_id_col="record_id", condition_map=None):
    """Parse a VCH JSON column from the wide df.

    Parameters
    ----------
    condition_map : dict or None
        Maps nominal condition proportion (0.0–0.75) → empirical_condition value
        written to the staging CSV.  If None, falls back to VCH_CONDITION_TO_EMPIRICAL
        (the current corrected non-hallucinator mapping).  Pass a value from
        VCH_STAGING_MAPPINGS to match MODEL_TYPE.

    Returns a long DataFrame with columns:
      record_id, subj_id, trial, component, response, rt, confidence,
      intensity, rank, condition (proportion 0–0.75), empirical_condition, modality
    """
    if condition_map is None:
        condition_map = VCH_CONDITION_TO_EMPIRICAL
    participant_dfs = []
    blocks = ["component_1", "component_2", "component_3", "component_4"]
    threshold_mapping = {0: 0, 1: 25, 2: 50, 3: 75}

    for _, row in df.iterrows():
        if not isinstance(row[task_col], str):
            continue
        try:
            data = _decompress_json(row[task_col])
            rid  = row[record_id_col]

            block_dfs = []
            for b_num, block in enumerate(blocks, start=1):
                if block not in data:
                    continue
                df1 = pd.DataFrame({
                    "response" : data[block]["response"],
                    "rt"       : data[block].get("responseTime", np.nan),
                    "confidence": data[block].get("ratings", np.nan),
                    "intensity": data[block]["contrasts"],
                    "component": b_num,
                    "record_id": rid,
                })
                # Rank contrasts within block (dense, 0-based) → condition %
                df1["rank"]      = df1["intensity"].rank(method="dense").astype(int) - 1
                df1["condition"] = df1["rank"].map(threshold_mapping)
                block_dfs.append(df1)

            participant_df = pd.concat(block_dfs, ignore_index=True)
            participant_dfs.append(participant_df)

        except Exception as e:
            print(f"  [VCH] Skipping record {row.get(record_id_col, '?')}: {e}")

    if not participant_dfs:
        return pd.DataFrame()

    df_out = pd.concat(participant_dfs, ignore_index=True)
    df_out["trial"]    = df_out.groupby("record_id").cumcount() + 1
    df_out["modality"] = "vch"
    df_out["record_id"]  = df_out["record_id"].astype(int)
    df_out["subj_id"]    = df_out["record_id"]
    df_out["condition"]  = df_out["condition"].astype(float)
    if df_out["condition"].max() > 1:
        df_out["condition"] = df_out["condition"] / 100
    df_out["empirical_condition"] = df_out["condition"].map(condition_map)
    return df_out


################################################################################
# ─── STEP 1: Parse task data ──────────────────────────────────────────────────
################################################################################

def step1_parse_data(cfg, df, tasks, test_mode=False):
    """Parse all requested task columns from the wide df.

    Parameters
    ----------
    cfg       : project config dict
    df        : wide DataFrame already loaded from disk
    tasks     : list of task keys to process, e.g. ["vch"]
    test_mode : if True, limit each task/timepoint to the first 3 record_ids

    Returns
    -------
    data_by_task_tp : dict  {task: {timepoint: long_df}}
    """
    print("\n" + "="*70)
    print("STEP 1: Parsing task data from wide df")
    print("="*70)

    data_by_task_tp = {}

    for task in tasks:
        if task not in cfg["task_cols"]:
            print(f"  [WARNING] Task '{task}' not configured for this project. Skipping.")
            continue

        task_cols  = cfg["task_cols"][task]
        timepoints = cfg["timepoints_vch"]

        if len(task_cols) != len(timepoints):
            raise ValueError(
                f"Mismatch: {len(task_cols)} task_cols but {len(timepoints)} "
                f"timepoints for task '{task}' in project '{ACTIVE_PROJECT}'"
            )

        data_by_task_tp[task] = {}

        for tp, col in zip(timepoints, task_cols):
            if col not in df.columns:
                print(f"  [WARNING] Column '{col}' not found in wide df. Skipping {task}/{tp}.")
                continue

            print(f"\n  Parsing {task.upper()} / {tp}  (column: {col}) ...")

            if task == "vch":
                vch_map = VCH_STAGING_MAPPINGS.get(MODEL_TYPE, VCH_CONDITION_TO_EMPIRICAL)
                long_df = parse_vch_data(df, col, condition_map=vch_map)
            else:
                print(f"  [WARNING] Unknown task '{task}'. Skipping.")
                continue

            if long_df.empty:
                print(f"  [WARNING] No data parsed for {task}/{tp}.")
                continue

            if test_mode:
                test_records = sorted(long_df["record_id"].dropna().unique())[:3]
                long_df = long_df[long_df["record_id"].isin(test_records)].copy()
                print(f"  [TEST MODE] Limited to 3 records: {test_records}")

            long_df["timepoint"] = tp
            data_by_task_tp[task][tp] = long_df

            n_subj = long_df["record_id"].nunique()
            print(f"  → {n_subj} participants, {len(long_df)} rows")

    return data_by_task_tp


################################################################################
# ─── STEP 2: Export per-participant CSVs ──────────────────────────────────────
################################################################################

def step2_export_csvs(cfg, data_by_task_tp, today):
    """Export per-participant staging CSVs and save master long DataFrames.

    CSVs → {LOCAL_JULIA_CH_DIR}/data_n_cmnds/{proj_key}/{tp}/vch_data/{rid}.csv

    Master long CSV → {julia_data_dir}/vch_master_{proj_key}_{tp}_{today}.csv
    Combined CH master → {long_df_dir}/ch_master_{proj_key}_{today}.csv

    Returns
    -------
    exported : dict  {task: {timepoint: sorted_list_of_record_ids}}
    """
    print("\n" + "="*70)
    print("STEP 2: Exporting per-participant staging CSVs")
    print("="*70)

    proj_key = cfg["project_hpc_key"]
    exported = {}
    ch_long_dfs = []  # collect vch dfs for combined ch_master

    # Julia needs these columns in the per-subject VCH CSVs
    VCH_JULIA_COLS = ["subj_id", "condition", "empirical_condition", "response"]

    for task, tp_data in data_by_task_tp.items():
        exported[task] = {}

        for tp, long_df in tp_data.items():

            # ── Stage per-participant VCH CSVs ────────────────────────────────
            stage_dir = os.path.join(
                LOCAL_JULIA_CH_DIR, "data_n_cmnds", proj_key, tp, f"{task}_data"
            )
            os.makedirs(stage_dir, exist_ok=True)

            # Save master long df to julia_data_dir for reference
            os.makedirs(cfg["julia_data_dir"], exist_ok=True)
            master_path = os.path.join(
                cfg["julia_data_dir"],
                f"{task}_master_{proj_key}_{tp}_{today}.csv"
            )
            long_df.to_csv(master_path, index=False)

            # Write per-participant CSVs
            avail_cols = [c for c in VCH_JULIA_COLS if c in long_df.columns]
            records = sorted(long_df["record_id"].unique())
            for rid in records:
                sub_df = long_df[long_df["record_id"] == rid][avail_cols]
                sub_df.to_csv(os.path.join(stage_dir, f"{rid}.csv"), index=False)

            # ── Remove stale CSVs not in current subject set ──────────────────
            # Julia sorts staging files lexicographically; any extra file
            # shifts indices and causes subjects (especially 9xx IDs) to
            # fall beyond the job-array range.
            valid_names = {f"{rid}.csv" for rid in records}
            stale = [f for f in os.listdir(stage_dir)
                     if f.endswith(".csv") and f not in valid_names]
            for sf in stale:
                os.remove(os.path.join(stage_dir, sf))
            if stale:
                print(f"  [CLEANUP] Removed {len(stale)} stale CSV(s): "
                      f"{sorted(stale)}")

            exported[task][tp] = records
            ch_long_dfs.append(long_df)
            print(f"  Staged {len(records)} {task.upper()} CSVs → {stage_dir}")
            print(f"  Master CSV → {master_path}")

    # Save combined VCH ch_master
    if ch_long_dfs:
        os.makedirs(cfg["long_df_dir"], exist_ok=True)
        ch_master = pd.concat(ch_long_dfs, ignore_index=True)
        ch_master_path = os.path.join(
            cfg["long_df_dir"], f"ch_master_{proj_key}_{today}.csv"
        )
        ch_master.to_csv(ch_master_path, index=False)
        print(f"\n  Combined CH master → {ch_master_path}")

    return exported


################################################################################
# ─── STEP 3: Generate job arrays + print transfer/submission commands ──────────
################################################################################

def _make_vch_job_line(hpc_dir, file_indices, modality, skip_or_nah,
                        timepoint, proj_key, model_type):
    """Build one dSQ job line (two Julia calls bundled together).

    file_indices is a list of 1 or 2 1-based file indices.
    """
    def julia_cmd(idx):
        return (
            f"export GKSwstype=nul; "
            f"julia --project=. --threads 4 model_fitting_singleagent_forarray.jl "
            f"{idx} {modality} {skip_or_nah} {timepoint} {proj_key} {model_type}"
        )
    cmds = "; ".join(julia_cmd(i) for i in file_indices)
    return f"module load {JULIA_MODULE_VCH}; cd {hpc_dir}; {cmds}"


def step3_generate_job_arrays(cfg, exported, skip_or_nah, model_type):
    """Generate dSQ job-array .txt files and print all transfer + submission commands.

    Each .txt file has one SLURM-task line per two subjects (bundled for efficiency).
    Printed commands are ready to copy-paste into your terminal in order.
    """
    print("\n" + "="*70)
    print("STEP 3: Generating job arrays + transfer commands")
    print("="*70)

    proj_key       = cfg["project_hpc_key"]
    array_dir      = cfg["array_scripts_dir"]
    timelimit_vch  = VCH_TIMELIMITS.get(model_type, "120:00")
    os.makedirs(array_dir, exist_ok=True)

    generated_vch_arrays = []  # (local_path, filename, timepoint, task)

    for task, tp_records in exported.items():
        for tp, records in tp_records.items():
            if not records:
                print(f"  [WARNING] No records for {task}/{tp}. Skipping job array.")
                continue

            n = len(records)

            # ── VCH job array ──────────────────────────────────────────────────
            fname    = f"{proj_key}_{tp}_run_{task}_hgf_{model_type}.txt"
            fpath    = os.path.join(array_dir, fname)
            dsq_cmd  = (
                f"module load dSQ; dsq --job-file {fname} "
                f"--mem-per-cpu 4g -t {timelimit_vch} --mail-type ALL"
            )
            with open(fpath, "w") as f:
                f.write(f"# to run: {dsq_cmd}\n")
                for i in range(1, n + 1, 2):
                    indices = [i, i + 1] if i + 1 <= n else [i]
                    f.write(_make_vch_job_line(
                        HPC_JULIA_CH_DIR, indices, task, skip_or_nah,
                        tp, proj_key, model_type
                    ) + "\n")

            generated_vch_arrays.append((fpath, fname, tp, task))
            print(f"  Wrote {fname}  ({n} subjects, {(n + 1) // 2} SLURM tasks)")

    # ── Print transfer commands ────────────────────────────────────────────────
    print("\n" + "-"*70)
    print("TRANSFER COMMANDS — copy-paste these into your terminal in order:")
    print("-"*70)

    if generated_vch_arrays:
        print("\n# ── Transfer VCH staging data ────────────────────────────────────")
        print(f"scp -r {LOCAL_JULIA_CH_DIR}/data_n_cmnds/{proj_key} "
              f"{HPC_USER}@{HPC_TRANSFER_HOST}:{HPC_JULIA_CH_DIR}/data_n_cmnds/")
        print(f"\n# ── Transfer VCH job-array scripts ───────────────────────────────")
        for _, fname, tp, task in generated_vch_arrays:
            print(f"scp {LOCAL_JULIA_CH_DIR}/data_n_cmnds/array_scripts/{fname} "
                  f"{HPC_USER}@{HPC_TRANSFER_HOST}:{HPC_JULIA_CH_DIR}/data_n_cmnds/array_scripts/")

    # ── Print submission commands ──────────────────────────────────────────────
    print("\n" + "-"*70)
    print("SUBMISSION COMMANDS — SSH into cluster, then run dSQ:")
    print("-"*70)
    print(f"\nssh {HPC_USER}@{HPC_LOGIN_HOST}")

    if generated_vch_arrays:
        print(f"\n# ── VCH submissions ──────────────────────────────────────────────")
        print(f"cd {HPC_JULIA_CH_DIR}/data_n_cmnds/array_scripts")
        for _, fname, tp, task in generated_vch_arrays:
            dsq = (f"module load dSQ; dsq --job-file {fname} "
                   f"--mem-per-cpu 4g -t {timelimit_vch} --mail-type ALL")
            print(f"# {task.upper()} / {tp}")
            print(dsq)


################################################################################
# ─── STEP 4: Print pull commands ──────────────────────────────────────────────
################################################################################

def step4_print_pull_commands(cfg, exported, model_type):
    """Print rsync commands to pull HPC results back to local directories.

    Results on HPC: {HPC_JULIA_CH_DIR}/results/{model_type}/{proj_key}/{tp}/vch/
    """
    print("\n" + "="*70)
    print("STEP 4: Pull commands — run AFTER jobs complete on HPC")
    print("="*70)

    proj_key = cfg["project_hpc_key"]

    for task, tp_records in exported.items():
        for tp in tp_records:
            hpc_src   = f"{HPC_JULIA_CH_DIR}/results/{model_type}/{proj_key}/{tp}/{task}/"
            local_dst = os.path.join(
                cfg["results_local"][task], model_type, tp, task
            )
            os.makedirs(local_dst, exist_ok=True)
            print(f"\n# {task.upper()} / {tp}")
            print(f"rsync -avz {HPC_USER}@{HPC_TRANSFER_HOST}:{hpc_src} {local_dst}/")


################################################################################
# ─── MAIN ─────────────────────────────────────────────────────────────────────
################################################################################

if __name__ == "__main__":
    today = datetime.today().strftime("%m-%d-%Y")

    cfg = PROJECT_CONFIGS[ACTIVE_PROJECT]

    print("=" * 70)
    print(f"HGF PIPELINE  —  project: {ACTIVE_PROJECT}")
    print(f"  tasks      : {ACTIVE_TASKS}")
    print(f"  model_type : {MODEL_TYPE}")
    print(f"  skip_existing : {SKIP_EXISTING}")
    print(f"  test_mode  : {TEST_MODE}")
    print(f"  date stamp : {today}")
    print("=" * 70)

    # ── Load wide df ──────────────────────────────────────────────────────────
    print("\nLoading wide df ...")
    if WIDE_DF_PATH:
        df_path = WIDE_DF_PATH
        print(f"  Using specified path: {df_path}")
    else:
        df_path = load_most_recent_csv(
            cfg["long_df_dir"], f"df_{cfg['proj_suffix']}_"
        )
    df = pd.read_csv(df_path, low_memory=False)
    print(f"  Loaded {len(df)} rows × {len(df.columns)} columns")

    # ── Step 1 ────────────────────────────────────────────────────────────────
    data_by_task_tp = step1_parse_data(cfg, df, ACTIVE_TASKS, test_mode=TEST_MODE)

    # ── Step 2 ────────────────────────────────────────────────────────────────
    exported = step2_export_csvs(cfg, data_by_task_tp, today)

    # ── Step 3 ────────────────────────────────────────────────────────────────
    step3_generate_job_arrays(cfg, exported, SKIP_EXISTING, MODEL_TYPE)

    # ── Step 4 ────────────────────────────────────────────────────────────────
    step4_print_pull_commands(cfg, exported, MODEL_TYPE)

    print("\n" + "="*70)
    print("Pipeline complete.")
    print("="*70)
