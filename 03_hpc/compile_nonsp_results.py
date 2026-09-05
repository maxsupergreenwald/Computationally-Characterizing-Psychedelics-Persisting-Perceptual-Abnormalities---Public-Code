"""
compile_nonsp_results.py
========================
Pull summary_dfs CSVs from Bouchet nonsp_predictor_analyses, compile into
one dataframe matching (and extending) the format of
  results/sensitivity_analyses_single_paths/existingresults_manuscript.csv

Directory structure on Bouchet:
  nonsp_predictor_analyses/
    {spvar}/                      ← e.g. vch_nu_normalized
      {model_type}/               ← e.g. nice_covariates, true_univariate_spusers
        results/
          summary_dfs/
            {dv}.csv              ← e.g. persist_vis_yn.csv, caps_vision.csv

Output
------
  results/sensitivity_analyses_single_paths/existingresults_manuscript.csv
  results/sensitivity_analyses_single_paths/existingresults_manuscript_counterfactual.csv

Columns:
  Estimate, Est.Error, l-95% CI, u-95% CI, Rhat, Bulk_ESS, Tail_ESS,
  var, model, covariates, N, num_divergents, prob_below_0, prob_above_0,
  dv, spvar, cov

  dv     = DV name (from CSV filename, e.g. persist_vis_yn)
  spvar  = focal predictor directory name (e.g. vch_nu_normalized)
  cov    = model type directory name (e.g. nice_covariates) — same as
           covariates column; kept for schema parity with existingresults

Strategy: SSH into Bouchet, tar only the summary_dfs CSVs into a single
archive, scp the archive, then walk the extracted tree locally.
With 8000+ CSVs, one tarball transfer is far faster than individual scp calls.

Usage
-----
  python compile_nonsp_results.py
  python compile_nonsp_results.py --predictors vch_nu_normalized vch_threshold_normalized
  python compile_nonsp_results.py --model-types nice_covariates true_univariate
  python compile_nonsp_results.py --dry-run

Requires an active SSH ControlMaster socket (~/.ssh/bouchet.sock).
Establish with:  ssh -MNf bouchet
"""

import argparse
import os
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────────
HPC_BASE   = "/nfs/roberts/scratch/pi_arp29/msg74/aim1_baseline_final/nonsp_predictor_analyses"
HPC_HOST   = "bouchet"
HPC_TMPTAR = "/tmp/nonsp_summary_dfs.tar.gz"

LOCAL_PROJECT = Path(__file__).resolve().parents[1]   # hppd_manuscript_public/
LOCAL_OUT_DIR  = LOCAL_PROJECT / "results" / "sensitivity_analyses_single_paths"

# Directories inside HPC_BASE that are NOT predictor dirs (skip these)
NON_PREDICTOR_DIRS = {"job_summaries"}

# DHARMa columns pulled from diagnostics CSVs and merged into compiled output.
DHARMA_COLS = [
    "dharma_uniformity_pval",
    "dharma_dispersion_pval",
    "dharma_outlier_pval",
    "dharma_quantiles_pval",
    "dharma_heteroscedasticity_pval",
    "dharma_heteroscedasticity_q25_pval",
    "dharma_heteroscedasticity_q50_pval",
    "dharma_heteroscedasticity_q75_pval",
    "dharma_quantiles_q25_pval",
    "dharma_quantiles_q50_pval",
    "dharma_quantiles_q75_pval",
    "ks_test_pval",
]

RHAT_CENTER = 1.0
RHAT_TOLERANCE = 0.1
RHAT_MIN = RHAT_CENTER - RHAT_TOLERANCE
RHAT_MAX = RHAT_CENTER + RHAT_TOLERANCE


def _filter_models_by_rhat(compiled: pd.DataFrame) -> pd.DataFrame:
    """Keep only models where all coefficient rows have Rhat within [0.9, 1.1]."""
    if "Rhat" not in compiled.columns:
        print("WARNING: Rhat column missing; skipping Rhat filtering.")
        return compiled

    group_cols = ["dv", "spvar", "cov"]
    rhat_numeric = pd.to_numeric(compiled["Rhat"], errors="coerce")

    group_stats = (
        compiled.assign(Rhat_numeric=rhat_numeric)
        .groupby(group_cols, dropna=False)["Rhat_numeric"]
        .agg(
            min_rhat="min",
            max_rhat="max",
            n_missing=lambda s: s.isna().sum(),
        )
        .reset_index()
    )

    valid_groups = group_stats[
        (group_stats["n_missing"] == 0)
        & (group_stats["min_rhat"] >= RHAT_MIN)
        & (group_stats["max_rhat"] <= RHAT_MAX)
    ][group_cols]

    filtered = compiled.merge(valid_groups, on=group_cols, how="inner")
    dropped_rows = len(compiled) - len(filtered)
    dropped_models = len(group_stats) - len(valid_groups)

    print(
        "\nRhat filter applied: "
        f"kept {len(filtered):,}/{len(compiled):,} rows; "
        f"dropped {dropped_rows:,} rows across {dropped_models:,} model groups "
        f"(required all Rhats in [{RHAT_MIN:.1f}, {RHAT_MAX:.1f}])."
    )

    return filtered


def _confirm_controlmaster():
    """Require the user to confirm DUO auth before any remote calls run."""
    print("\nThis script pulls files from Bouchet via SSH + scp.")
    print("It requires an active SSH ControlMaster socket (~/.ssh/bouchet.sock).\n")
    print("Have you run:  ssh -MNf bouchet")
    print("and approved the DUO push on your phone?\n")
    answer = input('Type "yes" to continue: ').strip().lower()
    if answer != "yes":
        print("Aborting. Run `ssh -MNf bouchet`, approve DUO, then retry.")
        raise SystemExit(1)
    print()


def _run(cmd: list, check: bool = True) -> subprocess.CompletedProcess:
    result = subprocess.run(cmd, capture_output=True, text=True)
    if check and result.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\n{result.stderr.strip()}"
        )
    return result


def _build_find_cmd(predictors: list | None, model_types: list | None) -> str:
    """Build the remote find command, optionally filtered by predictor / model type."""
    base = f"find {HPC_BASE}"

    # Restrict to specific predictor subdirs if requested
    if predictors:
        # find ... \( -path "*/vch_nu_normalized/*" -o -path "*/vch_threshold_normalized/*" \)
        path_clauses = " -o ".join(
            f'-path "{HPC_BASE}/{p}/*"' for p in predictors
        )
        base += f" \\( {path_clauses} \\)"

    # Match summary_dfs CSVs OR diagnostics CSVs in one find pass
    base += (
        ' \\( -name "*.csv" -path "*/summary_dfs/*"'
        ' -o -name "*_diagnostics.csv" -path "*/diagnostics/*" \\)'
    )

    # Post-filter by model type if requested (done via grep, simpler than find's path matching)
    if model_types:
        pattern = "|".join(f"/{mt}/" for mt in model_types)
        base += f" | grep -E '{pattern}'"

    return base


def pull_tarball(
    predictors: list | None,
    model_types: list | None,
    dry_run: bool,
) -> Path | None:
    """Create tarball on Bouchet containing only summary_dfs CSVs, scp it locally."""

    find_cmd = _build_find_cmd(predictors, model_types)

    # Count files first so we know what to expect
    count_result = _run(
        ["ssh", HPC_HOST, f"({find_cmd}) | wc -l"]
    )
    n_files = int(count_result.stdout.strip())
    print(f"Found {n_files:,} CSVs to pull (summary_dfs + diagnostics).")

    if n_files == 0:
        print("Nothing to pull — check predictor / model-type filters.")
        return None

    if dry_run:
        print(f"[dry-run] Would tar {n_files:,} CSVs from {HPC_BASE} → {HPC_TMPTAR}")
        print(f"[dry-run] Would scp {HPC_HOST}:{HPC_TMPTAR} → local tmp dir")
        return None

    # Create tarball on Bouchet
    print(f"Creating tarball on Bouchet ({n_files:,} files) ...")
    tar_cmd = (
        f"cd {HPC_BASE} && "
        f"({find_cmd} | xargs -I{{}} realpath --relative-to={HPC_BASE} {{}}) "
        f"| tar czf {HPC_TMPTAR} -C {HPC_BASE} -T -"
    )
    _run(["ssh", HPC_HOST, tar_cmd])
    print("  Tarball created.")

    # scp tarball to a local tempfile
    local_tar = Path(tempfile.mktemp(suffix=".tar.gz"))
    print(f"Downloading tarball → {local_tar} ...")
    _run(["scp", f"{HPC_HOST}:{HPC_TMPTAR}", str(local_tar)])
    print(f"  Downloaded ({local_tar.stat().st_size / 1e6:.1f} MB).")

    # Clean up remote tarball
    _run(["ssh", HPC_HOST, f"rm -f {HPC_TMPTAR}"])

    return local_tar


def compile_from_tarball(
    local_tar: Path,
    predictors: list | None,
    model_types: list | None,
    dvs: list | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Extract tarball and compile summary_dfs + diagnostics into two dataframes.

    summary_dfs CSVs provide coefficient-level rows.
    diagnostics CSVs provide DHARMa p-values which are merged in at the
    (spvar, cov, dv) level and broadcast across all coefficient rows.

    Returns (compiled, compiled_counterfactual).
    compiled_counterfactual contains rows from {dv}_counterfactual.csv files;
    its "dv" column holds the base DV name (suffix stripped) so diagnostics
    merge correctly.  Returns an empty DataFrame for compiled_counterfactual
    if no counterfactual CSVs are found.
    """
    frames = []
    cf_frames = []   # rows from {dv}_counterfactual.csv files
    diag_rows = []   # one entry per (spvar, model_type, dv)
    skipped = 0

    with tempfile.TemporaryDirectory() as tmpdir:
        print("Extracting tarball ...")
        with tarfile.open(local_tar, "r:gz") as tf:
            tf.extractall(tmpdir)

        root = Path(tmpdir)
        all_csvs = list(root.rglob("*.csv"))

        summary_csvs     = [p for p in all_csvs if p.parent.name == "summary_dfs"]
        diagnostic_csvs  = [p for p in all_csvs
                            if p.parent.name == "diagnostics"
                            and p.name.endswith("_diagnostics.csv")]

        print(f"  Extracted {len(summary_csvs):,} summary CSVs, "
              f"{len(diagnostic_csvs):,} diagnostics CSVs. Compiling ...")

        # ── Summary rows ─────────────────────────────────────────────────────
        for csv_path in sorted(summary_csvs):
            # Expected: {spvar}/{model_type}/results/summary_dfs/{dv}.csv
            # or {spvar}/{model_type}/results/summary_dfs/{dv}_counterfactual.csv
            parts = csv_path.relative_to(root).parts
            if len(parts) < 4:
                skipped += 1
                continue

            spvar      = parts[0]
            model_type = parts[1]
            dv_name    = csv_path.stem
            is_cf      = dv_name.endswith("_counterfactual")
            base_dv    = dv_name[: -len("_counterfactual")] if is_cf else dv_name

            if spvar in NON_PREDICTOR_DIRS:
                skipped += 1
                continue
            if predictors and spvar not in predictors:
                skipped += 1
                continue
            if model_types and model_type not in model_types:
                skipped += 1
                continue
            # Filter on base DV name so counterfactual files pass the same filter
            # as their parent DV (e.g. "lshs_total" passes for "lshs_total_counterfactual")
            if dvs and base_dv not in dvs:
                skipped += 1
                continue

            try:
                df = pd.read_csv(csv_path)
            except Exception as e:
                print(f"  WARNING: could not read {csv_path}: {e}")
                skipped += 1
                continue

            df["dv"]    = base_dv   # always store base name so diagnostics merge works
            df["spvar"] = spvar
            df["cov"]   = model_type
            if is_cf:
                cf_frames.append(df)
            else:
                frames.append(df)

        # ── Diagnostic rows ──────────────────────────────────────────────────
        for csv_path in sorted(diagnostic_csvs):
            # Expected: {spvar}/{model_type}/results/diagnostics/{dv}_diagnostics.csv
            parts = csv_path.relative_to(root).parts
            if len(parts) < 4:
                continue

            spvar      = parts[0]
            model_type = parts[1]
            dv_name    = csv_path.stem.replace("_diagnostics", "")

            if spvar in NON_PREDICTOR_DIRS:
                continue
            if predictors and spvar not in predictors:
                continue
            if model_types and model_type not in model_types:
                continue
            if dvs and dv_name not in dvs:
                continue

            try:
                d = pd.read_csv(csv_path)
                row = {"spvar": spvar, "cov": model_type, "dv": dv_name}
                for col in DHARMA_COLS:
                    row[col] = float(d[col].iloc[0]) if col in d.columns and len(d) > 0 else float("nan")
                diag_rows.append(row)
            except Exception as e:
                print(f"  WARNING: could not read diagnostics {csv_path}: {e}")

    if not frames:
        print("No data compiled — check filters and tarball contents.")
        return pd.DataFrame(), pd.DataFrame()

    compiled = pd.concat(frames, ignore_index=True)

    # Build diagnostic lookup once — keyed on base DV name, reused for counterfactual too
    diag_df = pd.DataFrame(diag_rows) if diag_rows else None

    # Merge DHARMa columns in, broadcast to all coefficient rows of the same (spvar, cov, dv)
    if diag_df is not None:
        compiled = compiled.merge(diag_df, on=["spvar", "cov", "dv"], how="left")
        n_with_diag = compiled[DHARMA_COLS[0]].notna().sum()
        print(f"  Merged DHARMa diagnostics for {len(diag_rows)} (spvar, model, dv) combos "
              f"({n_with_diag:,} rows with uniformity p-value).")
    else:
        print("  No diagnostics CSVs found in tarball — DHARMa columns will be NaN.")
        for col in DHARMA_COLS:
            compiled[col] = float("nan")

    # Canonical column order
    base_cols = [
        "Estimate", "Est.Error", "l-95% CI", "u-95% CI",
        "Rhat", "Bulk_ESS", "Tail_ESS",
        "var", "model", "covariates", "N", "num_divergents",
        "prob_below_0", "prob_above_0",
        "dv", "spvar", "cov",
    ] + DHARMA_COLS
    extra = [c for c in compiled.columns if c not in base_cols]
    compiled = compiled[[c for c in base_cols if c in compiled.columns] + extra]

    compiled = _filter_models_by_rhat(compiled)

    # ── Counterfactual compilation ─────────────────────────────────────────────
    # Mirrors the main compilation exactly; diagnostics are keyed on base DV name
    # so they merge correctly even though the CSV filename had _counterfactual.
    if cf_frames:
        compiled_cf = pd.concat(cf_frames, ignore_index=True)
        if diag_df is not None:
            compiled_cf = compiled_cf.merge(diag_df, on=["spvar", "cov", "dv"], how="left")
        else:
            for col in DHARMA_COLS:
                compiled_cf[col] = float("nan")
        cf_extra = [c for c in compiled_cf.columns if c not in base_cols]
        compiled_cf = compiled_cf[
            [c for c in base_cols if c in compiled_cf.columns] + cf_extra
        ]
        compiled_cf = _filter_models_by_rhat(compiled_cf)
        print(f"  Counterfactual: compiled {len(compiled_cf):,} rows "
              f"from {len(cf_frames)} CSVs.")
    else:
        compiled_cf = pd.DataFrame()

    if skipped:
        print(f"  ({skipped} paths skipped — non-predictor dirs or parse errors)")

    return compiled, compiled_cf



# ══════════════════════════════════════════════════════════════════════════════
# SINGLE SOURCE OF TRUTH: generate_hpc_jobs.py
# ══════════════════════════════════════════════════════════════════════════════
def load_generator_config():
    """Return generate_hpc_jobs.py's CONFIG namespace, without generating anything.

    The file is exec'd with HPC_JOBS_CONFIG_ONLY=1, which skips its data load and
    makes it stop at the END OF CONFIG DERIVATION marker — before the first line
    that writes anything. What comes back is the same MODEL_VARIANTS / ALL_PREDS
    / HPPD_CAPS_DVS / MEDIATION_ANALYSES / DVS the generator itself uses.

    Why this exists: compile_nonsp_results_config.py and
    compile_mediation_results_config.py used to hold hand-maintained copies of
    these lists. They drifted — on 2026-08-31 the mediation config named 3 of the
    9 submitted model types and 4 of the 6 mediators, so a compile silently
    pulled a fraction of what was on the cluster and reported success. A missing
    directory is a no-op for the pull, so nothing surfaced the gap. Both files
    are deleted; there is now exactly one place model types are configured.
    """
    gen = Path(__file__).resolve().parent / 'generate_hpc_jobs.py'
    if not gen.exists():
        sys.exit(f'ERROR: {gen} not found — the compile scripts read their '
                 'settings from it and cannot run without it.')
    ns = {'__name__': '__generator_config__', '__file__': str(gen)}
    env_before = os.environ.get('HPC_JOBS_CONFIG_ONLY')
    os.environ['HPC_JOBS_CONFIG_ONLY'] = '1'
    try:
        # The CONFIG region prints a variant summary; swallow it so it does not
        # interleave with this script's own output.
        import contextlib, io
        with contextlib.redirect_stdout(io.StringIO()):
            exec(compile(gen.read_text(), str(gen), 'exec'), ns)
    except SystemExit:
        pass                      # expected: the marker stops execution
    finally:
        if env_before is None:
            os.environ.pop('HPC_JOBS_CONFIG_ONLY', None)
        else:
            os.environ['HPC_JOBS_CONFIG_ONLY'] = env_before

    missing = [k for k in ('MODEL_VARIANTS', 'ALL_PREDS', 'HPPD_CAPS_DVS',
                           'MEDIATION_ANALYSES', 'DVS') if k not in ns]
    if missing:
        sys.exit(f'ERROR: generate_hpc_jobs.py did not define {missing} before its '
                 'END OF CONFIG DERIVATION marker. If a setting moved below the '
                 'marker, move it back above.')
    return ns


def main():
    # ── Subset mode: filters come from generate_hpc_jobs.py's CONFIG ─────────
    use_subset = len(sys.argv) > 1 and sys.argv[1].lower() == "subset"
    if use_subset:
        cfg = load_generator_config()
        # The compile script addresses spvar DIRECTORIES, which are named with
        # the _normalized suffix (see the note in 03_hpc/README.md), so the raw
        # predictor names from CONFIG are suffixed here rather than there.
        cfg_predictors  = [f"{p}_normalized" for p in cfg["ALL_PREDS"]]
        cfg_model_types = list(cfg["MODEL_VARIANTS"])
        cfg_dvs         = list(cfg["HPPD_CAPS_DVS"])
        if cfg.get("INCLUDE_SP_VCH_JOBS"):
            cfg_dvs += [d for d in cfg["VCH_AS_DVS"] if d not in cfg_dvs]
        print(f"Subset mode: {len(cfg_predictors)} predictors, "
              f"{len(cfg_model_types)} model types, "
              f"{len(cfg_dvs)} DVs — from generate_hpc_jobs.py")

        # Minimal arg parse for --dry-run only in subset mode
        dry_run = "--dry-run" in sys.argv
        if not dry_run:
            _confirm_controlmaster()

        local_tar = pull_tarball(
            predictors=cfg_predictors,
            model_types=cfg_model_types,
            dry_run=dry_run,
        )
        if dry_run or local_tar is None:
            return

        compiled, compiled_cf = compile_from_tarball(
            local_tar=local_tar,
            predictors=cfg_predictors,
            model_types=cfg_model_types,
            dvs=cfg_dvs,
        )

    else:
        # ── Standard argparse mode ────────────────────────────────────────────
        parser = argparse.ArgumentParser(description=__doc__,
                                         formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument(
            "--predictors", nargs="+", metavar="SPVAR",
            help="Pull only these predictor directories (e.g. vch_nu_normalized). "
                 "Default: all.",
        )
        parser.add_argument(
            "--model-types", nargs="+", metavar="MODEL_TYPE",
            help="Pull only these model-type directories (e.g. nice_covariates). "
                 "Default: all.",
        )
        parser.add_argument(
            "--dvs", nargs="+", metavar="DV",
            help="Keep only these DVs in the compiled output. Default: all.",
        )
        parser.add_argument(
            "--dry-run", action="store_true",
            help="Print what would be pulled without executing.",
        )
        parser.add_argument(
            "--no-confirm", action="store_true",
            help="Skip the SSH ControlMaster confirmation prompt. "
                 "Use when called from generate_hpc_jobs.py (user already confirmed).",
        )
        args = parser.parse_args()

        if not args.dry_run and not args.no_confirm:
            _confirm_controlmaster()

        local_tar = pull_tarball(
            predictors=args.predictors,
            model_types=args.model_types,
            dry_run=args.dry_run,
        )
        if args.dry_run or local_tar is None:
            return

        compiled, compiled_cf = compile_from_tarball(
            local_tar=local_tar,
            predictors=args.predictors,
            model_types=args.model_types,
            dvs=args.dvs,
        )

    local_tar.unlink(missing_ok=True)

    if compiled.empty:
        return

    LOCAL_OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = LOCAL_OUT_DIR / "existingresults_manuscript.csv"
    compiled.to_csv(out_path, index=False)

    print(f"\nCompiled {len(compiled):,} rows")
    print(f"  DVs:        {sorted(compiled['dv'].unique())}")
    print(f"  Predictors: {compiled['spvar'].nunique()} unique spvars")
    print(f"  Model types:{compiled['cov'].nunique()} unique model types")
    print(f"\nSaved → {out_path}")

    if not compiled_cf.empty:
        cf_out_path = LOCAL_OUT_DIR / "existingresults_manuscript_counterfactual.csv"
        compiled_cf.to_csv(cf_out_path, index=False)
        print(f"\nCounterfactual compiled {len(compiled_cf):,} rows")
        print(f"  DVs:        {sorted(compiled_cf['dv'].unique())}")
        print(f"  Predictors: {compiled_cf['spvar'].nunique()} unique spvars")
        print(f"  Model types:{compiled_cf['cov'].nunique()} unique model types")
        print(f"\nSaved → {cf_out_path}")


if __name__ == "__main__":
    main()
