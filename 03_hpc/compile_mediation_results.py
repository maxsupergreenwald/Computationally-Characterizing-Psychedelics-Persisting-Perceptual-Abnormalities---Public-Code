"""
compile_mediation_results.py
============================
Pull mediation result CSVs from Bouchet and extract them directly into
  results/{dv}/mediation_models/{model_name}/
at the repo root — matching the paths that 0X_all_figures.py reads.

Directory structure on Bouchet:
  hppd_manuscript_mediation_analyses/mediation_analyses/results/
    {dv_col}/
      mediation_models/
        {model_name}/
          summary_{equation}_{model_name}.csv
          mediation_results_{predictor}_{model_name}.csv
          mc_mediation_summary.csv
          hu_paths_summary.csv
          convergence_diagnostics.csv
          dharma_comprehensive_*.png
          dharma_diagnostics_*.png
          pp_check_{dv}_{model_name}.png
          pp_check_{mediator}_{model_name}.png
          ...

CSVs, DHARMa diagnostic PNGs (dharma_*.png), posterior predictive check PNGs
(pp_check_*.png), and .RData fit objects are pulled.
Tarball is created relative to the remote results root so extracted paths
land correctly under local results/.

Usage
-----
  python compile_mediation_results.py           # prompted for SSH + pulls all configured models
  python compile_mediation_results.py subset    # read model list from generate_hpc_jobs.py
  python compile_mediation_results.py --dry-run

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

# ── Paths ─────────────────────────────────────────────────────────────────────
HPC_MED_RESULTS = (
    "/nfs/roberts/scratch/pi_arp29/msg74/aim1_baseline_final"
    "/hppd_manuscript_mediation_analyses/mediation_analyses/results"
)
HPC_HOST    = "bouchet"
HPC_TMPTAR  = "/tmp/mediation_results.tar.gz"

_SCRIPT_DIR   = Path(__file__).resolve().parent
LOCAL_PROJECT = _SCRIPT_DIR.parent                  # hppd_manuscript_public/
LOCAL_RESULTS = LOCAL_PROJECT / "results"           # where 0X_all_figures.py reads from


def _derive_model_names(mediation_analyses, dvs):
    """
    Compute (dv_col, model_name) pairs from a MEDIATION_ANALYSES config list.

    mediation_analyses : list of dicts with keys spvar, mediator, dv, cov_types
    dvs                : dict mapping dv shorthand → column name (e.g. DVS from generate_hpc_jobs.py)

    Returns a list of (dv_col, model_name) tuples — one per (analysis × cov_type).
    """
    result = []
    for entry in mediation_analyses:
        dv_short       = entry["dv"]
        spvar_short    = entry["spvar"]
        mediator_short = entry["mediator"]
        dv_col         = dvs.get(dv_short, dv_short)
        for cov_type in entry.get("cov_types", ["nice_covariates"]):
            model_name = f"{dv_col}_{spvar_short}_{mediator_short}_{cov_type}"
            result.append((dv_col, model_name))
    return result


def _build_find_cmd(model_pairs):
    """
    Build a remote command that returns all CSVs and DHARMa PNGs under the
    specified (dv_col, model_name) directories.

    model_pairs : list of (dv_col, model_name) tuples.
                  If empty, returns all matching files under HPC_MED_RESULTS.

    Uses grep -E to filter by model name rather than embedding 32 -path
    clauses in the find call (avoids ARG_MAX and readability issues).

    Files pulled: *.csv  +  dharma_*.png  +  pp_check_*.png  +  *.RData fit objects.
    """
    base = (
        f"find {HPC_MED_RESULTS} -type f"
        r" \( -name '*.csv' -o -name 'dharma_*.png' -o -name 'pp_check_*.png' -o -name '*.RData' \)"
    )
    if not model_pairs:
        return base
    # Match any path containing one of the model_name directory names.
    # A model_name like caps_vision_avgdose_vchnu_nice_covariates only appears
    # as a directory component under .../mediation_models/, so this is safe.
    pattern = "|".join(f"/{model_name}/" for _, model_name in model_pairs)
    return f"{base} | grep -E '{pattern}'"


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


def pull_tarball(model_pairs, dry_run: bool) -> "Path | None":
    """
    Create tarball on Bouchet containing CSVs for specified models, scp locally.

    model_pairs : list of (dv_col, model_name) tuples; empty = pull everything.
    """
    find_cmd = _build_find_cmd(model_pairs)

    if dry_run:
        print(f"[dry-run] find command that would run on {HPC_HOST}:")
        print(f"  {find_cmd}")
        print(f"[dry-run] Would tar CSVs from {HPC_MED_RESULTS} → {HPC_TMPTAR}")
        print(f"[dry-run] Would extract into {LOCAL_RESULTS}")
        return None

    # Count before committing
    count_result = _run(["ssh", HPC_HOST, f"({find_cmd}) | wc -l"])
    n_files = int(count_result.stdout.strip())
    print(f"Found {n_files:,} files to pull (CSVs + dharma_*.png + pp_check_*.png + .RData).")

    if n_files == 0:
        print("Nothing to pull — check model names or HPC paths.")
        return None

    # Create tarball on Bouchet relative to HPC_MED_RESULTS so extracted paths
    # land as {dv}/mediation_models/{model_name}/file.csv under local results/
    print(f"Creating tarball on Bouchet ({n_files:,} files) ...")
    tar_cmd = (
        f"({find_cmd}"
        f" | xargs -I{{}} realpath --relative-to={HPC_MED_RESULTS} {{}})"
        f" | tar czf {HPC_TMPTAR} -C {HPC_MED_RESULTS} -T -"
    )
    _run(["ssh", HPC_HOST, tar_cmd])
    print("  Tarball created.")

    local_tar = Path(tempfile.mktemp(suffix=".tar.gz"))
    print(f"Downloading tarball → {local_tar} ...")
    _run(["scp", f"{HPC_HOST}:{HPC_TMPTAR}", str(local_tar)])
    print(f"  Downloaded ({local_tar.stat().st_size / 1e6:.1f} MB).")

    _run(["ssh", HPC_HOST, f"rm -f {HPC_TMPTAR}"])

    return local_tar


def extract_to_results(local_tar: Path, model_pairs) -> int:
    """
    Extract tarball into LOCAL_RESULTS.

    Paths inside the tarball are relative to HPC_MED_RESULTS, so they land
    at LOCAL_RESULTS/{dv}/mediation_models/{model_name}/file.csv.

    Returns number of files extracted.
    """
    LOCAL_RESULTS.mkdir(parents=True, exist_ok=True)
    n = 0
    with tarfile.open(local_tar, "r:gz") as tf:
        members = tf.getmembers()

        # Safety check: reject any absolute paths or path traversal attempts
        for m in members:
            if m.name.startswith("/") or ".." in m.name:
                raise ValueError(f"Unsafe tarball member path: {m.name!r}")

        # Optionally filter to only the expected model paths
        if model_pairs:
            expected_prefixes = tuple(
                f"{dv_col}/mediation_models/{model_name}/"
                for dv_col, model_name in model_pairs
            )
            members = [
                m for m in members
                if any(m.name.startswith(p) for p in expected_prefixes)
            ]

        tf.extractall(str(LOCAL_RESULTS), members=members)
        n = len(members)

    return n


def main():
    use_subset = len(sys.argv) > 1 and sys.argv[1].lower() == "subset"
    dry_run    = "--dry-run" in sys.argv

    if use_subset:
        # Filters come from generate_hpc_jobs.py's CONFIG — the same lists it
        # generated the jobs from, so this can never pull a subset of what was
        # submitted. See load_generator_config() for why the old
        # compile_mediation_results_config.py was deleted.
        sys.path.insert(0, str(_SCRIPT_DIR))
        from compile_nonsp_results import load_generator_config
        cfg = load_generator_config()

        model_pairs = _derive_model_names(cfg["MEDIATION_ANALYSES"], cfg["DVS"])
        print(f"Subset mode: {len(model_pairs)} (dv, model_name) pairs "
              f"from generate_hpc_jobs.py.")
        for dv_col, model_name in model_pairs:
            print(f"  {dv_col}/mediation_models/{model_name}")

    else:
        parser = argparse.ArgumentParser(
            description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
        )
        parser.add_argument(
            "--model-names", nargs="+", metavar="MODEL_NAME",
            help=(
                "Pull only these model_name directories "
                "(e.g. caps_vision_avgdose_vchnu_nice_covariates). "
                "Must also pass --dvs so the script knows which DV subdirectory each belongs to. "
                "Default: pull all CSVs under mediation results."
            ),
        )
        parser.add_argument(
            "--dvs", nargs="+", metavar="DV_COL",
            help=(
                "DV column names paired with --model-names "
                "(e.g. caps_vision persist_vis_yn). "
                "Must be the same length as --model-names."
            ),
        )
        parser.add_argument(
            "--dry-run", action="store_true",
            help="Print what would be pulled without executing.",
        )
        args = parser.parse_args()
        dry_run = args.dry_run

        if args.model_names and args.dvs:
            if len(args.model_names) != len(args.dvs):
                sys.exit("ERROR: --model-names and --dvs must have the same number of entries.")
            model_pairs = list(zip(args.dvs, args.model_names))
        elif args.model_names or args.dvs:
            sys.exit("ERROR: --model-names and --dvs must be provided together.")
        else:
            model_pairs = []   # pull everything
            print("No model filter — pulling all CSVs under mediation results.")

    if not dry_run:
        _confirm_controlmaster()

    local_tar = pull_tarball(model_pairs, dry_run=dry_run)
    if dry_run or local_tar is None:
        return

    print(f"\nExtracting into {LOCAL_RESULTS} ...")
    n_extracted = extract_to_results(local_tar, model_pairs)
    local_tar.unlink(missing_ok=True)

    print(f"\nExtracted {n_extracted:,} files → {LOCAL_RESULTS}")
    if model_pairs:
        dvs_found = sorted({dv_col for dv_col, _ in model_pairs})
        print(f"  DVs: {dvs_found}")
        print(f"  Models: {len(model_pairs)}")


if __name__ == "__main__":
    main()
