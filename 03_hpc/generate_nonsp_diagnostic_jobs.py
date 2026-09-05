#!/usr/bin/env /usr/local/bin/python3.12
"""
generate_nonsp_diagnostic_jobs.py
==================================
Reads the most recent combined_all_analyses.txt (produced by generate_hpc_jobs.py)
and generates a new dSQ job-array file whose jobs call nonsp_diagnostic_worker.R
to produce 2×2 diagnostic compilation PNGs for each nonsp single-path model.

Workflow
--------
1. Parse every `Rscript nonsp_predictors.R ...` call in combined_all_analyses.txt.
2. Deduplicate by (predictor, model, dv) — bundled job lines contain multiple DVs
   on one line; each DV gets its own diagnostic job.
3. Optionally filter to a specific model type (--model-type).
4. Write nonsp_diagnostic_jobs.txt.
5. Print transfer + dSQ commands.

Usage
-----
    cd 03_hpc

    # 1. Generate job array + print transfer/HPC commands:
    python generate_nonsp_diagnostic_jobs.py
    python generate_nonsp_diagnostic_jobs.py --model-type nice_covariates_spusers
    python generate_nonsp_diagnostic_jobs.py --input /path/to/combined_all_analyses.txt

    # 2. After HPC jobs finish, retrieve compilation PNGs programmatically:
    python generate_nonsp_diagnostic_jobs.py --retrieve
    #    Requires an active SSH ControlMaster: ssh -MNf bouchet
    #    Pulls only *_diagnostic_compilation.png files (~5 MB total) via tarball,
    #    extracts into data/final/nonsp_predictor_analyses/hpc_mirror/, then
    #    prompts to run compile_single_path_diagnostic_pdfs.py.
"""

# ==============================================================================
# CONFIG — keep in sync with generate_hpc_jobs.py
# ==============================================================================

# Path to combined_all_analyses.txt on the HPC (used in printed transfer commands).
HPC_PARENT = '/nfs/roberts/scratch/pi_arp29/msg74/aim1_baseline_final'
HPC_BASE   = f'{HPC_PARENT}/nonsp_predictor_analyses'
HPC_LOGIN  = 'msg74@bouchet.ycrc.yale.edu'
HPC_TRANSFER = 'msg74@transfer-bouchet.ycrc.yale.edu'

# SSH alias for programmatic access — must match ~/.ssh/config Host entry.
# Uses ControlMaster socket so DUO is only needed once (ssh -MNf bouchet).
HPC_HOST   = 'bouchet'

# Remote temp path for the tarball (on bouchet login node; used by --retrieve).
HPC_TMPTAR = '/tmp/nonsp_diagnostics.tar.gz'

# Module load prefix — must match generate_hpc_jobs.py MODULE_PREFIX
MODULE_PREFIX = (
    'module purge && module load foss/2022b && '
    'module load R/4.4.1-foss-2022b && export R_LIBS_USER=$HOME/R/4.4'
)

# Walltime per diagnostic job.  Observed runtime ~5 min; 15 min provides margin.
# DHARMa regenerates posterior predictive samples and assembles a 2×2 PNG.
WALLTIME = '15:00'

# Output job array filename (written alongside this script, to be rsync'd to HPC_BASE)
OUTPUT_STEM = 'nonsp_diagnostic_jobs'

# Default model-type filter: only generate diagnostics for this covariate set.
# Matches the manuscript-relevant model (matches HPPD_MODEL_TYPE in 0X_all_figures.py).
DEFAULT_MODEL_TYPE = 'nice_covariates_spusers'

# ==============================================================================
# END CONFIG
# ==============================================================================

import sys
import os
import re
import shlex
import subprocess
import tarfile
import tempfile
from pathlib import Path

_SCRIPT_DIR   = Path(__file__).resolve().parent
_PROJECT_ROOT = _SCRIPT_DIR.parent

# Default input: local copy of combined_all_analyses.txt produced by generate_hpc_jobs.py
DEFAULT_INPUT = _PROJECT_ROOT / 'data' / 'final' / 'nonsp_predictor_analyses' / 'combined_all_analyses.txt'


def parse_nonsp_calls(line: str) -> list[dict]:
    """
    Extract all (settings, model, dv, predictor) tuples from one job line.

    A single bundled line may contain multiple `Rscript nonsp_predictors.R ...`
    calls separated by `&&`.  We split on the Rscript invocation boundary and
    parse each segment with shlex so that quoted formula strings (sp_and_covs)
    are handled correctly.

    Mediation lines (`hpc_mediation.R`) are silently ignored.

    nonsp_predictors.R positional args:
        [0] settings    — brms family (e.g. hurdle_negbinom_huvary)
        [1] model       — covariate set (e.g. nice_covariates_spusers)
        [2] dv          — DV column name
        [3] iterations
        [4] sp_and_covs — quoted formula RHS
        [5] predictor   — quoted normalized spvar column name
        [6] longornot
        [7] vartype
        [8] df_path
    """
    results = []
    # re.split splits on ALL occurrences, yielding one segment per Rscript call.
    segments = re.split(r'&&\s+Rscript\s+nonsp_predictors\.R\s+', line)
    for seg in segments[1:]:
        seg = seg.strip()
        if not seg:
            continue
        try:
            tokens = shlex.split(seg)
        except ValueError:
            continue
        if len(tokens) < 6:
            continue
        results.append({
            'settings':  tokens[0],
            'model':     tokens[1],
            'dv':        tokens[2],
            'predictor': tokens[5],
        })
    return results


def make_diagnostic_job(settings: str, model: str, dv: str, predictor: str) -> str:
    """
    Build one dSQ job line that calls nonsp_diagnostic_worker.R.

    The `cd {HPC_BASE}` ensures relative paths (predictor/model/...) resolve
    correctly regardless of the dSQ submission directory.
    """
    return (
        f'{MODULE_PREFIX} && '
        f'cd {HPC_BASE} && '
        f'Rscript nonsp_diagnostic_worker.R '
        f'{settings} {model} {dv} "{predictor}"'
    )


def generate_diagnostic_jobs(
    input_path: Path,
    model_type: str | None = DEFAULT_MODEL_TYPE,
    output_dir: Path       = _SCRIPT_DIR,
) -> tuple[Path, list[dict]]:
    """
    Parse combined_all_analyses.txt, deduplicate, optionally filter, and write
    the diagnostic job-array .txt file.

    Returns (txt_path, jobs) where jobs is the list of
    {settings, model, dv, predictor} dicts that were written.
    """
    if not input_path.exists():
        raise FileNotFoundError(
            f'combined_all_analyses.txt not found at {input_path}.\n'
            f'Run generate_hpc_jobs.py first to produce it, then re-run this script.'
        )

    print(f'Reading: {input_path}')

    # Parse all (settings, model, dv, predictor) tuples from all job lines
    seen  = set()
    jobs  = []
    skipped_filter = 0
    n_lines = 0

    with open(input_path) as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line or line.startswith('#'):
                continue
            n_lines += 1
            for call in parse_nonsp_calls(line):
                key = (call['predictor'], call['model'], call['dv'])
                if key in seen:
                    continue
                seen.add(key)

                if model_type and call['model'] != model_type:
                    skipped_filter += 1
                    continue

                jobs.append(call)

    print(f'  Parsed {n_lines} job lines')
    print(f'  Unique (predictor, model, dv) tuples: {len(seen)}')
    if model_type:
        print(f'  After model_type filter ("{model_type}"): {len(jobs)} jobs')
        print(f'  Skipped (different model type): {skipped_filter}')
    else:
        print(f'  Total jobs (no model_type filter): {len(jobs)}')

    if not jobs:
        print('WARNING: No jobs matched the filter — output file will be empty.')

    # Write job array
    out_path = output_dir / f'{OUTPUT_STEM}.txt'
    with open(out_path, 'w') as fh:
        fh.write(f'# nonsp diagnostic compilation job array\n')
        fh.write(f'# Generated by generate_nonsp_diagnostic_jobs.py\n')
        fh.write(f'# model_type filter: {model_type!r}\n')
        fh.write(f'# {len(jobs)} jobs\n')
        fh.write(
            f'# to run: module load dSQ; '
            f'dsq --job-file {OUTPUT_STEM}.txt --mem-per-cpu 8g -t {WALLTIME} --mail-type ALL\n'
        )
        for j in jobs:
            fh.write(make_diagnostic_job(j['settings'], j['model'], j['dv'], j['predictor']) + '\n')

    print(f'\n  Written: {out_path} ({len(jobs)} jobs)')
    return out_path, jobs


def build_tar_retrieve_commands(
    jobs: list[dict],
    local_mirror: Path,
    tarball: str = 'nonsp_diagnostics.tar.gz',
) -> str:
    """
    Build a 3-step retrieve recipe using a targeted tarball.

    The tarball is written to HPC_BASE (NFS-mounted, visible from both the
    login node and the transfer node) — NOT /tmp/, which is local to each
    machine and not shared between bouchet and transfer-bouchet.

    Only the specific {predictor}/{model}/results/diagnostics/ subdirectories
    that appear in `jobs` are included — no fit RData or posterior draws.
    Returns the command block as a printable string.
    """
    # Unique (predictor, model) pairs — these are the only directories we need
    pairs = sorted({(j['predictor'], j['model']) for j in jobs})

    # Build the glob paths for tar on the HPC side:
    #   predictor/model/results/diagnostics/*_diagnostic_compilation.png
    tar_paths = ' \\\n    '.join(
        f'{pred}/{model}/results/diagnostics/*_diagnostic_compilation.png'
        for pred, model in pairs
    )

    hpc_tarball = f'{HPC_BASE}/{tarball}'

    lines = [
        f'# Step 1 — SSH to login node and create a small targeted tarball',
        f'#          (written to HPC_BASE on NFS — visible from transfer node):',
        f'ssh {HPC_LOGIN}',
        f'  cd {HPC_BASE}',
        f'  tar czf {hpc_tarball} \\',
        f'    {tar_paths}',
        f'  exit',
        f'',
        f'# Step 2 — Transfer the tarball via the transfer node (from local machine):',
        f'scp {HPC_TRANSFER}:{hpc_tarball} /tmp/{tarball}',
        f'',
        f'# Step 3 — Extract into local HPC mirror:',
        f'mkdir -p {local_mirror}',
        f'tar xzf /tmp/{tarball} -C {local_mirror}',
        f'rm /tmp/{tarball}',
        f'',
        f'# Step 4 — Clean up tarball on HPC (optional):',
        f'ssh {HPC_LOGIN} "rm {hpc_tarball}"',
        f'',
        f'# Then assemble PDFs:',
        f'cd {_PROJECT_ROOT / "04_visualizations"}',
        f'python compile_single_path_diagnostic_pdfs.py',
    ]
    return '\n'.join(lines)


# ==============================================================================
# Programmatic retrieval (--retrieve mode)
# ==============================================================================

def _run(cmd: list[str], check: bool = True) -> subprocess.CompletedProcess:
    """Run a subprocess command, raising RuntimeError on non-zero exit."""
    result = subprocess.run(cmd, capture_output=True, text=True)
    if check and result.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\n"
            f"stdout: {result.stdout.strip()}\n"
            f"stderr: {result.stderr.strip()}"
        )
    return result


def _confirm_controlmaster() -> None:
    """Require the user to confirm DUO auth before any remote calls run."""
    print('\nThis step pulls files from Bouchet via SSH + scp.')
    print('It requires an active SSH ControlMaster socket.')
    print('\nHave you run:  ssh -MNf bouchet')
    print('and approved the DUO push on your phone?\n')
    answer = input('Type "yes" to continue: ').strip().lower()
    if answer != 'yes':
        print('Skipping retrieval. Run `ssh -MNf bouchet`, approve DUO, then rerun.')
        raise SystemExit(0)
    print()


def retrieve_diagnostic_compilations(
    jobs:         list[dict],
    local_mirror: Path,
) -> None:
    """
    Pull *_diagnostic_compilation.png files from HPC via a targeted tarball.

    Steps:
      1. SSH to bouchet (via ControlMaster) and create /tmp/nonsp_diagnostics.tar.gz
         containing only the specific compilation PNGs for the given jobs.
      2. scp the tarball to a local tempfile.
      3. Extract into local_mirror, preserving the {predictor}/{model}/... path.
      4. Clean up the remote tarball.

    Requires an active ControlMaster: ssh -MNf bouchet
    """
    _confirm_controlmaster()

    pairs     = sorted({(j['predictor'], j['model']) for j in jobs})
    n_pairs   = len(pairs)
    n_jobs    = len(jobs)

    print(f'Retrieving compilation PNGs for {n_jobs} jobs ({n_pairs} predictor×model pairs)...')

    # ── Step 1: create tarball on bouchet ─────────────────────────────────────
    # Paths are relative to HPC_BASE so the extract lands at
    # {local_mirror}/{predictor}/{model}/results/diagnostics/*.png
    # Paths must be UNQUOTED so the remote bash -c shell can expand the * glob.
    # Double-quoting the paths (e.g. "pred/model/.../*...") prevents glob
    # expansion — bash does not expand wildcards inside double quotes — causing
    # ls/tar to find 0 files. HPC paths never contain spaces so unquoted is safe.
    glob_args = ' '.join(
        f'{pred}/{model}/results/diagnostics/*_diagnostic_compilation.png'
        for pred, model in pairs
    )
    # Use bash -c so the shell expands the globs; cd first so paths are relative.
    # `ls` check lets us report how many files matched before creating the archive.
    check_cmd = f'cd {HPC_BASE} && ls {glob_args} 2>/dev/null | wc -l'
    tar_cmd   = (
        f'cd {HPC_BASE} && '
        f'tar czf {HPC_TMPTAR} {glob_args} 2>/dev/null; '
        f'echo "exit:$?"'
    )

    print(f'  Counting available PNGs on {HPC_HOST}...')
    count_res = _run(['ssh', HPC_HOST, f'bash -c \'{check_cmd}\''])
    n_found   = int(count_res.stdout.strip() or 0)
    print(f'  Found {n_found} compilation PNGs (expected up to {n_jobs}).')

    if n_found == 0:
        print('  Nothing to retrieve — HPC diagnostic jobs may not have finished yet.')
        return

    print(f'  Creating tarball on {HPC_HOST}: {HPC_TMPTAR}')
    _run(['ssh', HPC_HOST, f'bash -c \'{tar_cmd}\''])
    print('  Tarball created.')

    # ── Step 2: scp to local tempfile ─────────────────────────────────────────
    local_tar = Path(tempfile.mktemp(suffix='.tar.gz'))
    print(f'  Downloading → {local_tar} ...')
    _run(['scp', f'{HPC_HOST}:{HPC_TMPTAR}', str(local_tar)])
    size_mb = local_tar.stat().st_size / 1e6
    print(f'  Downloaded ({size_mb:.1f} MB).')

    # ── Step 3: extract into local_mirror ─────────────────────────────────────
    local_mirror.mkdir(parents=True, exist_ok=True)
    print(f'  Extracting → {local_mirror} ...')
    with tarfile.open(local_tar, 'r:gz') as tf:
        tf.extractall(local_mirror)
    local_tar.unlink()

    # Count what was extracted
    n_extracted = sum(1 for _ in local_mirror.rglob('*_diagnostic_compilation.png'))
    print(f'  Extracted {n_extracted} compilation PNGs.')

    # ── Step 4: clean up remote tarball ───────────────────────────────────────
    _run(['ssh', HPC_HOST, f'rm -f {HPC_TMPTAR}'], check=False)
    print('  Remote tarball cleaned up.')

    # ── Prompt to assemble PDFs ────────────────────────────────────────────────
    pdf_script = _PROJECT_ROOT / '04_visualizations' / 'compile_single_path_diagnostic_pdfs.py'
    print(f'\nRun compile_single_path_diagnostic_pdfs.py now to assemble per-DV PDFs?')
    answer = input('  [y/N] ').strip().lower()
    if answer == 'y':
        subprocess.run([sys.executable, str(pdf_script)], check=False)
    else:
        print(f'  Skipped. Run manually: python {pdf_script}')


def main() -> None:
    args       = sys.argv[1:]
    model_type = DEFAULT_MODEL_TYPE
    input_path = DEFAULT_INPUT
    do_retrieve = False

    for i, a in enumerate(args):
        if a == '--model-type' and i + 1 < len(args):
            model_type = args[i + 1]
        if a == '--all-models':
            model_type = None   # include all model types
        if a == '--input' and i + 1 < len(args):
            input_path = Path(args[i + 1])
        if a == '--retrieve':
            do_retrieve = True

    local_mirror = _PROJECT_ROOT / 'data' / 'final' / 'nonsp_predictor_analyses' / 'hpc_mirror'

    # ── --retrieve mode: skip job generation, just pull compilation PNGs ──────
    if do_retrieve:
        print('=' * 70)
        print('generate_nonsp_diagnostic_jobs.py  [--retrieve mode]')
        print(f'  model_type  : {model_type!r}')
        print(f'  local_mirror: {local_mirror}')
        print('=' * 70)
        # Re-parse jobs from combined_all_analyses.txt so we know which
        # (predictor, model) pairs to include in the tarball.
        _, jobs = generate_diagnostic_jobs(
            input_path=Path(input_path),
            model_type=model_type,
            output_dir=_SCRIPT_DIR,
        )
        retrieve_diagnostic_compilations(jobs, local_mirror)
        return

    # ── Normal mode: generate job array + print instructions ──────────────────
    print('=' * 70)
    print('generate_nonsp_diagnostic_jobs.py')
    print(f'  input      : {input_path}')
    print(f'  model_type : {model_type!r}  (None = all)')
    print(f'  output stem: {OUTPUT_STEM}')
    print('=' * 70)

    out_path, jobs = generate_diagnostic_jobs(
        input_path=Path(input_path),
        model_type=model_type,
        output_dir=_SCRIPT_DIR,
    )

    # ── Print transfer + HPC commands ─────────────────────────────────────────
    out_abs = os.path.abspath(out_path)
    worker  = os.path.abspath(_SCRIPT_DIR / 'nonsp_diagnostic_worker.R')

    print(f'\n{"=" * 70}')
    print('TRANSFER (run locally — copy both files to HPC):')
    print(f'{"=" * 70}')
    print(f'scp {worker} {HPC_TRANSFER}:{HPC_BASE}/')
    print(f'scp {out_abs} {HPC_TRANSFER}:{HPC_BASE}/')

    print(f'\n{"=" * 70}')
    print(f'HPC COMMANDS (run on {HPC_LOGIN}):')
    print(f'{"=" * 70}')
    print(f'ssh {HPC_LOGIN}')
    print(f'  cd {HPC_BASE}')
    print(f'  module load dSQ')
    print(f'  dsq --job-file {OUTPUT_STEM}.txt --mem-per-cpu 8g -t {WALLTIME} --mail-type ALL')
    print(f'  # → run the sbatch command that dSQ prints')

    print(f'\n{"=" * 70}')
    print('RETRIEVE compilations when jobs finish:')
    print(f'  Option A (programmatic — requires ssh -MNf bouchet):')
    print(f'    python {os.path.basename(__file__)} --retrieve')
    print(f'  Option B (manual tarball commands):')
    print(f'{"=" * 70}')
    print(build_tar_retrieve_commands(jobs, local_mirror))
    print(f'{"=" * 70}')


if __name__ == '__main__':
    main()
