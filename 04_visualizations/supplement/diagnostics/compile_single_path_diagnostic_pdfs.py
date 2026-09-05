#!/usr/bin/env /usr/local/bin/python3.12
"""
compile_single_path_diagnostic_pdfs.py
=======================================
Assemble per-DV diagnostic PDFs for the nonsp single-path brms models
whose forest plots appear in the manuscript figures.

Reads 2×2 diagnostic compilation PNGs produced on HPC by
nonsp_diagnostic_worker.R and assembles them into two multi-page PDFs
(one per DV), organized by IV type column:
  col 1 — sp_predictors
  col 2 — vch_behavior
  col 3 — vch_computations

Output:
    results/supplement/diagnostics/single_path_diagnostics_{dv}.pdf

CONFIG (top of this file) mirrors the relevant config from 0X_all_figures.py.
When changing the manuscript model type or IV groups, update both files.

Usage
-----
    cd 04_visualizations/supplement/diagnostics
    python compile_single_path_diagnostic_pdfs.py

Prerequisites
-------------
Diagnostic compilation PNGs must be available locally in HPC_RESULTS_MIRROR.
After running nonsp_diagnostic_worker.R jobs on the HPC, retrieve with:

    cd 03_hpc
    python generate_nonsp_diagnostic_jobs.py --retrieve

This pulls only *_diagnostic_compilation.png files via a targeted tarball
into data/final/nonsp_predictor_analyses/hpc_mirror/. Requires an active
SSH ControlMaster socket (ssh -MNf bouchet). If the mirror directory is
missing, this script prints manual tarball retrieval commands instead.

See 04_visualizations/supplement/README.md § "Pipeline: regenerating
single-path diagnostic PDFs" for the full step-by-step workflow.

The expected file path for each model is:
    {HPC_RESULTS_MIRROR}/{predictor_normalized}/{model_type}/results/diagnostics/{dv}_diagnostic_compilation.png

Public API (called from 0X_all_figures.py)
------------------------------------------
    from compile_single_path_diagnostic_pdfs import compile_single_path_pdfs
    compile_single_path_pdfs()
"""

# ==============================================================================
# CONFIG — keep in sync with 0X_all_figures.py and generate_hpc_jobs.py
# ==============================================================================

from pathlib import Path
import sys
import os

_SCRIPT_DIR   = Path(__file__).parent.resolve()
_VIZ_DIR      = _SCRIPT_DIR.parent.parent       # 04_visualizations/
_PROJECT_ROOT = _VIZ_DIR.parent                  # hppd_manuscript_public/
sys.path.insert(0, str(_PROJECT_ROOT / 'modules'))

# ── IV type groups — mirrors IVTYPES in 0X_all_figures.py ────────────────────
IVTYPES = ["sp_predictors", "vch_behavior", "vch_computations"]

# ── DVs to assemble PDFs for ─────────────────────────────────────────────────
DVS = ["hppd_binary", "caps_vision"]

# ── Model type — mirrors HPPD_MODEL_TYPE in 0X_all_figures.py ────────────────
# Both hppd_binary and caps_vision use the spusers sample in the manuscript
# Change here if the manuscript model type changes.
MODEL_TYPE = "nice_covariates_spusers"

# ── Local HPC results mirror ──────────────────────────────────────────────────
# After rsync'ing from HPC, diagnostic compilations are at:
#   {HPC_RESULTS_MIRROR}/{predictor_normalized}/{MODEL_TYPE}/results/diagnostics/
#              {dv}_diagnostic_compilation.png
#
# Default: data/final/nonsp_predictor_analyses/hpc_mirror/ (one level below
# data/final/nonsp_predictor_analyses/ so that rsync of the entire HPC dir
# lands here without polluting the Python / CSV outputs).
HPC_RESULTS_MIRROR = _PROJECT_ROOT / 'data' / 'final' / 'nonsp_predictor_analyses' / 'hpc_mirror'

# ── Output directory ──────────────────────────────────────────────────────────
FIGURES_DIR = _PROJECT_ROOT / 'results' / 'supplement' / 'diagnostics'

# DPI used when embedding pages into the PDF
_PDF_DPI = 100

# ==============================================================================
# END CONFIG
# ==============================================================================

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.backends.backend_pdf import PdfPages

plt.rcParams['font.family'] = 'Arial'

# ── Predictor lookup ──────────────────────────────────────────────────────────
# Import iv_type_dict and VARIABLE_REGISTRY from master_config.
# iv_type_dict: single source of truth for predictor lists per IVTYPE.
# VARIABLE_REGISTRY: plot_label lookups for figure titles.
try:
    from master_config import iv_type_dict, VARIABLE_REGISTRY
except ImportError as e:
    raise ImportError(
        f'Cannot import master_config: {e}.\n'
        f'Run from 04_visualizations/ so that ../modules/ is on sys.path.'
    )


def _predictor_to_normalized(pred: str) -> str:
    """Return the _normalized column name used as the HPC directory name."""
    if pred.endswith('_normalized'):
        return pred
    return f'{pred}_normalized'


def _build_predictor_index() -> dict[str, list[str]]:
    """
    Return a dict mapping ivtype → list of normalized predictor column names,
    in the order they appear in iv_type_dict.
    """
    result = {}
    for ivtype in IVTYPES:
        preds = iv_type_dict.get(ivtype, [])
        result[ivtype] = [_predictor_to_normalized(p) for p in preds]
    return result


# ═════════════════════════════════════════════════════════════════════════════
# Model discovery
# ═════════════════════════════════════════════════════════════════════════════

def collect_models(dv: str, model_type: str = MODEL_TYPE) -> list[dict]:
    """
    Scan HPC_RESULTS_MIRROR for diagnostic_compilation.png files for the
    given DV and model type.  Returns one entry per available PNG, preserving
    the canonical IVTYPE ordering from iv_type_dict.

    Each entry:
        {ivtype: str, predictor: str, dv: str, png: Path}

    Predictors not in any IVTYPE are skipped with a warning.
    """
    pred_index = _build_predictor_index()

    # Flat ordered list: (ivtype, predictor_normalized) in display order
    ordered = []
    seen_preds = set()
    for ivtype in IVTYPES:
        for pred in pred_index[ivtype]:
            if pred not in seen_preds:
                ordered.append((ivtype, pred))
                seen_preds.add(pred)

    models = []
    missing = []
    for ivtype, pred in ordered:
        png_path = (HPC_RESULTS_MIRROR / pred / model_type
                    / 'results' / 'diagnostics'
                    / f'{dv}_diagnostic_compilation.png')
        if png_path.exists():
            models.append({
                'ivtype':    ivtype,
                'predictor': pred,
                'dv':        dv,
                'png':       png_path,
            })
        else:
            missing.append(f'  {ivtype:<22} {pred}')

    if missing:
        print(f'  [{dv}] {len(missing)} compilation PNGs not found (HPC jobs may not have run):')
        for m in missing[:15]:
            print(m)
        if len(missing) > 15:
            print(f'  ... and {len(missing) - 15} more')

    print(f'  [{dv}] {len(models)} compilations found')
    for ivtype in IVTYPES:
        n = sum(1 for m in models if m['ivtype'] == ivtype)
        print(f'    {ivtype:<22} {n}')

    return models


# ═════════════════════════════════════════════════════════════════════════════
# PDF assembly
# ═════════════════════════════════════════════════════════════════════════════

def _plot_label(varname: str) -> str:
    """Return plot_label from VARIABLE_REGISTRY; fallback to raw name."""
    entry = VARIABLE_REGISTRY.get(varname, {})
    lbl = entry.get('plot_label', None)
    if lbl and str(lbl) != 'nan':
        return str(lbl)
    # Strip _normalized suffix for a readable fallback
    return varname.replace('_normalized', '').replace('_', ' ')


def _write_image_page(pdf: PdfPages, png_path: Path,
                      predictor: str, dv: str) -> None:
    """Write one diagnostic_compilation.png as a full-page PDF page with title."""
    img   = mpimg.imread(str(png_path))
    h, w  = img.shape[:2]
    fig_w = w / _PDF_DPI
    fig_h = h / _PDF_DPI

    # Extra vertical space for the title
    title_h = 0.4
    fig, ax = plt.subplots(figsize=(fig_w, fig_h + title_h))
    fig.patch.set_facecolor('white')
    ax.set_position([0, 0, 1, fig_h / (fig_h + title_h)])
    ax.imshow(img, aspect='auto')
    ax.axis('off')

    # Title: {predictor plot_label} → {dv plot_label}
    pred_raw = predictor.replace('_normalized', '')
    pred_lbl = _plot_label(pred_raw)
    dv_lbl   = _plot_label(dv)
    fig.text(0.5, 1 - 0.12 / (fig_h + title_h),
             f'{pred_lbl}  \u2192  {dv_lbl}',
             ha='center', va='center',
             fontsize=13, fontweight='bold', color='#111111')

    pdf.savefig(fig, dpi=_PDF_DPI, bbox_inches='tight',
                pad_inches=0.05, facecolor='white')
    plt.close(fig)


def assemble_pdf(models: list[dict], dv: str, output_path: Path) -> None:
    """
    Collect compilation PNGs for all models of a given DV and write a
    multi-page PDF.  One page per predictor in canonical IVTYPE order.
    """
    dv_models = [m for m in models if m['dv'] == dv]
    if not dv_models:
        print(f'  [PDF {dv}] No models — skipping.')
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(f'\nAssembling {len(dv_models)}-model PDF → {output_path.name}')

    with PdfPages(str(output_path)) as pdf:
        d = pdf.infodict()
        d['Title']   = f'Single-path nonsp diagnostics — {dv}'
        d['Subject'] = (
            f'2x2 diagnostic compilation panels (PP check, trace, DHARMa, '
            f'DHARMa vs predictor) for all {dv} nonsp single-path models'
        )

        for i, m in enumerate(dv_models):
            _write_image_page(pdf, m['png'], m['predictor'], dv)

            if (i + 1) % 25 == 0:
                print(f'  {i + 1}/{len(dv_models)} pages written...')

    size_mb = output_path.stat().st_size / 1e6
    print(f'  Saved ({size_mb:.1f} MB): {output_path}')


# ═════════════════════════════════════════════════════════════════════════════
# Public API
# ═════════════════════════════════════════════════════════════════════════════

def compile_single_path_pdfs(
    dvs:        list[str] = DVS,
    model_type: str       = MODEL_TYPE,
    ivtypes:    list[str] = IVTYPES,
) -> None:
    """
    Scan HPC_RESULTS_MIRROR and assemble one PDF per DV.

    Parameters
    ----------
    dvs        : list of DV column names to include
    model_type : covariate-set suffix that matches the HPC directory names
    ivtypes    : IV type groups in display order (must match iv_type_dict keys)
    """
    print('=' * 70)
    print('compile_single_path_diagnostic_pdfs')
    print(f'  model_type        : {model_type}')
    print(f'  dvs               : {dvs}')
    print(f'  ivtypes           : {ivtypes}')
    print(f'  HPC_RESULTS_MIRROR: {HPC_RESULTS_MIRROR}')
    print('=' * 70)

    if not HPC_RESULTS_MIRROR.exists():
        # Build the list of (predictor, model_type) pairs we expect so the
        # tarball command only includes those exact directories.
        pred_index  = _build_predictor_index()
        pairs       = sorted({
            (pred, model_type)
            for ivtype in IVTYPES
            for pred in pred_index[ivtype]
        })
        tar_paths = ' \\\n    '.join(
            f'{pred}/{mt}/results/diagnostics/*_diagnostic_compilation.png'
            for pred, mt in pairs
        )
        hpc_base   = '/nfs/roberts/scratch/pi_arp29/msg74/aim1_baseline_final/nonsp_predictor_analyses'
        hpc_login  = 'msg74@bouchet.ycrc.yale.edu'
        hpc_xfer   = 'msg74@transfer-bouchet.ycrc.yale.edu'
        # Write tarball to HPC_BASE (NFS), NOT /tmp/ — /tmp/ is local to each
        # machine and is NOT shared between bouchet (login) and transfer-bouchet.
        tarball_name = 'nonsp_diagnostics.tar.gz'
        hpc_tarball  = f'{hpc_base}/{tarball_name}'
        print(
            f'\nWARNING: HPC_RESULTS_MIRROR does not exist:\n'
            f'  {HPC_RESULTS_MIRROR}\n\n'
            f'Run nonsp_diagnostic_worker.R jobs on the HPC, then retrieve with a\n'
            f'targeted tarball (only diagnostic PNGs — no fit RData / draws):\n\n'
            f'  # Step 1 — SSH to login node and create tarball (written to NFS):\n'
            f'  ssh {hpc_login}\n'
            f'    cd {hpc_base}\n'
            f'    tar czf {hpc_tarball} \\\n'
            f'      {tar_paths}\n'
            f'    exit\n\n'
            f'  # Step 2 — Transfer tarball via transfer node (from local machine):\n'
            f'  scp {hpc_xfer}:{hpc_tarball} /tmp/{tarball_name}\n\n'
            f'  # Step 3 — Extract:\n'
            f'  mkdir -p {HPC_RESULTS_MIRROR}\n'
            f'  tar xzf /tmp/{tarball_name} -C {HPC_RESULTS_MIRROR}\n'
            f'  rm /tmp/{tarball_name}\n\n'
            f'  # Step 4 — Clean up tarball on HPC (optional):\n'
            f'  ssh {hpc_login} "rm {hpc_tarball}"\n'
        )
        return

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    all_models = []
    for dv in dvs:
        print(f'\nCollecting models for dv={dv}:')
        models = collect_models(dv=dv, model_type=model_type)
        all_models.extend(models)

    for dv in dvs:
        out = FIGURES_DIR / f'single_path_diagnostics_{dv}.pdf'
        assemble_pdf(all_models, dv, out)

    print('\nDone.')


# ═════════════════════════════════════════════════════════════════════════════
# CLI entry point
# ═════════════════════════════════════════════════════════════════════════════

def main() -> None:
    compile_single_path_pdfs()


if __name__ == '__main__':
    main()
