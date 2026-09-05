#!/usr/bin/env /usr/local/bin/python3.12
"""
compile_mediation_diagnostic_pdfs.py

Batch generates diagnostic_compilation.png for parseable brms mediation model
directories, filtered to a specific covariate set (model_type), then assembles
per-DV PDFs.

Only models whose directory name ends with the model_type suffix are included.
Default: 'nice_covariates_spusers' — the covariate set used in the manuscript.

Workflow
--------
1. Scan results/{caps_vision,hppd_binary}/mediation_models/ for all directories
   that (a) contain a fit_*.RData file, (b) are parseable by parse_model_info(),
   and (c) whose name ends with model_type.
2. Run _compile_diagnostics_helper.R for each model in parallel (--workers N)
   to generate cached compiled_*.png files.  Already-cached models are skipped
   unless --force is passed.
3. Run matplotlib assembly sequentially (fast; R step already cached) to
   produce diagnostic_compilation.png in each model directory.
4. Assemble per-DV multi-page PDFs, one page per model, ordered by
   (spvar, mediator, model_name).

Output
------
    results/supplement/diagnostics/mediation_diagnostics_caps_vision.pdf
    results/supplement/diagnostics/mediation_diagnostics_hppd_binary.pdf

Usage (CLI)
-----------
    cd 04_visualizations/supplement/diagnostics
    python compile_mediation_diagnostic_pdfs.py                            # default: nice_covariates_spusers
    python compile_mediation_diagnostic_pdfs.py --model-type nice_covariates_spusers
    python compile_mediation_diagnostic_pdfs.py --workers 8  # faster
    python compile_mediation_diagnostic_pdfs.py --force      # re-run R even if cached
    python compile_mediation_diagnostic_pdfs.py --pdf-only   # skip R; assemble from
                                                              # existing PNGs only

Public API (called from 0X_all_figures.py)
------------------------------------------
    from compile_mediation_diagnostic_pdfs import compile_diagnostic_pdfs
    compile_diagnostic_pdfs(model_type='nice_covariates_spusers', workers=4)

    (0X_all_figures.py adds supplement/diagnostics/ to sys.path so this import works.)
"""

import sys
import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.backends.backend_pdf import PdfPages

plt.rcParams['font.family'] = 'Arial'

# ── Paths ─────────────────────────────────────────────────────────────────────
_SCRIPT_DIR = Path(__file__).parent.resolve()
_VIZ_DIR    = _SCRIPT_DIR.parent.parent     # 04_visualizations/
sys.path.insert(0, str(_SCRIPT_DIR))

from create_mediation_diagnostic_compilation import (   # noqa: E402
    parse_model_info,
    generate_compiled_figures,
    make_diagnostic_compilation,
    _compiled_paths,
)

RESULTS     = _VIZ_DIR.parent / 'results'
FIGURES_DIR = RESULTS / 'supplement' / 'diagnostics'
DV_DIRS     = ['caps_vision', 'hppd_binary']

# DPI used when embedding pages into the PDF (lower = smaller file, still readable)
_PDF_DPI = 100


# ═════════════════════════════════════════════════════════════════════════════
# Model discovery
# ═════════════════════════════════════════════════════════════════════════════

def collect_models(model_type: str | None = None) -> list[dict]:
    """
    Return a sorted list of parseable mediation model directories that
    contain a fit_*.RData file.

    Parameters
    ----------
    model_type : str | None
        If given, only include directories whose name ends with this suffix
        (exact match, case-sensitive).  E.g. 'nice_covariates_spusers'.
        Pass None to include all parseable models.

    Each entry: {dir: Path, info: dict, rdata: Path}

    Sorted by (dv, spvar, mediator, model_name) so that mediator variants for
    the same spvar × dv combination appear on consecutive PDF pages.
    """
    models = []
    for dv in DV_DIRS:
        med_dir = RESULTS / dv / 'mediation_models'
        if not med_dir.exists():
            continue
        for d in sorted(med_dir.iterdir()):
            if not d.is_dir():
                continue
            # Exclude interaction models (x_age, x_something naming convention)
            if '_x_' in d.name:
                continue
            # Apply covariate-set filter before any expensive ops
            if model_type and not d.name.endswith(model_type):
                continue
            rdata_files = list(d.glob('fit_*.RData'))
            if not rdata_files:
                continue
            try:
                info = parse_model_info(d)
            except ValueError:
                continue
            models.append({'dir': d, 'info': info, 'rdata': rdata_files[0]})

    models.sort(key=lambda m: (
        m['info']['dv'],
        m['info']['spvar'],
        m['info']['mediator'],
        m['dir'].name,
    ))
    return models


# ═════════════════════════════════════════════════════════════════════════════
# Parallel R generation
# ═════════════════════════════════════════════════════════════════════════════

def _generate_worker(model: dict, force: bool) -> tuple[Path, bool, str | None]:
    """Thread worker: run R helper for one model."""
    try:
        ok = generate_compiled_figures(model['dir'], model['info'], force=force)
        return model['dir'], ok, None
    except Exception as e:
        return model['dir'], False, str(e)


def batch_generate_r(models: list[dict], workers: int = 4,
                     force: bool = False) -> dict[Path, bool]:
    """
    Run _compile_diagnostics_helper.R for all models in parallel.

    Returns a dict mapping model_dir → success (bool).
    Already-cached models (all compiled_* files present) are skipped by
    generate_compiled_figures() unless force=True.
    """
    # Pre-check which need R vs are already cached
    needs_r = [m for m in models
               if force or not all(p.exists() for p in _compiled_paths(m['dir']).values())]
    cached  = [m for m in models if m not in needs_r]

    print(f"\nR step: {len(needs_r)} models need R  |  {len(cached)} already cached")
    if not needs_r:
        return {m['dir']: True for m in models}

    print(f"Running {workers} parallel Rscript workers...")
    results = {m['dir']: True for m in cached}

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_generate_worker, m, force): m for m in needs_r}
        done = 0
        for fut in as_completed(futures):
            mdir, ok, err = fut.result()
            done += 1
            tag = 'OK  ' if ok else f'FAIL'
            err_note = f'  ({err})' if err else ''
            print(f"  [{done:>3}/{len(needs_r)}] {tag}  {mdir.name}{err_note}")
            results[mdir] = ok

    n_ok   = sum(1 for v in results.values() if v)
    n_fail = len(results) - n_ok
    print(f"\nR generation complete: {n_ok} OK, {n_fail} failed.")
    return results


# ═════════════════════════════════════════════════════════════════════════════
# Sequential matplotlib assembly
# ═════════════════════════════════════════════════════════════════════════════

def batch_assemble_pngs(models: list[dict], r_results: dict[Path, bool]) -> None:
    """
    Run the matplotlib assembly step (fast) sequentially for every model
    where R succeeded.  Skips models whose diagnostic_compilation.png is
    already up-to-date (all compiled_* files present and PNG exists).
    """
    to_assemble = [
        m for m in models
        if r_results.get(m['dir'], False)
        and all(p.exists() for p in _compiled_paths(m['dir']).values())
    ]
    print(f"\nMatplotlib assembly: {len(to_assemble)} models...")
    for i, m in enumerate(to_assemble, 1):
        try:
            make_diagnostic_compilation(m['dir'], force_regenerate=False)
            if (i % 20 == 0) or i == len(to_assemble):
                print(f"  {i}/{len(to_assemble)} assembled")
        except Exception as e:
            print(f"  [WARN] Assembly failed for {m['dir'].name}: {e}")


# ═════════════════════════════════════════════════════════════════════════════
# PDF assembly
# ═════════════════════════════════════════════════════════════════════════════

def assemble_pdf(models: list[dict], dv: str, output_path: Path) -> None:
    """
    Collect diagnostic_compilation.png for all models of a given DV and
    write a multi-page PDF.  Pages are in (spvar, mediator, model_name) order.
    """
    dv_models = [m for m in models if m['info']['dv'] == dv]

    pages = []
    skipped = []
    for m in dv_models:
        png = m['dir'] / 'diagnostic_compilation.png'
        if png.exists():
            pages.append((png, m['dir'].name, m['info']))
        else:
            skipped.append(m['dir'].name)

    if skipped:
        print(f"  [PDF {dv}] Skipping {len(skipped)} models with no PNG:")
        for s in skipped[:10]:
            print(f"    {s}")
        if len(skipped) > 10:
            print(f"    ... and {len(skipped) - 10} more")

    if not pages:
        print(f"  [PDF {dv}] No pages — skipping.")
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"\nAssembling {len(pages)}-page PDF → {output_path.name}")

    with PdfPages(str(output_path)) as pdf:
        # Write PDF metadata
        d = pdf.infodict()
        d['Title']   = f'Mediation model diagnostics — {dv}'
        d['Subject'] = (f'MCMC traces, posterior predictive checks, and '
                        f'DHARMa diagnostics for all {dv} mediation models')

        for i, (png_path, model_name, info) in enumerate(pages):
            _write_image_page(pdf, png_path, model_name)

            if (i + 1) % 25 == 0:
                print(f"  {i + 1}/{len(pages)} pages written...")

    size_mb = output_path.stat().st_size / 1e6
    print(f"  Saved ({size_mb:.1f} MB): {output_path}")



def _write_image_page(pdf: PdfPages, png_path: Path, model_name: str) -> None:
    """Write one diagnostic_compilation.png as a PDF page."""
    img    = mpimg.imread(str(png_path))
    h, w   = img.shape[:2]
    fig_w  = w / _PDF_DPI
    fig_h  = h / _PDF_DPI

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.patch.set_facecolor('white')
    ax.imshow(img, aspect='auto')
    ax.axis('off')
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

    pdf.savefig(fig, dpi=_PDF_DPI, bbox_inches='tight',
                pad_inches=0, facecolor='white')
    plt.close(fig)


# ═════════════════════════════════════════════════════════════════════════════
# Public API
# ═════════════════════════════════════════════════════════════════════════════

def compile_diagnostic_pdfs(
    model_type: str = 'nice_covariates_spusers',
    workers:    int  = 4,
    force:      bool = False,
    pdf_only:   bool = False,
) -> None:
    """
    Generate diagnostic compilation PNGs and assemble per-DV PDFs.

    Parameters
    ----------
    model_type : str
        Only process model directories whose name ends with this suffix.
        Matches the covariate-set string used in directory names, e.g.
        'nice_covariates_spusers' (the manuscript default).
    workers : int
        Number of parallel Rscript processes for the R generation step.
    force : bool
        Re-run R helper even if compiled_*.png files already exist.
    pdf_only : bool
        Skip R generation and matplotlib assembly; only assemble PDFs from
        existing diagnostic_compilation.png files.
    """
    print("=" * 70)
    print("compile_mediation_diagnostic_pdfs")
    print(f"  model_type : {model_type}")
    print(f"  workers    : {workers}")
    print(f"  force      : {force}")
    print(f"  pdf-only   : {pdf_only}")
    print("=" * 70)

    models = collect_models(model_type=model_type)
    print(f"\nFound {len(models)} models matching '{model_type}'")
    for dv in DV_DIRS:
        n = sum(1 for m in models if m['info']['dv'] == dv)
        print(f"  {dv}: {n}")

    if not pdf_only:
        r_results = batch_generate_r(models, workers=workers, force=force)
        batch_assemble_pngs(models, r_results)
    else:
        print("\n[pdf-only] Skipping R generation and matplotlib assembly.")

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for dv in DV_DIRS:
        out = FIGURES_DIR / f'mediation_diagnostics_{dv}.pdf'
        assemble_pdf(models, dv, out)

    print("\nDone.")


# ═════════════════════════════════════════════════════════════════════════════
# CLI entry point
# ═════════════════════════════════════════════════════════════════════════════

def main() -> None:
    args       = sys.argv[1:]
    force      = '--force'    in args
    pdf_only   = '--pdf-only' in args
    workers    = 4
    model_type = 'nice_covariates_spusers'
    for i, a in enumerate(args):
        if a == '--workers'    and i + 1 < len(args):
            workers = int(args[i + 1])
        if a == '--model-type' and i + 1 < len(args):
            model_type = args[i + 1]

    compile_diagnostic_pdfs(
        model_type=model_type,
        workers=workers,
        force=force,
        pdf_only=pdf_only,
    )


if __name__ == '__main__':
    main()
