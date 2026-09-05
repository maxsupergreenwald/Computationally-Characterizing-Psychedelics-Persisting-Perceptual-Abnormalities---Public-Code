#!/usr/bin/env /usr/local/bin/python3.12
"""
create_mediation_diagnostic_compilation.py

Generates a multi-panel MCMC and model diagnostics compilation figure for
a single brms mediation model results directory.

ALL diagnostic figures are (re)generated from the .RData fit object by
calling the companion R script _compile_diagnostics_helper.R.  This ensures
white/transparent backgrounds and properly labeled axes.  Generated figures
are cached in the model directory with a "compiled_" prefix; subsequent calls
skip re-generation unless force_regenerate=True.

Layout — dynamic based on DV:

  caps_vision (5 rows × 2–3 columns):
    Row 0  (2 cols):  Posterior predictive checks — left = DV,  right = mediator
    Row 1  (3 cols):  MCMC traces, mu submodel only — spvar→DV | med→DV | spvar→med
    Row 2  (3 cols):  MCMC traces, hu submodel only — spvar→DV | med→DV | spvar→med
                        (col 3 blank — mediator submodel has no hu component)
    Row 3  (2 cols):  DHARMa comprehensive — left = DV,  right = mediator
    Row 4  (3 cols):  DHARMa residuals vs. predictors:
                        col 1: DV resid. vs. spvar
                        col 2: DV resid. vs. mediator
                        col 3: mediator resid. vs. spvar

  hppd_binary and others (4 rows × 2–3 columns):
    Row 0  (2 cols):  Posterior predictive checks
    Row 1  (3 cols):  MCMC traces, mu submodel — spvar→DV | med→DV | spvar→med
    Row 2  (2 cols):  DHARMa comprehensive
    Row 3  (3 cols):  DHARMa residuals vs. predictors

The mu trace row is generated at half height (H3/2) for caps_vision (to match the
split layout) and full height (H3) for all other DVs.  hu traces are always H3/2.

Column headers above row 0 use DAG notation from VARIABLE_REGISTRY plot_labels:
    Left:  "<spvar_label>  →  <dv_label>"
    Right: "<mediator_label>  →  <dv_label>"

Public API:
    make_diagnostic_compilation(model_dir, output_path=None, dpi=150,
                                 force_regenerate=False)

Integration:
    Called from the RUN_DIAGNOSTICS section of 0X_all_figures.py:

        from create_mediation_diagnostic_compilation import make_diagnostic_compilation
        make_diagnostic_compilation('../results/caps_vision/mediation_models/'
                                    'caps_vision_avgdose_vchrate_nice_covariates_spusers')

CLI usage:
    python create_mediation_diagnostic_compilation.py [model_dir] [--force]

    --force : regenerate compiled_ figures even if they already exist.
"""

import sys
import os
import subprocess
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from PIL import Image
import string

plt.rcParams['font.family'] = 'Arial'

# ─── Paths ────────────────────────────────────────────────────────────────────
SCRIPT_DIR     = Path(__file__).parent.resolve()
_PROJECT_ROOT  = SCRIPT_DIR.parent.parent.parent   # hppd_manuscript_public/
_R_HELPER      = SCRIPT_DIR / '_compile_diagnostics_helper.R'

sys.path.insert(0, str(_PROJECT_ROOT / 'modules'))
from master_config import VARIABLE_REGISTRY  # noqa: E402

# ─── Shorthand → full variable name (from 0X_all_figures.py / MEMORY.md) ─────
SPVAR_TO_FULL = {
    'spage':     'psychedelic_age',
    'avgdose':   'avg_life_dose',
    'lifenomic': 'psycheduse_life_nomic',
}
MEDIATOR_TO_FULL = {
    'vchrate':      'vch_bl_yes_0',
    'vchthreshold': 'vch_threshold',
    'vchnu':        'vch_nu',
    'vchbeta':      'vch_beta',
}

# ─── Default model (used by CLI) ─────────────────────────────────────────────
DEFAULT_MODEL_DIR = (
    _PROJECT_ROOT
    / 'results/caps_vision/mediation_models'
    / 'caps_vision_avgdose_vchrate_nice_covariates_spusers'
)

# ─── R figure dimensions (must match gridspec row heights below) ──────────────
# 2-column panels (PP checks, DHARMa comprehensive):  8.5" × 5.2"
# 3-column trace / residual panels at full height:    5.5" × 4.1"
# 3-column trace panels at half height (caps_vision): 5.5" × 2.05"
_W2, _H2  = 8.5, 5.2
_W3, _H3  = 5.5, 4.1
_H3_HALF  = _H3 / 2   # 2.05" — mu trace height for caps_vision; always used for hu traces
_R_DPI    = 150   # DPI for R-generated figures


# ═════════════════════════════════════════════════════════════════════════════
# Model metadata helpers
# ═════════════════════════════════════════════════════════════════════════════

def _get_label(varname: str) -> str:
    """plot_label from VARIABLE_REGISTRY; fallback to the raw column name."""
    entry = VARIABLE_REGISTRY.get(varname, {})
    lbl = entry.get('plot_label', None)
    if lbl and str(lbl) != 'nan':
        return str(lbl)
    return varname


def _mediator_in_dv(mediator: str) -> str:
    """
    Return the column name used as *predictor* in the DV formula for a mediator.
    For need_non_normalized=True variables, the DV formula uses the _normalized version.
    For inplace_normalized=True, the same column is used in both formulas.
    """
    entry = VARIABLE_REGISTRY.get(mediator, {})
    if entry.get('need_non_normalized', False):
        return mediator + '_normalized'
    return mediator          # inplace_normalized=True or student_t renorm-within-sample


def parse_model_info(model_dir) -> dict:
    """
    Infer DV, spvar, and mediator from the model directory name.
    Convention:  {dv}_{spvar_short}_{mediator_short}_{covariate_set}[_flags]
    Returns dict: {dv, spvar, mediator, mediator_in_dv, model_name}
    """
    name = Path(model_dir).name
    DV_OPTIONS = ['caps_vision', 'hppd_binary', 'persist_vis_yn']
    dv = next((d for d in DV_OPTIONS if name.startswith(d + '_')), None)
    if dv is None:
        raise ValueError(f"Cannot infer DV from model name: '{name}'")

    remainder = name[len(dv) + 1:]

    spvar_full = None
    for short, full in SPVAR_TO_FULL.items():
        if remainder.startswith(short + '_'):
            spvar_full = full
            remainder = remainder[len(short) + 1:]
            break

    mediator_full = None
    for short, full in MEDIATOR_TO_FULL.items():
        if remainder.startswith(short + '_') or remainder == short:
            mediator_full = full
            break

    if spvar_full is None or mediator_full is None:
        raise ValueError(
            f"Cannot parse spvar/mediator from '{name}'.\n"
            f"  spvar → {spvar_full!r},  mediator → {mediator_full!r}\n"
            f"  Known spvar shorthands:    {list(SPVAR_TO_FULL)}\n"
            f"  Known mediator shorthands: {list(MEDIATOR_TO_FULL)}"
        )

    return {
        'dv':            dv,
        'spvar':         spvar_full,
        'mediator':      mediator_full,
        'mediator_in_dv': _mediator_in_dv(mediator_full),
        'model_name':    name,
    }


def _compiled_paths(model_dir) -> dict:
    """
    Return canonical paths for the R-generated compiled_ figures.

    Trace files are split into mu-only and hu-only variants:
    - compiled_traces_mu_*: mu submodel coefficients.
      Generated at H3 (full height) for non-caps_vision DVs;
      at H3_HALF for caps_vision (which has a separate hu row).
    - compiled_traces_hu_*: hu/zi submodel coefficients, always H3_HALF.
      For DVs without a hu submodel (e.g. hppd_binary), R saves blank
      placeholder PNGs so this cache check always passes.
    """
    d = Path(model_dir)
    return {
        'pp_dv':               d / 'compiled_pp_check_dv.png',
        'pp_med':              d / 'compiled_pp_check_med.png',
        'dharma_dv':           d / 'compiled_dharma_dv.png',
        'dharma_med':          d / 'compiled_dharma_med.png',
        'trace_mu_spvar_dv':   d / 'compiled_traces_mu_spvar_dv.png',
        'trace_mu_med_dv':     d / 'compiled_traces_mu_med_dv.png',
        'trace_mu_spvar_med':  d / 'compiled_traces_mu_spvar_med.png',
        'trace_hu_spvar_dv':   d / 'compiled_traces_hu_spvar_dv.png',
        'trace_hu_med_dv':     d / 'compiled_traces_hu_med_dv.png',
        'trace_hu_spvar_med':  d / 'compiled_traces_hu_spvar_med.png',
        'resid_dv_sp':         d / 'compiled_resid_dv_sp.png',
        'resid_dv_med':        d / 'compiled_resid_dv_med.png',
        'resid_med_sp':        d / 'compiled_resid_med_sp.png',
    }


# ═════════════════════════════════════════════════════════════════════════════
# R figure generation
# ═════════════════════════════════════════════════════════════════════════════

def _r_str(s) -> str:
    """Escape a Python value as a shell-safe string for passing to Rscript."""
    return str(s)


def generate_compiled_figures(model_dir, info: dict,
                               force: bool = False, timeout: int = 600) -> bool:
    """
    Call _compile_diagnostics_helper.R to generate all compiled_ figures.

    Returns True if the R script completed without error.
    Skips execution if all compiled_ figures already exist, the fit .RData
    is not newer than the oldest compiled PNG, and force=False.
    """
    paths = _compiled_paths(model_dir)
    rdata = Path(model_dir) / f"fit_{info['model_name']}.RData"

    if not force and all(p.exists() for p in paths.values()):
        # Staleness check: regenerate if the fit file is newer than the
        # oldest compiled PNG (e.g. after re-syncing updated HPC results).
        if rdata.exists():
            fit_mtime = rdata.stat().st_mtime
            oldest_compiled = min(p.stat().st_mtime for p in paths.values())
            if fit_mtime > oldest_compiled:
                print("  [R] Fit file is newer than compiled PNGs — regenerating.")
            else:
                print("  [R] All compiled_ figures exist and up to date; skipping.")
                return True
        else:
            print("  [R] All compiled_ figures exist (fit not found); skipping.")
            return True

    if not rdata.exists():
        print(f"  [R] .RData fit not found at {rdata}; cannot generate figures.")
        return False

    if not _R_HELPER.exists():
        print(f"  [R] R helper script not found at {_R_HELPER}")
        return False

    cmd = [
        'Rscript', '--vanilla', str(_R_HELPER),
        _r_str(rdata),
        _r_str(model_dir),
        _r_str(info['spvar']),
        _r_str(_get_label(info['spvar'])),
        _r_str(info['mediator']),
        _r_str(_get_label(info['mediator'])),
        _r_str(info['mediator_in_dv']),
        _r_str(info['dv']),
        _r_str(_get_label(info['dv'])),
        _r_str(_W2), _r_str(_H2),
        _r_str(_W3), _r_str(_H3),
        _r_str(_R_DPI),
        _r_str(_H3_HALF),   # arg 15: half-height for split trace rows
    ]

    print(f"  [R] Generating compiled_ figures (timeout={timeout}s)…")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        # Print R stdout (trimmed to key lines)
        for line in result.stdout.splitlines():
            if any(k in line for k in ('Saved', 'WARN', 'Error', 'Done', '===',
                                        'fit', 'Family', 'Loaded', 'residuals')):
                print(f"    R: {line}")
        if result.returncode != 0:
            print(f"  [R] Rscript exited {result.returncode}:")
            for line in result.stderr.splitlines()[-25:]:
                print(f"    {line}")
            return False
        return True
    except subprocess.TimeoutExpired:
        print(f"  [R] Timed out after {timeout}s.")
        return False
    except FileNotFoundError:
        print("  [R] Rscript not found in PATH.")
        return False


# ═════════════════════════════════════════════════════════════════════════════
# Image loading and panel rendering
# ═════════════════════════════════════════════════════════════════════════════

def _load_img(path) -> np.ndarray | None:
    """Load image as uint8 RGB array; return None if missing."""
    p = Path(path) if path else None
    if p is None or not p.exists():
        return None
    return np.asarray(Image.open(p).convert('RGB'))


def _show_panel(ax, img, title: str | None = None,
                missing_note: str = "Not available",
                label: str | None = None) -> None:
    """
    Render an image (or placeholder) into a matplotlib axes panel.

    All compiled_ figures are generated at dimensions matching the gridspec
    panel sizes, so aspect='auto' (fill panel) introduces negligible distortion
    (< 3 %) compared to the native image aspect ratio.
    """
    ax.set_facecolor('#f8f8f8')
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    if img is not None:
        ax.imshow(img, aspect='auto')
    else:
        ax.text(0.5, 0.5, missing_note,
                transform=ax.transAxes,
                ha='center', va='center',
                fontsize=8, color='#999999', style='italic',
                multialignment='center')

    if title is not None:
        ax.set_title(title, fontsize=9, pad=3, loc='left',
                     fontweight='bold', color='#333333')

    if label is not None:
        ax.text(0.01, 0.98, label,
                transform=ax.transAxes, ha='left', va='top',
                fontsize=10, fontweight='bold', color='#111111',
                bbox=dict(facecolor='white', edgecolor='none',
                          alpha=0.7, pad=2))


# ═════════════════════════════════════════════════════════════════════════════
# Main compilation function
# ═════════════════════════════════════════════════════════════════════════════

def make_diagnostic_compilation(model_dir, output_path=None,
                                 dpi: int = 150,
                                 force_regenerate: bool = False) -> Path:
    """
    Create a diagnostic compilation figure for a brms mediation model directory.

    Parameters
    ----------
    model_dir        : str | Path  — mediation model results directory.
    output_path      : str | Path | None — output path; defaults to
                        model_dir/diagnostic_compilation.png.
    dpi              : int — output DPI for the final compilation PNG.
    force_regenerate : bool — if True, re-run R helper even if compiled_ files exist.

    Returns
    -------
    Path to the saved figure.
    """
    model_dir = Path(model_dir).resolve()
    if output_path is None:
        output_path = model_dir / 'diagnostic_compilation.png'

    # ── Parse model metadata ──────────────────────────────────────────────────
    info = parse_model_info(model_dir)

    dv_label  = _get_label(info['dv'])
    sp_label  = _get_label(info['spvar'])
    med_label = _get_label(info['mediator'])

    print(f"\n{'='*70}")
    print(f"Diagnostic compilation: {info['model_name']}")
    print(f"  DV        : {info['dv']}  ({dv_label})")
    print(f"  spvar     : {info['spvar']}  ({sp_label})")
    print(f"  mediator  : {info['mediator']}  ({med_label})")
    print(f"  med_in_dv : {info['mediator_in_dv']}")
    print(f"{'='*70}")

    # ── Generate compiled_ figures ────────────────────────────────────────────
    generate_compiled_figures(model_dir, info, force=force_regenerate)

    # ── Load images ───────────────────────────────────────────────────────────
    cpaths = _compiled_paths(model_dir)
    imgs   = {k: _load_img(v) for k, v in cpaths.items()}

    for k, img in imgs.items():
        status = "OK" if img is not None else "MISSING"
        flag   = "  ←" if img is None else ""
        print(f"  [{status}] {k:20s}  {cpaths[k].name}{flag}")

    # ── Dynamic layout: 5 rows for caps_vision (split mu/hu trace), 4 for others
    is_caps = info['dv'] == 'caps_vision'

    # R generates mu traces at H3_HALF for caps_vision (split layout) and H3
    # for all other DVs (single trace row).  hu traces are always H3_HALF.
    if is_caps:
        ROW_HEIGHTS = [_H2, _H3_HALF, _H3_HALF, _H2, _H3]
    else:
        ROW_HEIGHTS = [_H2, _H3, _H2, _H3]

    FIG_W = 18.0
    FIG_H = sum(ROW_HEIGHTS) + 3.2   # 3.2" for headers, inter-row gaps, margins

    fig = plt.figure(figsize=(FIG_W, FIG_H))
    fig.patch.set_facecolor('white')

    gs = gridspec.GridSpec(
        len(ROW_HEIGHTS), 6,
        figure        = fig,
        height_ratios = ROW_HEIGHTS,
        hspace        = 0.30,
        wspace        = 0.05,
        left          = 0.05,
        right         = 0.98,
        top           = 0.88,
        bottom        = 0.02,
    )

    # ── Panel label iterator (a, b, c, ...) ────────────────────────────────────
    _lbl_iter = iter(string.ascii_lowercase)

    # ── Row 0: PP checks (2-col, always) ─────────────────────────────────────
    ax_pp_dv  = fig.add_subplot(gs[0, 0:3])
    ax_pp_med = fig.add_subplot(gs[0, 3:6])
    _show_panel(ax_pp_dv,  imgs['pp_dv'],  title='Posterior predictive check',
                label=next(_lbl_iter))
    _show_panel(ax_pp_med, imgs['pp_med'], title='Posterior predictive check',
                label=next(_lbl_iter))

    # ── Trace rows: split (rows 1+2) for caps_vision, single (row 1) otherwise
    ax_mu1 = fig.add_subplot(gs[1, 0:2])
    ax_mu2 = fig.add_subplot(gs[1, 2:4])
    ax_mu3 = fig.add_subplot(gs[1, 4:6])
    _show_panel(ax_mu1, imgs['trace_mu_spvar_dv'],  label=next(_lbl_iter))
    _show_panel(ax_mu2, imgs['trace_mu_med_dv'],    label=next(_lbl_iter))
    _show_panel(ax_mu3, imgs['trace_mu_spvar_med'], label=next(_lbl_iter))

    if is_caps:
        ax_hu1 = fig.add_subplot(gs[2, 0:2])
        ax_hu2 = fig.add_subplot(gs[2, 2:4])
        ax_hu3 = fig.add_subplot(gs[2, 4:6])
        _show_panel(ax_hu1, imgs['trace_hu_spvar_dv'],  label=next(_lbl_iter))
        _show_panel(ax_hu2, imgs['trace_hu_med_dv'],    label=next(_lbl_iter))
        # col 3 (spvar→mediator): mediator has no hu submodel — show placeholder
        _show_panel(ax_hu3, imgs['trace_hu_spvar_med'],
                    missing_note='No hu component\nin mediator submodel',
                    label=next(_lbl_iter))
        dh_row = 3
        rd_row = 4
    else:
        dh_row = 2
        rd_row = 3

    # ── DHARMa comprehensive (2-col) ─────────────────────────────────────────
    ax_dh_dv  = fig.add_subplot(gs[dh_row, 0:3])
    ax_dh_med = fig.add_subplot(gs[dh_row, 3:6])
    _show_panel(ax_dh_dv,  imgs['dharma_dv'],
                title=f'DHARMa  (response: {dv_label})',
                label=next(_lbl_iter))
    _show_panel(ax_dh_med, imgs['dharma_med'],
                title=f'DHARMa  (response: {med_label})',
                label=next(_lbl_iter))

    # ── DHARMa residuals vs. specific predictors (3-col) ─────────────────────
    ax_r1 = fig.add_subplot(gs[rd_row, 0:2])
    ax_r2 = fig.add_subplot(gs[rd_row, 2:4])
    ax_r3 = fig.add_subplot(gs[rd_row, 4:6])
    _show_panel(ax_r1, imgs['resid_dv_sp'],  label=next(_lbl_iter))
    _show_panel(ax_r2, imgs['resid_dv_med'], label=next(_lbl_iter))
    _show_panel(ax_r3, imgs['resid_med_sp'], label=next(_lbl_iter))

    # ── Top column headers (above row 0 PP panels) ───────────────────────────
    for ax, txt in [
        (ax_pp_dv,  f'{sp_label}  →  {dv_label}'),
        (ax_pp_med, f'{med_label}  →  {dv_label}'),
    ]:
        ax.annotate(
            txt,
            xy=(0.5, 1.0), xycoords='axes fraction',
            xytext=(0, 52), textcoords='offset points',
            ha='center', va='bottom',
            fontsize=12, fontweight='bold', color='#111111',
            annotation_clip=False,
        )

    # ── Per-path headers for mu trace row (row 1) ─────────────────────────────
    for ax, txt in [
        (ax_mu1, f'{sp_label}  →  {dv_label}'),
        (ax_mu2, f'{med_label}  →  {dv_label}'),
        (ax_mu3, f'{sp_label}  →  {med_label}'),
    ]:
        ax.annotate(
            txt,
            xy=(0.5, 1.0), xycoords='axes fraction',
            xytext=(0, 28), textcoords='offset points',
            ha='center', va='bottom',
            fontsize=9.5, fontweight='bold', color='#333333',
            annotation_clip=False,
        )

    # ── Per-path headers for DHARMa residuals row ─────────────────────────────
    for ax, txt in [
        (ax_r1, f'{sp_label}  →  {dv_label}'),
        (ax_r2, f'{med_label}  →  {dv_label}'),
        (ax_r3, f'{sp_label}  →  {med_label}'),
    ]:
        ax.annotate(
            txt,
            xy=(0.5, 1.0), xycoords='axes fraction',
            xytext=(0, 28), textcoords='offset points',
            ha='center', va='bottom',
            fontsize=9.5, fontweight='bold', color='#333333',
            annotation_clip=False,
        )

    # ── Row labels (left margin) ───────────────────────────────────────────────
    mu_lbl = 'MCMC\nTraces (μ)' if is_caps else 'MCMC\nTraces'
    row_label_axes = [
        (ax_pp_dv, 'PP Check'),
        (ax_mu1,   mu_lbl),
    ]
    if is_caps:
        row_label_axes.append((ax_hu1, 'MCMC\nTraces (hu)'))
    row_label_axes += [
        (ax_dh_dv, 'DHARMa'),
        (ax_r1,    'Resid.\nvs. Pred.'),
    ]
    for ax, lbl in row_label_axes:
        ax.annotate(
            lbl,
            xy=(0.0, 0.5), xycoords='axes fraction',
            xytext=(-24, 0), textcoords='offset points',
            ha='right', va='center',
            fontsize=9, fontweight='bold', color='#555555',
            rotation=90, annotation_clip=False,
        )

    # ── Figure title: {predictor} → {mediator} → {dv} using plot_labels ─────
    fig.text(0.5, 0.988,
             f'{sp_label}  \u2192  {med_label}  \u2192  {dv_label}',
             ha='center', va='top',
             fontsize=13, fontweight='bold', color='#111111')

    # ── Save ──────────────────────────────────────────────────────────────────
    plt.savefig(str(output_path), dpi=dpi, bbox_inches='tight',
                facecolor='white', edgecolor='none', pad_inches=0.12)
    plt.close(fig)
    print(f"\nSaved: {output_path}")
    return output_path


# ─── CLI entry point ──────────────────────────────────────────────────────────
if __name__ == '__main__':
    args  = [a for a in sys.argv[1:] if not a.startswith('--')]
    flags = {a for a in sys.argv[1:] if a.startswith('--')}

    mdir  = Path(args[0]) if args else DEFAULT_MODEL_DIR
    force = '--force' in flags

    make_diagnostic_compilation(mdir, force_regenerate=force)
