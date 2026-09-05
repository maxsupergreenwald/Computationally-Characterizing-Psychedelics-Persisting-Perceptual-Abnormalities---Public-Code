"""
state_trajectories.py
=====================
Publication-quality block-level HGF state trajectory plots for mediation figures.

Plots two vertically stacked panels:
  Top:    xprob    — Contingency Belief (P(stimulus | cue) at end of trial)
  Bottom: xbin_pred — Contingency Prediction (before-trial prediction)

Group lines are split by a DV-specific grouping scheme:

  DV             Grouping scheme     Colors
  -----------    -----------------   ------
  hppd_binary    ppa_history         binary_palette (2-group: PPA (-) vs PPA (+))
  caps_vision    ppa_caps_3group     electric_blue_palette spaced evenly (3-group)

Other DVs can be registered in DV_GROUPING_MAP below.

Aggregation — MEANS, not medians:
  Trial-level states → (subject, block) MEAN → (group, block) MEAN + bootstrap CI
  of the mean (percentile bootstrap, 1000 resamples, 95% CI by default).

  NOTE: The wide-format df_bl columns vch_xprob_block_{n} and vch_xprob_median are
  block MEDIANS (and median-of-medians) from an earlier model fit. They do NOT match
  what this module plots. Do not use vch_xprob_median as a proxy for the trajectory
  line value — load the states CSV and compute means if you need a per-subject summary
  consistent with these figures.

  The single NaN in xbin_pred (trial 1 of block 1) is handled by pandas mean() which
  skips NaN by default; the block-level mean is unaffected.

Inputs:
  Most recent vch_master_withstates*.csv from data_dir (excludes "3lev" and "nominal"
  variants — these are sensitivity analyses, not the primary states).
  Subject-level group labels are derived from df_bl (must contain persist_vis_yn,
  hppd_current, caps_vision; group-label columns are added if not already present).

Output:
  results/{dv}/vch_computations/trajectories/{var_name}.png  — one file per variable
  (relative to project root, inferred as data_dir.parent.parent, or supplied via save_dir).

Key parameters:
  dv            : str  — outcome variable name (determines grouping scheme)
  df_bl         : pd.DataFrame — subject-level data (the shipped analysis dataframe)
  data_dir      : str | Path — directory containing vch_master_withstates*.csv
  sp_only       : bool — if True, restrict to psycheduse_yn == "Yes" (default False)
  savepath      : str | Path | None — explicit save path; overrides default
  show_plot     : bool — whether to call plt.show() (default True)
"""

from __future__ import annotations

import re
import textwrap
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import pandas as pd


def _wrap_ylabel(label: str, max_chars: int = 14) -> str:
    """Wrap a y-axis label string if it exceeds max_chars per line.

    Uses textwrap to split at word boundaries. This prevents matplotlib from
    clipping long rotated y-axis labels that extend beyond the figure canvas
    at large fontsizes (e.g. fontsize >= 24 on a 4.5" tall figure).
    """
    if len(label) <= max_chars:
        return label
    lines = textwrap.wrap(label, width=max_chars, break_long_words=False)
    return '\n'.join(lines) if lines else label


# ---------------------------------------------------------------------------
# DV → grouping configuration
# ---------------------------------------------------------------------------
# Each entry maps a DV name to the grouping scheme used to split trajectory lines.
# 'col'    — column name added by _add_group_labels() to the subject-level df
# 'groups' — ordered list of group values (controls plot order and color assignment)
# 'colors' — resolved at import time from master_config; see _build_palettes()
# ---------------------------------------------------------------------------

def _build_palettes():
    """Import binary_palette and electric_blue_palette from master_config."""
    try:
        import sys as _sys, os as _os
        _here = _os.path.dirname(_os.path.abspath(__file__))
        _modules = _os.path.dirname(_here)
        if _modules not in _sys.path:
            _sys.path.insert(0, _modules)
        from master_config import binary_palette as _bp, electric_blue_palette as _eb
        _three = [_eb[0], _eb[3], _eb[6]]
        return _bp, _three
    except Exception:
        # Fallback: dark vs. light grey for binary; three evenly spaced blues
        _bp = ['#cccccc', '#333333']
        _three = ['#aec6cf', '#4682b4', '#00008b']
        return _bp, _three


_BINARY_PAL, _THREE_GROUP_PAL = _build_palettes()


DV_GROUPING_MAP: dict[str, dict] = {
    'hppd_binary': {
        'col':    'ppa_history_label',
        'groups': ['PPA (-)', 'PPA (+)'],
        'colors': _BINARY_PAL,
    },
    'caps_vision': {
        # Grouped by current CAPS status (CAPS- vs CAPS+), ignoring PPA history.
        # CAPS- = caps_vision == 0; CAPS+ = caps_vision > 0.
        'col':    'caps_binary_label',
        'groups': ['CAPS (-)', 'CAPS (+)'],
        'colors': _BINARY_PAL,
    },
}

# ---------------------------------------------------------------------------
# Y-axis label lookup — dynamic, registry-backed, works for any variable
# ---------------------------------------------------------------------------

def _load_variable_registry():
    """Import VARIABLE_REGISTRY from master_config; returns None on failure."""
    try:
        import sys as _sys, os as _os
        _here = _os.path.dirname(_os.path.abspath(__file__))
        _modules = _os.path.dirname(_here)
        if _modules not in _sys.path:
            _sys.path.insert(0, _modules)
        from master_config import VARIABLE_REGISTRY as _VR
        return _VR
    except Exception:
        return None


_VARIABLE_REGISTRY = _load_variable_registry()


def get_state_label(var: str) -> str:
    """Return the display label for any state variable from VARIABLE_REGISTRY.

    Lookup order: plot_label_verbose → plot_label → raw variable name.
    Works for any variable registered in master_config.VARIABLE_REGISTRY —
    no hardcoding required.  Adding a new variable to VARIABLE_REGISTRY is
    all that is needed to get a correct label here.

    Handles both dict-style entries (current master_config format) and
    namedtuple/object-style entries via getattr fallback.
    """
    if _VARIABLE_REGISTRY is not None and var in _VARIABLE_REGISTRY:
        entry = _VARIABLE_REGISTRY[var]
        for key in ('plot_label',):
            # dict-style entry (current master_config format)
            if isinstance(entry, dict):
                val = entry.get(key)
            else:
                val = getattr(entry, key, None)
            if val is not None and not (isinstance(val, float) and np.isnan(val)):
                return str(val)
    return var   # final fallback: raw column name

# ---------------------------------------------------------------------------
# Helper: find most recent primary vch_master_withstates CSV
# ---------------------------------------------------------------------------

def _find_most_recent_states_csv(data_dir: Path) -> Path:
    """Return the most recent vch_master_withstates CSV that is NOT a 3lev or
    nominal sensitivity variant.

    Two filename forms are accepted:

        vch_master_withstates_MM-DD-YYYY.csv   dated form, as produced upstream
        vch_master_public.csv                  the form shipped in this repo

    Sensitivity variants carry extra tokens between "withstates" and the date
    ("3lev", "nominal", "avg_avg") and are deliberately excluded: they come from
    different HGF fits than the vch_nu / vch_beta / vch_omega parameters the
    manuscript reports, so plotting their trajectories would mix two fits.  The
    dated regex is what enforces that.

    vch_master_public.csv is this repo's copy of the primary dated file, verified
    identical to vch_master_withstates_06-22-2026.csv (77,760 rows, 24 columns,
    all state columns matching).  It is matched by exact name rather than by
    pattern, so no sensitivity variant can ever satisfy it.

    Raises FileNotFoundError if no matching file exists.
    """
    _primary_re = re.compile(
        r'^vch_master_withstates_\d{2}-\d{2}-\d{4}\.csv$'
    )
    _public_name = 'vch_master_public.csv'
    candidates = [
        p for p in data_dir.glob('vch_master_withstates*.csv')
        if _primary_re.match(p.name)
    ]
    _public = data_dir / _public_name
    if _public.is_file():
        candidates.append(_public)
    if not candidates:
        raise FileNotFoundError(
            f'No primary states CSV found in {data_dir}. Expected either '
            f'vch_master_withstates_MM-DD-YYYY.csv or {_public_name}. '
            f'Files present: {[p.name for p in data_dir.glob("vch_master*.csv")]}'
        )
    # Sort by modification time descending; take the most recent
    candidates.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return candidates[0]


# ---------------------------------------------------------------------------
# Helper: add group-label columns to subject-level df
# ---------------------------------------------------------------------------

def _add_group_labels(df_bl: pd.DataFrame) -> pd.DataFrame:
    """Add group-label columns to df_bl.

    Adds:
      ppa_history_label  — 'PPA (-)' / 'PPA (+)'  (from persist_vis_yn)
      caps_binary_label  — 'CAPS (-)' / 'CAPS (+)' (caps_vision == 0 vs > 0);
                           NaN when caps_vision is NaN (non-CAPS completers excluded)

    Required columns in df_bl: persist_vis_yn, caps_vision.
    Returns a copy with additional label columns; existing columns are preserved.
    """
    data = df_bl.copy()
    data['ppa_history_label'] = data['persist_vis_yn'].map({0: 'PPA (-)', 1: 'PPA (+)'})
    # NaN caps_vision → NaN label (participant did not complete CAPS;
    # excluded from trajectory groups rather than folded into CAPS (+)).
    _caps_label = pd.Series(np.nan, index=data.index, dtype=object)
    _caps_valid = data['caps_vision'].notna()
    _caps_label[_caps_valid & (data['caps_vision'] == 0)] = 'CAPS (-)'
    _caps_label[_caps_valid & (data['caps_vision'] >  0)] = 'CAPS (+)'
    data['caps_binary_label'] = _caps_label
    return data


# ---------------------------------------------------------------------------
# Helper: percentile bootstrap CI
# ---------------------------------------------------------------------------

def _bootstrap_ci(
    arr: np.ndarray,
    n_boot: int = 1000,
    ci_level: float = 0.95,
    seed: int = 42,
) -> tuple[float, float, float]:
    """Percentile bootstrap CI for the MEAN (not median).

    The center of each trajectory line is np.mean(arr); bootstrap resamples
    also take the mean of each resample. This is intentional — trajectory
    figures use means throughout. See module docstring for details.

    Returns (mean, lo, hi).  Returns (NaN, NaN, NaN) when n == 0 and
    (mean, NaN, NaN) when n < 2.
    """
    arr = np.asarray(arr, dtype=float)
    arr = arr[~np.isnan(arr)]
    n = len(arr)
    if n == 0:
        return np.nan, np.nan, np.nan
    center = float(np.mean(arr))
    if n < 2:
        return center, np.nan, np.nan
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(n_boot, n))
    boot_means = arr[idx].mean(axis=1)
    alpha = (1.0 - ci_level) / 2.0
    lo = float(np.percentile(boot_means, 100.0 * alpha))
    hi = float(np.percentile(boot_means, 100.0 * (1.0 - alpha)))
    return center, lo, hi


# ---------------------------------------------------------------------------
# Stats computation (separated from plotting so callers can pre-compute
# shared y-axis limits across multiple DVs before generating figures)
# ---------------------------------------------------------------------------

def compute_state_stats(
    dv: str,
    df_bl: pd.DataFrame,
    data_dir: str | Path,
    variables: list[str] | None = None,
    sp_only: bool = False,
    n_boot: int = 1000,
    ci_level: float = 0.95,
    rng_seed: int = 42,
) -> tuple[dict, list]:
    """Compute group-level mean + bootstrap CI for state variables without plotting.

    Performs the data loading, aggregation, and bootstrap CI computation that
    plot_state_trajectories uses internally.  Call this to pre-compute statistics
    for multiple DVs before plotting — e.g. to derive shared y-axis limits so
    panels for different DVs share the same scale.

    Parameters
    ----------
    dv : str
        Outcome variable — must be in DV_GROUPING_MAP.
    df_bl : pd.DataFrame
        The shipped analysis dataframe.
    data_dir : str | Path
        Directory containing vch_master_withstates*.csv.
    variables : list[str] | None
        State variables to compute statistics for. Default ['xprob', 'xbin_pred'].
    sp_only : bool
        If True, restrict to SP users. Default False.
    n_boot, ci_level, rng_seed : int, float, int
        Bootstrap CI parameters. Defaults match plot_state_trajectories.

    Returns
    -------
    var_stats : dict[str, dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]]
        Nested dict: var_stats[variable][group] = (centers, lower_ci, upper_ci).
        Each array has length == number of blocks.
    blocks : list[int]
        Sorted list of block numbers (x-axis positions).
    """
    if variables is None:
        variables = ['xprob', 'xbin_pred']

    data_dir = Path(data_dir)

    if dv not in DV_GROUPING_MAP:
        warnings.warn(
            f"compute_state_stats: DV '{dv}' not in DV_GROUPING_MAP "
            f"({sorted(DV_GROUPING_MAP)}). Returning empty.",
            stacklevel=2,
        )
        return {}, []

    gcfg = DV_GROUPING_MAP[dv]

    states_csv = _find_most_recent_states_csv(data_dir)
    print(f"[state_trajectories] States CSV: {states_csv.name}")
    states = pd.read_csv(states_csv)

    if gcfg['col'] not in df_bl.columns:
        df_bl = _add_group_labels(df_bl)

    required_subj_cols = ['record_id', gcfg['col']]
    if sp_only:
        if 'psycheduse_yn' not in df_bl.columns:
            raise KeyError("sp_only=True requires 'psycheduse_yn' column in df_bl.")
        subj = df_bl[df_bl['psycheduse_yn'] == 'Yes'][required_subj_cols].copy()
    else:
        subj = df_bl[required_subj_cols].copy()

    subj = subj[subj[gcfg['col']].notna()].copy()
    valid_ids = set(subj['record_id'])

    states = states[states['record_id'].isin(valid_ids)].copy()
    if len(states) == 0:
        warnings.warn(
            f"compute_state_stats: No states rows after filtering for DV '{dv}'.",
            stacklevel=2,
        )
        return {}, []

    print(
        f"[state_trajectories] {len(valid_ids)} subjects, "
        f"{states['record_id'].nunique()} with states data, "
        f"{states['block'].nunique()} blocks"
    )

    agg = (
        states
        .groupby(['record_id', 'block'], as_index=False)[variables]
        .mean()
    )
    agg = agg.merge(subj, on='record_id', how='inner')

    blocks = sorted(agg['block'].unique())
    var_stats: dict[str, dict[str, tuple]] = {}
    for var in variables:
        var_stats[var] = {}
        for group in gcfg['groups']:
            gdata = agg[agg[gcfg['col']] == group]
            centers, los, his = [], [], []
            for blk in blocks:
                vals = gdata[gdata['block'] == blk][var].dropna().values
                c, lo, hi = _bootstrap_ci(vals, n_boot=n_boot, ci_level=ci_level, seed=rng_seed)
                centers.append(c)
                los.append(lo)
                his.append(hi)
            var_stats[var][group] = (
                np.array(centers, dtype=float),
                np.array(los,     dtype=float),
                np.array(his,     dtype=float),
            )

    return var_stats, blocks


# ---------------------------------------------------------------------------
# Main public function
# ---------------------------------------------------------------------------

def plot_state_trajectories(
    dv: str,
    df_bl: pd.DataFrame,
    data_dir: str | Path,
    variables: list[str] | None = None,
    sp_only: bool = False,
    save_dir: str | Path | None = None,
    n_boot: int = 1000,
    ci_level: float = 0.95,
    rng_seed: int = 42,
    figsize: tuple[float, float] = (6.5, 4.5),
    linewidth: float = 3.0,
    alpha_fill: float = 0.18,
    tick_fontsize: float = 18,
    label_fontsize: float = 20,
    legend_fontsize: float = 16,
    dpi: int = 200,
    show_plot: bool = True,
    xlabel_variables: list[str] | None = None,
    ylim: dict[str, tuple[float, float]] | None = None,
) -> dict[str, plt.Figure]:
    """Plot one single-panel block trajectory figure per state variable.

    Each variable gets its own figure, saved as:
        {save_dir}/{var_name}.png

    By default, save_dir = {project_root}/results/{dv}/vch_computations/trajectories/
    where project_root = data_dir.parent.parent.

    This modular layout makes it easy to add new state variables later: simply
    extend the `variables` list and each gets its own file without disrupting
    the others.

    Parameters
    ----------
    dv : str
        Outcome variable name — determines the grouping scheme.
        Must be a key in DV_GROUPING_MAP (e.g. 'hppd_binary', 'caps_vision').
    df_bl : pd.DataFrame
        The shipped analysis dataframe.  Must contain:
        record_id, persist_vis_yn, hppd_current, caps_vision.
        Group-label columns are added automatically if absent.
    data_dir : str | Path
        Directory containing vch_master_withstates*.csv.  The most recent
        primary (non-3lev, non-nominal) file is selected automatically.
    variables : list[str] | None
        State variables to plot.  Each produces one figure.
        Default: ['xprob', 'xbin_pred'].
    sp_only : bool
        If True, restrict to psycheduse_yn == "Yes" subjects.  Default False.
    save_dir : str | Path | None
        Output directory for all variable figures.  If None, defaults to
        {project_root}/results/{dv}/vch_computations/trajectories/.
    n_boot, ci_level, rng_seed : int, float, int
        Bootstrap CI parameters.  Default: 1000 resamples, 95% CI, seed=42.
    figsize : tuple[float, float]
        Figure size per panel in inches (width, height).  Default (6.5, 4.5).
    linewidth : float
        Trajectory line width.  Default 3.0.
    alpha_fill : float
        CI band opacity.  Default 0.18.
    tick_fontsize, label_fontsize, legend_fontsize : float
        Font sizes for tick labels, axis labels, and legend.
    dpi : int
        Output resolution.  Default 200.
    show_plot : bool
        Call plt.show() after each figure.  Default True.
    xlabel_variables : list[str] | None
        Subset of `variables` that should display the 'Block' x-axis label.
        If None (default), all variables show the label.  Pass e.g.
        ['xbin_pred'] to show the label only on the bottom panel when
        compositing into a multi-row figure layout.
    ylim : dict[str, tuple[float, float]] | None
        Optional forced y-axis limits per variable, e.g.
        {'xprob': (-2.5, 1.0), 'xbin_pred': (-5.0, -0.5)}.
        When provided, overrides matplotlib's automatic scaling for the
        named variables.  Variables absent from the dict are auto-scaled.
        Use compute_state_stats() across all DVs to derive shared limits
        before calling this function.

    Returns
    -------
    dict[str, matplotlib.figure.Figure]
        Mapping of {variable_name: figure}.  Empty dict if DV not supported.
    """
    if variables is None:
        variables = ['xprob', 'xbin_pred']

    # ── 1–5. Load data, aggregate, compute bootstrap CI ──────────────────────
    var_stats, blocks = compute_state_stats(
        dv=dv, df_bl=df_bl, data_dir=data_dir, variables=variables,
        sp_only=sp_only, n_boot=n_boot, ci_level=ci_level, rng_seed=rng_seed,
    )
    if not var_stats:
        return {}

    gcfg = DV_GROUPING_MAP[dv]

    # ── 6. Resolve output directory ───────────────────────────────────────────
    if save_dir is None:
        project_root = data_dir.parent.parent
        save_dir = project_root / 'results' / dv / 'vch_computations' / 'trajectories'
    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    # ── 7. One figure per variable ────────────────────────────────────────────
    figures: dict[str, plt.Figure] = {}
    for var in variables:
        # constrained_layout=True uses a more reliable margin algorithm than
        # tight_layout() for long rotated y-axis labels at large fontsizes.
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

        for g_idx, group in enumerate(gcfg['groups']):
            color = gcfg['colors'][g_idx]
            centers, los, his = var_stats[var][group]
            valid = ~np.isnan(centers)
            if not valid.any():
                continue
            bx = np.array(blocks)[valid]
            ax.plot(bx, centers[valid], color=color, lw=linewidth, zorder=3)
            ci_mask = valid & ~np.isnan(los) & ~np.isnan(his)
            if ci_mask.any():
                ax.fill_between(
                    np.array(blocks)[ci_mask],
                    los[ci_mask], his[ci_mask],
                    color=color, alpha=alpha_fill, zorder=2,
                )

        ax.set_ylabel(get_state_label(var), fontsize=label_fontsize)
        _show_xlabel = (xlabel_variables is None) or (var in xlabel_variables)
        ax.set_xlabel('Block' if _show_xlabel else '', fontsize=label_fontsize)

        # X-ticks: bottom panel (xlabel_variables) shows only 1 and 12;
        # all other panels show no ticks or labels to reduce clutter.
        if _show_xlabel:
            ax.set_xticks([blocks[0], blocks[-1]])
            ax.set_xticklabels([str(blocks[0]), str(blocks[-1])])
        else:
            ax.set_xticks([])
            ax.set_xticklabels([])

        ax.tick_params(axis='y', labelsize=tick_fontsize, length=6, width=1.5)
        ax.tick_params(
            axis='x',
            labelsize=tick_fontsize,
            length=6 if _show_xlabel else 0,
            width=1.5,
        )
        for spine in ax.spines.values():
            spine.set_visible(False)
        # Legend omitted — caller assembles panels into a figure with a shared legend

        # Apply forced y-axis limits when provided (e.g. to match sibling panels)
        if ylim is not None and var in ylim:
            ax.set_ylim(*ylim[var])

        # Do NOT call tight_layout() — constrained_layout handles margins automatically
        out_path = save_dir / f'{var}.png'
        fig.savefig(out_path, dpi=dpi, bbox_inches='tight', pad_inches=0.1)
        # Also save vector SVG for true-vector figure assembly
        fig.savefig(save_dir / f'{var}.svg', format='svg', bbox_inches='tight', pad_inches=0.1)
        print(f"[state_trajectories] Saved → {out_path}")

        if show_plot:
            plt.show()
        else:
            plt.close(fig)

        figures[var] = fig

    return figures


__all__ = ['plot_state_trajectories', 'compute_state_stats', 'get_state_label', 'DV_GROUPING_MAP']
