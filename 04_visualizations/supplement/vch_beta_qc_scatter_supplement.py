#!/usr/bin/env /usr/local/bin/python3.12
"""
vch_beta_qc_scatter_supplement.py
==================================
Publication-quality supplementary scatter figure: vch_beta (decision
precision) vs. a curated selection of task quality-control metrics.

Scientific purpose
------------------
Tests whether vch_beta correlates continuously with behavioral indicators
of task disengagement.  If low-beta participants were simply button-mashing
or otherwise disengaged, QC metrics should be worse for them.  The figure
presents this test at publication resolution for the supplement.

Sample
------
SP users only (psycheduse_yn == "Yes").

Panels (6 active; 2 × 3 grid, no empty slots)
----------------------------------------------
Row 0:  effort_qc            Self-report: effort (Likert y-ticks)
        distraction_qc       Self-report: non-distraction (Likert y-ticks)
        n_timeout_trials     Trials with ≥1 Trial Repeat
Row 1:  threshold_empiric    |vch_bl_yes_75 − 0.75|  (calibration)
        z_composite_rt       Streak × rectified RT drop [multiplicative, experimental]
        z_composite_accuracy Streak × rectified accuracy drop [multiplicative, experimental]

Note: n_timeouts_total (total trial repeats) is omitted — n_timeout_trials
(proportion of trials affected) is the more interpretable metric for this figure.
Composite scores use the multiplicative + directional formula from
an earlier exploratory version of this scatter.

Design choices
--------------
- No panel titles (per supplement convention).
- Y-axis label = plot_label from VARIABLE_REGISTRY.
- Long y-labels wrapped at 28 characters using textwrap.fill().
- Likert y-ticks ({1: "Strongly Disagree" … 5: "Strongly Agree"}) on
  effort_qc and distraction_qc.
- Points colored by vch_beta via electric_blue_palette (reversed):
  lower beta → electric blue, higher beta → dark gray.
- Spearman ρ and p annotated in upper-right corner of each panel.
- Global colorbar placed outside the grid on the right edge.
- No grid lines; no axis spines.

Reads
-----
  data/final/df_public_*.csv — the shipped analysis dataframe (most recent wins)
  task_data_vch_short_psychedelic_bl (JSON embedded in df)
  VARIABLE_REGISTRY, electric_blue_palette (modules/master_config.py)

Output
------
  results/supplement/vch_beta_qc_scatter_supplement/supplementary_figure_s7.png

Run from repo root:
  /usr/local/bin/python3.12 04_visualizations/supplement/vch_beta_qc_scatter_supplement.py
"""

import sys
import base64
import gzip
import json
import textwrap
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.stats import spearmanr

warnings.filterwarnings('ignore')

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT  = SCRIPT_DIR.parents[1]   # supplement/ → 04_visualizations/ → repo root

sys.path.insert(0, str(REPO_ROOT / 'modules'))
from master_config import (
    VARIABLE_REGISTRY, electric_blue_palette,
    SP_USER_COL, SP_USER_VALUE,
)
from data_prep import most_recent_public_df

# ══════════════════════════════════════════════════════════════════════════════
# CONFIG
# ══════════════════════════════════════════════════════════════════════════════

DATA_DIR = str(REPO_ROOT / 'data' / 'final')
OUT_DIR  = REPO_ROOT / 'results' / 'supplement' / 'vch_beta_qc_scatter_supplement'
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Continuous colormap: reversed electric_blue_palette so that
#   low vch_beta  → electric blue (clinically relevant direction)
#   high vch_beta → dark gray
BETA_CMAP = LinearSegmentedColormap.from_list(
    'beta_cmap',
    list(reversed(electric_blue_palette)),
)

# Likert tick labels per the user specification, with "Neither Agree nor
# Disagree" wrapped to two lines to prevent overlap on the y-axis.
LIKERT_TICKS = {
    1: "Strongly Disagree",
    2: "Disagree",
    3: "Neither Agree\nnor Disagree",
    4: "Agree",
    5: "Strongly Agree",
}
LIKERT_COLS = {'effort_qc', 'distraction_qc'}   # columns that get Likert y-ticks

# Y-label line-wrap width (characters)
# Set wide enough that quoted Likert statements ("I did the best I could on the tasks",
# "I was not distracted during the tasks") stay on a single line.
YLABEL_WRAP = 45

matplotlib.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

# ── JSON condition → integer % ────────────────────────────────────────────────
_THRESHOLD_MAP = {0: 0, 1: 25, 2: 50, 3: 75}


# ══════════════════════════════════════════════════════════════════════════════
# JSON DECODING HELPERS  (mirrors 02_beta_qc_scatter.py)
# ══════════════════════════════════════════════════════════════════════════════

def _decompress_vch(row, col='task_data_vch_short_psychedelic_bl'):
    raw = row.get(col) if isinstance(row, dict) else getattr(row, col, None)
    if not isinstance(raw, str):
        return None
    try:
        return json.loads(gzip.decompress(base64.b64decode(raw)).decode('utf-8'))
    except Exception:
        return None


def _find_longest_streak(responses):
    n = len(responses)
    if n == 0:
        return 0, 0, 0, -1
    max_len = best_start = best_end = 0
    best_val = -1
    for target in (0, 1):
        i = 0
        while i < n:
            if responses[i] == target:
                j = i
                while j < n and responses[j] == target:
                    j += 1
                if j - i > max_len:
                    max_len, best_start, best_end, best_val = j - i, i, j, target
                i = j
            else:
                i += 1
    return max_len, best_start, best_end, best_val


def compute_qc_metrics(row):
    """
    Decode task JSON and return dict with:
      n_timeouts_total, n_timeout_trials,
      longest_streak_length, delta_median_rt, delta_accuracy

    Returns {} if JSON is unavailable.
    """
    data = _decompress_vch(row)
    if data is None:
        return {}

    result = {}

    # ── Timeout counts ────────────────────────────────────────────────────────
    total_timeouts = 0
    total_timeout_trials = 0
    for comp in ['component_1', 'component_2', 'component_3', 'component_4']:
        if comp not in data or 'timestamps' not in data[comp]:
            continue
        ts     = data[comp]['timestamps']
        n_stim = sum(1 for e in ts if e[0] == 'stim')
        n_resp = sum(1 for e in ts if e[0] == 'resp')
        total_timeouts += n_stim - n_resp
        stim_count = 0
        for event in ts:
            if event[0] == 'stim':
                stim_count += 1
            elif event[0] == 'resp':
                if stim_count >= 2:
                    total_timeout_trials += 1
                stim_count = 0
    result['n_timeouts_total'] = total_timeouts
    result['n_timeout_trials'] = total_timeout_trials

    # ── Trial-level reconstruction ────────────────────────────────────────────
    all_responses  = []
    all_rts        = []
    all_conditions = []
    for comp in ['component_1', 'component_2', 'component_3', 'component_4']:
        if comp not in data:
            continue
        responses   = data[comp].get('response', [])
        rts         = data[comp].get('responseTime', [])
        intensities = data[comp].get('contrasts', [])
        if not responses or not rts or not intensities:
            continue
        ranks = (pd.Series(intensities, dtype=float).rank(method='dense').astype(int) - 1)
        all_responses.extend(responses)
        all_rts.extend(rts)
        all_conditions.extend(ranks.map(_THRESHOLD_MAP).tolist())

    if not all_responses:
        result['longest_streak_length'] = np.nan
        result['delta_median_rt']       = np.nan
        result['delta_accuracy']        = np.nan
        return result

    resp_arr = np.array(all_responses, dtype=int)
    rt_arr   = np.array(all_rts,       dtype=float)
    cond_arr = np.array(all_conditions, dtype=int)

    # ── Longest same-response streak ──────────────────────────────────────────
    slen, s_start, s_end, _ = _find_longest_streak(resp_arr.tolist())
    result['longest_streak_length'] = slen

    if slen > 0:
        within_mask  = np.zeros(len(resp_arr), dtype=bool)
        within_mask[s_start:s_end] = True
        outside_mask = ~within_mask

        rt_w, rt_o = rt_arr[within_mask], rt_arr[outside_mask]
        result['delta_median_rt'] = (
            float(np.median(rt_w) - np.median(rt_o))
            if len(rt_w) and len(rt_o) else np.nan
        )

        correct_arr = (cond_arr > 0).astype(int)
        acc_w = float((resp_arr[within_mask] == correct_arr[within_mask]).mean())
        result['delta_accuracy'] = (
            acc_w - float((resp_arr[outside_mask] == correct_arr[outside_mask]).mean())
            if outside_mask.sum() > 0 else np.nan
        )
    else:
        result['delta_median_rt'] = np.nan
        result['delta_accuracy']  = np.nan

    return result


# ══════════════════════════════════════════════════════════════════════════════
# LOAD DATA
# ══════════════════════════════════════════════════════════════════════════════
print('=' * 65)
print('vch_beta_qc_scatter_supplement.py')
print('=' * 65)

# The shipped analysis dataframe is already fully prepared; read it directly.
# most_recent_public_df() picks the newest data/final/df_public_*.csv.
print('\n[1] Loading analysis dataframe...')
df = pd.read_csv(most_recent_public_df(DATA_DIR), low_memory=False)
df_sp = df[df[SP_USER_COL] == SP_USER_VALUE].copy()   # canonical SP-user filter
print(f'  Full df N={len(df)}  |  SP users: N={len(df_sp)}')

# ══════════════════════════════════════════════════════════════════════════════
# COMPUTE QC METRICS
# ══════════════════════════════════════════════════════════════════════════════
print('\n[2] Decoding task JSON...')
sp = df_sp.copy()

json_rows = []
for _, row in sp.iterrows():
    m = compute_qc_metrics(row)
    m['record_id'] = int(row['record_id'])
    json_rows.append(m)
json_df = pd.DataFrame(json_rows)
sp = sp.merge(json_df, on='record_id', how='left')
print(f'  Decoded {json_df["n_timeouts_total"].notna().sum()}/{len(sp)} participants')

# ── Derived metrics ────────────────────────────────────────────────────────────
sp['threshold_empiric_v_nominal'] = (sp['vch_bl_yes_75'] - 0.75).abs()

# ── Z-score streak metrics ─────────────────────────────────────────────────────
print('\n[3] Z-scoring streak metrics...')
for src, z_col in [
    ('longest_streak_length', 'z_streak_length'),
    ('delta_median_rt',       'z_delta_rt'),
    ('delta_accuracy',        'z_delta_accuracy'),
]:
    vals = sp[src].dropna()
    sp[z_col] = (sp[src] - vals.mean()) / vals.std()
    print(f'  {src}: n={len(vals)}, mean={vals.mean():.3f}, sd={vals.std():.3f}')

# ── Composite scores (multiplicative + directional, experimental) ──────────────
# Formula mirrors 02_beta_qc_scatter.py:
#   Step 1: delta components already z-scored above
#   Step 2: rectify — retain only the concerning direction
#   Step 3: multiply by raw streak length (not z-scored)
#   Step 4: z-score the product
sp['_rt_drop_z']  = np.maximum(0, -sp['z_delta_rt'])
sp['_acc_drop_z'] = np.maximum(0, -sp['z_delta_accuracy'])
sp['_crt_raw']    = sp['longest_streak_length'] * sp['_rt_drop_z']
sp['_cacc_raw']   = sp['longest_streak_length'] * sp['_acc_drop_z']
for src, z_col in [('_crt_raw', 'z_composite_rt'), ('_cacc_raw', 'z_composite_accuracy')]:
    vals = sp[src].dropna()
    sp[z_col] = (sp[src] - vals.mean()) / vals.std()


# ══════════════════════════════════════════════════════════════════════════════
# PANEL DEFINITIONS
# 6 panels in a 2×3 grid — no empty slots.
# plot_label pulled from VARIABLE_REGISTRY.
# ══════════════════════════════════════════════════════════════════════════════
PANELS = [
    # (column_name, row_idx, col_idx)
    # Row 0: self-report QC + trial repeats
    ('effort_qc',                    0, 0),
    ('distraction_qc',               0, 1),
    ('n_timeout_trials',             0, 2),
    # Row 1: behavioral calibration + composite engagement
    ('threshold_empiric_v_nominal',  1, 0),
    ('z_composite_rt',               1, 1),
    ('z_composite_accuracy',         1, 2),
]

N_ROWS, N_COLS = 2, 3


# ══════════════════════════════════════════════════════════════════════════════
# SPEARMAN SUMMARY (terminal)
# ══════════════════════════════════════════════════════════════════════════════
print('\n' + '=' * 65)
print('SPEARMAN CORRELATIONS: vch_beta vs. QC metrics (SP users)')
print('=' * 65)

beta_arr = sp['vch_beta'].values.astype(float)

spearman_results = {}
for col, _, _ in PANELS:
    if col not in sp.columns:
        print(f'  {col:<45} NOT FOUND')
        spearman_results[col] = (np.nan, np.nan, np.nan, 'MISSING')
        continue
    mv = sp[col].values.astype(float)
    mask = ~np.isnan(beta_arr) & ~np.isnan(mv)
    n = mask.sum()
    if n < 3:
        print(f'  {col:<45} n={n} — insufficient')
        spearman_results[col] = (n, np.nan, np.nan, 'INSUFF')
        continue
    rho, p = spearmanr(beta_arr[mask], mv[mask])
    sig = ('***' if p < 0.001 else '**' if p < 0.01
           else '*' if p < 0.05 else '~' if p < 0.10 else 'ns')
    p_str = f'p={p:.4f}'
    print(f'  {col:<45} n={n:<4} ρ={rho:+.3f}  {p_str}  {sig}')
    spearman_results[col] = (n, rho, p, sig)

print('=' * 65)


# ══════════════════════════════════════════════════════════════════════════════
# SCATTER FIGURE
# ══════════════════════════════════════════════════════════════════════════════
print('\n[4] Generating supplement scatter figure...')

beta_all  = sp['vch_beta'].dropna().values
beta_norm = Normalize(vmin=beta_all.min(), vmax=beta_all.max())

fig, axes = plt.subplots(N_ROWS, N_COLS, figsize=(N_COLS * 4.8, N_ROWS * 4.2))
fig.subplots_adjust(hspace=0.50, wspace=0.55)

for panel_idx, (col, row, col_idx) in enumerate(PANELS):
    ax = axes[row, col_idx]
    panel_label = chr(ord('a') + panel_idx)
    ax.text(-0.12, 1.0, panel_label,
            transform=ax.transAxes, ha='left', va='bottom',
            fontsize=36, fontweight='bold', fontfamily='Arial',
            clip_on=False)

    # Y-axis label from VARIABLE_REGISTRY, wrapped for readability
    plot_label = VARIABLE_REGISTRY[col]['plot_label']
    ylabel     = textwrap.fill(str(plot_label), width=YLABEL_WRAP)

    res = spearman_results.get(col, (np.nan, np.nan, np.nan, ''))
    n_valid, rho, p, sig = res

    if col not in sp.columns or np.isnan(n_valid) or n_valid < 3:
        ax.text(0.5, 0.5, f'{col}\nnot available',
                ha='center', va='center', transform=ax.transAxes, color='gray')
        ax.set_ylabel(ylabel, fontsize=8.5)
        continue

    mv   = sp[col].values.astype(float)
    mask = ~np.isnan(beta_arr) & ~np.isnan(mv)
    x    = beta_arr[mask]
    y    = mv[mask]

    colors = BETA_CMAP(beta_norm(x))
    ax.scatter(x, y, c=colors, s=24, alpha=0.65,
               linewidths=0.0, zorder=3, rasterized=True)

    # Line of best fit (OLS, no CI band)
    slope, intercept = np.polyfit(x, y, 1)
    x_line = np.linspace(x.min(), x.max(), 200)
    ax.plot(x_line, slope * x_line + intercept,
            color='#333333', linewidth=1.2, linestyle='-', zorder=4)

    # Spearman annotation (upper-right corner)
    p_str  = f'p={p:.3f}' if p >= 0.001 else 'p<0.001'
    annot  = f'ρ={rho:+.3f}, {p_str} {sig}\n(n={int(n_valid)})'
    ax.text(0.97, 0.97, annot,
            transform=ax.transAxes, ha='right', va='top',
            fontsize=7.5, color='#222222',
            bbox=dict(facecolor='white', alpha=0.78, edgecolor='none', pad=2))

    # Likert y-ticks for effort/distraction items
    if col in LIKERT_COLS:
        ax.set_yticks(list(LIKERT_TICKS.keys()))
        ax.set_yticklabels(list(LIKERT_TICKS.values()), fontsize=6.5)
        ax.set_ylim(0.5, 5.5)

    ax.set_xlabel(
        VARIABLE_REGISTRY['vch_beta']['plot_label'],   # "Decision Precision (β)"
        fontsize=8.5,
    )
    ax.set_ylabel(ylabel, fontsize=8.5)
    ax.tick_params(labelsize=7.5)

    # Remove all spines and grid
    for spine in ax.spines.values():
        spine.set_visible(False)

# ── Shared colorbar (right edge) ──────────────────────────────────────────────
sm = plt.cm.ScalarMappable(cmap=BETA_CMAP, norm=beta_norm)
sm.set_array([])
cbar_ax = fig.add_axes([0.93, 0.12, 0.015, 0.76])
cbar    = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label(
    f'{VARIABLE_REGISTRY["vch_beta"]["plot_label"]}\n(low = electric blue)',
    fontsize=8, labelpad=6,
)
cbar.ax.tick_params(labelsize=7)

out_path = OUT_DIR / 'supplementary_figure_s7.png'
fig.savefig(out_path, dpi=200, bbox_inches='tight')
tiff_path = OUT_DIR / 'supplementary_figure_s7.tiff'
fig.savefig(tiff_path, dpi=200, bbox_inches='tight')
plt.close(fig)
print(f'  Saved → {out_path}')
print(f'  Saved → {tiff_path}')
print('\nDone.')
