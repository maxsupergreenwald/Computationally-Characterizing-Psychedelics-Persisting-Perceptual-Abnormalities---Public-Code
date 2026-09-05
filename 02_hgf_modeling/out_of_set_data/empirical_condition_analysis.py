#!/usr/bin/env python3
"""
empirical_condition_analysis.py
================================
Interrogates the empirical_condition likelihoods derived from the out-of-set cohort,
VCH (visual contrast hallucination) modality only.

This cohort (132 subjects, sex_matched == False) was used to compute the
empirical_condition lookup values that the production HGF fits consume.

Step 1 — Verification
    Checks which computation reproduces the stored ``empirical_condition`` column
    in the VCH data.  Compares three sources per condition:
      (a) unfiltered grand mean  (all subjects with VCH data, no QC)
      (b) 6-flag QC-filtered grand mean  (n=114 vis-passing)
      (c) hgf_pipeline.py hardcoded VCH_CONDITION_TO_EMPIRICAL dict
          — what the production HGF actually consumed

    Key finding:
      • stored CSV column   = unfiltered grand mean
      • pipeline dict       = QC-filtered grand mean  ← used in production
      • VCH condition=0 is hardcoded to 0.0 in the pipeline (not the empirical
        ~11 % false-alarm rate); comment in hgf_pipeline.py: "veridical contrast is 0"

Step 2 — Figure: empirical detection probability by condition
    Two overlaid series, VCH only:
      (A) Full 6-flag QC-passing sample  (n=114)
      (B) Non-hallucinator subset        (n=29)
          Defined as: max_vh_freq == 0 AND max_ah_freq == 0,
          treating NaN as 0 (see behavioral_data_OUT_OF_SET_README.txt:
          "Missing values are subjects who didn't complete the relevant CHAT items
          (typically pure controls who skipped the symptom batteries).")

    For each series and each condition (0, 25, 50, 75 %):
      • Per-subject mean response rates from QC-filtered trials
      • Group mean ± 95 % CI (t-distribution across subjects)
      • One-sample t-test vs. the nominal condition proportion (0, .25, .50, .75)
    Dotted reference lines at the nominal proportions.

Reads:
    behavioral_data_OUT_OF_SET_with_metadata.csv  (same directory)

Writes:
    empirical_condition_verification.csv            Step 1 comparison table (VCH)
    figures/empirical_condition_by_condition.png    Step 2 combined figure

Usage:
    /usr/local/bin/python3.12 empirical_condition_analysis.py

Author: Max Greenwald
Date:   2026-06-17
"""


# ══════════════════════════════════════════════════════════════════════════════
# DATA AVAILABILITY — READ BEFORE RUNNING
#
# This script requires two inputs that are NOT distributed with this
# repository, because both contain individual participant data:
#
#   behavioral_data_OUT_OF_SET_with_metadata.csv
#       Trial-level data from the out-of-set COPE normative cohort
#       (132 subjects, 86,040 trials).  Collected for a separate study;
#       not ours to redistribute.  Its column dictionary and QC-flag
#       definitions ARE included, in
#       behavioral_data_OUT_OF_SET_README.txt.
#
#   <repo>/data/final/df_public_*.csv
#       The main HPPD sample, used to draw the electric-blue comparison series.
#       Read directly via load_public_wide_df(); load_and_prepare_data() is NOT
#       run on it (see below).
#
# The DERIVED outputs of this script are included, so the numbers behind
# the mapping can be checked without the raw data:
#   empirical_condition_verification.csv
#   empirical_condition_full_vs_nonhall.csv
#   figures/empirical_condition_by_condition{,_spusers}.png
#
# WHY THIS SCRIPT MATTERS: it is the derivation of the four numbers that
# every production HGF fit consumes as stimulus intensities.  Those values
# are hardcoded in hgf_pipeline.py, model_fitting_singleagent_forarray.jl,
# ppc_classic_vch.jl, and prior_recovery_vch_mcmc.jl; this is the record of
# where they came from.
# ══════════════════════════════════════════════════════════════════════════════

# ══════════════════════════════════════════════════════════════════════════════
# IMPORTS
# ══════════════════════════════════════════════════════════════════════════════
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from scipy import stats


# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION: PATHS, DATASETS, AND THE MAPPINGS UNDER TEST
#
# Three things are set up here:
#   1. Paths to the out-of-set trial data and the figure output directory.
#   2. The MAIN HPPD sample, read from the public wide df, which is ALREADY the
#      output of load_and_prepare_data() and must not be passed through it again.
#   3. Two frozen copies of the condition → intensity lookup: the values
#      used before 2026-06-17 and the current non-hallucinator values.
#      Keeping both lets Step 1 verify which one production actually ran,
#      rather than assuming.
# ══════════════════════════════════════════════════════════════════════════════
# ── Paths ──────────────────────────────────────────────────────────────────────
THIS_DIR  = os.path.dirname(os.path.abspath(__file__))
DATA_FILE = os.path.join(THIS_DIR, 'behavioral_data_OUT_OF_SET_with_metadata.csv')
FIG_DIR   = os.path.join(THIS_DIR, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

# ── Main dataset (public, PII-free wide df) ───────────────────────────────────
# Read directly.  load_and_prepare_data() is NOT called: the public file is
# already that function's output, so its derived columns (including the
# psycheduse_yn SP-user filter used below) are present, and re-deriving them
# from source fields that were removed for release is not possible.
HGF_MODULE_DIR = os.path.dirname(THIS_DIR)          # 02_hgf_modeling/
sys.path.insert(0, HGF_MODULE_DIR)
from hgf_pipeline import load_public_wide_df
df_main, _MAIN_DF_PATH = load_public_wide_df(
    project='hppd_manuscript', require=['psycheduse_yn']
)
print(f"  Main sample: {len(df_main)} participants "
      f"({os.path.basename(_MAIN_DF_PATH)})")

CONDITIONS = [0, 25, 50, 75]

# ── Pipeline lookup — PRE-2026-06-17 values (frozen for historical comparison) ─
# These were the values hardcoded in hgf_pipeline.py BEFORE this analysis.
# VCH condition=0 was deliberately set to 0.0 ("veridical contrast is 0").
# Conditions 25/50/75 were the QC-filtered full-sample (n=114) grand means.
# On 2026-06-17, the pipeline was switched to VCH_NONHALL (see below).
VCH_PIPELINE_PRE_20260617 = {
    0:  0.0,
    25: 0.44184090362893536,
    50: 0.7010130961205832,
    75: 0.8690302144249513,
}

# ── Pipeline lookup — current values (non-hallucinator, post-2026-06-17) ───────
# Non-hallucinator subset: max_vh_freq == 0 AND max_ah_freq == 0 (NaN → 0), n=29.
# IMPORTANT — 0 % condition is ALWAYS 0.0, not the empirical false-alarm rate.
# There is literally no stimulus at 0 % contrast; ground-truth probability = 0.0.
VCH_NONHALL = {
    0:  0.0,   # ground-truth: no stimulus → 0.0 (NOT the empirical FA rate)
    25: 0.4180444024563061,
    50: 0.7115104419621175,
    75: 0.8994252873563219,
}

# Alias used by Step 1 verification (compare against old pipeline values)
VCH_PIPELINE = VCH_PIPELINE_PRE_20260617


# ══════════════════════════════════════════════════════════════════════════════
# LOAD THE OUT-OF-SET COHORT AND APPLY QC
#
# Two nested subsets are built and both are kept:
#
#   df_vch_all  every VCH trial, no QC.  Step 1 needs this because the
#               empirical_condition column stored in the CSV turns out to
#               be the UNFILTERED grand mean.
#   df_vch      VCH trials passing all six QC flags and with a positive
#               response time.  This is the basis of the reported values.
#
# The non-hallucinator subset is then carved out of df_vch.  Its NaN policy
# is a real analytic decision, not a convenience: see the comment on the
# criterion below for why NaN is read as 0 and what the strict alternative
# would have cost (n=7 instead of n=29).
# ══════════════════════════════════════════════════════════════════════════════
# ── Load data ──────────────────────────────────────────────────────────────────
print("Loading data …")
df = pd.read_csv(DATA_FILE)
print(f"  All rows: {len(df):,}  |  subjects: {df['sudo_rec'].nunique()}")

# Base QC: 6-flag pass + impossible RTs removed (see README)
df_qc = df[df['six_flag_qc_pass'] & (df['responseTime'] > 0)].copy()

# VCH-only slices
df_vch_all = df[df['modality'] == 'v'].copy()          # unfiltered VCH (Step 1)
df_vch     = df_qc[df_qc['modality'] == 'v'].copy()    # QC-filtered VCH (Step 2)

n_full = df_vch['sudo_rec'].nunique()
print(f"  VCH QC-passing: {n_full} subjects, {len(df_vch):,} trials")

# ── Non-hallucinator subset ────────────────────────────────────────────────────
# Criterion: max_vh_freq == 0 AND max_ah_freq == 0.
# NaN is treated as 0 because, per behavioral_data_OUT_OF_SET_README.txt:
#   "Missing values are subjects who didn't complete the relevant CHAT items
#    (typically pure controls who skipped the symptom batteries)."
# Strict == 0 (excluding NaN) yields only n=7; treating NaN as 0 gives n=29,
# which is sufficient for CI estimation.
subj_meta = df_vch[['sudo_rec', 'max_vh_freq', 'max_ah_freq']].drop_duplicates('sudo_rec')
nonhall_mask = (
    (subj_meta['max_vh_freq'].fillna(0) == 0) &
    (subj_meta['max_ah_freq'].fillna(0) == 0)
)
nonhall_ids  = set(subj_meta.loc[nonhall_mask, 'sudo_rec'])
df_nonhall   = df_vch[df_vch['sudo_rec'].isin(nonhall_ids)].copy()
n_nonhall    = df_nonhall['sudo_rec'].nunique()
print(f"  Non-hallucinator subset (VCH, QC-passing, max_v/ah_freq == 0 or NaN): "
      f"{n_nonhall} subjects, {len(df_nonhall):,} trials")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Verification (VCH only)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n─── STEP 1: Verification (VCH) ───")

rows = []
gm_all    = df_vch_all.groupby('condition')['response'].mean()
gm_qc     = df_vch.groupby('condition')['response'].mean()
stored_ec = df_vch_all.groupby('condition')['empirical_condition'].first()

for cond in CONDITIONS:
    rows.append({
        'condition':              cond,
        'stored_csv':             stored_ec[cond],
        'unfiltered_grand_mean':  gm_all[cond],
        'qc_filtered_grand_mean': gm_qc[cond],
        'pipeline_hardcoded':     VCH_PIPELINE[cond],
        'stored_minus_unfiltered':  stored_ec[cond] - gm_all[cond],
        'stored_minus_qc':          stored_ec[cond] - gm_qc[cond],
        'pipeline_minus_qc':        VCH_PIPELINE[cond] - gm_qc[cond],
    })

verif = pd.DataFrame(rows)
out_csv = os.path.join(THIS_DIR, 'empirical_condition_verification.csv')
verif.to_csv(out_csv, index=False)
print(f"  Saved: {out_csv}\n")
print(verif[[
    'condition', 'stored_csv', 'unfiltered_grand_mean',
    'qc_filtered_grand_mean', 'pipeline_hardcoded',
    'stored_minus_unfiltered', 'pipeline_minus_qc',
]].to_string(index=False, float_format='{:.6f}'.format))

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Figure
# ═══════════════════════════════════════════════════════════════════════════════
print("\n─── STEP 2: Figure ───")


# ── Inference: bootstrap over SUBJECTS, not trials ──────────────────────────
# Resampling subjects respects the nesting of trials within participants; a
# trial-level bootstrap would treat correlated trials as independent and
# understate the intervals.  The seed is fixed so the published CIs reproduce.
N_BOOT = 10_000
RNG    = np.random.default_rng(seed=5124)


def bootstrap_test(vals, null, n_boot=N_BOOT, rng=RNG):
    """
    Two-tailed bootstrap test of whether the population mean of `vals`
    equals `null`.  Returns (observed_mean, ci_lo, ci_hi, p_value).
    CI is the percentile (3rd–97th) of bootstrap means (94% CI).
    p-value: proportion of bootstrap means (under the null-shifted
    distribution) at least as extreme as the observed mean.
    """
    vals = np.asarray(vals)
    n = len(vals)
    obs_mean = vals.mean()

    # Bootstrap CIs (resample subjects with replacement)
    boot_means = np.array([rng.choice(vals, size=n, replace=True).mean()
                           for _ in range(n_boot)])
    ci_lo = np.percentile(boot_means, 3.0)
    ci_hi = np.percentile(boot_means, 97.0)

    # Two-tailed p-value: shift data to null, count how often bootstrap
    # means are as extreme as the observed mean
    shifted = vals - obs_mean + null
    boot_null = np.array([rng.choice(shifted, size=n, replace=True).mean()
                          for _ in range(n_boot)])
    obs_diff = abs(obs_mean - null)
    p = np.mean(np.abs(boot_null - null) >= obs_diff)

    return obs_mean, ci_lo, ci_hi, p


# ── Per-condition summary: subject means → group mean + CI + test ───────────
def compute_stats(trial_df, label):
    """
    Given a trial-level df (VCH, QC-filtered), compute per-subject means,
    group mean ± 95 % bootstrap CI, and two-tailed bootstrap test vs.
    nominal proportion.
    """
    subs = trial_df.groupby(['sudo_rec', 'condition'])['response'].mean()
    n_subj = trial_df['sudo_rec'].nunique()
    means, ci_lo, ci_hi, p_vals = [], [], [], []

    for cond in CONDITIONS:
        rates = subs.xs(cond, level='condition').values
        null  = cond / 100.0
        m, lo, hi, p = bootstrap_test(rates, null)
        means.append(m)
        ci_lo.append(lo)
        ci_hi.append(hi)
        p_vals.append(p)
        print(f"  [{label}]  cond={cond:2d}%:  "
              f"mean={m:.3f} (94%CI [{lo:.3f}, {hi:.3f}])  "
              f"nominal={null:.2f}  p={p:.4g}  (n={len(rates)})")

    return {
        'means':   np.array(means),
        'ci_lo':   np.array(ci_lo),
        'ci_hi':   np.array(ci_hi),
        'p_vals':  p_vals,
        'n':       n_subj,
        'label':   label,
    }


def sig_label(p):
    if p < 0.001: return '***'
    if p < 0.01:  return '**'
    if p < 0.05:  return '*'
    return 'ns'


r_full    = compute_stats(df_vch,     f'Full QC-passing (n={n_full})')
r_nonhall = compute_stats(df_nonhall, f'Non-hallucinator (n={n_nonhall})')


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE: EMPIRICAL DETECTION PROBABILITY BY CONDITION
#
# Supplementary Fig. S2a.  Three series are overlaid against dotted
# reference lines at the nominal QUEST-targeted proportions:
#
#   green         full out-of-set QC-passing cohort (n=114)
#   red           out-of-set non-hallucinator subset (n=29) — the source
#                 of the production mapping
#   electric blue the current HPPD sample
#
# The figure's argument is that all three sit well above the nominal
# proportions in the target-present conditions, and agree with one another
# — i.e. QUEST systematically under-targets here, consistently across
# independent samples.  That is the empirical case for feeding measured
# detection rates to the HGF instead of the nominal proportions.
#
# Two variants are emitted: one using the full main sample, one restricted
# to SP users (psycheduse_yn == 'Yes', the canonical filter).
# ══════════════════════════════════════════════════════════════════════════════
# ── Plot ───────────────────────────────────────────────────────────────────────
GREEN  = '#2CA02C'
RED    = '#D7191C'
ELECTRIC_BLUE = '#00BFF0'
GREY   = '#777777'

# ── Compute main-dataset stats (reusable for both variants) ───────────────────
def compute_main_stats(df_sub, label):
    """Bootstrap stats for main dataset (or subset) using vch_bl_yes_* columns."""
    n = len(df_sub)
    means, ci_lo, ci_hi, pvals = [], [], [], []
    for cond in CONDITIONS:
        col = f'vch_bl_yes_{cond}'
        vals = df_sub[col].dropna().values
        null = cond / 100.0
        m, lo, hi, p = bootstrap_test(vals, null)
        means.append(m)
        ci_lo.append(lo)
        ci_hi.append(hi)
        pvals.append(p)
        print(f"  [{label}]  cond={cond:2d}%:  "
              f"mean={m:.3f} (94%CI [{lo:.3f}, {hi:.3f}])  "
              f"nominal={null:.2f}  p={p:.4g}  (n={len(vals)})")
    return n, np.array(means), np.array(ci_lo), np.array(ci_hi), pvals

n_main, main_means, main_ci_lo, main_ci_hi, main_pvals = \
    compute_main_stats(df_main, f'Main dataset (n={len(df_main)})')

# SP-user subset: psycheduse_yn == "Yes" (canonical filter)
df_main_sp = df_main[df_main['psycheduse_yn'] == 'Yes'].copy()
n_main_sp, main_sp_means, main_sp_ci_lo, main_sp_ci_hi, main_sp_pvals = \
    compute_main_stats(df_main_sp, f'Main dataset SP users (n={len(df_main_sp)})')


# ── Plotting function ─────────────────────────────────────────────────────────
def make_figure(main_n, main_m, main_lo, main_hi, main_p, main_label, out_name):
    """Plot the 3-series empirical condition figure."""
    x       = np.array(CONDITIONS, dtype=float)
    x_full  = x - 3.5   # left offset
    x_nh    = x          # center
    x_main  = x + 3.5    # right offset

    fig, ax = plt.subplots(figsize=(8, 5.5))

    # Reference lines at nominal proportions
    for nominal in [0.0, 0.25, 0.50, 0.75]:
        ax.axhline(nominal, color=GREY, linestyle=':', linewidth=1.0, alpha=0.7, zorder=1)

    # ── Full out-of-sample series (green) ──
    ax.errorbar(
        x_full, r_full['means'],
        yerr=[r_full['means'] - r_full['ci_lo'], r_full['ci_hi'] - r_full['means']],
        fmt='o', color=GREEN, markersize=7, linewidth=1.8,
        capsize=4, capthick=1.5, zorder=3,
    )
    for xi, ci_top, m, p in zip(x_full, r_full['ci_hi'], r_full['means'], r_full['p_vals']):
        sl = sig_label(p)
        sig_y = ci_top + 0.030
        ax.text(xi, sig_y, sl, ha='center', va='bottom',
                fontsize=9, color='#1a6b1a' if sl != 'ns' else '#999999', fontweight='bold')
        ax.text(xi, sig_y + 0.050, f'{m:.0%}', ha='center', va='bottom',
                fontsize=6.5, fontfamily='Arial', color='#1a6b1a')
    for xi, m, cond in zip(x_full, r_full['means'], CONDITIONS):
        nominal = cond / 100.0
        ax.plot([xi, xi], [m, nominal], color='k', linestyle=':', linewidth=0.8, alpha=0.5, zorder=2)

    # ── Main dataset series (electric blue) ──
    ax.errorbar(
        x_main, main_m,
        yerr=[main_m - main_lo, main_hi - main_m],
        fmt='D', color=ELECTRIC_BLUE, markersize=7, linewidth=1.8,
        capsize=4, capthick=1.5, zorder=3,
    )
    for xi, ci_top, m, p in zip(x_main, main_hi, main_m, main_p):
        sl = sig_label(p)
        sig_y = ci_top + 0.030
        ax.text(xi, sig_y, sl, ha='center', va='bottom',
                fontsize=9, color='#0080A0' if sl != 'ns' else '#999999', fontweight='bold')
        ax.text(xi, sig_y + 0.050, f'{m:.0%}', ha='center', va='bottom',
                fontsize=6.5, fontfamily='Arial', color='#0080A0')
    for xi, m, cond in zip(x_main, main_m, CONDITIONS):
        nominal = cond / 100.0
        ax.plot([xi, xi], [m, nominal], color='k', linestyle=':', linewidth=0.8, alpha=0.5, zorder=2)

    # ── Out-of-sample no-hallucination-history series (red) ──
    ax.errorbar(
        x_nh, r_nonhall['means'],
        yerr=[r_nonhall['means'] - r_nonhall['ci_lo'], r_nonhall['ci_hi'] - r_nonhall['means']],
        fmt='s', color=RED, markersize=7, linewidth=1.8,
        capsize=4, capthick=1.5, zorder=3,
    )
    for xi, ci_top, m, p in zip(x_nh, r_nonhall['ci_hi'], r_nonhall['means'], r_nonhall['p_vals']):
        sl = sig_label(p)
        sig_y = ci_top + 0.030
        ax.text(xi, sig_y, sl, ha='center', va='bottom',
                fontsize=9, color='#7f0000' if sl != 'ns' else '#999999', fontweight='bold')
        ax.text(xi, sig_y + 0.050, f'{m:.0%}', ha='center', va='bottom',
                fontsize=6.5, fontfamily='Arial', color='#7f0000')
    for xi, m, cond in zip(x_nh, r_nonhall['means'], CONDITIONS):
        nominal = cond / 100.0
        ax.plot([xi, xi], [m, nominal], color='k', linestyle=':', linewidth=0.8, alpha=0.5, zorder=2)

    # ── Axes ──
    ax.set_xlim(-15, 85)
    ax.set_xticks(CONDITIONS)
    ax.set_xticklabels([f'{c}%\nDetection\nProbability' for c in CONDITIONS], fontsize=8)
    ax.set_xlabel('QUEST-Derived Stimulus Intensity (% Detection Probability)', fontsize=18)
    ax.set_ylabel('Empiric Detection Probability', fontsize=18)
    ax.set_ylim(-0.02, 1.12)
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0))
    ax.yaxis.set_major_locator(mtick.MultipleLocator(0.25))
    ax.tick_params(axis='both', labelsize=9)
    ax.spines[['top', 'right']].set_visible(False)

    legend_handles = [
        plt.Line2D([0], [0], marker='o', color=GREEN, markerfacecolor=GREEN,
                   markersize=7, linewidth=1.8, label=f'Full Out-of-Sample  (n={n_full})'),
        plt.Line2D([0], [0], marker='s', color=RED, markerfacecolor=RED,
                   markersize=7, linewidth=1.8,
                   label=f'Out-of-Sample No-Hallucination History  (n={n_nonhall})'),
        plt.Line2D([0], [0], marker='D', color=ELECTRIC_BLUE, markerfacecolor=ELECTRIC_BLUE,
                   markersize=7, linewidth=1.8, label=f'{main_label}  (n={main_n})'),
        plt.Line2D([0], [0], color=GREY, linestyle=':', linewidth=1.0,
                   label='Nominal proportion (0%, 25%, 50%, 75%)'),
    ]
    ax.legend(handles=legend_handles, loc='upper left', fontsize=8.5, frameon=False)

    plt.tight_layout()
    fig_path = os.path.join(FIG_DIR, out_name)
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  Saved: {fig_path}")


# ── Generate both figure variants ─────────────────────────────────────────────
make_figure(n_main, main_means, main_ci_lo, main_ci_hi, main_pvals,
            'Main Dataset', 'empirical_condition_by_condition.png')

make_figure(n_main_sp, main_sp_means, main_sp_ci_lo, main_sp_ci_hi, main_sp_pvals,
            'Main Dataset SP Users', 'empirical_condition_by_condition_spusers.png')


# ══════════════════════════════════════════════════════════════════════════════
# OUTPUT: FULL-SAMPLE vs NON-HALLUCINATOR COMPARISON TABLE
#
# Writes empirical_condition_full_vs_nonhall.csv — the per-condition
# statistics behind the choice of the non-hallucinator values, including
# the difference between the two candidate mappings.  This table is what
# supports the claim that the two curves are near-identical, so restricting
# to non-hallucinators costs nothing in precision.
# ══════════════════════════════════════════════════════════════════════════════
# ── Comparison CSV ─────────────────────────────────────────────────────────────
comp_rows = []
for cond, nominal in zip(CONDITIONS, [c / 100.0 for c in CONDITIONS]):
    i = CONDITIONS.index(cond)
    comp_rows.append({
        'condition':                cond,
        'nominal_proportion':       nominal,
        # Full out-of-sample
        'full_oos_n':               r_full['n'],
        'full_oos_mean':            r_full['means'][i],
        'full_oos_boot_ci_lo_94':   r_full['ci_lo'][i],
        'full_oos_boot_ci_hi_94':   r_full['ci_hi'][i],
        'full_oos_boot_p_vs_nominal': r_full['p_vals'][i],
        # Non-hallucinator subset
        'nonhall_oos_n':            r_nonhall['n'],
        'nonhall_oos_mean':         r_nonhall['means'][i],
        'nonhall_oos_boot_ci_lo_94': r_nonhall['ci_lo'][i],
        'nonhall_oos_boot_ci_hi_94': r_nonhall['ci_hi'][i],
        'nonhall_oos_boot_p_vs_nominal': r_nonhall['p_vals'][i],
        # Main dataset (full)
        'main_n':                   n_main,
        'main_mean':                main_means[i],
        'main_boot_ci_lo_94':       main_ci_lo[i],
        'main_boot_ci_hi_94':       main_ci_hi[i],
        'main_boot_p_vs_nominal':   main_pvals[i],
        # Main dataset (SP users only)
        'main_sp_n':                n_main_sp,
        'main_sp_mean':             main_sp_means[i],
        'main_sp_boot_ci_lo_94':    main_sp_ci_lo[i],
        'main_sp_boot_ci_hi_94':    main_sp_ci_hi[i],
        'main_sp_boot_p_vs_nominal': main_sp_pvals[i],
        # Difference (non-hallucinator minus full OOS)
        'mean_diff_nonhall_minus_full_oos': r_nonhall['means'][i] - r_full['means'][i],
    })

comp_df = pd.DataFrame(comp_rows)
comp_csv = os.path.join(THIS_DIR, 'empirical_condition_full_vs_nonhall.csv')
comp_df.to_csv(comp_csv, index=False, float_format='%.6f')
print(f"  Saved: {comp_csv}")
print()
print(comp_df.to_string(index=False, float_format='{:.4f}'.format))
print("Done.")
