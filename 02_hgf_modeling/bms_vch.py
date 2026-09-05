"""
bms_vch.py  ─  RFX-BMS comparing N VCH HGF model variants
============================================================

Implements Random-Effects Bayesian Model Selection (Stephan et al., 2009,
NeuroImage 46:1004–1017) with Protected Exceedance Probabilities (Rigoux
et al., 2014, NeuroImage 84:971–985) in pure Python/NumPy.

Algorithm overview
------------------
Given per-subject log model evidence LME_nm (subject n, model m):

1. Variational E/M updates on a Dirichlet distribution over model frequencies r:
     E-step:  u_nm  ∝ exp( LME_nm + ψ(α_m) − ψ(Σα) )   [normalised across m]
     M-step:  α_m   =  α0 + Σ_n u_nm
   Iterate to convergence.  Starting α0 = 1 (uniform Dirichlet prior).

2. Expected frequency:  Ef_m = α_m / Σα

3. Exceedance probability XP_m = P(r_m > r_j ∀j≠m | data)
   Estimated by drawing 10⁶ samples from Dirichlet(α).

4. Bayesian Omnibus Risk (BOR, Rigoux et al. 2014):
     log p(y|H0) = Σ_n log[ (1/K) Σ_m exp(LME_nm) ]   [null: all models equal]
     log p(y|H1) ≈ F_VB                                  [VB free energy]
     BOR = sigmoid( log p(H0|y) )   [equal priors on H0, H1]

5. Protected XP:  PXP_m = (1 − BOR)·XP_m + BOR/K

Sections controlled by toggle flags:
  DO_RUN   — launch bms_vch.jl locally (skips participants with existing output)
  DO_BMS   — load LME CSVs and run RFX-BMS
  DO_PLOT  — bar charts of Ef and PXP
"""

################################################################################
# ─── ROLE IN THE MANUSCRIPT ───────────────────────────────────────────────────
#
# MODEL SELECTION, part 2 of 2.  Produces Supplementary Fig. S3c.
#
#     bms_vch.jl (per-subject evidence) ──▶ THIS SCRIPT (group-level RFX-BMS)
#
# Why random-effects BMS rather than a fixed-effects comparison
#     Fixed-effects comparison assumes one model generated everyone's data, so a
#     handful of subjects with extreme evidence can dominate the group result.
#     RFX-BMS instead treats model identity as a random variable across subjects
#     and estimates the population frequency of each model, which is the right
#     assumption for a clinical sample that may well be heterogeneous.
#
#     Protected exceedance probability (PXP) guards the headline number against
#     the possibility that no model is genuinely better than any other: it mixes
#     the raw exceedance probability with the uniform null in proportion to the
#     Bayesian Omnibus Risk.  A low BOR means the models really do differ.
#
# Reproducibility note
#     bms_vch.jl does not set a random seed for its MAP random-restart fallback,
#     so its per-subject evidence output — and therefore the Ef/PXP/BOR this
#     script computes from it — varies slightly between runs. The XP-sampling
#     step below (n_xp_samples draws from Dirichlet(alpha)) is itself seeded
#     (default seed=42), but that only fixes the sampling noise on top of
#     whatever per-subject evidence bms_vch.jl produced that run.
#
# Output
#     model_comparison/bms/bms_summary.csv  and  figures/bms_ef_pxp.png.
#     The published panel was re-rendered from that CSV, with only cosmetic
#     changes, by hppd_manuscript_master/03_visualization/supplement/
#     hgf_bms_modified.py.
################################################################################

import os
import sys
import glob
import subprocess
import numpy as np
import pandas as pd
from scipy.special import digamma, gammaln, logsumexp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

################################################################################
# ─── FLAGS ────────────────────────────────────────────────────────────────────
################################################################################

DO_RUN  = True    # run bms_vch.jl (skips subjects already done)
DO_BMS  = True    # run RFX-BMS on collected LME CSVs
DO_PLOT = True    # generate Ef + PXP bar charts

################################################################################
# ─── ACTIVE CONFIG ────────────────────────────────────────────────────────────
# MODEL_KEYS: which models to include in the BMS comparison.
#   Any subset of model types that have been computed by bms_vch.jl.
#   Supported strings (must have corresponding lme_{key} column in LME CSVs):
#     "2level_empiric"             — 2-level HGF, current empirical conditions
#     "2level_nominal"             — 2-level HGF, raw condition proportions (0/0.25/0.5/0.75)
#     "3level_empiric"             — 3-level HGF, current empirical conditions
#     "3level_nominal"             — 3-level HGF, raw condition proportions
################################################################################

MODEL_KEYS = ["2level_empiric", "3level_empiric", "2level_nominal", "3level_nominal"]

# Display labels (edit to taste)
MODEL_LABELS = {
    "2level_empiric"             : "2-level\n(empirical)",
    "2level_nominal"             : "2-level\n(nominal)",
    "3level_empiric"             : "3-level\n(empirical)",
    "3level_nominal"             : "3-level\n(nominal)",
}

# Bar colours — cycles if more models than colours
_PALETTE = ["#4472C4", "#ED7D31", "#70AD47", "#FFC000",
            "#5B9BD5", "#C55A11", "#375623"]

################################################################################
# ─── PATHS ────────────────────────────────────────────────────────────────────
################################################################################

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
BMS_DIR     = os.path.join(SCRIPT_DIR, "model_comparison", "bms")
RESULTS_DIR = os.path.join(BMS_DIR, "results")
FIGURES_DIR = os.path.join(SCRIPT_DIR, "model_comparison", "figures")
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

################################################################################
# ─── RFX-BMS IMPLEMENTATION ───────────────────────────────────────────────────
################################################################################

def rfx_bms(log_evidence, alpha0=1.0, max_iter=1000, tol=1e-6,
            n_xp_samples=1_000_000, seed=42):
    """
    Random-Effects Bayesian Model Selection (Stephan et al. 2009) with
    Protected Exceedance Probabilities (Rigoux et al. 2014).

    Parameters
    ----------
    log_evidence : ndarray, shape (n_subjects, n_models)
        Laplace log marginal evidence per subject × model.
    alpha0 : float
        Dirichlet prior concentration parameter (1 = uniform).
    max_iter, tol : VB convergence settings.
    n_xp_samples  : Monte-Carlo samples for XP estimation.
    seed           : RNG seed for reproducibility.

    Returns
    -------
    alpha : ndarray (n_models,) — posterior Dirichlet parameters
    Ef    : ndarray (n_models,) — expected model frequencies
    xp    : ndarray (n_models,) — exceedance probabilities
    pxp   : ndarray (n_models,) — protected exceedance probabilities
    bor   : float               — Bayesian Omnibus Risk
    F     : float               — variational free energy
    """
    n, k = log_evidence.shape
    alpha0_vec = np.ones(k) * alpha0
    alpha = alpha0_vec + n / k

    u = np.zeros((n, k))
    for _ in range(max_iter):
        alpha_old = alpha.copy()
        log_u  = log_evidence + (digamma(alpha) - digamma(alpha.sum()))
        log_u -= log_u.max(axis=1, keepdims=True)
        u      = np.exp(log_u)
        u     /= u.sum(axis=1, keepdims=True)
        alpha = alpha0_vec + u.sum(axis=0)
        if np.max(np.abs(alpha - alpha_old)) < tol:
            break

    Ef = alpha / alpha.sum()

    rng     = np.random.default_rng(seed)
    samples = rng.dirichlet(alpha, size=n_xp_samples)
    winners = samples.argmax(axis=1)
    xp      = np.array([(winners == m).mean() for m in range(k)])

    F_data = (u * log_evidence).sum()
    kl = (gammaln(alpha.sum())    - gammaln(alpha0_vec.sum())
          - gammaln(alpha).sum()  + gammaln(alpha0_vec).sum()
          + ((alpha - alpha0_vec) * (digamma(alpha) - digamma(alpha.sum()))).sum())
    F = F_data - kl

    F0      = (logsumexp(log_evidence, axis=1) - np.log(k)).sum()
    log_bor = F0 - np.logaddexp(F0, F)
    bor     = float(np.clip(np.exp(log_bor), 0.0, 1.0))
    pxp     = (1 - bor) * xp + bor / k

    return alpha, Ef, xp, pxp, bor, F

################################################################################
# ─── DO_RUN ───────────────────────────────────────────────────────────────────
################################################################################

if DO_RUN:
    cmd = ["julia", "--project=.", "--threads", "1", "bms_vch.jl"]
    print(f"\nRunning: {' '.join(cmd)}")
    subprocess.run(cmd, cwd=SCRIPT_DIR, check=True)

################################################################################
# ─── DO_BMS ───────────────────────────────────────────────────────────────────
################################################################################

if DO_BMS:
    lme_files = sorted(glob.glob(os.path.join(RESULTS_DIR, "*_lme.csv")))
    if not lme_files:
        print(f"\n[BMS] No LME files found in {RESULTS_DIR}. Run DO_RUN first.")
    else:
        dfs    = [pd.read_csv(f) for f in lme_files]
        lme_df = pd.concat(dfs, ignore_index=True)

        # Validate that all requested model keys have LME columns
        missing_cols = [k for k in MODEL_KEYS
                        if f"lme_{k}" not in lme_df.columns]
        if missing_cols:
            raise ValueError(
                f"LME columns missing for models: {missing_cols}. "
                f"Re-run bms_vch.jl with these models in MODEL_TYPES."
            )

        # Include only participants who converged for ALL requested models
        conv_cols  = [f"converged_{k}" for k in MODEL_KEYS]
        lme_cols   = [f"lme_{k}"       for k in MODEL_KEYS]
        valid      = lme_df[conv_cols].all(axis=1)
        lme_valid  = lme_df[valid].copy()
        n_excluded = (~valid).sum()

        print(f"\n[BMS] {valid.sum()}/{len(lme_df)} participants with valid LME "
              f"for all {len(MODEL_KEYS)} models ({n_excluded} excluded)")
        if n_excluded > 0:
            print(f"  Excluded IDs: {lme_df[~valid]['record_id'].tolist()}")

        log_ev = lme_valid[lme_cols].values
        alpha, Ef, xp, pxp, bor, F = rfx_bms(log_ev)

        print(f"\n{'─'*60}")
        print(f"  RFX-BMS results  (n = {len(lme_valid)} subjects, "
              f"{len(MODEL_KEYS)} models)")
        print(f"  BOR (Bayesian Omnibus Risk) = {bor:.4f}")
        print(f"  {'Model':<22s}  {'α':>7s}  {'Ef':>7s}  {'XP':>7s}  {'PXP':>7s}")
        for i, key in enumerate(MODEL_KEYS):
            print(f"  {key:<22s}  {alpha[i]:>7.2f}  {Ef[i]:>7.3f}  "
                  f"{xp[i]:>7.3f}  {pxp[i]:>7.3f}")
        print(f"{'─'*60}\n")

        n_models = len(MODEL_KEYS)
        summary = pd.DataFrame({
            "model"      : MODEL_KEYS,
            "alpha"      : alpha,
            "Ef"         : Ef,
            "XP"         : xp,
            "PXP"        : pxp,
            "BOR"        : [bor]      * n_models,
            "F_VB"       : [F]        * n_models,
            "n_subjects" : [len(lme_valid)] * n_models,
        })
        summary_path = os.path.join(BMS_DIR, "bms_summary.csv")
        summary.to_csv(summary_path, index=False)
        print(f"Saved: {summary_path}")

        lme_valid.to_csv(os.path.join(BMS_DIR, "lme_per_subject.csv"), index=False)
        print(f"Saved: {os.path.join(BMS_DIR, 'lme_per_subject.csv')}")

################################################################################
# ─── DO_PLOT ──────────────────────────────────────────────────────────────────
################################################################################

if DO_PLOT:
    matplotlib.rcParams.update({
        "font.family":     "Arial",
        "axes.labelsize":  20,
        "xtick.labelsize": 13,
        "ytick.labelsize": 14,
    })

    summary_path = os.path.join(BMS_DIR, "bms_summary.csv")
    if not os.path.exists(summary_path):
        print("\n[PLOT] No BMS summary found. Run DO_BMS first.")
    else:
        summary = pd.read_csv(summary_path)
        n_models = len(summary)
        bor      = summary["BOR"].iloc[0]
        n_sub    = summary["n_subjects"].iloc[0]
        x        = np.arange(n_models)
        colors   = [_PALETTE[i % len(_PALETTE)] for i in range(n_models)]
        labels   = [MODEL_LABELS.get(k, k) for k in summary["model"]]
        bar_w    = max(0.3, min(0.6, 2.4 / n_models))
        fig_w    = max(9, 2.5 * n_models)

        fig, axes = plt.subplots(1, 2, figsize=(fig_w, 5))

        # ── Ef ────────────────────────────────────────────────────────────────
        axes[0].bar(x, summary["Ef"], color=colors, width=bar_w, edgecolor="none")
        axes[0].axhline(1 / n_models, color="gray", lw=1.2,
                        linestyle="--", alpha=0.7)
        axes[0].set_xticks(x)
        axes[0].set_xticklabels(labels)
        axes[0].set_ylabel("Ef")
        axes[0].set_ylim(0, 1)
        for xi, ef in zip(x, summary["Ef"]):
            axes[0].text(xi, ef + 0.02, f"{ef:.3f}", ha="center", va="bottom",
                         fontsize=11, fontweight="bold")
        sns.despine(ax=axes[0])

        # ── PXP ───────────────────────────────────────────────────────────────
        axes[1].bar(x, summary["PXP"], color=colors, width=bar_w, edgecolor="none")
        axes[1].axhline(1 / n_models, color="gray", lw=1.2,
                        linestyle="--", alpha=0.7)
        axes[1].set_xticks(x)
        axes[1].set_xticklabels(labels)
        axes[1].set_ylabel("PXP")
        axes[1].set_ylim(0, 1)
        for xi, pxp_v in zip(x, summary["PXP"]):
            axes[1].text(xi, pxp_v + 0.02, f"{pxp_v:.3f}", ha="center",
                         va="bottom", fontsize=11, fontweight="bold")
        sns.despine(ax=axes[1])

        fig.suptitle(f"RFX-BMS  |  n = {n_sub}  |  BOR = {bor:.3f}",
                     fontsize=14, y=1.02)
        fig.tight_layout()
        out = os.path.join(FIGURES_DIR, "bms_ef_pxp.png")
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"\n[PLOT] Saved: {out}")
