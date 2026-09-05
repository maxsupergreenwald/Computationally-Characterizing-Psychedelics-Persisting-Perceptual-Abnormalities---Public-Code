#!/usr/bin/env /usr/local/bin/python3.12
"""
fdr_correction.py
=================
Benjamini-Hochberg FDR correction of every frequentist test reported in the
manuscript for the three canonical predictor blocks (`sp_predictors`,
`vch_behavior`, `vch_computations`).

Two independent test families, one per outcome:

  Family 1 — caps_vision   : 9 Spearman correlations
  Family 2 — hppd_binary   : 9 Mann-Whitney U tests

Each family is corrected two ways:

  p_FDR+ (m = 9)   across all nine frequentist tests for that outcome.
  p_FDR- (m = 3)   within one iv_type_dict block only, run separately for each
                   of the three blocks.

BH is run with statsmodels `multipletests(method="fdr_bh")`, which returns the
step-up adjusted p-values (monotonicity enforced), so a value can be compared
directly against alpha = 0.05.

Note that p_FDR- is NOT uniformly smaller than p_FDR+. BH adjusts by m/i --
family size over rank -- not by m, so a smaller family does not imply a smaller
adjusted p; that intuition holds for Bonferroni and fails for a step-up
procedure. `fdr_correction_verify.py` documents the two mechanisms and verifies
both columns against an independent reimplementation and against R's p.adjust.

-------------------------------------------------------------------------------
Which p-value is used for each test
-------------------------------------------------------------------------------
caps_vision x sp_predictors uses the AGE-CONTROLLED partial Spearman
(`correlation_grid_age_control.csv`) because that is the panel shown in the
manuscript figure (Fig. 1b) and the value quoted by
`05_results_narrative/results_narrative.py`. The zero-order
`correlation_grid.csv` is NOT used for that block. The two vch blocks have no
age-controlled variant; they use `correlation_grid.csv`.

Reads
  results/caps_vision/sp_predictors/data_visualization/summary_results/
      correlation_grid_age_control.csv
  results/caps_vision/{vch_behavior,vch_computations}/data_visualization/
      summary_results/correlation_grid.csv
  results/hppd_binary/{ivtype}/data_visualization/summary_results/
      boxplot_grid.csv                                  (test == "mann_whitney")

Writes
  results/supplement/fdr_correction/fdr_correction.csv

Imported by
  mann_whitney_table_hppd_binary.py   (p_FDR+ / p_FDR- columns)
  spearman_table_caps_vision.py       (p_FDR+ / p_FDR- columns)
  fdr_correction_verify.py            (verification harness)

Those scripts call `family_frame(...)` rather than re-deriving anything, so the
family definitions, the age-controlled-grid mapping, and the BH settings live in
exactly one place. Importing this module runs `os.chdir(_HERE)`; the callers sit
in this same directory and chdir there themselves, so the effect is a no-op.

Usage
  cd 04_visualizations
  /usr/local/bin/python3.12 supplement/fdr_correction.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
os.chdir(_HERE)
sys.path.insert(0, str(_HERE.parent.parent / "modules"))

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from master_config import iv_type_dict, iv_type_group_labels, dv_to_lab_short

# ── CONFIG ───────────────────────────────────────────────────────────────────

RESULTS_BASE = Path("../../results")
OUTPUT_DIR   = RESULTS_BASE / "supplement" / "fdr_correction"
OUTPUT_CSV   = OUTPUT_DIR / "fdr_correction.csv"

# Predictor blocks, in manuscript order. Nine predictors total per outcome.
IVTYPES = ["sp_predictors", "vch_behavior", "vch_computations"]

# sp_predictors is the only block with an age-controlled correlation grid, and
# it is the one the manuscript figure shows. See module docstring.
CAPS_GRID_BY_IVTYPE = {
    "sp_predictors":    "correlation_grid_age_control.csv",
    "vch_behavior":     "correlation_grid.csv",
    "vch_computations": "correlation_grid.csv",
}

ALPHA = 0.05

FDR_METHOD = "fdr_bh"   # Benjamini-Hochberg

# Column names the supplement tables render as p_FDR+ / p_FDR-. Imported by
# mann_whitney_table_hppd_binary.py and spearman_table_caps_vision.py so the
# two tables can never drift onto different passes.
FDR_PLUS_COL  = "p_fdr_m9"   # BH across all 9 frequentist tests for the outcome
FDR_MINUS_COL = "p_fdr_m3"   # BH within the predictor's own iv_type_dict block


# ══════════════════════════════════════════════════════════════════════════════
# LOADERS
# ══════════════════════════════════════════════════════════════════════════════

def _require(path: Path) -> Path:
    """Fail loudly rather than silently substituting a different results file."""
    if not path.exists():
        raise FileNotFoundError(
            f"Expected results file not found: {path}\n"
            "Re-run 0X_all_figures.py (RUN_HPPD_CAPS_FIGS=True) to regenerate the "
            "data_visualization summary CSVs, or re-pull the single-path compile."
        )
    return path


def _one_row(df: pd.DataFrame, mask: pd.Series, context: str) -> pd.Series:
    """Return the single matching row, or raise. Never picks the first of many."""
    hits = df[mask]
    if len(hits) != 1:
        raise ValueError(f"Expected exactly 1 row for {context}; found {len(hits)}.")
    return hits.iloc[0]


def load_caps_spearman() -> pd.DataFrame:
    """9 Spearman rows for caps_vision, one per predictor."""
    out = []
    for ivtype in IVTYPES:
        path = _require(
            RESULTS_BASE / "caps_vision" / ivtype / "data_visualization"
            / "summary_results" / CAPS_GRID_BY_IVTYPE[ivtype]
        )
        grid = pd.read_csv(path)
        for var in iv_type_dict[ivtype]:
            r = _one_row(
                grid,
                (grid["row_var"] == var) & (grid["column_var"] == "caps_vision"),
                f"Spearman {var} x caps_vision in {path.name}",
            )
            out.append({
                "family":          "caps_vision",
                "test_type":       "spearman",
                "test_detail":     ("partial Spearman | age_v2"
                                    if ivtype == "sp_predictors" else "Spearman"),
                "source_file":     str(path.relative_to(RESULTS_BASE)),
                "ivtype":          ivtype,
                "ivtype_label":    iv_type_group_labels[ivtype],
                "predictor":       var,
                "predictor_label": dv_to_lab_short.get(var, var),
                "n":               r["n"],
                "effect_name":     "rho",
                "effect":          r["rho"],
                "p_raw":           r["p_value"],
                "p_tails":         "two-tailed",
            })
    return pd.DataFrame(out)


def load_hppd_mwu() -> pd.DataFrame:
    """9 Mann-Whitney U rows for hppd_binary, one per predictor."""
    out = []
    for ivtype in IVTYPES:
        path = _require(
            RESULTS_BASE / "hppd_binary" / ivtype / "data_visualization"
            / "summary_results" / "boxplot_grid.csv"
        )
        grid = pd.read_csv(path)
        for var in iv_type_dict[ivtype]:
            r = _one_row(
                grid,
                (grid["dv"] == var) & (grid["test"] == "mann_whitney"),
                f"Mann-Whitney {var} x hppd_binary in {path.name}",
            )
            out.append({
                "family":          "hppd_binary",
                "test_type":       "mann_whitney",
                "test_detail":     "Mann-Whitney U",
                "source_file":     str(path.relative_to(RESULTS_BASE)),
                "ivtype":          ivtype,
                "ivtype_label":    iv_type_group_labels[ivtype],
                "predictor":       var,
                "predictor_label": dv_to_lab_short.get(var, var),
                "n":               r["sample_size"],
                "effect_name":     "U",
                "effect":          r["statistic"],
                "p_raw":           r["p_value"],
                "p_tails":         "two-tailed",
            })
    return pd.DataFrame(out)


# ══════════════════════════════════════════════════════════════════════════════
# BH CORRECTION
# ══════════════════════════════════════════════════════════════════════════════

def _bh(pvals: pd.Series, suffix: str) -> pd.DataFrame:
    """Run BH on `pvals` and return rank / adjusted-p / significance columns."""
    reject, p_adj, _, _ = multipletests(pvals.to_numpy(), alpha=ALPHA, method=FDR_METHOD)
    return pd.DataFrame(
        {
            f"m_{suffix}":       len(pvals),
            f"rank_{suffix}":    pvals.rank(method="min").astype(int),
            f"p_fdr_{suffix}":   p_adj,
            f"sig_fdr_{suffix}": reject,
        },
        index=pvals.index,
    )


def _bh_within_block(freq: pd.DataFrame) -> pd.DataFrame:
    """BH applied separately inside each iv_type_dict block (m = 3 per block).

    Each block is treated as its own family, so `m_m3` is 3 and the ranks
    restart at 1 in every block.
    """
    parts = [
        _bh(freq.loc[freq["ivtype"] == ivtype, "p_raw"], "m3")
        for ivtype in IVTYPES
    ]
    return pd.concat(parts).reindex(freq.index)


def correct(freq: pd.DataFrame) -> pd.DataFrame:
    """Apply both passes to one family's nine tests and return the frame."""
    pass_plus  = _bh(freq["p_raw"], "m9")
    pass_minus = _bh_within_block(freq)
    return freq.join(pass_plus).join(pass_minus)


def family_frame(family: str) -> pd.DataFrame:
    """Fully corrected nine-row frame for one outcome family.

    This is the entry point the supplement table scripts import.
    """
    if family == "caps_vision":
        return correct(load_caps_spearman())
    if family == "hppd_binary":
        return correct(load_hppd_mwu())
    raise ValueError(
        f"Unknown family {family!r}; expected 'caps_vision' or 'hppd_binary'."
    )


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main() -> None:
    out = pd.concat(
        [family_frame("caps_vision"), family_frame("hppd_binary")],
        ignore_index=True,
    )

    out["sig_raw"] = out["p_raw"] < ALPHA

    col_order = [
        "family", "test_type", "test_detail", "ivtype", "ivtype_label",
        "predictor", "predictor_label", "n",
        "effect_name", "effect",
        "p_raw", "p_tails", "sig_raw",
        "m_m9", "rank_m9", "p_fdr_m9", "sig_fdr_m9",
        "m_m3", "rank_m3", "p_fdr_m3", "sig_fdr_m3",
        "source_file",
    ]
    missing = [c for c in col_order if c not in out.columns]
    if missing:
        raise KeyError(f"Column(s) missing from assembled frame: {missing}")
    out = out[col_order]

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    # utf-8-sig so Excel / Google Sheets read the Greek and arrow glyphs.
    out.to_csv(OUTPUT_CSV, index=False, encoding="utf-8-sig")

    # ── console summary ──────────────────────────────────────────────────────
    pd.set_option("display.width", 200)
    for fam in ["caps_vision", "hppd_binary"]:
        sub = out[out["family"] == fam]
        print(f"\n{'=' * 100}\nFAMILY: {fam}\n{'=' * 100}")
        show = sub[[
            "test_type", "predictor", "n", "effect_name", "effect",
            "p_raw", "p_fdr_m9", "sig_fdr_m9", "p_fdr_m3", "sig_fdr_m3",
        ]].copy()
        for c in ["effect", "p_raw", "p_fdr_m9", "p_fdr_m3"]:
            show[c] = show[c].map(lambda v: "" if pd.isna(v) else f"{v:.4f}")
        print(show.to_string(index=False))

        print(f"\n  tests surviving alpha={ALPHA}: "
              f"raw {int(sub['sig_raw'].sum())}/9 | "
              f"p_FDR+ (m=9) {int(sub['sig_fdr_m9'].sum())}/9 | "
              f"p_FDR- (m=3) {int(sub['sig_fdr_m3'].sum())}/9")

    print(f"\nWrote {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
