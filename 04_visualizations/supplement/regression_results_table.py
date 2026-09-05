#!/usr/bin/env /usr/local/bin/python3.12
"""
regression_results_table.py
================================
Same Bayesian regression data as the original v1 script, rendered in the
publication-table style produced by generate_publication_table.py:

  - ax.table() renderer with Arial (matches clinical_table_*.png style)
  - Bold-italic section-header rows with top+bottom structural lines only
  - No per-row borders on data rows; thick bottom line closes each panel
  - Amber highlights on probability-of-effect cells where P ≥ 0.90
  - No super-header rows — column names encode the section structure
  - No footnote

Panel a — hppd_binary (PPA History):
  Predictor | β | β 94% HDI | P(β≠0) | Δ | Δ 94% HDI | P(Δ≠0) | N

Panel b — caps_vision (adds hurdle component):
  Predictor | β | β 94% HDI | P(β≠0) | hu β | hu β 94% HDI | P(hu_β≠0) | Δ | Δ 94% HDI | P(Δ≠0) | N

Data sources
  β, HDI94, P  →  results/sensitivity_analyses_single_paths/existingresults_manuscript.csv
  Δ, Δ HDI, P  →  results/sensitivity_analyses_single_paths/existingresults_manuscript_counterfactual.csv
  Covariate set: nice_covariates_spusers

Output
  results/supplement/tables/supplementary_table_s4.png
  results/supplement/tables/supplementary_table_s4.csv

Usage
  cd 04_visualizations
  /usr/local/bin/python3.12 supplement/regression_results_table.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
os.chdir(_HERE)
sys.path.insert(0, str(_HERE.parent.parent / "modules"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from master_config import iv_type_dict, iv_type_group_labels, dv_to_lab_short
from master_config import point_estimate

# ── CONFIG ──────────────────────────────────────────────────────────────────────

RESULTS_BASE       = Path("../../results")
SENSITIVITY_CSV    = RESULTS_BASE / "sensitivity_analyses_single_paths" / "existingresults_manuscript.csv"
COUNTERFACTUAL_CSV = RESULTS_BASE / "sensitivity_analyses_single_paths" / "existingresults_manuscript_counterfactual.csv"
OUTPUT_DIR         = RESULTS_BASE / "supplement" / "tables"
OUTPUT_FILE        = OUTPUT_DIR / "supplementary_table_s4.png"
OUTPUT_CSV         = OUTPUT_DIR / "supplementary_table_s4.csv"
OUTPUT_DOCX        = OUTPUT_DIR / "supplementary_table_s4.docx"

IVTYPES       = ["sp_predictors", "vch_behavior", "vch_computations"]
CANONICAL_COV = "nice_covariates_spusers"
PROB_THRESHOLD = 0.90
OUTPUT_DPI     = 400

FONT_FAMILY = "Arial"
C_WHITE     = "#FFFFFF"
C_HIGHLIGHT = "#FFF3CD"   # amber: P ≥ PROB_THRESHOLD

# ── COLUMN HEADER STRINGS ───────────────────────────────────────────────────────

BETA_COLS  = ["\u03b2", "\u03b2 94% HDI", "P(\u03b2\u22600)"]
HU_COLS    = ["hu \u03b2", "hu \u03b2 94% HDI", "P(hu_\u03b2\u22600)"]
DELTA_COLS = ["\u0394", "\u0394 94% HDI", "P(\u0394\u22600)"]


# ══════════════════════════════════════════════════════════════════════════════
# DATA LOADING
# ══════════════════════════════════════════════════════════════════════════════

def load_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load and pre-filter both CSVs to the canonical covariate set."""
    sens = pd.read_csv(SENSITIVITY_CSV)
    cf   = pd.read_csv(COUNTERFACTUAL_CSV)
    sens = sens[sens["covariates"] == CANONICAL_COV].copy()
    cf   = cf[cf["cov"]           == CANONICAL_COV].copy()
    sens["prob"] = sens[["prob_above_0", "prob_below_0"]].max(axis=1)
    cf["prob"]   = cf[["prob_above_0",   "prob_below_0"]].max(axis=1)
    return sens, cf


def lookup_beta(sens: pd.DataFrame, var: str, dv: str) -> dict:
    row = sens[(sens["var"] == var) & (sens["dv"] == dv)]
    if len(row) == 0:
        return dict(estimate=np.nan, hdi_lo=np.nan, hdi_hi=np.nan, prob=np.nan, N=np.nan)
    r = row.iloc[0]
    return dict(estimate=r["Estimate"], hdi_lo=r["hdi_lower_94"],
                hdi_hi=r["hdi_upper_94"], prob=r["prob"], N=r["N"])


def lookup_delta(cf: pd.DataFrame, spvar: str, dv: str) -> dict:
    row = cf[(cf["spvar"] == spvar) & (cf["dv"] == dv)]
    if len(row) == 0:
        return dict(estimate=np.nan, hdi_lo=np.nan, hdi_hi=np.nan, prob=np.nan)
    r = row.iloc[0]
    # point_estimate() reads POINT_ESTIMATE_COL ('mean') — the same summary of
    # the posterior that lookup_beta()'s "Estimate" already is, so the two
    # halves of each table row are now comparable. Raises if the CSV predates
    # the mean columns rather than falling back to 'median'.
    return dict(estimate=point_estimate(r, source=str(COUNTERFACTUAL_CSV)),
                hdi_lo=r["hdi_lower_94"],
                hdi_hi=r["hdi_upper_94"], prob=r["prob"])


def fmt_num(v: float) -> str:
    if pd.isna(v):
        return "\u2014"
    return f"{v:+.3f}"


def fmt_prob(v: float) -> str:
    if pd.isna(v):
        return "\u2014"
    if v >= 0.999:
        return ">.999"
    return f"{v:.3f}"


def fmt_hdi(lo: float, hi: float) -> str:
    """Format HDI bounds as '[lo, hi]' with explicit signs."""
    if pd.isna(lo) or pd.isna(hi):
        return "\u2014"
    return f"[{lo:+.3f}, {hi:+.3f}]"


# ══════════════════════════════════════════════════════════════════════════════
# ROW BUILDERS
# ══════════════════════════════════════════════════════════════════════════════

def build_rows_panel_a(sens: pd.DataFrame, cf: pd.DataFrame) -> list[dict]:
    """Panel a — hppd_binary. 7 data cells: β | β HDI | P(β) | Δ | Δ HDI | P(Δ) | N."""
    DV = "hppd_binary"
    rows: list[dict] = []
    for ivtype in IVTYPES:
        rows.append({"type": "group_header", "label": iv_type_group_labels[ivtype]})
        for pred in iv_type_dict[ivtype]:
            pred_norm = pred + "_normalized"
            b = lookup_beta(sens, pred_norm, DV)
            d = lookup_delta(cf, pred_norm, DV)
            n_val = int(b["N"]) if not pd.isna(b["N"]) else None
            cells = [
                fmt_num(b["estimate"]), fmt_hdi(b["hdi_lo"], b["hdi_hi"]), fmt_prob(b["prob"]),
                fmt_num(d["estimate"]), fmt_hdi(d["hdi_lo"], d["hdi_hi"]), fmt_prob(d["prob"]),
                str(n_val) if n_val is not None else "\u2014",
            ]
            probs = [None, None, b["prob"], None, None, d["prob"], None]
            rows.append({"type": "data", "label": dv_to_lab_short.get(pred, pred),
                         "cells": cells, "probs": probs})
    return rows


def build_rows_panel_b(sens: pd.DataFrame, cf: pd.DataFrame) -> list[dict]:
    """Panel b — caps_vision. 10 data cells: β | β HDI | P(β) | hu β | hu HDI | hu P | Δ | Δ HDI | P(Δ) | N."""
    DV = "caps_vision"
    rows: list[dict] = []
    for ivtype in IVTYPES:
        rows.append({"type": "group_header", "label": iv_type_group_labels[ivtype]})
        for pred in iv_type_dict[ivtype]:
            pred_norm = pred + "_normalized"
            hu_var    = "hu_" + pred_norm
            b  = lookup_beta(sens, pred_norm, DV)
            hu = lookup_beta(sens, hu_var,    DV)
            d  = lookup_delta(cf, pred_norm,  DV)
            n_val = int(b["N"]) if not pd.isna(b["N"]) else None
            cells = [
                fmt_num(b["estimate"]),  fmt_hdi(b["hdi_lo"],  b["hdi_hi"]),  fmt_prob(b["prob"]),
                fmt_num(hu["estimate"]), fmt_hdi(hu["hdi_lo"], hu["hdi_hi"]), fmt_prob(hu["prob"]),
                fmt_num(d["estimate"]),  fmt_hdi(d["hdi_lo"],  d["hdi_hi"]),  fmt_prob(d["prob"]),
                str(n_val) if n_val is not None else "\u2014",
            ]
            probs = [None, None, b["prob"], None, None, hu["prob"], None, None, d["prob"], None]
            rows.append({"type": "data", "label": dv_to_lab_short.get(pred, pred),
                         "cells": cells, "probs": probs})
    return rows


# ══════════════════════════════════════════════════════════════════════════════
# TABLE HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def rows_to_display(
    rows: list[dict], include_hu: bool
) -> tuple[pd.DataFrame, pd.DataFrame, list[list]]:
    """Convert list-of-row-dicts → (display_df, row_meta, prob_rows)."""
    if include_hu:
        col_headers = ["Predictor"] + BETA_COLS + HU_COLS + DELTA_COLS + ["N"]
    else:
        col_headers = ["Predictor"] + BETA_COLS + DELTA_COLS + ["N"]

    n_data = len(col_headers) - 1
    display_rows: list[list] = []
    meta_rows:    list[dict] = []
    prob_rows:    list[list] = []

    for row in rows:
        if row["type"] == "group_header":
            display_rows.append([row["label"]] + [""] * n_data)
            meta_rows.append({"row_type": "section_header"})
            prob_rows.append([None] * n_data)
        else:
            display_rows.append([row["label"]] + row["cells"])
            meta_rows.append({"row_type": "data"})
            prob_rows.append(row["probs"])

    display_df = pd.DataFrame(display_rows, columns=col_headers)
    row_meta   = pd.DataFrame(meta_rows)
    return display_df, row_meta, prob_rows


# ══════════════════════════════════════════════════════════════════════════════
# PANEL RENDERER
# ══════════════════════════════════════════════════════════════════════════════

def draw_panel(
    ax: plt.Axes,
    display_df: pd.DataFrame,
    row_meta: pd.DataFrame,
    prob_rows: list[list],
    panel_label: str,
) -> None:
    """Render one panel table on *ax* in the reference publication style.

    Style mirrors _save_png_table_thickonly from generate_publication_table.py:
      - Header row: bold, visible_edges="TB", linewidth=1.5
      - Section-header rows: italic+bold, visible_edges="TB", linewidth=1.0
      - Data rows: no visible borders (visible_edges="")
      - Last data row: thick bottom border to close the table
      - Amber cell background when P ≥ PROB_THRESHOLD
    """
    ax.axis("off")
    plt.rcParams["font.family"] = FONT_FAMILY

    n_rows = len(display_df)
    n_cols = len(display_df.columns)   # includes Predictor

    # ── panel label ────────────────────────────────────────────────────────────
    ax.text(
        -0.005, 1.0, panel_label,
        transform=ax.transAxes,
        fontsize=13, fontweight="bold", fontfamily=FONT_FAMILY,
        ha="right", va="top",
    )

    # ── build table ────────────────────────────────────────────────────────────
    table = ax.table(
        cellText=display_df.values,
        colLabels=display_df.columns.tolist(),
        loc="upper left",
        cellLoc="center",
        colLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9.5)

    # ── column widths ─────────────────────────────────────────────────────────
    # Predictor column gets 22%; remaining data cols split the rest equally.
    pred_frac = 0.22
    data_frac = (1.0 - pred_frac) / (n_cols - 1)

    # ── row-height allocation ──────────────────────────────────────────────────
    HEADER_W  = 1.4
    SECTION_W = 1.3
    DATA_W    = 1.0
    total_w   = HEADER_W + sum(
        SECTION_W if row_meta.iloc[i]["row_type"] == "section_header" else DATA_W
        for i in range(n_rows)
    )
    per       = 0.96 / total_w
    header_h  = HEADER_W  * per
    section_h = SECTION_W * per
    data_h    = DATA_W    * per

    # ── style each cell ────────────────────────────────────────────────────────
    for (row_i, col_i), cell in table.get_celld().items():
        cell.set_facecolor(C_WHITE)
        cell.set_edgecolor("#222222")
        cell.set_width(pred_frac if col_i == 0 else data_frac)
        cell.PAD = 0.05

        # ─── column header row ───
        if row_i == 0:
            cell.visible_edges = "TBR" if col_i == 0 else "TB"
            cell.set_linewidth(1.5)
            cell.set_height(header_h)
            cell.set_text_props(weight="bold", ha="center", va="center")
            if col_i == 0:
                cell.set_text_props(weight="bold", ha="left", va="center")
            continue

        # ─── data / section rows ───
        row_type = row_meta.iloc[row_i - 1]["row_type"]
        is_last  = (row_i == n_rows)

        if row_type == "section_header":
            cell.set_height(section_h)
            cell.visible_edges = "TBR" if col_i == 0 else "TB"
            cell.set_linewidth(1.0)
            if col_i == 0:
                cell.set_text_props(style="italic", weight="bold", ha="left", va="center")
            else:
                cell.set_text_props(ha="center", va="center")

        else:
            cell.set_height(data_h)
            if col_i == 0:
                cell.visible_edges = "BR" if is_last else "R"
                cell.set_linewidth(1.5 if is_last else 0.7)
                cell.set_text_props(ha="left", va="center")
            else:
                cell.visible_edges = "B" if is_last else ""
                cell.set_linewidth(1.5 if is_last else 0.4)
                p_idx = col_i - 1
                probs = prob_rows[row_i - 1]
                prob  = probs[p_idx] if (probs is not None and p_idx < len(probs)) else None
                if prob is not None and not pd.isna(prob) and prob >= PROB_THRESHOLD:
                    cell.set_facecolor(C_HIGHLIGHT)
                    cell.set_text_props(weight="bold", ha="center", va="center")
                else:
                    cell.set_text_props(ha="center", va="center")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def build_csv_data(sens: pd.DataFrame, cf: pd.DataFrame) -> pd.DataFrame:
    """Build a structured CSV with raw numeric values for both panels.

    Column names use the same Unicode symbols as the PNG table headers
    (β, Δ, ≠) so the CSV is publication-ready for Scientific Reports.
    """
    records: list[dict] = []
    for ivtype in IVTYPES:
        for pred in iv_type_dict[ivtype]:
            pred_norm = pred + "_normalized"
            for dv in ["hppd_binary", "caps_vision"]:
                b = lookup_beta(sens, pred_norm, dv)
                d = lookup_delta(cf, pred_norm, dv)
                rec = {
                    "Panel": "a" if dv == "hppd_binary" else "b",
                    "DV": dv,
                    "Section": iv_type_group_labels[ivtype],
                    "Predictor": dv_to_lab_short.get(pred, pred),
                    "Predictor Variable": pred,
                    "\u03b2": b["estimate"],
                    "\u03b2 94% HDI Lower": b["hdi_lo"],
                    "\u03b2 94% HDI Upper": b["hdi_hi"],
                    "P(\u03b2\u22600)": b["prob"],
                    "\u0394": d["estimate"],
                    "\u0394 94% HDI Lower": d["hdi_lo"],
                    "\u0394 94% HDI Upper": d["hdi_hi"],
                    "P(\u0394\u22600)": d["prob"],
                    "N": int(b["N"]) if not pd.isna(b["N"]) else np.nan,
                }
                # Add hurdle component for caps_vision
                if dv == "caps_vision":
                    hu_var = "hu_" + pred_norm
                    hu = lookup_beta(sens, hu_var, dv)
                    rec["hu \u03b2"] = hu["estimate"]
                    rec["hu \u03b2 94% HDI Lower"] = hu["hdi_lo"]
                    rec["hu \u03b2 94% HDI Upper"] = hu["hdi_hi"]
                    rec["P(hu_\u03b2\u22600)"] = hu["prob"]
                records.append(rec)
    return pd.DataFrame(records)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading data …")
    sens, cf = load_data()

    print("Building rows …")
    rows_a = build_rows_panel_a(sens, cf)
    rows_b = build_rows_panel_b(sens, cf)

    df_a, meta_a, probs_a = rows_to_display(rows_a, include_hu=False)
    df_b, meta_b, probs_b = rows_to_display(rows_b, include_hu=True)

    # ── save CSV (UTF-8 BOM for Excel/Google Sheets compatibility) ───────────
    csv_df = build_csv_data(sens, cf)
    csv_df.to_csv(OUTPUT_CSV, index=False, encoding="utf-8-sig")
    print(f"Saved CSV → {OUTPUT_CSV}")

    # ── save Word table (primary submission artifact) ─────────────────────────
    from _docx_helper import save_docx_tables
    save_docx_tables(
        panels=[
            {"display_df": df_a, "row_meta": meta_a, "prob_rows": probs_a, "panel_label": "a"},
            {"display_df": df_b, "row_meta": meta_b, "prob_rows": probs_b, "panel_label": "b"},
        ],
        output_path=OUTPUT_DOCX,
        prob_threshold=PROB_THRESHOLD,
        landscape=True,
    )
    print(f"Saved DOCX → {OUTPUT_DOCX}")

    # ── figure sizing ──────────────────────────────────────────────────────────
    n_groups = len(IVTYPES)
    n_preds  = sum(len(iv_type_dict[ivt]) for ivt in IVTYPES)

    panel_h = 0.38 + n_groups * 0.30 + n_preds * 0.25
    GAP_H   = 0.45
    total_h = 2 * panel_h + GAP_H

    # Width driven by panel b: 3 (β) + 3 (hu) + 3 (Δ) + 1 (N) = 10 data cols.
    n_max_data = len(BETA_COLS) + len(HU_COLS) + len(DELTA_COLS) + 1   # 10
    fig_w = 2.8 + n_max_data * 1.35   # ≈ 16.3 in

    print(f"Figure: {fig_w:.1f} × {total_h:.1f} in")
    fig = plt.figure(figsize=(fig_w, total_h))

    panel_frac = panel_h / total_h
    gap_frac   = GAP_H   / total_h

    ax_a = fig.add_axes([0.03, panel_frac + gap_frac, 0.97, panel_frac])
    draw_panel(ax_a, df_a, meta_a, probs_a, panel_label="a")

    ax_b = fig.add_axes([0.03, 0.0, 0.97, panel_frac])
    draw_panel(ax_b, df_b, meta_b, probs_b, panel_label="b")

    print(f"Saving → {OUTPUT_FILE}")
    fig.savefig(OUTPUT_FILE, dpi=OUTPUT_DPI, bbox_inches="tight")
    plt.close(fig)
    print("Done.")


if __name__ == "__main__":
    main()
