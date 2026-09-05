#!/usr/bin/env /usr/local/bin/python3.12
"""
mediation_results_table.py
===========================
Publication-style table summarising all mediation analysis results for the
HPPD manuscript (nice_covariates_spusers covariate set).

Panels
------
  a — hppd_binary   (PPA History;      SP predictor = psychedelic_age / spage)
  b — caps_vision   (CAPS Vision score; SP predictor = avg_life_dose  / avgdose)

Mediators (4 per panel):  vchbeta, vchnu, vchrate, vchthreshold

Column structure (same for both panels):
  Path | β | β 94% HDI | P(β≠0) | Δ | Δ 94% HDI | P(Δ≠0)

Rows per mediator section:
  hppd_binary : a, b, c′, NIE, NDE, Total, PMed
  caps_vision : a, b, c′, b (hu), c′ (hu), NIE, NDE, Total, PMed

Data sources (per model directory results/{dv}/mediation_models/{dv}_{spvar}_{med}_{cov}/)
  path_coefficients_summary.csv  → β, hdi_lower_94, hdi_upper_94, prob  (paths a/b/c_prime/b_hu/c_hu)
  path_counterfactual_summary.csv → Δ, hdi_low, hdi_high, p  (effects "A path"/"B path"/"C' path")
  mc_mediation_summary.csv       → Δ, hdi_low, hdi_high, p  (effects NIE/NDE/TE/PMed)

Two posterior summaries in one Δ column
---------------------------------------
The a/b/c′ rows report the posterior MEAN of the counterfactual contrast
(POINT_ESTIMATE_COL in modules/master_config.py). The NIE/NDE/Total/PMed rows
report the posterior MEDIAN (MC_EFFECT_POINT_ESTIMATE_COL): those four are
produced by Monte-Carlo integration over the mediator's posterior predictive,
whose heavy tails make the mean unusable — in the caps_vision × vch_beta model
5 draws out of 16,000 supplied ~100% of the mean, putting it at 6.4e7 against a
median of 0.12. The 94% HDI and the direction probabilities are quantile-based
and identical under either summary, which is why only the point estimate differs.
The column header names both ("Δ/Δ_med"), and the CSV carries a "Δ summary"
column recording which one each row used.

Styling mirrors regression_results_table.py / generate_publication_table.py:
  Arial font, bold-italic section headers, TB structural lines, amber P-cells ≥ 0.90.

Output
------
  results/supplement/tables/supplementary_table_s5.png
  results/supplement/tables/supplementary_table_s5.csv

Usage
-----
  cd 04_visualizations
  /usr/local/bin/python3.12 supplement/mediation_results_table.py
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

from master_config import point_estimate


# ── CONFIG ─────────────────────────────────────────────────────────────────────

RESULTS_BASE = Path("../../results")
OUTPUT_DIR   = RESULTS_BASE / "supplement" / "tables"
OUTPUT_FILE  = OUTPUT_DIR / "supplementary_table_s5.png"
OUTPUT_CSV   = OUTPUT_DIR / "supplementary_table_s5.csv"
OUTPUT_DOCX  = OUTPUT_DIR / "supplementary_table_s5.docx"

COV_SET = "nice_covariates_spusers"

PROB_THRESHOLD = 0.90
OUTPUT_DPI     = 400
FONT_FAMILY    = "Arial"
C_WHITE        = "#FFFFFF"
C_HIGHLIGHT    = "#FFF3CD"   # amber: P ≥ PROB_THRESHOLD

# ── OUTCOME / MEDIATOR METADATA ────────────────────────────────────────────────

PANEL_CONFIGS = [
    {
        "dv":        "hppd_binary",
        "spvar":     "spage",
        "dv_label":  "PPA History (hppd_binary)",
        "has_hu":    False,
        "panel_id":  "a",
    },
    {
        "dv":        "caps_vision",
        "spvar":     "avgdose",
        "dv_label":  "CAPS Vision Score (caps_vision)",
        "has_hu":    True,
        "panel_id":  "b",
    },
]

MEDIATORS = ["vchbeta", "vchnu", "vchrate", "vchthreshold"]

MEDIATOR_LABELS = {
    "vchbeta":      "VCH Decision Noise (\u03b2)",       # β
    "vchnu":        "VCH Top-down Bias (\u03bd)",         # ν
    "vchrate":      "VCH Response Rate",
    "vchthreshold": "VCH 75% Detection Threshold",
}

# Column header strings
BETA_COLS  = ["\u03b2", "\u03b2 94% HDI", "P(\u03b2\u22600)"]   # β, ≠
# The Δ point-estimate column carries TWO different posterior summaries, so its
# header names both: the a/b/c′ rows report the posterior mean (POINT_ESTIMATE_COL)
# while the NIE/NDE/Total/PMed rows report the posterior median
# (MC_EFFECT_POINT_ESTIMATE_COL) — see modules/master_config.py for why.
# Only this column is marked: the 94% HDI and the direction probabilities are
# quantile-based and are identical under either summary, so they need no caveat.
# "$_{med}$" is matplotlib mathtext for the PNG; DELTA_SUBSCRIPT_MAP renders the
# same string as a real Word subscript in the .docx.
DELTA_EST_HEADER = "\u0394/\u0394$_{med}$"
DELTA_SUBSCRIPT_MAP = {
    DELTA_EST_HEADER: [
        {"text": "\u0394/\u0394"},
        {"text": "med", "subscript": True},
    ],
}
DELTA_COLS = [DELTA_EST_HEADER, "\u0394 94% HDI", "P(\u0394\u22600)"]    # Δ, ≠


# ══════════════════════════════════════════════════════════════════════════════
# DATA LOADING
# ══════════════════════════════════════════════════════════════════════════════

def model_dir(dv: str, spvar: str, mediator: str) -> Path:
    return RESULTS_BASE / dv / "mediation_models" / f"{dv}_{spvar}_{mediator}_{COV_SET}"


def load_model_data(dv: str, spvar: str, mediator: str) -> dict[str, pd.DataFrame]:
    """Load the three result CSVs for one mediation model."""
    d = model_dir(dv, spvar, mediator)
    return {
        "path_coef": pd.read_csv(d / "path_coefficients_summary.csv"),
        "path_cf":   pd.read_csv(d / "path_counterfactual_summary.csv"),
        "mc":        pd.read_csv(d / "mc_mediation_summary.csv"),
    }


_NAN_RESULT = dict(estimate=np.nan, hdi_lo=np.nan, hdi_hi=np.nan, prob=np.nan)


def get_path_coef(df: pd.DataFrame, path: str) -> dict:
    """Extract one path's β, 94% HDI, and max-direction P from path_coefficients_summary."""
    row = df[df["path"] == path]
    if len(row) == 0:
        return _NAN_RESULT.copy()
    r = row.iloc[0]
    # Reported point estimate is the posterior mean (POINT_ESTIMATE_COL); see
    # modules/master_config.py. Raises rather than falling back to 'median' if
    # this CSV predates the mean columns, in which case the model needs refitting.
    return dict(
        estimate = point_estimate(r, source="path_coefficients_summary.csv"),
        hdi_lo = r["hdi_lower_94"],
        hdi_hi = r["hdi_upper_94"],
        prob   = max(r["prob_above_0"], r["prob_below_0"]),
    )


def get_cf_effect(df: pd.DataFrame, effect_prefix: str, *,
                  mc_integrated: bool = False) -> dict:
    """Extract one effect's Δ, HDI, and max-direction P from counterfactual/MC CSV.

    Matches by str.startswith() because the 'effect' column carries verbose labels
    (e.g. "NIE  Indirect (via vch_beta)").
    """
    row = df[df["effect"].str.startswith(effect_prefix)]
    if len(row) == 0:
        return _NAN_RESULT.copy()
    r = row.iloc[0]
    # Which posterior summary is reported depends on where the row came from,
    # and is decided by the caller via mc_integrated (see modules/master_config.py):
    #   path_counterfactual_summary.csv (A/B/C' paths)  -> POINT_ESTIMATE_COL ('mean')
    #   mc_mediation_summary.csv (NIE/NDE/TE/PMed)      -> MC_EFFECT_POINT_ESTIMATE_COL ('median')
    # The MC effects are integrated over the mediator's posterior predictive, whose
    # heavy tails make the mean unusable (a handful of draws can supply ~100% of it);
    # the median, the 94% HDI and the direction probabilities are all quantile-based
    # and stable, which is why only the point estimate differs between the two.
    return dict(
        estimate = point_estimate(r,
                                  source="path_counterfactual_summary.csv / mc_mediation_summary.csv",
                                  mc_integrated=mc_integrated),
        hdi_lo = r["hdi_low"],
        hdi_hi = r["hdi_high"],
        prob   = max(r["p_above_0"], r["p_below_0"]),
    )


# ══════════════════════════════════════════════════════════════════════════════
# FORMATTING HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def fmt_num(v: float) -> str:
    if pd.isna(v):
        return "\u2014"   # em dash
    return f"{v:+.3f}"


def fmt_prob(v: float) -> str:
    if pd.isna(v):
        return "\u2014"
    if v >= 0.999:
        return ">.999"
    return f"{v:.3f}"


def fmt_hdi(lo: float, hi: float) -> str:
    if pd.isna(lo) or pd.isna(hi):
        return "\u2014"
    return f"[{lo:+.3f}, {hi:+.3f}]"


# ══════════════════════════════════════════════════════════════════════════════
# ROW BUILDERS
# ══════════════════════════════════════════════════════════════════════════════

def _make_data_row(label: str, beta: dict, delta: dict) -> dict:
    """Build one display row given pre-fetched β and Δ result dicts."""
    cells = [
        fmt_num(beta["estimate"]),
        fmt_hdi(beta["hdi_lo"], beta["hdi_hi"]),
        fmt_prob(beta["prob"]),
        fmt_num(delta["estimate"]),
        fmt_hdi(delta["hdi_lo"], delta["hdi_hi"]),
        fmt_prob(delta["prob"]),
    ]
    # probs list aligns 1-to-1 with the 6 data cells (None = no highlighting)
    probs = [None, None, beta["prob"], None, None, delta["prob"]]
    return {"type": "data", "label": label, "cells": cells, "probs": probs}


def build_rows(dv: str, spvar: str, has_hu: bool) -> list[dict]:
    """Build all display rows for one panel (one DV)."""
    rows: list[dict] = []
    _nan = _NAN_RESULT.copy()

    for mediator in MEDIATORS:
        rows.append({"type": "section_header", "label": MEDIATOR_LABELS[mediator]})

        data = load_model_data(dv, spvar, mediator)
        pc   = data["path_coef"]
        cf   = data["path_cf"]
        mc   = data["mc"]

        # ── path coefficients — hu rows immediately follow their cognate path ──
        rows.append(_make_data_row("a",  get_path_coef(pc, "a"),       get_cf_effect(cf, "A path")))
        rows.append(_make_data_row("b",  get_path_coef(pc, "b"),       get_cf_effect(cf, "B path")))
        if has_hu:
            rows.append(_make_data_row("b  (hu)", get_path_coef(pc, "b_hu"), _nan))
        rows.append(_make_data_row("c\u2032", get_path_coef(pc, "c_prime"), get_cf_effect(cf, "C' path")))
        if has_hu:
            rows.append(_make_data_row("c\u2032 (hu)", get_path_coef(pc, "c_hu"), _nan))

        # ── MC mediation effects (Δ only; no normalised β) ───────────────────
        mc_specs = [
            ("NIE",    "NIE"),
            ("NDE",    "NDE"),
            ("Total",  "TE"),
            ("PMed",   "PMed"),   # proportion mediated (0–1 scale)
        ]
        for label, mc_prefix in mc_specs:
            rows.append(_make_data_row(label, _nan,
                                       get_cf_effect(mc, mc_prefix, mc_integrated=True)))

    return rows


# ══════════════════════════════════════════════════════════════════════════════
# TABLE CONVERSION
# ══════════════════════════════════════════════════════════════════════════════

def rows_to_display(rows: list[dict]) -> tuple[pd.DataFrame, pd.DataFrame, list[list]]:
    """Convert list-of-row-dicts → (display_df, row_meta, prob_rows)."""
    col_headers = ["Path"] + BETA_COLS + DELTA_COLS   # 7 columns

    n_data = len(col_headers) - 1   # 6 data cells
    display_rows: list[list] = []
    meta_rows:    list[dict] = []
    prob_rows:    list[list] = []

    for row in rows:
        if row["type"] == "section_header":
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
# PANEL RENDERER  (mirrors regression_results_table.py)
# ══════════════════════════════════════════════════════════════════════════════

def draw_panel(
    ax: plt.Axes,
    display_df: pd.DataFrame,
    row_meta: pd.DataFrame,
    prob_rows: list[list],
    panel_label: str,
) -> None:
    """Render one panel on *ax* in the publication style from generate_publication_table.py."""
    ax.axis("off")
    plt.rcParams["font.family"] = FONT_FAMILY

    n_rows = len(display_df)
    n_cols = len(display_df.columns)   # includes Path predictor column

    # ── panel label (bold letter, left of table) ───────────────────────────────
    ax.text(
        -0.005, 1.0, panel_label,
        transform=ax.transAxes,
        fontsize=13, fontweight="bold", fontfamily=FONT_FAMILY,
        ha="right", va="top",
    )

    # ── build table ────────────────────────────────────────────────────────────
    table = ax.table(
        cellText  = display_df.values,
        colLabels = display_df.columns.tolist(),
        loc       = "upper left",
        cellLoc   = "center",
        colLoc    = "center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9.0)

    # ── column widths ──────────────────────────────────────────────────────────
    # "Path" column must fit section-header labels (longest: "VCH 75% Detection
    # Threshold", ~29 chars) as well as short path labels ("a", "b  (hu)", etc.).
    # HDI cells need room for "[+0.000, +0.000]" (≈16 chars).
    pred_frac = 0.22
    data_frac = (1.0 - pred_frac) / (n_cols - 1)

    # ── row-height allocation ──────────────────────────────────────────────────
    HEADER_W  = 1.4   # column header row
    SECTION_W = 1.25  # bold-italic section header rows
    DATA_W    = 1.0   # ordinary data rows
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

        # ─── column header row (row_i == 0) ───
        if row_i == 0:
            cell.visible_edges = "TBR" if col_i == 0 else "TB"
            cell.set_linewidth(1.5)
            cell.set_height(header_h)
            if col_i == 0:
                cell.set_text_props(weight="bold", ha="left",   va="center")
            else:
                cell.set_text_props(weight="bold", ha="center", va="center")
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

        else:  # data row
            cell.set_height(data_h)
            if col_i == 0:
                # Thin right divider always; thick bottom only on last row
                cell.visible_edges = "BR" if is_last else "R"
                cell.set_linewidth(1.5 if is_last else 0.7)
                cell.set_text_props(ha="left", va="center")
            else:
                cell.visible_edges = "B" if is_last else ""
                cell.set_linewidth(1.5 if is_last else 0.4)
                # Amber highlight for notable probability cells
                p_idx = col_i - 1   # 0-based index into probs list
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

def build_csv_data() -> pd.DataFrame:
    """Build a structured CSV with raw numeric values for all mediation results.

    Column names use the same Unicode symbols as the PNG table headers
    (β, Δ, ≠) so the CSV is publication-ready for Scientific Reports.
    """
    records: list[dict] = []
    for cfg in PANEL_CONFIGS:
        dv, spvar, has_hu = cfg["dv"], cfg["spvar"], cfg["has_hu"]
        for mediator in MEDIATORS:
            data = load_model_data(dv, spvar, mediator)
            pc = data["path_coef"]
            cf = data["path_cf"]
            mc = data["mc"]

            # Path coefficients (a, b, c') — each has both beta and delta
            path_specs = [
                ("a",  "a",       "A path"),
                ("b",  "b",       "B path"),
                ("c\u2032", "c_prime", "C' path"),
            ]
            for label, path_key, cf_prefix in path_specs:
                beta  = get_path_coef(pc, path_key)
                delta = get_cf_effect(cf, cf_prefix)
                records.append({
                    "Panel": cfg["panel_id"], "DV": dv, "SP Predictor": spvar,
                    "Mediator": MEDIATOR_LABELS[mediator], "Path": label,
                    "\u03b2": beta["estimate"],
                    "\u03b2 94% HDI Lower": beta["hdi_lo"],
                    "\u03b2 94% HDI Upper": beta["hdi_hi"],
                    "P(\u03b2\u22600)": beta["prob"],
                    "\u0394": delta["estimate"],
                    "\u0394 summary": "mean",
                    "\u0394 94% HDI Lower": delta["hdi_lo"],
                    "\u0394 94% HDI Upper": delta["hdi_hi"],
                    "P(\u0394\u22600)": delta["prob"],
                })

            # Hurdle paths (caps_vision only) — beta only, no delta
            if has_hu:
                for label, path_key in [("b (hu)", "b_hu"),
                                        ("c\u2032 (hu)", "c_hu")]:
                    beta = get_path_coef(pc, path_key)
                    records.append({
                        "Panel": cfg["panel_id"], "DV": dv, "SP Predictor": spvar,
                        "Mediator": MEDIATOR_LABELS[mediator], "Path": label,
                        "\u03b2": beta["estimate"],
                        "\u03b2 94% HDI Lower": beta["hdi_lo"],
                        "\u03b2 94% HDI Upper": beta["hdi_hi"],
                        "P(\u03b2\u22600)": beta["prob"],
                        "\u0394": np.nan,
                        "\u0394 summary": "",
                        "\u0394 94% HDI Lower": np.nan,
                        "\u0394 94% HDI Upper": np.nan,
                        "P(\u0394\u22600)": np.nan,
                    })

            # MC mediation effects (NIE/NDE/Total/PMed) — delta only, no beta
            for label, mc_prefix in [("NIE", "NIE"), ("NDE", "NDE"),
                                     ("Total", "TE"), ("PMed", "PMed")]:
                delta = get_cf_effect(mc, mc_prefix, mc_integrated=True)
                records.append({
                    "Panel": cfg["panel_id"], "DV": dv, "SP Predictor": spvar,
                    "Mediator": MEDIATOR_LABELS[mediator], "Path": label,
                    "\u03b2": np.nan,
                    "\u03b2 94% HDI Lower": np.nan,
                    "\u03b2 94% HDI Upper": np.nan,
                    "P(\u03b2\u22600)": np.nan,
                    "\u0394": delta["estimate"],
                    # median, not mean — see MC_EFFECT_POINT_ESTIMATE_COL
                    "\u0394 summary": "median",
                    "\u0394 94% HDI Lower": delta["hdi_lo"],
                    "\u0394 94% HDI Upper": delta["hdi_hi"],
                    "P(\u0394\u22600)": delta["prob"],
                })
    return pd.DataFrame(records)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    n_med = len(MEDIATORS)   # 4

    # Rows per mediator section in each panel
    rows_per_med_a = 3 + 4   # a/b/c'  +  NIE/NDE/Total/PMed
    rows_per_med_b = 5 + 4   # a/b/c'/b_hu/c'_hu  +  NIE/NDE/Total/PMed

    # Approximate panel heights (inches): col-header + section headers + data rows
    COL_H = 0.30
    SEC_H = 0.26
    DAT_H = 0.21
    GAP_H = 0.40

    panel_a_h = COL_H + n_med * SEC_H + n_med * rows_per_med_a * DAT_H
    panel_b_h = COL_H + n_med * SEC_H + n_med * rows_per_med_b * DAT_H
    total_h   = panel_a_h + panel_b_h + GAP_H

    # Width: wider pred column (0.22) + 6 data cols at 1.35 in each
    n_data_cols = len(BETA_COLS) + len(DELTA_COLS)   # 6
    fig_w = 3.2 + n_data_cols * 1.35   # ≈11.3 in

    print(f"Figure: {fig_w:.1f} × {total_h:.1f} in")

    # ── save CSV (UTF-8 BOM for Excel/Google Sheets compatibility) ───────────
    csv_df = build_csv_data()
    csv_df.to_csv(OUTPUT_CSV, index=False, encoding="utf-8-sig")
    print(f"Saved CSV → {OUTPUT_CSV}")

    # ── build display tables ───────────────────────────────────────────────────
    panels: list[dict] = []
    for cfg in PANEL_CONFIGS:
        print(f"  Building rows for {cfg['dv']} …")
        rows = build_rows(cfg["dv"], cfg["spvar"], cfg["has_hu"])
        df, meta, probs = rows_to_display(rows)
        panels.append({
            "df":     df,
            "meta":   meta,
            "probs":  probs,
            "label":  cfg["panel_id"],
            "height": panel_a_h if not cfg["has_hu"] else panel_b_h,
        })

    # ── save Word table (primary submission artifact) ─────────────────────────
    from _docx_helper import save_docx_tables
    save_docx_tables(
        panels=[
            {"display_df": p["df"], "row_meta": p["meta"],
             "prob_rows": p["probs"], "panel_label": p["label"]}
            for p in panels
        ],
        output_path=OUTPUT_DOCX,
        prob_threshold=PROB_THRESHOLD,
        subscript_map=DELTA_SUBSCRIPT_MAP,
    )
    print(f"Saved DOCX → {OUTPUT_DOCX}")

    # ── layout figure (panel a on top, panel b on bottom) ─────────────────────
    fig = plt.figure(figsize=(fig_w, total_h))

    panel_frac_a = panels[0]["height"] / total_h
    panel_frac_b = panels[1]["height"] / total_h
    gap_frac     = GAP_H / total_h

    # Panel a: top; bottom edge = panel_frac_b + gap_frac
    ax_a = fig.add_axes([0.03, panel_frac_b + gap_frac, 0.97, panel_frac_a])
    draw_panel(ax_a, panels[0]["df"], panels[0]["meta"], panels[0]["probs"],
               panel_label=panels[0]["label"])

    # Panel b: bottom
    ax_b = fig.add_axes([0.03, 0.0, 0.97, panel_frac_b])
    draw_panel(ax_b, panels[1]["df"], panels[1]["meta"], panels[1]["probs"],
               panel_label=panels[1]["label"])

    print(f"Saving → {OUTPUT_FILE}")
    fig.savefig(OUTPUT_FILE, dpi=OUTPUT_DPI, bbox_inches="tight")
    plt.close(fig)
    print("Done.")


if __name__ == "__main__":
    main()
