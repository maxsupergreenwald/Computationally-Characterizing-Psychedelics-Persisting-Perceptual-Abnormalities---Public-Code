"""
Spearman correlations between caps_vision and potential confounding variables.

This script is a supplementary confound check that complements the group-comparison
tables (clinical_table_hppd_split_caps_split.png, demographics_table_hppd_split_caps_split.png).
Rather than dichotomising caps_vision into CAPS(-)/CAPS(+), each variable is correlated
with continuous caps_vision using Spearman's rho, preserving power.

Sample
------
SP users with valid caps_vision (psycheduse_yn == "Yes" AND caps_vision is not null),
replicating the right-side group of the combined HPPD+CAPS split tables (n ≈ 130).
Same QC filter as the main tables: raven_total > 0 AND qc_passed > 0.

Variables included
------------------
All variables from DEMOGRAPHICS_SECTIONS and CLINICAL_SECTIONS in 0X_all_figures.py,
plus four additional 6-month substance-use count variables:
  mj_6month, atypicals_6month, stimulants_6month, sedatives_6month

Multi-level unordered categoricals (sex_v2, race_v2, location_summary) are recoded
as one binary dummy per level. Levels absent from the analysis sample (n_positive == 0)
are skipped. Levels with n_positive ≤ 5 are included but flagged with ‡.

Reads
-----
  data/final/df_public_*.csv  (the shipped analysis dataframe; most recent wins)

Outputs (all in results/supplement/caps_vision_confounds/)
-----------------------------------------------------------
  supplementary_table_s1.png   — publication-quality matplotlib table
  supplementary_table_s1.csv   — raw results (for downstream use)

Run from any directory:
  /usr/local/bin/python3.12 04_visualizations/supplement/caps_vision_confounds_spearman.py
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

# ── Path setup ─────────────────────────────────────────────────────────────────
SCRIPT_DIR  = Path(__file__).resolve().parent
REPO_ROOT   = SCRIPT_DIR.parent.parent
MODULES_DIR = REPO_ROOT / "modules"
sys.path.insert(0, str(MODULES_DIR))

from data_prep import most_recent_public_df


OUT_DIR = REPO_ROOT / "results" / "supplement" / "tables"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Data loading ──────────────────────────────────────────────────────────────
# This repository ships the analysis dataframe already prepared — every derived
# column is present in the CSV, so there is no preparation step to run.
# most_recent_public_df() picks the newest data/final/df_public_*.csv and prints
# which one it chose.
df = pd.read_csv(most_recent_public_df(REPO_ROOT / "data" / "final"), low_memory=False)

# ── Sample selection — identical to 0X_all_figures.py ─────────────────────────
# df_wraven: all QC-passed participants with valid Raven scores
df_wraven = df[(df["raven_total"] > 0) & (df["qc_passed"] > 0)].copy()

# Derive any-substance indicators — replicates 0X_all_figures.py lines 737–753.
# NaN substance-use flags are treated as 0 (not reported = not used), consistent
# with the existing pipeline that generates the clinical tables.
_month_candidates    = ["alc_month_yn", "ghb_month_yn", "opioids_month_yn",
                        "sedatives_month_yn", "mj_month_yn", "atypicals_month_yn",
                        "stimulants_month_yn"]
_lifetime_candidates = ["alc_lifetime", "ghb_lifetime", "opioids_lifetime",
                        "sedatives_life_yn", "mj_lifetime", "atypicals_life_yn",
                        "stimulants_life_yn"]
_month_cols    = [c for c in _month_candidates    if c in df_wraven.columns]
_lifetime_cols = [c for c in _lifetime_candidates if c in df_wraven.columns]
if not _month_cols:
    raise KeyError("No valid past-month substance columns found for any_substance_month_yn")
if not _lifetime_cols:
    raise KeyError("No valid lifetime substance columns found for any_substance_lifetime_yn")
df_wraven["any_substance_month_yn"] = (
    df_wraven[_month_cols].apply(pd.to_numeric, errors="coerce").fillna(0).gt(0).any(axis=1)
).astype(int)
df_wraven["any_substance_lifetime_yn"] = (
    df_wraven[_lifetime_cols].apply(pd.to_numeric, errors="coerce").fillna(0).gt(0).any(axis=1)
).astype(int)

# Location summary — same coding as tables_and_descriptives.py / 0X_all_figures.py.
df_wraven.loc[df_wraven["us_loc_v2"] == 1, "outside_us_v2"] = 0
df_wraven["location_summary"] = "Other"
df_wraven.loc[df_wraven["outside_us_v2"].isin([0, 32]),
              "location_summary"] = "United States"
df_wraven.loc[df_wraven["outside_us_v2"].isin(
              [64, 82, 154, 159, 127, 136, 17, 45, 9, 74, 60, 137, 139, 180]),
              "location_summary"] = "Europe"
df_wraven.loc[df_wraven["outside_us_v2"].isin([122, 9]),
              "location_summary"] = "Oceania"
df_wraven.loc[df_wraven["outside_us_v2"].isin([134, 110]),
              "location_summary"] = "Latin America"

# Analysis sample: SP users with valid caps_vision
df_spusers = df_wraven[df_wraven["psycheduse_yn"] == "Yes"].copy()
df_caps    = df_spusers[df_spusers["caps_vision"].notna()].copy()
N          = len(df_caps)
print(f"Analysis sample: n = {N} (SP users with valid caps_vision)")

# ── Binary dummies for multi-level categoricals ────────────────────────────────
# Each level becomes a binary 0/1 column.  Levels with n_positive == 0 in the
# analysis sample (df_caps) yield constant columns — these are detected later
# and skipped during the Spearman computation.
_sex_spec = {
    "sex_binary_male":          (df_caps["sex_v2"] == 1.0, "Sex: Male"),
    "sex_binary_female":        (df_caps["sex_v2"] == 2.0, "Sex: Female"),
    "sex_binary_other":         (df_caps["sex_v2"] == 3.0, "Sex: Other gender"),
}
_race_spec = {
    "race_binary_aian":         (df_caps["race_v2"] == 1.0, "Race: Am. Indian/Alaska Native"),
    "race_binary_asian":        (df_caps["race_v2"] == 2.0, "Race: Asian"),
    "race_binary_nhpi":         (df_caps["race_v2"] == 3.0, "Race: Native Hawaiian/Pacific Islander"),
    "race_binary_black":        (df_caps["race_v2"] == 4.0, "Race: Black/African American"),
    "race_binary_white":        (df_caps["race_v2"] == 5.0, "Race: White"),
    "race_binary_multiracial":  (df_caps["race_v2"] == 6.0, "Race: Multiracial"),
    "race_binary_unknown":      (df_caps["race_v2"] == 7.0, "Race: Unknown/Prefer not to say"),
    "race_binary_latino":       (df_caps["race_v2"] == 8.0, "Race: Latino/a"),
}
_loc_spec = {
    "loc_binary_us":            (df_caps["location_summary"] == "United States", "Location: United States"),
    "loc_binary_eur":           (df_caps["location_summary"] == "Europe",        "Location: Europe"),
    "loc_binary_oce":           (df_caps["location_summary"] == "Oceania",       "Location: Oceania"),
    "loc_binary_latam":         (df_caps["location_summary"] == "Latin America", "Location: Latin America"),
    "loc_binary_other":         (df_caps["location_summary"] == "Other",         "Location: Other"),
}
for _spec in [_sex_spec, _race_spec, _loc_spec]:
    for _col, (_mask, _) in _spec.items():
        df_caps[_col] = _mask.astype(int)

# ── Balanced education variable (canonical version now in data_prep.py) ───────
# Collapses only the two sparse tails: levels 2+3 → 1, levels 8+9 → 6.
# All resulting cells have n ≥ 11 in this sample.
# Level map: 1=HS or less, 2=Some college, 3=2-yr degree,
#            4=College student (enrolled), 5=Bachelor's, 6=Master's or higher.
_ed_balanced_map = {2: 1, 3: 1, 4: 2, 5: 3, 6: 4, 7: 5, 8: 6, 9: 6}
df_caps['highest_education_balanced'] = df_caps['highest_education'].map(_ed_balanced_map)

# ── Variable specification ─────────────────────────────────────────────────────
# Each entry: (section_label, row_label, column_in_df_caps)
# Ordered to mirror the layout of the original publication tables.
VARIABLE_SPEC = [
    # ── Demographic Characteristics ───────────────────────────────────────────
    ("Demographic Characteristics", "Age",                                   "age_v2"),
    ("Demographic Characteristics", "Sex: Male",                             "sex_binary_male"),
    ("Demographic Characteristics", "Sex: Female",                           "sex_binary_female"),
    ("Demographic Characteristics", "Sex: Other gender",                     "sex_binary_other"),
    ("Demographic Characteristics", "Race: White",                           "race_binary_white"),
    ("Demographic Characteristics", "Race: Latino/a",                        "race_binary_latino"),
    ("Demographic Characteristics", "Race: Asian",                           "race_binary_asian"),
    ("Demographic Characteristics", "Race: Multiracial",                     "race_binary_multiracial"),
    ("Demographic Characteristics", "Race: Am. Indian/Alaska Native",        "race_binary_aian"),
    ("Demographic Characteristics", "Race: Native Hawaiian/Pacific Islander","race_binary_nhpi"),
    ("Demographic Characteristics", "Race: Black/African American",          "race_binary_black"),
    ("Demographic Characteristics", "Race: Unknown/Prefer not to say",       "race_binary_unknown"),
    ("Demographic Characteristics", "Location: United States",               "loc_binary_us"),
    ("Demographic Characteristics", "Location: Europe",                      "loc_binary_eur"),
    ("Demographic Characteristics", "Location: Oceania",                     "loc_binary_oce"),
    ("Demographic Characteristics", "Location: Latin America",               "loc_binary_latam"),
    ("Demographic Characteristics", "Location: Other",                       "loc_binary_other"),
    ("Demographic Characteristics", "Education (9-level)",                   "highest_education"),
    ("Demographic Characteristics", "Education (3-level)",                   "highest_education_3level_final"),
    ("Demographic Characteristics", "Education (6-level)",                   "highest_education_balanced"),
    ("Demographic Characteristics", "Raven's Progressive Matrices (9-item)", "raven_total"),
    # ── Psychiatric History ───────────────────────────────────────────────────
    ("Psychiatric History",         "Any Diagnosis",                         "mental_illness2_v2"),
    ("Psychiatric History",         "Psychotic Spectrum",                    "psych_spectrum_v2"),
    ("Psychiatric History",         "Schizophrenia Spectrum",                "schizophrenia_spectrum"),
    ("Psychiatric History",         "Unipolar Mood Disorder",                "mood_disorder"),
    ("Psychiatric History",         "Bipolar Disorder (no psychosis)",       "mental_illness_dx_current_9"),
    ("Psychiatric History",         "Substance Use Disorder",                "addiction"),
    ("Psychiatric History",         "Anxiety Disorder",                      "anxiety_disorder"),
    ("Psychiatric History",         "OCD",                                   "mental_illness_dx_current_30"),
    ("Psychiatric History",         "Trauma-Related Disorder (PTSD/c-PTSD)","trauma_disorder"),
    ("Psychiatric History",         "Eating Disorder",                       "eating_disorder"),
    ("Psychiatric History",         "Personality Disorder",                  "personality_disorder"),
    ("Psychiatric History",         "Autism Spectrum Disorder",              "mental_illness_dx_current_6"),
    ("Psychiatric History",         "ADHD",                                  "mental_illness_dx_current_5"),
    ("Psychiatric History",         "Sleep Disorder",                        "sleep_disorder"),
    # ── Psychiatric Medications ───────────────────────────────────────────────
    ("Psychiatric Medications",     "Any Medication",                        "medication_current_v2"),
    ("Psychiatric Medications",     "Antipsychotic",                         "antipsychotic"),
    ("Psychiatric Medications",     "Antidepressant",                        "antidepressants"),
    ("Psychiatric Medications",     "Stimulant Medication",                  "simulant_medication"),
    ("Psychiatric Medications",     "Benzodiazepine",                        "benzos"),
    ("Psychiatric Medications",     "Anxiolytic (Non-Benzodiazepine)",       "nonbenzo_anxiolytics"),
    ("Psychiatric Medications",     "Sedative (Non-Benzodiazepine)",         "nonbenzo_sedatives"),
    ("Psychiatric Medications",     "Opioid Antagonist",                     "opioid_antagonists"),
    # ── Other Substance Use — Past Month ─────────────────────────────────────
    ("Other Substance Use (Past Month)", "Any Substance",                    "any_substance_month_yn"),
    ("Other Substance Use (Past Month)", "Alcohol",                          "alc_month_yn"),
    ("Other Substance Use (Past Month)", "Sedative-Hypnotic",                "ghb_month_yn"),
    ("Other Substance Use (Past Month)", "Opioid",                           "opioids_month_yn"),
    ("Other Substance Use (Past Month)", "Cannabis (binary)",                "mj_month_yn"),
    ("Other Substance Use (Past Month)", "Atypical Psychedelics (binary)",   "atypicals_month_yn"),
    ("Other Substance Use (Past Month)", "Stimulants (binary)",              "stimulants_month_yn"),
    ("Other Substance Use (Past Month)", "Cannabis (count, past 6 months)",  "mj_6month"),
    ("Other Substance Use (Past Month)", "Atypical Psychedelics (count, past 6 months)", "atypicals_6month"),
    ("Other Substance Use (Past Month)", "Stimulants (count, past 6 months)","stimulants_6month"),
    ("Other Substance Use (Past Month)", "Sedative-Hypnotic (count, past 6 months)", "sedatives_6month"),
    # ── Other Substance Use — Lifetime ────────────────────────────────────────
    ("Other Substance Use (Lifetime)", "Any Substance",                      "any_substance_lifetime_yn"),
    ("Other Substance Use (Lifetime)", "Alcohol",                            "alc_lifetime"),
    ("Other Substance Use (Lifetime)", "Sedative-Hypnotic",                  "ghb_lifetime"),
    ("Other Substance Use (Lifetime)", "Opioid",                             "opioids_lifetime"),
    ("Other Substance Use (Lifetime)", "Cannabis",                           "mj_lifetime"),
    ("Other Substance Use (Lifetime)", "Atypical Psychedelics",              "atypicals_life_yn"),
    ("Other Substance Use (Lifetime)", "Stimulants",                         "stimulants_life_yn"),
]

# ── Compute Spearman correlations ──────────────────────────────────────────────
# Rows with NaN in either variable are dropped by scipy automatically.
# Constant columns (e.g. absent race categories with n_positive = 0) yield
# undefined rho; these are detected and skipped before calling spearmanr.

SPARSE_THRESHOLD = 5   # n_positive ≤ this → include but mark with ‡

dv = df_caps["caps_vision"]
results = []

for section, label, col in VARIABLE_SPEC:
    if col not in df_caps.columns:
        print(f"  WARNING: '{col}' not in df_caps — skipping '{label}'")
        continue

    x = df_caps[col]

    # Count positives for binary dummies; for continuous use n_valid
    n_pos = int(x.sum()) if set(x.dropna().unique()).issubset({0, 1, 0.0, 1.0}) else None

    # Skip all-zero (constant) columns — spearmanr would return NaN with a warning
    if x.nunique(dropna=True) <= 1:
        print(f"  Skipping '{label}' ({col}): constant column"
              f" (n_positive = {int(x.sum())})")
        continue

    # Valid pairs for the correlation
    valid = x.notna() & dv.notna()
    n_valid = int(valid.sum())

    rho, p = spearmanr(x[valid], dv[valid])

    display_label = label

    results.append({
        "section": section,
        "label":   display_label,
        "col":     col,
        "n":       n_valid,
        "n_pos":   n_pos,   # None for continuous variables
        "rho":     float(rho),
        "p":       float(p),
    })

results_df = pd.DataFrame(results)
print(f"\nComputed {len(results_df)} Spearman correlations.")

# ── Save CSV ───────────────────────────────────────────────────────────────────
csv_path = OUT_DIR / "supplementary_table_s1.csv"
results_df.to_csv(csv_path, index=False)
print(f"Saved: {csv_path}")

# ── Formatting helpers ─────────────────────────────────────────────────────────

def _fmt_rho(rho: float) -> str:
    """Format Spearman rho: signed, 3 decimal places (e.g. '+0.123')."""
    return f"{rho:+.3f}"


def _fmt_p(p: float) -> str:
    """Format p-value with leading-zero stripped for values < 1."""
    if p < 0.001:
        return "< .001"
    s = f"{p:.3f}"        # e.g. "0.045"
    return s[1:] if s.startswith("0.") else s   # ".045"


# ── Build rows for matplotlib table ───────────────────────────────────────────
# Each entry in cell_text is [label_str, rho_str, p_str].
# Section header rows have rho_str = "" and p_str = "".
# row_types is a parallel list tracking "section_header" vs "data".

INDENT = "    "   # 4-space indent for variable rows

cell_text  = []
row_types  = []

current_section = None
for _, row in results_df.iterrows():
    if row["section"] != current_section:
        current_section = row["section"]
        cell_text.append([row["section"], "", ""])
        row_types.append("section_header")
    cell_text.append([INDENT + row["label"], _fmt_rho(row["rho"]), _fmt_p(row["p"])])
    row_types.append("data")

col_labels = ["Variable", "ρ", "p"]

# ── Render PNG — matching _save_png_table_thickonly style ─────────────────────
plt.rcParams["font.family"] = "Arial"

n_rows = len(cell_text)
n_cols = len(col_labels)

FIG_W = 8.5
FIG_H = max(3.8, 0.5 + (n_rows * 0.16))

fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
ax.axis("off")

tbl = ax.table(
    cellText=cell_text,
    colLabels=col_labels,
    loc="upper left",
    cellLoc="left",
    colLoc="center",
)
tbl.auto_set_font_size(False)
fontsize = 10.5 if n_rows <= 35 else 9.5
tbl.set_fontsize(fontsize)

# Column widths — first column wider to accommodate long variable names
FIRST_COL_W = 0.62
OTHER_COL_W = (1.0 - FIRST_COL_W) / (n_cols - 1)   # 0.19 each

# Row height weights — section headers get a little extra space
row_weights = [1.28 if rt == "section_header" else 1.0 for rt in row_types]
header_weight = 1.35
table_fill = 0.95
per_weight = table_fill / (header_weight + sum(row_weights))
header_height = header_weight * per_weight

for (row_i, col_i), cell in tbl.get_celld().items():
    cell.set_facecolor("white")
    cell.set_edgecolor("#222222")
    cell.set_width(FIRST_COL_W if col_i == 0 else OTHER_COL_W)

    # Column-header row
    if row_i == 0:
        cell.visible_edges = "TB"
        cell.set_text_props(weight="bold", ha="center", va="center")
        cell.set_linewidth(1.5)
        cell.set_height(header_height)
        continue

    row_type = row_types[row_i - 1]
    cell.set_height(row_weights[row_i - 1] * per_weight)

    if col_i == 0:
        cell.set_text_props(ha="left")
    else:
        cell.set_text_props(ha="center")

    if row_type == "section_header":
        cell.visible_edges = "TB"
        if col_i == 0:
            cell.set_text_props(style="italic", weight="bold", ha="left")
        cell.set_linewidth(1.0)
    else:
        cell.visible_edges = ""

fig.subplots_adjust(left=0.02, right=0.98, top=0.95, bottom=0.02)

png_path = OUT_DIR / "supplementary_table_s1.png"
fig.savefig(png_path, dpi=400, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {png_path}")

# ── Save DOCX via shared _docx_helper ─────────────────────────────────────────
sys.path.insert(0, str(SCRIPT_DIR))
from _docx_helper import save_docx_tables

# Build display_df and row_meta for save_docx_tables()
docx_rows = []
docx_meta = []

current_section = None
for _, row in results_df.iterrows():
    if row["section"] != current_section:
        current_section = row["section"]
        docx_rows.append({"Variable": row["section"], "ρ": "", "p": ""})
        docx_meta.append({"row_type": "section_header"})
    docx_rows.append({
        "Variable": INDENT + row["label"],
        "ρ": _fmt_rho(row["rho"]),
        "p": _fmt_p(row["p"]),
    })
    docx_meta.append({"row_type": "data"})

display_df = pd.DataFrame(docx_rows)
row_meta   = pd.DataFrame(docx_meta)

docx_path = OUT_DIR / "supplementary_table_s1.docx"
save_docx_tables(
    panels=[{
        "display_df":  display_df,
        "row_meta":    row_meta,
        "prob_rows":   None,
        "panel_label": None,
    }],
    output_path=docx_path,
)
print(f"Saved: {docx_path}")

# ── Console summary ────────────────────────────────────────────────────────────
sig_rows = results_df[results_df["p"] < 0.05].copy()
print(f"\n{'='*60}")
print(f"Significant (p < .05) correlations with caps_vision (n = {N}):")
print(f"{'='*60}")
if sig_rows.empty:
    print("  None")
else:
    for _, r in sig_rows.iterrows():
        print(f"  {r['label']:<50s}  ρ = {r['rho']:+.3f}  p = {_fmt_p(r['p'])}")
print()
