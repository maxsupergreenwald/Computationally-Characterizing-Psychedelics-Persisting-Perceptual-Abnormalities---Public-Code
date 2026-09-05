#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Supplement: CONSORT-style Participant Flow Diagram

Renders the full participant flow from screening through to the analytic
subsamples, so that every N reported in the manuscript can be traced to an
explicit inclusion/exclusion step. Written in response to a reviewer request:

    "The participant counts fluctuate considerably across analyses (N = 228
     total, n = 186 SP users, n = 160 for VCH data, n = 130 for CAPS data,
     n = 113 for combined models). A comprehensive CONSORT-style participant
     flow diagram should be added to the Supplementary Materials to clarify
     exclusions at each stage."

All four of those subsample Ns are reproduced by this script directly from the
shipped analysis dataframe; none is hardcoded.


ONE STAGE PER PARTICIPANT
-------------------------
Every participant is counted at exactly one stage, so each arrow subtracts. That
is what makes this a flow rather than a set of independent tallies, and it is the
property the reviewer asked for. The chain closes with no residual:

    1249 screened
     -757 ineligible at screening                    -> Supp. ineligibility table
    = 492 eligible
     -195 failed to complete minimum study measures
              68 screening battery never scored
             118 no quality-control record; did not complete
               9 failed behavioural tasks or incomplete record
    = 297 completed minimum study measures
     - 69 failed post-completion quality control     -> Supp. QC failure table
    = 228 final analytic sample  (186 SP users, 42 SP-naive)

The 27 salvage records - participants who timed out but had completed at least
the full first questionnaire and passed all QC - are inside the final 228
throughout. They never leave the spine, so the count is annotated on the final
box rather than drawn as a re-entry arrow.

The script asserts this closure at runtime and raises if it ever stops holding
(e.g. after a data refresh), rather than silently rendering a broken figure.

The quality-control box carries only the eligible subset of the QC failures.
Fourteen participants were screened ineligible and subsequently also failed QC;
they leave the flow at the screening step. The si_2_v2 restriction applied below
is what keeps the box and the QC table reporting the same 69 people; this script
derives that population itself and renders those four counts as a table.

`results_narrative.py` also computes *marginal* counts over everyone in the
recruitment CSV rather than within a stage - RQ (347), RFQ (83), RTS (27). Each is
individually correct but they reach outside their own box and cannot be
differenced. RS (1249), RE (492) and RG (228) are identical under both definitions.


THE THREE "MINIMUM STUDY MEASURES" BUCKETS
------------------------------------------
Assigned by fixed precedence so the parts sum exactly to 195. These are exclusions
for incomplete or unscorable data, not for data quality.

Do not relabel any of them as a fraud exclusion. The 68 whose screening battery was
never scored are subject-pool participants: all have `student_yn == 1`, none has a
fraud-associated phone prefix, none is flagged ineligible, none has `qc_bad_data`
set, and those carrying `qc_notes` read "credit granted" / "credit was assigned" /
"couldn't assign her credit because the sona id was not found". Fraud exclusions in
this study occur at two *other* stages, both of which appear on this figure with
their own counts: screening (545 of 757) and post-completion QC (6 of 69).


Data prep is copied verbatim from `05_results_narrative/results_narrative.py`
so that the df_recruit population is identical. Do not diverge the two without
updating both.

Both sets of exclusion reasons are counted here, off the masks this script
already builds — no sibling script and no intermediate CSV. Keeping the
screening loop and the reason tallies in one file is what stops the diagram and
the counts it cites from drifting apart.

Reads
-----
  data/final/df_public_*.csv             - shipped analysis df
  RECRUIT_CSV (modules/master_config.py) - REDCap recruitment tracking CSV

Outputs (both in results/supplement/consort_diagram/)
------------------------------------------------------
  supplementary_figure_s9.png    - manuscript drafting / Google Docs preview
  supplementary_figure_s9.tiff   - journal submission (same DPI)
  supplementary_figure_s9.svg    - vector version for editing in Illustrator/Inkscape

Run from any directory:
  /usr/local/bin/python3.12 04_visualizations/supplement/consort_diagram.py
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

# Arial font — set before importing pyplot so rcParams take effect.
# Matches the convention used by every other figure script in this directory.
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

warnings.filterwarnings("ignore")

# ── Path setup ─────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
MODULES_DIR = REPO_ROOT / "modules"
sys.path.insert(0, str(MODULES_DIR))

from data_prep import most_recent_public_df
from master_config import RECRUIT_CSV
# QC-failure categories were pattern-matched from qc_notes (RA free text, since
# withheld) and frozen as record_id sets. See the module docstring.
from qc_redacted_categories import (
    QCCAT_INCONSISTENT, QCCAT_ATTN_CHECK, QCCAT_CHALLENGE,
    QCCAT_FRAUD, QCCAT_FRAUD_SIGNAL,
)
from collections import OrderedDict

# Recruitment CSV — the same file results_narrative.py reads. Update the
# RECRUIT_CSV constant in modules/master_config.py if the file moves.

OUT_DIR = REPO_ROOT / "results" / "supplement" / "consort_diagram"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Output resolution. 600 matches beta_sigmoid_creator.py; drop to 150 while
# iterating on layout, put it back to 600 before regenerating for submission.
DPI = 600

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: DATA PREP  (copied verbatim from results_narrative.py Section 2)
# ══════════════════════════════════════════════════════════════════════════════

_df_path = most_recent_public_df(REPO_ROOT / "data" / "final")
df = pd.read_csv(_df_path, low_memory=False)
print(f"Analysis df ({_df_path.name}): N = {len(df)}")

# Canonical SP-user filter — psycheduse_yn == "Yes", never != "No".
df_sp = df[df["psycheduse_yn"] == "Yes"].copy()

df_recruit_raw = pd.read_csv(RECRUIT_CSV, low_memory=False).copy()

df_recruit = df_recruit_raw.copy()
# Exclude test records (record_id < 203 AND not a student).
df_recruit = df_recruit[
    ~((df_recruit["record_id"] < 203) & (df_recruit["student_yn"] < 1))
].copy()
# Restrict to records within the study cutoff.
df_recruit = df_recruit[df_recruit["record_id"] <= df["record_id"].max()].copy()
# Fix for record 1858 — applied identically in results_narrative.py.
df_recruit.loc[df_recruit["record_id"] == 1858, ["salvage_yn", "qc_passed"]] = 0

df_recruit["raven_total"] = (
    df_recruit[[f"correct_answer{x}_v2" for x in [2, 3]]].sum(axis=1, min_count=1)
    + df_recruit["correct_answer_v2"]
    + df_recruit["raven_total_score_v2"]
)
df_recruit.loc[
    (df_recruit["raven_total"] < 1) & (df_recruit["student_yn"] > 0), "raven_total"
] = np.nan

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: STAGE MASKS
# ══════════════════════════════════════════════════════════════════════════════

# Membership in the shipped analysis df is the ground truth for "retained".
IN_DF = df_recruit["record_id"].isin(set(df["record_id"]))

# Completed screening. screening_survey_complete == 2 is REDCap's "Complete".
SCREENED = ~(df_recruit["screening_survey_complete"] < 2)

# Completed all study measures (results_narrative.py RQ).
STUDENTS_WHO_COMPLETED = (
    (df_recruit["student_yn"] == 1)
    & (df_recruit["qc_passed"] > 0)
    & (df_recruit["task_data_prltask_present"] == 1)
)
COMPLETED = df_recruit["honesty_qc"].notna() | STUDENTS_WHO_COMPLETED

# Student task-fail records (results_narrative.py RSTF_MASK): non-SP students who
# failed the behavioural tasks. Credit is awarded per institutional requirement but
# the data cannot be used. Not a data-quality failure.
RSTF = (
    (df_recruit["student_yn"] == 1)
    & (df_recruit["psycheduse_yn"] == 2)
    & (df_recruit["qc_passed"] < 1)
    & (df_recruit["raven_total_score_v2"] >= 1)
    & (df_recruit["raven_total_score_v2"].notna())
    & (df_recruit["qc_bad_data"].fillna(0) < 1)
)

# Documented post-completion QC failure (results_narrative.py RFQ, N=83 marginal).
QC_FAILED = (
    (df_recruit["qc_passed"] < 1)
    & (df_recruit["honesty_qc"].notna() | (df_recruit["student_yn"] > 0))
    & (df_recruit["raven_total_score_v2"] >= 1)
    & (df_recruit["raven_total_score_v2"].notna())
    & ~RSTF
)

# RAVEN component of the screening battery never scored.
NO_RAVEN = (df_recruit["raven_total"] < 1) | (df_recruit["raven_total"].isna())

# Salvage: timed out before finishing, but completed at least the full first
# questionnaire and passed all QC checks.
SALVAGED = df_recruit["salvage_yn"] == 1

# ── Screening eligibility ─────────────────────────────────────────────────────
# Hierarchical screening loop, copied verbatim from results_narrative.py. Tiers
# are checked in order; a participant only reaches a lower tier if the one above
# passed. Reason strings are concatenated, so a participant may carry several.
_reason = pd.Series("", index=df_recruit.index)
for index, row in df_recruit[
    ~(df_recruit["qc_passed"] == 1) & (df_recruit["screening_survey_complete"] > 0)
].iterrows():
    if row["screening_pass"] < 1 or pd.isna(row["geo_crit"]):
        _reason[index] += "Fraud-associated phone # or IP address"
    else:
        nonnegotiables = False
        if row["no_computer"] > 0:
            nonnegotiables = True
            _reason[index] += "No computer"
        if row["english_fluency"] < 1:
            nonnegotiables = True
            _reason[index] += "Non-English speaking"
        if not nonnegotiables:
            majorcriteria = False
            if row["age_v2"] > 65:
                majorcriteria = True
                _reason[index] += ">65 years old"
            if row["age_v2"] < 18:
                majorcriteria = True
                _reason[index] += "<18 years old"
            if row["cognition_screener_v2"] > 0:
                majorcriteria = True
                _reason[index] += "Neurocognitive Impairment"
            if row["seizure_hx_v2"] > 0:
                majorcriteria = True
                _reason[index] += "Epilepsy"
            if row["intox_screen_v2"] > 0:
                majorcriteria = True
                _reason[index] += "Active intoxication "
            if row["raven_total_score_v2"] < 1:
                majorcriteria = True
                _reason[index] += "0 RAVEN score"
            if not majorcriteria:
                if (
                    (row["activecannabisuse_lastuse"] < 28 and row["cannabis_frequency"] > 9)
                    or (row["activecannabisuse_lastuse"] < 14 and row["cannabis_frequency"] in [8, 9])
                    or (row["activecannabisuse_lastuse"] < 7 and row["cannabis_frequency"] == 7)
                    or row["activecannabisuse_lastuse"] < 3
                ):
                    _reason[index] += "Recent heavy cannabis use"
                if row["atypical_since_sp"] > 0:
                    _reason[index] += "Atypical psychedelic more recent than SP"
                if (row["psycheduse_yn"] > 1) or (row["sp_naiive"] < 1):
                    _reason[index] += "No SP Use "

# The loop also tags ~148 *eligible* early-cohort participants as fraud because
# geo_crit did not exist as a field when they enrolled. Restricting to si_2_v2
# null removes them, yielding the canonical ineligible count of 757
# (screened 1249 - eligible 492).
INELIGIBLE = (_reason.str.strip() != "") & df_recruit["si_2_v2"].isna()
ELIGIBLE = SCREENED & ~INELIGIBLE


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: STAGE-CONDITIONED COUNTS
# ══════════════════════════════════════════════════════════════════════════════
# Each participant is counted at exactly one stage, so every arrow subtracts.
#
#   screened − ineligible                      = eligible
#   eligible − failed minimum study measures   = completed minimum measures
#   completed minimum measures − failed QC     = final analytic sample
#
# The quality-control box carries only the eligible subset of the supplementary
# quality-control table. Participants who were screened ineligible and then also
# failed QC leave the flow at the screening step; the same si_2_v2 restriction is
# applied to the QC categories, so both describe the same population.

N_SCREENED = int(SCREENED.sum())
N_INELIGIBLE = int(INELIGIBLE.sum())
N_ELIGIBLE = int(ELIGIBLE.sum())

# Eligible participants who never reached the final analysis df.
LOST = ELIGIBLE & ~IN_DF

# Quality-control failures among the eligible — the population of Table S2.
QC_ELIGIBLE = QC_FAILED & ELIGIBLE
N_QC = int(QC_ELIGIBLE.sum())

# Everyone else who was lost: they never produced the minimum set of measures
# required for analysis. Bucketed by fixed precedence so the parts sum exactly.
MIN_MEASURES = LOST & ~QC_ELIGIBLE
B_NO_RAVEN = MIN_MEASURES & NO_RAVEN
B_NO_RECORD = MIN_MEASURES & df_recruit["qc_passed"].isna() & ~B_NO_RAVEN
B_TASK = MIN_MEASURES & ~B_NO_RAVEN & ~B_NO_RECORD

N_MIN_MEASURES = int(MIN_MEASURES.sum())
N_NO_RAVEN = int(B_NO_RAVEN.sum())
N_NO_RECORD = int(B_NO_RECORD.sum())
N_TASK = int(B_TASK.sum())

N_AFTER_MIN = N_ELIGIBLE - N_MIN_MEASURES
N_FINAL = int(IN_DF.sum())
N_SP = int((df["psycheduse_yn"] == "Yes").sum())
N_NAIVE = N_FINAL - N_SP

# Salvage records are inside the final sample; they never leave the spine. The
# count is annotated on the final box rather than drawn as a re-entry arrow.
N_SALVAGED = int((ELIGIBLE & SALVAGED & IN_DF).sum())

# ── Analytic subsamples (the four Ns the reviewer asked about) ────────────────
# vch_threshold is non-null only for participants whose VCH data passed task QC;
# task_data_vch_* is non-null for everyone who attempted the task.
N_SP_VCH_ATTEMPTED = int(df_sp["task_data_vch_short_psychedelic_bl"].notna().sum())
N_SP_VCH = int(df_sp["vch_threshold"].notna().sum())
N_SP_VCH_QC_FAIL = N_SP_VCH_ATTEMPTED - N_SP_VCH
N_SP_VCH_NOT_DONE = N_SP - N_SP_VCH_ATTEMPTED
N_SP_CAPS = int(df_sp["caps_vision"].notna().sum())
N_SP_BOTH = int((df_sp["vch_threshold"].notna() & df_sp["caps_vision"].notna()).sum())

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4: CLOSURE ASSERTIONS
# ══════════════════════════════════════════════════════════════════════════════
# The whole point of this figure is that the arithmetic closes. If a data refresh
# breaks any of these, stop and surface it rather than render a misleading figure.

_checks = [
    ("screened - ineligible = eligible",
     N_SCREENED - N_INELIGIBLE, N_ELIGIBLE),
    ("eligible - failed minimum measures = completed minimum measures",
     N_ELIGIBLE - N_MIN_MEASURES, N_AFTER_MIN),
    ("completed minimum measures - QC failures = final sample",
     N_AFTER_MIN - N_QC, N_FINAL),
    ("minimum-measure buckets sum to total",
     N_NO_RAVEN + N_NO_RECORD + N_TASK, N_MIN_MEASURES),
    ("SP users + SP-naive = final",
     N_SP + N_NAIVE, N_FINAL),
    ("final analytic sample matches shipped df",
     N_FINAL, len(df)),
]
_failures = [(label, got, want) for label, got, want in _checks if got != want]
if _failures:
    msg = "\n".join(f"  {label}: got {got}, expected {want}" for label, got, want in _failures)
    raise AssertionError(
        "Participant flow no longer closes — the data or the stage masks have "
        f"changed. Resolve before regenerating the figure:\n{msg}"
    )

# No participant may be excluded at one stage and still appear in the final df.
for _label, _mask in [("ineligible", INELIGIBLE), ("QC-failed", QC_ELIGIBLE)]:
    _overlap = int((_mask & IN_DF).sum())
    if _overlap:
        raise AssertionError(
            f"DATA INTEGRITY ERROR: {_overlap} {_label} record(s) appear in the "
            "final analysis df."
        )

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5: EXCLUSION-REASON BULLETS
# ══════════════════════════════════════════════════════════════════════════════
# Screening-exclusion bullets, counted off the `_reason` strings this script
# already builds in Section 2 — not read from a separate table.
#
# Counting them here rather than in a sibling script keeps a single copy of the
# screening loop — a second copy is exactly the drift risk this directory's
# README warns about.
#
# Reasons at tiers 3 and 4 can co-occur for one participant (e.g. recent
# cannabis AND atypical psychedelics), so each bullet asks whether that reason
# appears anywhere in the participant's concatenated reason string.  Percentages
# are out of the total ineligible N and therefore do not sum to 100%.
#
# Each entry is (display label, substring to match).  The substrings are the
# discriminating fragments of the reason strings written by the loop above —
# not the full strings, because two of them ("Active intoxication ",
# "No SP Use ") are written with a trailing space that the display label drops.
INELIG_REASONS = [
    ("Fraud-associated phone # or IP address",   "Fraud"),
    ("No computer",                              "No computer"),
    ("Non-English speaking",                     "Non-English"),
    (">65 years old",                            ">65"),
    ("<18 years old",                            "<18"),
    ("Neurocognitive Impairment",                "Neurocognitive"),
    ("Epilepsy",                                 "Epilepsy"),
    ("Active intoxication",                      "Active intoxication"),
    ("0 RAVEN score",                            "0 RAVEN"),
    ("Recent heavy cannabis use",                "Recent heavy cannabis"),
    ("Atypical psychedelic more recent than SP", "Atypical psychedelic"),
    ("No SP Use",                                "No SP Use"),
]

_inelig_reasons = _reason[INELIGIBLE]
_n_inelig = len(_inelig_reasons)

# Every ineligible participant must be captured by at least one bullet; if the
# loop ever writes a reason string no substring matches, the bullets would
# silently under-count the box they sit in.
_captured = pd.Series(False, index=_inelig_reasons.index)
for _, _sub in INELIG_REASONS:
    _captured |= _inelig_reasons.str.contains(_sub, na=False, regex=False)
if not _captured.all():
    raise ValueError(
        f"{(~_captured).sum()} ineligible participants match no INELIG_REASONS "
        f"substring. Unmatched reason strings:\n"
        f"{_inelig_reasons[~_captured].value_counts().to_string()}"
    )

INELIG_BULLETS = []
for _label, _sub in INELIG_REASONS:
    _n = int(_inelig_reasons.str.contains(_sub, na=False, regex=False).sum())
    INELIG_BULLETS.append(f"{_label} — {_n} ({100 * _n / _n_inelig:.1f}%)")

# ── Quality-control failure reasons ───────────────────────────────────────────
# The four QC-failure counts that fill the quality-control box of the diagram,
# computed off the QC_ELIGIBLE mask this script already derives.  qc_notes (RA
# free text) and phone_number do not survive de-identification, so the category
# memberships are frozen record_id sets in modules/qc_redacted_categories.py.
# No participant-level QC artifact is written.
QC_CATEGORIES = OrderedDict((
    ("Inconsistent Answers",        QCCAT_INCONSISTENT),
    ("Failed Attention Checks",     QCCAT_ATTN_CHECK),
    ("Failed Challenge Questions",  QCCAT_CHALLENGE),
))

df_qc_failed = df_recruit[QC_ELIGIBLE].copy()
_rid = df_qc_failed["record_id"].astype(int)
for _label, _ids in QC_CATEGORIES.items():
    df_qc_failed[_label] = _rid.isin(_ids).astype(int)

# Category 4 is the residual: a fraud signal with no independent Cat 1-3 hit.
# Recomputed rather than looked up so the residual logic stays visible.
_fraud_signal = _rid.isin(QCCAT_FRAUD_SIGNAL)
df_qc_failed["Fraud-associated Phone / IP"] = (
    _fraud_signal
    & ~df_qc_failed[list(QC_CATEGORIES)].any(axis=1)
).astype(int)
if set(_rid[df_qc_failed["Fraud-associated Phone / IP"] == 1]) != set(QCCAT_FRAUD):
    raise AssertionError(
        "Frozen Category 4 membership no longer matches the residual computed "
        "from Categories 1-3 and the fraud signal — regenerate "
        "modules/qc_redacted_categories.py from the private recruitment export."
    )

QC_LABELS = list(QC_CATEGORIES) + ["Fraud-associated Phone / IP"]
_uncaptured = int((df_qc_failed[QC_LABELS].sum(axis=1) == 0).sum())
if _uncaptured:
    raise AssertionError(
        f"{_uncaptured} of {N_QC} QC-failed participants fall into no category; "
        "the frozen sets no longer describe this population."
    )

QC_BULLETS = [
    f"{lab} — {int(df_qc_failed[lab].sum())} "
    f"({100 * df_qc_failed[lab].sum() / N_QC:.1f}%)"
    for lab in QC_LABELS
]

print("\nParticipant flow (all closure checks passed):")
print(f"  {N_SCREENED} screened - {N_INELIGIBLE} ineligible = {N_ELIGIBLE} eligible")
print(f"  {N_ELIGIBLE} eligible - {N_MIN_MEASURES} failed minimum study measures "
      f"({N_NO_RAVEN} unscored battery / {N_NO_RECORD} no QC record / {N_TASK} task) "
      f"= {N_AFTER_MIN}")
print(f"  {N_AFTER_MIN} - {N_QC} QC failures = {N_FINAL} final "
      f"({N_SP} SP users, {N_NAIVE} SP-naive; incl. {N_SALVAGED} salvaged)")
print(f"  SP subsamples: VCH {N_SP_VCH}, CAPS {N_SP_CAPS}, both {N_SP_BOTH}")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6: FIGURE
# ══════════════════════════════════════════════════════════════════════════════
# Half-page width: 91.5 mm, half of the 183 mm Scientific Reports double-column
# measure recorded in modules/figure_assembly.py. This matches FIG_WIDTH_IN in
# hardware_keydown_check.py (3.602 in), the repo's existing half-page figures.
#
# At this width there is no room for an exclusion column beside the spine, so
# the flow is a single vertical stack: white boxes are stages, grey boxes are
# exclusions leaving the stage above them. Layout is cursor-driven rather than
# hardcoded — each element is appended below the previous one and the figure
# height is derived from the content, so changing a font size or adding a bullet
# reflows the figure instead of silently overlapping.

MM_PER_IN = 25.4
FIG_W = (183.0 / 2) / MM_PER_IN      # 91.5 mm — half of the 183 mm page width
X_SPAN = 100.0                       # x coordinate space
UNIT_PT = FIG_W * 72.0 / X_SPAN      # points per coordinate unit (~2.49)

def pt(points):
    """Convert points to coordinate units, so type sizes drive the layout."""
    return points / UNIT_PT

MAIN_FC = "white"                    # stage box fill
SIDE_FC = "#F2F2F2"                  # exclusion box fill
EDGE = "black"
LW = 0.6

# Type sizes in points, as they will appear on the printed page.
FS_MAIN = 7.0                        # stage box
FS_FINAL = 8.0                       # final-sample box
FS_SIDE_HEAD = 6.5                   # exclusion box headline
FS_BULLET = 5.6                      # exclusion bullets
FS_SUB = 6.0                         # subsample and SP-naive boxes
FS_NOTE = 5.2                        # annotation inside the final box

LINE = 1.45                          # line spacing multiplier
PAD_Y = pt(4.0)                      # vertical padding inside a box
ARROW = pt(9.0)                      # arrow length between elements
GAP = pt(3.0)                        # gap between paired boxes

L, R = 1.0, 99.0                     # full-width box edges
W_FULL = R - L
MID = (L + R) / 2
W_HALF = (W_FULL - GAP) / 2
X_LEFT = L + W_HALF / 2              # centre of the left half-width box
X_RIGHT = R - W_HALF / 2             # centre of the right half-width box


def text_h(n_lines, fs):
    """Height of an n-line text block at font size fs, in coordinate units."""
    return pt(n_lines * fs * LINE)


# ── Pass 1: measure ───────────────────────────────────────────────────────────
# Heights are computed before anything is drawn so the figure can be sized to
# its content instead of being trimmed by bbox_inches afterwards.
H_SCREEN = text_h(2, FS_MAIN) + 2 * PAD_Y
H_INELIG = text_h(1, FS_SIDE_HEAD) + text_h(len(INELIG_BULLETS), FS_BULLET) + 2.6 * PAD_Y
H_ELIG = text_h(2, FS_MAIN) + 2 * PAD_Y
H_MINM = text_h(1, FS_SIDE_HEAD) + 2 * PAD_Y
H_AFTER = text_h(2, FS_MAIN) + 2 * PAD_Y
H_QC = text_h(1, FS_SIDE_HEAD) + text_h(len(QC_BULLETS), FS_BULLET) + 2.6 * PAD_Y
H_FINAL = text_h(2, FS_FINAL) + text_h(1, FS_NOTE) + 2 * PAD_Y
H_SP = text_h(2, FS_SUB) + 2 * PAD_Y
H_SUB = text_h(4, FS_SUB) + 2 * PAD_Y
H_COMB = text_h(2, FS_SUB) + 2 * PAD_Y

TOTAL_H = (
    H_SCREEN + ARROW + H_INELIG + ARROW + H_ELIG + ARROW + H_MINM + ARROW
    + H_AFTER + ARROW + H_QC + ARROW + H_FINAL + ARROW + H_SP + ARROW
    + H_SUB + ARROW + H_COMB + 2 * PAD_Y
)

FIG_H = FIG_W * TOTAL_H / X_SPAN     # 1 x-unit and 1 y-unit are the same length

fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
ax.set_xlim(0, X_SPAN)
ax.set_ylim(0, TOTAL_H)
ax.axis("off")

cursor = TOTAL_H - PAD_Y             # top of the next element, walking downward


def draw_box(h, *, xc=MID, w=W_FULL, fc=MAIN_FC, top=None):
    """Draw a box of height h whose top edge is at `top` (default: the cursor)."""
    y_top = cursor if top is None else top
    ax.add_patch(FancyBboxPatch(
        (xc - w / 2, y_top - h), w, h,
        boxstyle="square,pad=0", linewidth=LW,
        edgecolor=EDGE, facecolor=fc, zorder=2))
    return y_top - h


def stage(text, h, *, fs=FS_MAIN, weight="normal", xc=MID, w=W_FULL,
          fc=MAIN_FC, top=None):
    """White stage box with centred text."""
    y_top = cursor if top is None else top
    y_bot = draw_box(h, xc=xc, w=w, fc=fc, top=y_top)
    ax.text(xc, (y_top + y_bot) / 2, text, ha="center", va="center",
            fontsize=fs, fontweight=weight, zorder=3, linespacing=LINE)
    return y_bot


def exclusion(headline, bullets, h, *, top=None):
    """Grey exclusion box: bold headline, then one bulleted line per reason."""
    y_top = cursor if top is None else top
    y_bot = draw_box(h, fc=SIDE_FC, top=y_top)
    if not bullets:
        ax.text(MID, (y_top + y_bot) / 2, headline, ha="center", va="center",
                fontsize=FS_SIDE_HEAD, fontweight="bold", zorder=3)
        return y_bot
    y = y_top - PAD_Y
    ax.text(L + pt(4), y, headline, ha="left", va="top",
            fontsize=FS_SIDE_HEAD, fontweight="bold", zorder=3)
    y -= text_h(1, FS_SIDE_HEAD) + pt(1.5)
    for b in bullets:
        ax.text(L + pt(7), y, "• " + b, ha="left", va="top",
                fontsize=FS_BULLET, zorder=3)
        y -= text_h(1, FS_BULLET)
    return y_bot


def arrow(y_from, *, x=MID, length=ARROW):
    """Downward arrow of the standard length; returns the new cursor."""
    ax.add_patch(FancyArrowPatch(
        (x, y_from), (x, y_from - length), arrowstyle="-|>", mutation_scale=5,
        linewidth=LW, color=EDGE, shrinkA=0, shrinkB=0, zorder=1))
    return y_from - length


# ── Pass 2: draw, walking down the page ───────────────────────────────────────
# Thousands separator must be brace-escaped inside mathtext, otherwise the comma
# is typeset as a punctuation operator and renders as "1, 249".
cursor = stage(
    f"Completed screening\n$\\bf{{n = {N_SCREENED // 1000}{{,}}{N_SCREENED % 1000:03d}}}$",
    H_SCREEN)
cursor = arrow(cursor)
cursor = exclusion(f"Ineligible at screening (n = {N_INELIGIBLE})",
                   INELIG_BULLETS, H_INELIG)
cursor = arrow(cursor)
cursor = stage(f"Eligible\n$\\bf{{n = {N_ELIGIBLE}}}$", H_ELIG)
cursor = arrow(cursor)
cursor = exclusion(f"Failed to Complete Minimum Measures (n = {N_MIN_MEASURES})",
                   [], H_MINM)
cursor = arrow(cursor)
cursor = stage(f"Completed minimum study measures\n$\\bf{{n = {N_AFTER_MIN}}}$", H_AFTER)
cursor = arrow(cursor)
cursor = exclusion(f"Failed quality control (n = {N_QC})", QC_BULLETS, H_QC)
cursor = arrow(cursor)

# Final sample, with the salvage records annotated rather than drawn as an inflow.
_top = cursor
_bot = draw_box(H_FINAL, top=_top)
ax.text(MID, _top - PAD_Y - text_h(2, FS_FINAL) / 2,
        f"FINAL ANALYTIC SAMPLE\n$\\bf{{n = {N_FINAL}}}$",
        ha="center", va="center", fontsize=FS_FINAL, fontweight="bold",
        zorder=3, linespacing=LINE)
ax.text(MID, _bot + PAD_Y + text_h(1, FS_NOTE) / 2,
        f"includes {N_SALVAGED} who timed out with sufficient data",
        ha="center", va="center", fontsize=FS_NOTE, zorder=3)
cursor = _bot - ARROW
ax.plot([MID, MID], [_bot, _bot + ARROW / 2 - ARROW], color=EDGE, lw=LW, zorder=1)

# SP users and SP-naive split the width; the flow continues from SP users.
_top = cursor
stage(f"SP users\n$\\bf{{n = {N_SP}}}$", H_SP, fs=FS_SUB,
      xc=X_LEFT, w=W_HALF, top=_top)
stage(f"SP-naïve\n$\\bf{{n = {N_NAIVE}}}$", H_SP, fs=FS_SUB,
      xc=X_RIGHT, w=W_HALF, fc=SIDE_FC, top=_top)
ax.plot([X_LEFT, X_RIGHT], [_top + ARROW / 2, _top + ARROW / 2],
        color=EDGE, lw=LW, zorder=1)
for _x in (X_LEFT, X_RIGHT):
    ax.add_patch(FancyArrowPatch(
        (_x, _top + ARROW / 2), (_x, _top), arrowstyle="-|>", mutation_scale=5,
        linewidth=LW, color=EDGE, shrinkA=0, shrinkB=0, zorder=1))
cursor = _top - H_SP

# VCH and CAPS side by side, fed from SP users.
_bus = cursor - ARROW / 2
ax.plot([X_LEFT, X_LEFT], [cursor, _bus], color=EDGE, lw=LW, zorder=1)
ax.plot([X_LEFT, X_RIGHT], [_bus, _bus], color=EDGE, lw=LW, zorder=1)
_top = cursor - ARROW
for _x in (X_LEFT, X_RIGHT):
    ax.add_patch(FancyArrowPatch(
        (_x, _bus), (_x, _top), arrowstyle="-|>", mutation_scale=5,
        linewidth=LW, color=EDGE, shrinkA=0, shrinkB=0, zorder=1))
stage((f"VCH task data\nFailed to complete: −{N_SP_VCH_NOT_DONE}\n"
       f"Excluded at task QC: −{N_SP_VCH_QC_FAIL}\n$\\bf{{n = {N_SP_VCH}}}$ analysed"),
      H_SUB, fs=FS_SUB, xc=X_LEFT, w=W_HALF, top=_top)
stage((f"CAPS data\nCompleted CAPS: {N_SP_CAPS}\n\n"
       f"$\\bf{{n = {N_SP_CAPS}}}$ analysed"),
      H_SUB, fs=FS_SUB, xc=X_RIGHT, w=W_HALF, top=_top)
cursor = _top - H_SUB

# Converge both subsamples into the combined-model box.
_conv = cursor - ARROW / 2
for _x in (X_LEFT, X_RIGHT):
    ax.plot([_x, _x], [cursor, _conv], color=EDGE, lw=LW, zorder=1)
ax.plot([X_LEFT, X_RIGHT], [_conv, _conv], color=EDGE, lw=LW, zorder=1)
ax.add_patch(FancyArrowPatch(
    (MID, _conv), (MID, cursor - ARROW), arrowstyle="-|>", mutation_scale=5,
    linewidth=LW, color=EDGE, shrinkA=0, shrinkB=0, zorder=1))
stage(f"VCH and CAPS both present\n$\\bf{{n = {N_SP_BOTH}}}$ analysed",
      H_COMB, fs=FS_SUB, top=cursor - ARROW)

fig.savefig(OUT_DIR / "supplementary_figure_s9.png", dpi=DPI, bbox_inches="tight",
            facecolor="white")
fig.savefig(OUT_DIR / "supplementary_figure_s9.tiff", dpi=DPI, bbox_inches="tight",
            facecolor="white")
fig.savefig(OUT_DIR / "supplementary_figure_s9.svg", bbox_inches="tight",
            facecolor="white")
plt.close(fig)

print(f"\nFigure: {FIG_W * MM_PER_IN:.0f} mm x {FIG_H * MM_PER_IN:.0f} mm "
      f"({FIG_W:.2f} x {FIG_H:.2f} in) at {DPI} dpi")
print(f"Figure saved to: {OUT_DIR}/supplementary_figure_s9.{{png,tiff,svg}}")
