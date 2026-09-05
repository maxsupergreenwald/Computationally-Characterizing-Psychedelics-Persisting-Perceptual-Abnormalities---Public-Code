# Out-of-Set Cohort — Empirical Condition Analysis (VCH Only)

> **Scope:** All analyses in this directory are restricted to the **VCH (visual
> contrast hallucination)** modality. ACH (auditory) data are present in the
> source CSV but are not analyzed here.

## Background

The out-of-set cohort (132 subjects with `sex_matched == False`) was used to
derive the **empirical condition lookup** that maps each nominal contrast level
(0 %, 25 %, 50 %, 75 %) to an empirically-measured hit rate. Those four scalars
are fed into the production VCH HGF fits as the trial-by-trial stimulus strength
instead of the nominal contrast values.

This directory contains analyses interrogating those lookup values: how they
were computed, whether the pipeline consumed the correct values, and whether the
empirical hit rates differ significantly from the nominal proportions — both for
the full QC-passing sample and for a "non-hallucinator" subset.

---

## Analysis: `empirical_condition_analysis.py`

**Run with:** `/usr/local/bin/python3.12 empirical_condition_analysis.py`

---

### Step 1 — Verification: did the pipeline use QC-filtered subjects?

We checked which computation reproduces the stored `empirical_condition` column
(broadcast as a constant per condition cell into every VCH trial row).
Three sources were compared:

| Source | Computation |
|--------|-------------|
| **stored CSV** | `empirical_condition` column already in the file |
| **unfiltered grand mean** | `mean(response)` over all VCH subjects, no QC |
| **QC-filtered grand mean** | `mean(response)` over 6-flag QC-passing VCH subjects (n=114) |
| **pipeline hardcoded** | `VCH_CONDITION_TO_EMPIRICAL` dict in `hgf_pipeline.py` |

**Finding:** The stored CSV column reproduces exactly from the **unfiltered**
grand mean. The pipeline dict reproduces exactly from the **QC-filtered** grand
mean. These differ by up to ~2 percentage points per condition. The
`empirical_condition` column in the CSV is a preliminary, pre-QC artifact and
does **not** reflect what the HGF consumed.

**Conclusion: the production pipeline correctly used the QC-filtered
subjects (n=114 VCH-passing).** The hardcoded pipeline values match the
6-flag-QC-passing grand means to machine precision.

**VCH condition=0 note:** The pipeline hardcodes condition=0 to `0.0`.
The empirical false-alarm rate at condition=0 is ~11 %, significantly above
zero, but it is deliberately NOT used — see *Resolved Decisions* below.

Full comparison table: `empirical_condition_verification.csv`.

---

### Step 2 — Empirical detection probability by condition, full vs. non-hallucinator

Two overlaid series are plotted for VCH:

#### Series A — Full QC-passing sample (n=114)

All VCH subjects passing the 6-flag QC (`six_flag_qc_pass == True`,
`responseTime > 0`).

#### Series B — Non-hallucinator subset (n=29)

Subjects satisfying `max_vh_freq == 0 AND max_ah_freq == 0`, where `NaN` is
treated as 0.

**Rationale for treating NaN as 0:** Alexandria's `behavioral_data_OUT_OF_SET_README.txt`
states:

> *"Missing values are subjects who didn't complete the relevant CHAT items
> (typically pure controls who skipped the symptom batteries)."*

Subjects who skipped the symptom battery are presumed to have no hallucination
history (i.e., they are controls who were never asked about symptoms because
they did not endorse them). Strict `== 0` (excluding NaN) yields only n=7,
which is insufficient for CI estimation. Including NaN-as-0 yields n=29.

**Statistical test:** one-sample t-test per condition vs. the nominal proportion
(0, .25, .50, .75). Group mean ± 95 % CI computed from per-subject means
(t-distribution across subjects).

**Figure:** `figures/empirical_condition_by_condition.png`

---

## Key Findings and Decision (2026-06-17)

1. **Pipeline integrity confirmed.** The production HGF fits (through 2026-06-16)
   consumed the QC-filtered full-sample grand means (n=114 VCH-passing), not the
   unfiltered values stored in the CSV column.

2. **Empirical rates are systematically and significantly higher than nominal
   proportions at every condition.** At 25 % contrast, VCH detection ≈ 44 %
   (vs. nominal 25 %, p < 0.001). This motivates using empirical rather than
   nominal inputs to the HGF.

3. **Non-hallucinator and full-sample curves are nearly identical** (all 95 % CIs
   overlap; Δ ≤ 3 pp at any condition). See
   `figures/empirical_condition_by_condition.png`.

4. **Decision: switched to non-hallucinator values on 2026-06-17.** The HGF
   detection likelihood should be grounded in the perceptual behavior of
   participants without hallucination histories. All scripts updated simultaneously;
   see `../README.md §Stimulus intensity: empirical vs. nominal` for the current
   values and the list of scripts that hardcode them.

   New values (non-hallucinator, n=29):

   | Condition | New value | Old value | Note |
   |-----------|-----------|-----------|------|
   | 0 % | 0.0 | 0.0 | **Unchanged** — ground-truth (no stimulus); do NOT use empirical FA rate |
   | 25 % | 0.4180444024563061 | 0.44184090362893536 | |
   | 50 % | 0.7115104419621175 | 0.7010130961205832 | |
   | 75 % | 0.8994252873563219 | 0.8690302144249513 | |

   The 0 % condition was *incorrectly* set to 0.08620689655172414 (the
   non-hallucinator false-alarm rate) in the original 2026-06-17 commit.
   This was corrected: the 0 % value is always 0.0 because there is no
   stimulus on those trials and the ground-truth detection probability is
   exactly 0.0.

---

## Resolved Decisions

- **VCH condition=0 policy — RESOLVED:** The 0 % condition value is fixed at
  **0.0** and is never replaced with the empirical false-alarm rate. At 0 %
  contrast there is literally no stimulus; the ground-truth detection probability
  is therefore exactly 0.0. Using the empirical FA rate here would misrepresent
  the task structure to the HGF.
- **Non-hallucinator-derived lookup — RESOLVED:** Switched to non-hallucinator
  values (n=29) for 25/50/75 % conditions on 2026-06-17. The full-sample and
  non-hallucinator curves are nearly identical (all 95 % CIs overlap; Δ ≤ 3 pp),
  confirming robustness. All production scripts updated simultaneously.

---

## Files

| File | Description |
|------|-------------|
| `behavioral_data_OUT_OF_SET_with_metadata.csv` | **NOT distributed** — trial-level data for the out-of-set cohort. Belongs to a separate study; the script cannot be re-run without it, but every derived output below is included. |
| `behavioral_data_OUT_OF_SET_README.txt` | Column-level documentation for that file, incl. QC-flag definitions |
| `empirical_condition_analysis.py` | Analysis script (Steps 1 & 2) |
| `empirical_condition_verification.csv` | Step 1: VCH comparison table (all sources × all conditions) |
| `empirical_condition_full_vs_nonhall.csv` | Step 2: per-condition statistics, full sample vs. non-hallucinator subset |
| `figures/empirical_condition_by_condition.png` | Step 2: full vs. non-hallucinator figure |
| `figures/empirical_condition_by_condition_spusers.png` | Step 2, SP-user variant (`psycheduse_yn == "Yes"`) — Supplementary Fig. S2a |
| `README_empirical_condition_redetermination.md` | This file |
