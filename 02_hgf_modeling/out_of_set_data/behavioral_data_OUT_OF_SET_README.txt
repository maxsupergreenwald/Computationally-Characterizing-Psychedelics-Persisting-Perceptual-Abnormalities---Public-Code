NOTE (public release): the CSV described below, behavioral_data_OUT_OF_SET_with_metadata.csv,
is NOT distributed with this repository - it is participant-level data from a
separate study. This column dictionary is included so the derivation of the
empirical_condition values can be understood without the raw file. The document
describes the source dataset as collected, which covers both the visual (v) and
auditory (a) modalities; only the visual data were used in this manuscript.
------------------------------------------------------------------------------

behavioral_data_OUT_OF_SET_with_metadata.csv
============================================

Trial-level behavioral data for the OUT-OF-SET cohort (sex_matched == False).
This cohort is independent of the sex-matched analysis sample and was used to
compute the empirical_condition lookup that the production HGF fits consume.

Subjects: 132 subjects (A+V+ 88, A-V- 44), 86,040 trials total.
After 6-flag QC: 109 aud-passing cells, 114 vis-passing cells, 127 unique subjects
with ≥1 passing modality, 96 paired (both modalities pass).
Groups represented in this cohort: A+V+, A-V- only.
A+V- and A-V+ have 0 subjects here by sex-match construction — those small groups
were retained entirely in the sex-matched main cohort.

Cohort selection:
  - From behavioral_data_all_subjects_with_qc.csv
  - Filter: sex_matched == False

Per-(subject, modality) QC flags (the same 6-flag set used in the production fits, from
comprehensive_qc_with_redcap.py:142-149):
  fail_neg  : logistic regression of response ~ condition has β < 0
  fail_ns   : logistic regression β not significant at p < 0.05
  fail_lv   : variance(response) < 0.01
  fail_hzr  : response rate at condition=0 > 0.50  (high false-alarm)
  fail_flat : mean(response) < 0.05 OR > 0.95
  fail_nie  : (response@75 − response@0) < 0.10   (no intensity effect)

A cell passes 6-flag QC iff NONE of the above is True.

Subject counts:
  Total subjects in file:         132   (A+V+ = 88, A-V- = 44)
  6-flag QC pass, ≥1 modality:    127
    aud-passing:                  109
    vis-passing:                  114
  6-flag QC pass, paired (both):   96

Symptom-frequency coverage (from REDCap CHAT items):
  Subjects with non-null max_ah_freq:  107 / 132
  Subjects with non-null max_vh_freq:  112 / 132
  Missing values are subjects who didn't complete the relevant CHAT items
  (typically pure controls who skipped the symptom batteries).

Columns:
  sudo_rec            int, anonymized subject id
  project             str, data source: 'cope', 'nocb', or 'cb'
  group               str, phenotype: 'A+V+', 'A+V-', 'A-V+', 'A-V-'
  sex                 str, 'F' or 'M' (corrected from earlier boolean-encoding bug)
  sex_matched         bool, always False in this file
  qc_pass             bool, pre-existing subject-level QC flag from source
  six_flag_qc_pass    bool, our recomputed per-cell 6-flag QC (use this)
  max_ah_freq         float, max(chat_ah_freq1..7) from REDCap CHAT baseline items.
                       Auditory-hallucination frequency on the highest-rated item.
                       NaN for subjects who didn't complete CHAT (typically pure
                       controls — see counts above).
  max_vh_freq         float, max(chat_vh_freq1..4) from REDCap CHAT baseline items.
                       Visual-hallucination frequency on the highest-rated item.
                       Follow-up (fu_6mo) items are NOT included in the max.
  trial               int, 1..360 within modality
  modality            str, 'a' or 'v'
  condition           int, nominal % signal: 0, 25, 50, or 75
  response            int, binary yes(1)/no(0)
  confidence_binned   int, 1..5
  responseTime        float, RT in ms (some trials have negative or 0 values - drop those)
  empirical_condition float, per-(modality, condition) mean response rate computed
                       FROM THIS COHORT (out-of-set + 6-flag QC). Constant within
                       any (modality, condition) cell — the same 8 values
                       (4 conditions × 2 modalities) stored row-wise across all trials.
                       Identical to the rows in empirical_condition_qc_pass.csv.
                       That lookup CSV is the per-cell aggregation OF THIS dataset
                       (after 6-flag QC), then broadcast back into the trial-level
                       data as a constant lookup for downstream HGF fits.
  fail_neg, fail_ns, fail_lv, fail_hzr, fail_flat, fail_nie  bool, individual QC flags

To use with downstream analyses:
  df = pd.read_csv('behavioral_data_OUT_OF_SET_with_metadata.csv')
  df_qc = df[df['six_flag_qc_pass']]              # apply 6-flag QC
  df_qc = df_qc[df_qc['responseTime'] > 0]        # drop impossible RTs (3 trials)
