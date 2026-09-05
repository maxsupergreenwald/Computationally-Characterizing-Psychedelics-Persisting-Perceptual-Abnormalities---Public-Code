"""
Frozen QC-failure category assignments — replaces the redacted `qc_notes` column.

WHY THIS FILE EXISTS
--------------------
`qc_notes` is the research assistant's free-text note explaining why a
participant failed quality control.  It is genuinely load-bearing: the QC
failure table and the fraud-overlap audit both classify participants by
pattern-matching that text.  It is also personally identifying — the notes name
a participant's country and describe fraud investigations in prose — so it is
NOT shipped in `data/final/df_recruit_public_<date>.csv`.

The classification was therefore run once against the original notes and the
resulting record_id memberships frozen here.  Every published count is
reproduced exactly: Categories 1-4 hold 39 / 12 / 24 / 6 of the 69 eligible
QC-failed participants, matching Supplementary Table S2.

WHAT THIS MEANS FOR REPRODUCIBILITY
-----------------------------------
`consort_diagram.py` consumes these sets to draw the quality-control box and
to render Supplementary Table S2.  That classification is therefore a replay of
a frozen result, not a live computation; every other count in the participant
flow still computes from survey columns that DO ship and will change if the data
change.  A data refresh therefore requires these constants to be regenerated
from the private recruitment export alongside
`06_submission/deidentify_recruit_csv.py` in the master repository.

The fraud signal in these lists is a phone/IP-based flag recorded during data
collection; the majority of flagged numbers carried the +234 country code.  The
underlying phone numbers are not shipped — `df_recruit_public_<date>.csv`
carries only the derived `phone_cc_234` boolean.

Categories are NON-EXCLUSIVE: a participant may appear in more than one.
"""

# QC-failed participants whose qc_notes cell was empty; their exclusion reason
# was reconstructed from survey columns rather than from a note.
NAN_NOTES_RECORD_IDS = frozenset({
    225, 230, 232, 233, 316, 377, 423, 440, 442, 447, 458, 461, 485,
    492, 1989, 1997, 2028, 2214, 2345, 2353, 2382, 2418, 2430, 2443,
    2448, 2475, 2578, 2666
})   # n = 28

# Category 1 - Inconsistent Answers (eligible QC-failed, n=69 population)
QCCAT_INCONSISTENT = frozenset({
    209, 233, 238, 280, 316, 325, 340, 341, 344, 345, 349, 350, 353,
    377, 381, 382, 385, 400, 420, 422, 426, 433, 439, 440, 442, 444,
    447, 458, 481, 483, 484, 485, 486, 487, 489, 490, 624, 1858, 2028
})   # n = 39

# Category 2 - Failed Attention Checks (same population)
QCCAT_ATTN_CHECK = frozenset({
    209, 211, 221, 354, 360, 385, 386, 388, 389, 430, 431, 469
})   # n = 12

# Category 3 - Failed Challenge Questions (same population)
QCCAT_CHALLENGE = frozenset({
    221, 232, 233, 316, 377, 392, 440, 446, 447, 458, 469, 485, 512,
    978, 1055, 1583, 1673, 1684, 1773, 1792, 2028, 2448, 2475, 2614
})   # n = 24

# Category 4 - Fraud-associated phone / IP, residual of Cats 1-3
QCCAT_FRAUD = frozenset({
    225, 230, 435, 664, 2018, 2456
})   # n = 6

# Any fraud signal, before the Cat 1-3 residual is applied
QCCAT_FRAUD_SIGNAL = frozenset({
    211, 221, 225, 230, 232, 233, 316, 325, 345, 354, 360, 377, 386,
    388, 389, 400, 430, 431, 435, 440, 442, 447, 458, 469, 485, 664,
    2018, 2456
})   # n = 28

# Section 2 flag `_qccat_fraud` over all 87 QC-failed participants
ALL87_FRAUD = frozenset({
    211, 221, 325, 345, 354, 360, 386, 388, 389, 400, 430, 431, 435,
    469, 664, 2018, 2199, 2344, 2359, 2365, 2456
})   # n = 21

# Section 2 flag `_qccat_inconsistent_sp_use` over all 87 QC-failed participants
ALL87_INCONSISTENT_SP_USE = frozenset({
    209, 238, 280, 340, 344, 345, 349, 350, 353, 381, 382, 400, 420,
    426, 433, 439, 444, 483, 484, 486, 487, 489, 490
})   # n = 23

# Section 2 flag `_qccat_inconsistent_demographics` over all 87 QC-failed participants
ALL87_INCONSISTENT_DEMOGRAPHICS = frozenset({
    325, 341, 345, 382, 385, 400, 422, 426, 481, 1858
})   # n = 10

# Section 2 flag `_qccat_implausible_psych` over all 87 QC-failed participants
ALL87_IMPLAUSIBLE_PSYCH = frozenset({
    221, 392, 469, 512, 624, 978, 1055, 1583, 1673, 1684, 1773, 1792,
    2614
})   # n = 13

# Section 2 flag `_qccat_attn_check` over all 87 QC-failed participants
ALL87_ATTN_CHECK = frozenset({
    209, 211, 221, 354, 360, 385, 386, 388, 389, 430, 431, 469
})   # n = 12

# Section 2 flag `_qccat_self_admission` over all 87 QC-failed participants
ALL87_SELF_ADMISSION = frozenset({
    446
})   # n = 1

# Section 2 flag `_qccat_bad_data_no_note` over all 87 QC-failed participants
ALL87_BAD_DATA_NO_NOTE = frozenset({
    233, 316, 377, 423, 440, 442, 447, 458, 461, 485, 2028
})   # n = 11

# Section 2 flag `_has_fraud_signal` over all 87 QC-failed participants
ALL87_HAS_FRAUD_SIGNAL = frozenset({
    211, 221, 225, 230, 232, 233, 316, 325, 345, 354, 360, 377, 386,
    388, 389, 400, 423, 430, 431, 435, 440, 442, 447, 458, 469, 485,
    492, 664, 2018, 2199, 2344, 2359, 2365, 2456
})   # n = 34
