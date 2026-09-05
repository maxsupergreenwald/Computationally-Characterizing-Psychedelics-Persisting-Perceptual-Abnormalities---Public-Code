"""
master_config.py — Single source of truth for all variable metadata,
pipeline configuration, and label dictionaries.

Replaces: hpc_config.py, variable_labels.py, dv_to_lab_supershort.py,
          dv_to_lab_organized.py, dv_to_lab_dict_full.py,
          dv_to_lab_dict_gen_audience.py, dv_to_lab_dict_concise.py,
          dv_to_lab_supershort_sobp.py, dv_to_lab_organized_hppdmansucript.py

################################################################################
!! CRITICAL: VARIABLE_REGISTRY SCHEMA RULES !!
################################################################################
VARIABLE_REGISTRY maps every analysis variable name to a metadata dict with a
FIXED set of keys. ALL rules below are mandatory:

1. Every variable entry MUST have ALL keys listed in _REGISTRY_SCHEMA. Missing
   keys will cause KeyError in downstream code. Default to np.nan when unknown.

2. **When adding a NEW KEY to the schema:**
   - Add it to _REGISTRY_SCHEMA below
   - Add it to EVERY existing variable entry (np.nan as default)
   - Update this docstring
   - **AI AGENTS: YOU MUST ASK THE USER for the correct value. NEVER silently
     fill in a new key with a guess. Prompt the user explicitly.**

3. When adding a NEW VARIABLE:
   - Add ALL keys from _REGISTRY_SCHEMA (np.nan for unknowns)
   - Set distribution, data_type, need_non_normalized, is_categorical_factor
     correctly — these affect model fitting and MUST NOT be left as np.nan
     for any variable that enters a statistical model.

4. redcap_field_label is the human-readable question text from REDCap.
   redcap_choices is a dict {int: str} mapping response codes to labels (e.g.
   {1: 'Male', 2: 'Female', 3: 'Other'}). Only ~26 of 239 registry variables
   are REDCap source fields (the rest are derived/computed); leave both as
   np.nan for derived variables. redcap_data_dictionary.csv in this directory
   is the authoritative source for raw values.
################################################################################
"""

from pathlib import Path

import numpy as np
from matplotlib.colors import LinearSegmentedColormap, to_hex

# ── Shipped analysis dataframe ───────────────────────────────────────────────
# The wide, participant-level dataframe distributed with this repository lives in
# data/final/ and is named df_public_<date>.csv.  It is already fully prepared:
# every derived column (hppd_binary, subtle, baggot_total, the SDT/metacognition
# block, recalculated VCH hit rates, corrected avg_life_dose, ...) is present.
# Analysis scripts read it directly — there is no data-preparation step to run.
#
# Scripts resolve the actual file with data_prep.most_recent_public_df(), which
# picks the most recently modified match, mirroring how 02_hgf_modeling's
# load_public_wide_df() selects its input.  Only the prefix is declared here.
PUBLIC_DF_PREFIX = 'df_public_'

# Canonical SP-user filter, used wherever an analysis is restricted to
# psychedelic users:  df[df[SP_USER_COL] == SP_USER_VALUE]
# NOT `!= "No"` and NOT `psycheduse_life_nomic > 0` — both give a different N.
SP_USER_COL   = 'psycheduse_yn'
SP_USER_VALUE = 'Yes'

# ── Recruitment / screening dataframe ────────────────────────────────────────
# The consort diagram, ineligibility table, QC-failure tables and the results
# narrative all report recruitment-pipeline counts, which come from the
# screening export rather than from the analysis dataframe.
#
# This is the de-identified recruitment export that ships with the repository.
# It is produced from the private REDCap screening CSV by
# 06_submission/deidentify_recruit_csv.py in the master repository, which keeps
# only the 50 columns the published pipeline reads and replaces every
# identifying column with the derived value the code computes from it:
#   phone_number        -> phone_cc_234
#   continue_date       -> continue_date_passed
#   timestamp_survey_bl -> timesincesurveystart_gt14
#   task_data_prltask   -> task_data_prltask_present
#   qc_notes            -> withheld; categories frozen in
#                          modules/qc_redacted_categories.py
# Resolved relative to the repository root so it works on any machine.
RECRUIT_CSV = str(
    Path(__file__).resolve().parent.parent / 'data' / 'final'
    / 'df_recruit_public_09-03-2026.csv'
)

# ── Reported point estimates ─────────────────────────────────────────────────
# Every brms summary CSV stores both `median` and `mean`; these two constants
# name which one is reported, and callers must say which case they are in:
# point_estimate(row, source=..., mc_integrated=True) for the second.
#
#   MEAN   — everything: regression coefficients, a/b/c' path coefficients, and
#            single-path counterfactual deltas.  Matches brms' `Estimate`, which
#            is itself a posterior mean, and means are exactly additive.
#
#   MEDIAN — NIE / NDE / TE / PMed out of mc_mediation_summary.csv ONLY.  Those
#            posteriors are unstable and their means are uninformative: the MC
#            integration draws mediator values far outside the observed range
#            and the log link exponentiates them, so a handful of draws
#            dominate the mean (in caps_vision x vch_beta, 5 draws of 16,000
#            give a mean of 1.1e7 against a median of 1.28).  Quantile summaries
#            — median, 94% HDI, direction probabilities — are stable throughout.
#            Affects caps_vision models only; bernoulli's logit link bounds
#            E[Y] in [0, 1].
#
# Cost: medians are not additive, so reported NIE + NDE will not sum exactly to
# TE.  That is a property of the summary, not a mediated interaction, and the
# methods section must say so.
#
# Both columns are written to every CSV, so either choice can be revisited from
# results already on disk without refitting.
POINT_ESTIMATE_COL           = 'mean'
MC_EFFECT_POINT_ESTIMATE_COL = 'median'


class PointEstimateColumnMissing(KeyError):
    """Raised when a results CSV lacks the reported point-estimate column.

    Carries an explicit instruction rather than surfacing a bare pandas
    KeyError, because the fix is never "read whichever column is present" —
    substituting the other summary would silently report a different estimand
    than the manuscript claims. See point_estimate().

    Subclasses KeyError so that any existing `except KeyError` around a results
    read still catches it. __str__ is overridden because KeyError inherits
    LookupError's repr-based __str__, which would render the multi-line message
    as one line with literal \n escapes.
    """

    def __str__(self):
        return self.args[0] if self.args else ''


def point_estimate(row, source='', mc_integrated=False):
    """Return the reported point estimate from one results row.

    Parameters
    ----------
    row : pandas.Series | pandas.DataFrame | Mapping
        A single row of a brms summary CSV, or the whole frame (in which case
        only the presence of the column is checked and the column is returned).
    source : str
        Human-readable identifier of the file/model the row came from. It is
        quoted verbatim in the error message, so pass something that lets the
        reader find the file on disk.
    mc_integrated : bool
        True ONLY for rows out of ``mc_mediation_summary.csv`` — the NIE / NDE /
        TE / PMed effects that hpc_mediation.R produces by Monte Carlo
        integration over the mediator's posterior predictive. Those report
        MC_EFFECT_POINT_ESTIMATE_COL; everything else reports
        POINT_ESTIMATE_COL. See the comment on those two constants for why they
        differ. This is an explicit argument rather than something inferred from
        `source` so that a renamed file can never silently switch the estimand.

    Raises
    ------
    PointEstimateColumnMissing
        If the required column is absent. No fallback to the other summary
        column is attempted — that would report a different summary of the
        posterior than the manuscript describes.
    """
    col = MC_EFFECT_POINT_ESTIMATE_COL if mc_integrated else POINT_ESTIMATE_COL
    const_name = 'MC_EFFECT_POINT_ESTIMATE_COL' if mc_integrated else 'POINT_ESTIMATE_COL'

    try:
        keys = row.columns if hasattr(row, 'columns') else row.index
    except AttributeError:              # plain dict / Mapping
        keys = row.keys()

    if col not in keys:
        raise PointEstimateColumnMissing(
            f"column '{col}' not found"
            + (f" in {source}" if source else "")
            + f". Available: {sorted(keys)}.\n"
            "Every summary CSV written by hpc_mediation.R / nonsp_predictors.R "
            "carries both `median` and `mean`. A file missing the reported "
            "column was produced by an older version of those scripts and the "
            "model must be REFIT on the cluster before it can be reported.\n"
            "Not falling back to the other column: the reported estimate is "
            f"the posterior {col} ({const_name} in modules/master_config.py), "
            "and silently substituting the other would report a different "
            "summary of the posterior than the manuscript describes."
        )
    return row[col]


# ── Registry schema (canonical key list) ─────────────────────────────────────
# ALL variable entries in VARIABLE_REGISTRY must have EXACTLY these keys.
_REGISTRY_SCHEMA = [
    'plot_label',               # str | np.nan — short label (from dv_to_lab_supershort)
    'plot_label_verbose',       # str | np.nan — longer label (from dv_to_lab_organized / dv_to_lab_dict_full)
    'plot_label_lay_audience',  # str | np.nan — lay-audience label (from dv_to_lab_dict_gen_audience)
    'distribution',             # str | np.nan — brms family string passed to the R pipeline.
                               # Valid values (nonsp_predictors.R / hpc_mediation.R):
                               #   student_t, gaussian, beta, zero_inflated_beta, gamma,
                               #   zero_negbinomial, bernoulli, negbinomial, lognormal,
                               #   hurdle_negbinomial, lognormal_hierarchial,
                               #   student_t_hierarchial, binomial_hierarchial, hurdle_negbinom_huvary
                               # Special value:
                               #   "ordinal" — brms cumulative() family for ordered categorical
                               #   response variables (i.e., ordinal regression). Implements the
                               #   adjacent-category / proportional-odds model. Both R scripts
                               #   accept "ordinal" (canonical) and "cumulative" (legacy alias)
                               #   as equivalent strings; VARIABLE_REGISTRY always uses "ordinal".
                               #   Reference: Bürkner & Vuorre (2019) Advances in Methods and
                               #   Practices in Psychological Science, 2(1), 77–101.
                               #   https://journals.sagepub.com/doi/10.1177/2515245918823199
    'data_type',                # str | np.nan — 'continuous', 'categorical', 'binary', 'ordinal'
    'need_non_normalized',      # bool | np.nan — True = raw col for mediator, _normalized col for DV
    'inplace_normalized',       # bool | np.nan — True = normalized in-place; same col in both formulas
    'is_categorical_factor',    # bool | np.nan — True = must be as.factor() in R, never Gelman-normalized
    'predictor',                # bool — True if this variable is ever used as a RHS predictor in any model
    'dv',                       # bool — True if this variable is ever used as a model outcome
    'mediator',                 # bool — True if this variable is ever used as a mediator in a causal mediation model
    'redcap_field_label',       # str | np.nan — human-readable field label from REDCap data dictionary
    'redcap_choices',           # str | np.nan — response options from REDCap ("1, Label | 2, Label | ...")
]

# ── Helper to build a registry entry ─────────────────────────────────────────
def _r(plot_label=np.nan, plot_label_verbose=np.nan, plot_label_lay_audience=np.nan,
       distribution=np.nan, data_type=np.nan,
       need_non_normalized=np.nan, inplace_normalized=np.nan,
       is_categorical_factor=np.nan,
       predictor=False, dv=False, mediator=False,
       redcap_field_label=np.nan, redcap_choices=np.nan):
    return {
        'plot_label': plot_label,
        'plot_label_verbose': plot_label_verbose,
        'plot_label_lay_audience': plot_label_lay_audience,
        'distribution': distribution,
        'data_type': data_type,
        'need_non_normalized': need_non_normalized,
        'inplace_normalized': inplace_normalized,
        'is_categorical_factor': is_categorical_factor,
        'predictor': predictor,
        'dv': dv,
        'mediator': mediator,
        'redcap_field_label': redcap_field_label,
        'redcap_choices': redcap_choices,
    }


###############################################################################
# SECTION 1: VARIABLE_REGISTRY
###############################################################################

VARIABLE_REGISTRY = {

    # ── Clinical — High RPW ───────────────────────────────────────────────────
    'phq9_tot': _r(
        plot_label='PHQ9',
        plot_label_verbose='Depression (PHQ9)',
        plot_label_lay_audience='Depressive Symptoms',
        distribution='zero_negbinomial', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'aaq_tot_bl': _r(
        plot_label='AAQ-II',
        plot_label_verbose='INflexibility (AAQ-II)',
        plot_label_lay_audience='Cognitive Inflexibility',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Acceptance and Action Questionnaire Total Score',
        redcap_choices=np.nan,  # computed field (sum of items) — no response options
    ),
    'pdi_total': _r(
        plot_label='PDI Total',
        plot_label_verbose='Delusional Ideation (PDI)',
        plot_label_lay_audience='Delusional Ideation',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'absorption_tot_bl': _r(
        plot_label='ETAS',
        plot_label_verbose='Absorption (ETAS)',
        plot_label_lay_audience='Trait Absorption (Full Score)',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='ETAS total score',
        redcap_choices=np.nan,  # computed field (sum of items) — no response options
    ),
    'lshs_total': _r(
        plot_label='LSHS',
        plot_label_verbose='Trait Hallucinatory Predisposition (LSHS)',
        plot_label_lay_audience='Trait Hallucinatory Predisposition',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'baggot_total': _r(
        plot_label='# PPA Symptoms',
        plot_label_verbose='HPPD Sx Count',
        plot_label_lay_audience='HPPD Symptom Count',
        # hurdle_negbinom_huvary chosen over zero_negbinomial (June 2026) based on
        # DHARMa comparison across 3 SP predictors: zinb showed consistent q25/q50
        # quantile heteroscedasticity flags for psycheduse_life_nomic and avg_life_dose
        # (p < 0.005); hurdle model resolved life_nomic flags entirely (0/3 predictors
        # flagged) and substantially improved avg_life_dose (flags reduced). Estimates
        # and pd values were highly consistent across families. The hurdle parameterization
        # is also conceptually appropriate: HPPD symptom count has a true zero-generating
        # process (no symptoms at all) distinct from the count-generating process.
        distribution='hurdle_negbinom_huvary', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
        redcap_field_label='How many persisting visual symptoms were reported?',
        redcap_choices=np.nan,  # computed field (sum of binary symptom items) — no response options
    ),
    'hppd_sx_count': _r(
        plot_label='# PPA Symptoms',
        plot_label_verbose=np.nan,
        plot_label_lay_audience=np.nan,
        distribution='zero_negbinomial', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'persistvis_time': _r(
        plot_label='PPA Chronicity',
        plot_label_verbose='HPPD Chronicity',
        plot_label_lay_audience='HPPD Chronicity',
        distribution='zero_negbinomial', data_type='ordinal',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
        redcap_field_label='What is the longest you experienced any of these visual effects after a serotonergic psychedelic?',
        redcap_choices={
            1: "I've never experienced any visual effect the day after taking a psychedelic",
            2: '< 1 day',
            3: '1 - 3 days',
            4: '3 days - 1 week',
            5: '1 week - 1 month',
            6: '1 - 6 months',
            7: '6 months - 1 year',
            8: '>1 year',
            9: 'I still experience these effects',
        },
    ),
    'hppd_true_chronicity': _r(
        plot_label='PPA Chronicity',
        plot_label_verbose='HPPD Chronicity',
        plot_label_lay_audience=np.nan,
        distribution='gamma', data_type='ordinal',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'persistvis_duration': _r(
        plot_label='PPA Duration',
        plot_label_verbose='HPPD Sx Duration',
        plot_label_lay_audience='HPPD Symptom Duration',
        distribution='ordinal', data_type='ordinal',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
        redcap_field_label='Were these visual effects constant, present for brief spurts (minutes or seconds), or came and went in longer periods (hours to days)?',
        redcap_choices={
            0: 'I have never had visual effects starting the day after taking a psychedelic',
            1: 'Brief spurts (seconds to minutes)',
            2: 'Longer periods (hours to days)',
            3: 'Constant',
        },
    ),
    'persist_vis_yn': _r(
        plot_label='SP-related PPA History',
        plot_label_verbose='HPPD History (Y/N)',
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'baggot_binary': _r(
        plot_label='SP-related PPA History',
        plot_label_verbose='HPPD Risk',
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'hppd_binary': _r(
        plot_label='SP-related PPA History',
        plot_label_verbose='HPPD History (Y/N)',
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'hppd_6mo_binary': _r(
        plot_label='PPA Duration',
        plot_label_verbose='HPPD Sx Duration',
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'chr_pos_5plus_weighted': _r(
        plot_label='CHR (YALEPRIME)',
        plot_label_verbose='Clinical High Risk for Psychosis (YALEPRIME)',
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    # CAPS variables
    'caps_total': _r(
        plot_label='CAPS',
        plot_label_verbose='Anomalous Percepts in Past Month (CAPS)',
        plot_label_lay_audience='Total Anomalous Percepts',
        distribution='hurdle_negbinom_huvary', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'caps_distress': _r(
        plot_label='CAPS Distress',
        plot_label_verbose='Hallucinatory Distress (CAPS)',
        plot_label_lay_audience='Anomalous Percept-Related Distress',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'caps_frequency': _r(
        plot_label='CAPS Frequency',
        plot_label_verbose='Hallucinatory Frequency (CAPS)',
        plot_label_lay_audience='Anomalous Percepts Frequency',
        distribution='hurdle_negbinom_huvary', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'caps_temporal_lobe': _r(
        plot_label='CAPS Temporal Lobe',
        plot_label_verbose='CAPS Temporal Lobe Experiences',
        plot_label_lay_audience='Temporal Lobe Anomalous Percepts',
        distribution='hurdle_negbinom_huvary', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'caps_flooding': _r(
        plot_label='CAPS Flooding',
        plot_label_verbose='CAPS Flooding Experiences',
        plot_label_lay_audience='Flooding Anomalous Percepts',
        distribution='zero_negbinomial', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'caps_intensity': _r(
        plot_label='CAPS Intensity',
        plot_label_verbose='CAPS Intensity Experiences',
        plot_label_lay_audience='Anomalously Increased Intensity Percepts',
        distribution='hurdle_negbinom_huvary', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'caps_intensity_frequency': _r(
        plot_label='CAPS Intensity Frequency',
        plot_label_verbose='CAPS Intensity Frequency',
        plot_label_lay_audience=np.nan,
        distribution='zero_negbinomial', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'caps_intensity_maximum_frequency': _r(
        plot_label='CAPS Intensity Frequency (Max)',
        plot_label_verbose='CAPS Intensity Frequency (Maximum)',
        plot_label_lay_audience=np.nan,
        distribution='hurdle_negbinom_huvary', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'caps_vision': _r(
        plot_label='CAPS Visual Items Endorsed',
        # 'hurdle_negbinom_huvary' = hurdle_negbinomial() with a VARYING hu submodel:
        # the predictor and all covariates affect BOTH the hurdle probability (whether
        # any anomalous percepts occur) AND the count given occurrence.  The
        # nonsp-predictor pipeline uses the same family as the mediation pipeline.
        plot_label_verbose='CAPS Vision Experiences',
        plot_label_lay_audience='Anomalous Visual Percepts',
        distribution='hurdle_negbinom_huvary', data_type='continuous',
        need_non_normalized=True, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'caps_vision_formed': _r(
        plot_label='CAPS Vision (Halluc)',
        plot_label_verbose='CAPS Vision (Hallucinations)',
        plot_label_lay_audience=np.nan,
        distribution='hurdle_negbinom_huvary', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'caps_vision_bottomup': _r(
        plot_label='CAPS Vision (Elementary)',
        plot_label_verbose='CAPS Vision (Elementary)',
        plot_label_lay_audience=np.nan,
        distribution='hurdle_negbinom_huvary', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'caps_vision_frequency': _r(
        plot_label='CAPS Vision Frequency',
        plot_label_verbose='CAPS Vision Frequency',
        plot_label_lay_audience=np.nan,
        distribution='zero_negbinomial', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'caps_vision_maximum_frequency': _r(
        plot_label='CAPS Vision Frequency (Max)',
        plot_label_verbose='CAPS Vision Frequency (Maximum)',
        plot_label_lay_audience=np.nan,
        distribution='hurdle_negbinom_huvary', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'caps_vision_intrusiveness': _r(
        plot_label='CAPS Vision Intrusiveness',
        plot_label_verbose='CAPS Vision Intrusiveness',
        plot_label_lay_audience=np.nan,
        distribution='zero_negbinomial', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'caps_audition': _r(
        plot_label='CAPS Audition',
        plot_label_verbose='CAPS Audition Experiences',
        plot_label_lay_audience='Anomalous Auditory Percepts',
        distribution='zero_negbinomial', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'caps_tactile': _r(
        plot_label='CAPS Tactile',
        plot_label_verbose='CAPS Tactile Experiences',
        plot_label_lay_audience='Anomalous Tactile Percepts',
        distribution='zero_negbinomial', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'caps_taste': _r(
        plot_label='CAPS Taste',
        plot_label_verbose='CAPS Taste Experiences',
        plot_label_lay_audience='Anomalous Gustatory Percepts',
        distribution='zero_negbinomial', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'caps_maximum_frequency': _r(
        plot_label='CAPS Frequency (Max)',
        plot_label_verbose='CAPS Maximum Frequency',
        plot_label_lay_audience=np.nan,
        distribution='hurdle_negbinom_huvary', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'caps_vision_max_weighted_frequency': _r(
        plot_label='CAPS Vision Freq (Max Wt)',
        plot_label_verbose='CAPS Vision Max Weighted Frequency',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'caps_intensity_max_weighted_frequency': _r(
        plot_label='CAPS Intensity Freq (Max Wt)',
        plot_label_verbose='CAPS Intensity Max Weighted Frequency',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'caps_max_weighted_frequency': _r(
        plot_label='CAPS Freq (Max Wt)',
        plot_label_verbose='CAPS Max Weighted Frequency',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=True, mediator=False,
    ),
    'colors': _r(
        plot_label='Intense colors',
        plot_label_verbose='Brighter or more intense colors',
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Brighter or more intense colors',
        redcap_choices={0: 'No', 1: 'Yes'},
    ),
    'moving': _r(
        plot_label='Stationary things moving',
        plot_label_verbose='Stationary things appearing to move, breathe, grow, or shrink',
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Stationary things appearing to move, breathe, grow, or shrink',
        redcap_choices={0: 'No', 1: 'Yes'},
    ),
    'subtle': _r(
        plot_label='Subtle',
        plot_label_verbose='Subtle',
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),

    # ── PPA (Prior Psychedelic-related Anomalies) — item-level ───────────────
    # Verbal labels are exact Field Label values from REDCap data dictionary
    # (PsychedelicsPerceptionFollowup_DataDictionary_2026-06-05.csv).
    'grids': _r(
        plot_label='Grid patterns/distortion',
        plot_label_verbose='Distortion, movement, or patterns in grids, gratings, or closely spaced lines',
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Distortion, movement, or patterns in grids, gratings, or closely spaced lines',
        redcap_choices={0: 'No', 1: 'Yes'},
    ),
    'different': _r(
        plot_label='Things look different',
        plot_label_verbose='Things just looked/seemed different',
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Things just looked/seemed different',
        redcap_choices={0: 'No', 1: 'Yes'},
    ),
    'oscillating': _r(
        plot_label='Oscillating light intensity',
        plot_label_verbose='Increased intensity in oscillating or flashing lights (TV, light bulbs, etc.)',
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Increased intensity in oscillating or flashing lights (TV, light bulbs, etc.)',
        redcap_choices={0: 'No', 1: 'Yes'},
    ),
    'halos': _r(
        plot_label='Halos/auras',
        plot_label_verbose='Halos or auras around things',
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Halos or auras around things',
        redcap_choices={0: 'No', 1: 'Yes'},
    ),
    'still': _r(
        plot_label='Moving objects frozen',
        plot_label_verbose='Moving objects appear to not be moving',
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Moving objects appear to not be moving',
        redcap_choices={0: 'No', 1: 'Yes'},
    ),
    'trails': _r(
        plot_label='Object afterimage trails',
        plot_label_verbose='Afterimages left behind moving objects',
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Afterimages left behind moving objects',
        redcap_choices={0: 'No', 1: 'Yes'},
    ),
    'objects': _r(
        plot_label='Seeing absent objects',
        plot_label_verbose="Seeing objects that aren't really there",
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label="Seeing objects that aren't really there",
        redcap_choices={0: 'No', 1: 'Yes'},
    ),
    'pattern_open': _r(
        plot_label='Patterns (eyes open)',
        plot_label_verbose="Seeing patterns or textures that aren't there with eyes open",
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label="Seeing patterns or textures that aren't there with eyes open",
        redcap_choices={0: 'No', 1: 'Yes'},
    ),
    'pattern_closed': _r(
        plot_label='Patterns (eyes closed)',
        plot_label_verbose='Seeing patterns or textures with eyes closed',
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Seeing patterns or textures with eyes closed',
        redcap_choices={0: 'No', 1: 'Yes'},
    ),
    'attention': _r(
        plot_label='Heightened attention',
        plot_label_verbose='Noticing more things in your surrounding',
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Noticing more things in your surrounding',
        redcap_choices={0: 'No', 1: 'Yes'},
    ),

    # ── CAPS item-level endorsement (caps_bl_1 – caps_bl_32) ─────────────────
    # Verbal labels are exact Field Label values from REDCap data dictionary
    # (PsychedelicsPerceptionFollowup_DataDictionary_2026-06-05.csv).
    # Items 4, 5, 19, 20, 23 have no direct field-name entry in REDCap (only split
    # sub-items exist there); verbose labels below are composite summaries of the
    # split descriptions. redcap_field_label=np.nan for those five items.
    'caps_bl_1': _r(
        plot_label='Sounds louder',
        plot_label_verbose='Do you ever notice that sounds are much louder than they normally would be?',
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Do you ever notice that sounds are much louder than they normally would be?',
        redcap_choices={0: 'No', 1: 'Yes'},
    ),
    'caps_bl_4': _r(
        plot_label='Shapes, Lights, Colors',
        plot_label_verbose='Do you ever see shapes, lights, or colors even though there is nothing really there?',
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,  # composite of caps_bl_4_split_1/2/3 in REDCap
        redcap_choices={0: 'No', 1: 'Yes'},
        
    ),
    'caps_bl_19': _r(
        plot_label='Shape/Size/Color Changes',
        plot_label_verbose='Do you ever find that the appearance of things or people seems to change in a puzzling way (in shape, size, color, or other)?',
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,  # composite of caps_bl_19_split_shape/size/color/other in REDCap
        redcap_choices={0: 'No', 1: 'Yes'},
    ),
    'caps_bl_22': _r(
        plot_label='Face looks different',
        plot_label_verbose='Do you ever look in the mirror and think that your face seems different from usual?',
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Do you ever look in the mirror and think that your face seems different from usual?',
        redcap_choices={0: 'No', 1: 'Yes'},
    ),
    'caps_bl_23': _r(
        plot_label='Brighter Colors/Light',
        plot_label_verbose='Do you ever have days where lights seem brighter or colors seem more intense than usual?',
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,  # composite of caps_bl_23_split_1/2 in REDCap
        redcap_choices={0: 'No', 1: 'Yes'},
    ),
    'caps_bl_26': _r(
        plot_label='Things look abnormal',
        plot_label_verbose='Do you ever think that everyday things look abnormal to you?',
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Do you ever think that everyday things look abnormal to you?',
        redcap_choices={0: 'No', 1: 'Yes'},
    ),
    'caps_bl_31': _r(
        plot_label='See things others cannot',
        plot_label_verbose='Do you ever see things that other people cannot?',
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Do you ever see things that other people cannot?',
        redcap_choices={0: 'No', 1: 'Yes'},
    ),
    'caps_bl_32': _r(
        plot_label='Hear unshared sounds',
        plot_label_verbose="Do you ever hear sound or music that people near you don't hear?",
        plot_label_lay_audience=np.nan,
        distribution='bernoulli', data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label="Do you ever hear sound or music that people near you don't hear?",
        redcap_choices={0: 'No', 1: 'Yes'},
    ),

    # ── CAPS item-level frequency (caps_bl_1c – caps_bl_32c) ─────────────────
    # Frequency rating for endorsed items only. Values: 0 = not endorsed (structural
    # zero); 1–5 = frequency when endorsed (1=hardly, 2=not often, 3=sometimes,
    # 4=regularly, 5=all the time). Distribution: hurdle_negbinom_huvary.
    # Items 4c, 5c, 19c, 20c, 23c are composites matching their parent items above.
    'caps_bl_1c': _r(
        plot_label='Sounds louder (freq)',
        plot_label_verbose='Do you ever notice that sounds are much louder than they normally would be? — frequency',
        plot_label_lay_audience=np.nan,
        distribution='hurdle_negbinom_huvary', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='How often does this happen?',
        redcap_choices={0: 'Not endorsed', 1: 'Hardly happens at all', 2: 'Does not happen often', 3: 'Happens sometimes', 4: 'Happens regularly', 5: 'Happens all the time'},
    ),
    'caps_bl_32c': _r(
        plot_label='Hear unshared sounds (freq)',
        plot_label_verbose="Do you ever hear sound or music that people near you don't hear? — frequency",
        plot_label_lay_audience=np.nan,
        distribution='hurdle_negbinom_huvary', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='How often does this happen?',
        redcap_choices={0: 'Not endorsed', 1: 'Hardly happens at all', 2: 'Does not happen often', 3: 'Happens sometimes', 4: 'Happens regularly', 5: 'Happens all the time'},
    ),

    # ── Clinical — Low RPW ────────────────────────────────────────────────────
    'ffmq_total_bl': _r(
        plot_label='FFMQ',
        plot_label_verbose='Trait Mindfulness (FFMQ)',
        plot_label_lay_audience='Trait Mindfulness',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'asi_tot': _r(
        plot_label='ASI',
        plot_label_verbose='Aberrant Salience (ASI)',
        plot_label_lay_audience='Aberrant Salience',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),

    # ── PRL — High RPW ────────────────────────────────────────────────────────
    'perseverative_errors': _r(
        plot_label='Perseverative Errors',
        plot_label_verbose='PRL Perseverative Errors',
        plot_label_lay_audience='Rigid Learning (PRL)',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),

    # ── PRL — Low RPW ─────────────────────────────────────────────────────────
    'regressive_errors': _r(
        plot_label='JTC Errors',
        plot_label_verbose='PRL JTC Errors',
        plot_label_lay_audience='Hasty Learning (PRL)',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'beta': _r(
        plot_label='β',
        plot_label_verbose='PRL Decision Precision (β)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'omega': _r(
        plot_label='ω',
        plot_label_verbose='PRL Contingency Belief Evolution Rate (ω)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'pwPE_negative': _r(
        plot_label='pwPEs (-)',
        plot_label_verbose='PRL pwPEs (-)',
        plot_label_lay_audience=np.nan,
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'pwPE_positive': _r(
        plot_label='pwPEs (+)',
        plot_label_verbose='PRL pwPEs (+)',
        plot_label_lay_audience=np.nan,
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'pwPE_negative_0': _r(
        plot_label='pwPEs (-) No Reward',
        plot_label_verbose='PRL pwPEs (-) No Reward',
        plot_label_lay_audience=np.nan,
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'pwPE_positive_0': _r(
        plot_label='pwPEs (+) No Reward',
        plot_label_verbose='PRL pwPEs (+) No Reward',
        plot_label_lay_audience=np.nan,
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'pwPE_negative_1': _r(
        plot_label='pwPEs (-) Reward',
        plot_label_verbose='PRL pwPEs (-) Reward',
        plot_label_lay_audience=np.nan,
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'pwPE_positive_1': _r(
        plot_label='pwPEs (+) Reward',
        plot_label_verbose='PRL pwPEs (+) Reward',
        plot_label_lay_audience=np.nan,
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'pwPE_negative_correct': _r(
        plot_label='pwPEs (-) Correct',
        plot_label_verbose='PRL pwPEs (-) Correct Choice',
        plot_label_lay_audience=np.nan,
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'pwPE_positive_correct': _r(
        plot_label='pwPEs (+) Correct',
        plot_label_verbose='PRL pwPEs (+) Correct Choice',
        plot_label_lay_audience=np.nan,
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'pwPE_negative_incorrect': _r(
        plot_label='pwPEs (-) Incorrect',
        plot_label_verbose='PRL pwPEs (-) Incorrect Choice',
        plot_label_lay_audience=np.nan,
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'pwPE_positive_incorrect': _r(
        plot_label='pwPEs (+) Incorrect',
        plot_label_verbose='PRL pwPEs (+) Incorrect Choice',
        plot_label_lay_audience=np.nan,
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),

    # ── ACH — High RPW ────────────────────────────────────────────────────────
    'ach_nu': _r(
        plot_label='Prior Weighting (ν)',
        plot_label_verbose='ACH Prior:Input Weighting',
        plot_label_lay_audience='RPW (Auditory)',
        distribution='gamma', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'ach_xprob_0': _r(
        plot_label='xprob (No Target)',
        plot_label_verbose='Contingency belief (No Target)',
        plot_label_lay_audience='Pavlovian Conditioning Credence (No Target)',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'ach_xprob_75': _r(
        plot_label='xprob (75%)',
        plot_label_verbose='Contingency belief (75% Trials)',
        plot_label_lay_audience='Pavlovian Conditioning Credence (75% Trials)',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'ach_xprob_present': _r(
        plot_label='xprob (Target)',
        plot_label_verbose='Contingency belief (Target Present)',
        plot_label_lay_audience='Pavlovian Conditioning Credence (Target Present)',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'total_ch_trials': _r(
        plot_label='Trials',
        plot_label_verbose='ACH Trials',
        plot_label_lay_audience='Auditory Conditioned Hallucinations',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'ach_bl_yes_0': _r(
        plot_label='ACH Rate',
        plot_label_verbose='ACH Rate',
        plot_label_lay_audience=np.nan,
        distribution='zero_inflated_beta', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'hits_75': _r(
        plot_label='Hits (75%)',
        plot_label_verbose='ACH Hits (75%)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'ach_bl_yes_75': _r(
        plot_label='Hit Rate (75% Contrast Trials)',
        plot_label_verbose='Hit Rate (75% Contrast Trials)',
        plot_label_lay_audience=np.nan,
        distribution='beta', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'ach_bl_yes_25': _r(
        plot_label='Hit Rate (25% Contrast Trials)',
        plot_label_verbose='Hit Rate (25% Contrast Trials)',
        plot_label_lay_audience=np.nan,
        distribution='beta', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'ach_bl_yes_50': _r(
        plot_label='Hit Rate (50% Contrast Trials)',
        plot_label_verbose='Hit Rate (50% Contrast Trials)',
        plot_label_lay_audience=np.nan,
        distribution='beta', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),

    'ach_pwPE_negative': _r(
        plot_label='pwPEs (-)',
        plot_label_verbose='ACH pwPEs (-)',
        plot_label_lay_audience='Auditory Belief Updates (-)',
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'ach_pwPE_positive': _r(
        plot_label='pwPEs (+)',
        plot_label_verbose='ACH pwPEs (+)',
        plot_label_lay_audience='Auditory Belief Updates (+)',
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'ach_pwPE_negative_0': _r(
        plot_label='pwPEs (-) No Target',
        plot_label_verbose='ACH pwPEs (-) No Target',
        plot_label_lay_audience='Auditory Belief Updates (-) No Target',
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'ach_pwPE_positive_0': _r(
        plot_label='pwPEs (+) No Target',
        plot_label_verbose='ACH pwPEs (+) No Target',
        plot_label_lay_audience='Auditory Belief Updates (+) No Target',
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'ach_pwPE_negative_75': _r(
        plot_label='pwPEs (-) 75%',
        plot_label_verbose='ACH pwPEs (-) 75% Trials',
        plot_label_lay_audience='Auditory Belief Updates (-) 75% Trials',
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'ach_pwPE_positive_75': _r(
        plot_label='pwPEs (+) 75%',
        plot_label_verbose='ACH pwPEs (+) 75% Trials',
        plot_label_lay_audience='Auditory Belief Updates (+) 75% Trials',
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'ach_pwPE_negative_present': _r(
        plot_label='pwPEs (-) Target',
        plot_label_verbose='ACH pwPEs (-) Target Present',
        plot_label_lay_audience='Auditory Belief Updates (-) Target Present',
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'ach_pwPE_positive_present': _r(
        plot_label='pwPEs (+) Target',
        plot_label_verbose='ACH pwPEs (+) Target Present',
        plot_label_lay_audience='Auditory Belief Updates (+) Target Present',
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),


    # ── VCH — High Decision Noise ────────────────────────────────────────────────────────
    'd_prime_overall': _r(
        plot_label="Sensitivity (d')",
        plot_label_verbose="Sensitivity (d')",
        plot_label_lay_audience='Sensitivity',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=True, is_categorical_factor=False,
        predictor=True, dv=True, mediator=True,
    ),
    'criterion_overall': _r(
        plot_label="Criterion (c)",
        plot_label_verbose="Criterion (c)",
        plot_label_lay_audience='Criterion',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=True, is_categorical_factor=False,
        predictor=True, dv=True, mediator=True,
    ),
    'mean_conf_fas': _r(
        plot_label="VCH Confidence",
        plot_label_verbose="VCH Confidence",
        plot_label_lay_audience='VCH Confidence',
        distribution='gamma', data_type='continuous',
        need_non_normalized=True, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=True,
    ),
    'vch_hit_rate': _r(
        plot_label='VCH Hit Rate',
        plot_label_verbose='VCH Hit Rate (signal trials)',
        plot_label_lay_audience='VCH Hit Rate',
        distribution='zero_inflated_beta', data_type='continuous',
        need_non_normalized=True, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=True,
    ),

    # ── VCH — High RPW ────────────────────────────────────────────────────────
    'vch_nu': _r(
        plot_label='Prior Weighting (ν)',
        plot_label_verbose='VCH Prior Weighting (ν)',
        plot_label_lay_audience='VCH Prior:Input Weighting',
        distribution='gamma', data_type='continuous',
        need_non_normalized=True, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=True,
    ),
    'vch_nu_log': _r(
        plot_label=np.nan,
        plot_label_verbose='VCH Prior:Input (log)',
        plot_label_lay_audience='VCH Prior:Input (log)',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_xprob_median': _r(
        plot_label='Contingency Belief',
        plot_label_verbose='VCH Contingency Belief (xprob)',
        plot_label_lay_audience='VCH Pavlovian Conditioning Credence (xprob)',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_belief_median': _r(
        plot_label='Belief',
        plot_label_verbose='VCH Belief in Target Presence',
        plot_label_lay_audience='VCH Belief in Target Presence',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_xprob_0': _r(
        plot_label='Contingency Belief (No Target)',
        plot_label_verbose='Contingency belief (No Target)',
        plot_label_lay_audience='Pavlovian Conditioning Credence (No Target)',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_xprob_75': _r(
        plot_label='Contingency Belief (75%)',
        plot_label_verbose='Contingency belief (75% Trials)',
        plot_label_lay_audience='Pavlovian Conditioning Credence (75% Trials)',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_xprob_present': _r(
        plot_label='Contingency Belief (Target)',
        plot_label_verbose='Contingency belief (Target Present)',
        plot_label_lay_audience='Pavlovian Conditioning Credence (Target Present)',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'total_vch_trials': _r(
        plot_label='Trials',
        plot_label_verbose='VCH Trials',
        plot_label_lay_audience='Visual Conditioned Hallucinations',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_bl_yes_0': _r(
        plot_label='VCH Rate',
        plot_label_verbose='VCH Rate',
        plot_label_lay_audience=np.nan,
        distribution='zero_inflated_beta', data_type='continuous',
        need_non_normalized=True, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=True,
    ),
    'vch_bl_yes_0_normalized': _r(
        plot_label='VCH Rate',
        plot_label_verbose='VCH Rate',
        plot_label_lay_audience=np.nan,
        distribution='zero_inflated_beta', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=False, mediator=True,
    ),
    'vch_hits_75': _r(
        plot_label='Hits (75%)',
        plot_label_verbose='VCH Hits (75%)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_bl_yes_75': _r(
        plot_label='Hit Rate (75% Contrast Trials)',
        plot_label_verbose='Hit Rate (75% Contrast Trials)',
        plot_label_lay_audience=np.nan,
        distribution='beta', data_type='continuous',
        need_non_normalized=True, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=True,
    ),
    'vch_bl_yes_25': _r(
        plot_label='Hit Rate (25% Contrast Trials)',
        plot_label_verbose='Hit Rate (25% Contrast Trials)',
        plot_label_lay_audience=np.nan,
        distribution='beta', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_bl_yes_50': _r(
        plot_label='Hit Rate (50% Contrast Trials)',
        plot_label_verbose='Hit Rate (50% Contrast Trials)',
        plot_label_lay_audience=np.nan,
        distribution='beta', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_belief': _r(
        plot_label='Belief',
        plot_label_verbose='VCH Belief in Target Presence',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    # VCH state keys (top-level in dv_to_lab_supershort['vch'])
    'xprob': _r(
        plot_label='P(Target | Cue)',
        plot_label_verbose='Contingency Belief',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'xbin_pred': _r(
        plot_label='Trialwise P(Target | Cue)',
        plot_label_verbose="Contingency Prediction",
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'belief': _r(
        plot_label='Prior-influenced Belief',
        plot_label_verbose=np.nan,
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'pwPE': _r(
        plot_label='pwPE',
        plot_label_verbose=np.nan,
        plot_label_lay_audience=np.nan,
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),

    # ── VCH — Low RPW ─────────────────────────────────────────────────────────
    'vch_omega': _r(
        plot_label='Contingency Belief Evolution Rate (ω)',
        plot_label_verbose='VCH Belief Evolution Rate (ω)',
        plot_label_lay_audience='Visual Contingency Belief Evolution Rate',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=True, is_categorical_factor=False,
        predictor=True, dv=True, mediator=True,
    ),
    'vch_pwPE_median': _r(
        plot_label='pwPE',
        plot_label_verbose='pwPE (Visual)',
        plot_label_lay_audience='Visual Belief Updates',
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_xprob_change': _r(
        plot_label='Δxprob',
        plot_label_verbose='VCH Contingency Belief Change (Δxprob)',
        plot_label_lay_audience='Change in Pavlovian Conditioning Credence',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_belief_change': _r(
        plot_label='ΔBelief',
        plot_label_verbose='VCH Target Presence Belief Change (ΔBelief)',
        plot_label_lay_audience='Visual Detection Credence Change (ΔBelief)',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_beta': _r(
        plot_label='Decision Precision (β)',
        plot_label_verbose='VCH Decision Precision (β)',
        plot_label_lay_audience='Visual Error',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=True, is_categorical_factor=False,
        predictor=True, dv=True, mediator=True,
    ),
    'vch_threshold': _r(
        plot_label='75% Threshold',
        plot_label_verbose='Visual 75% Threshold',
        plot_label_lay_audience='Visual Threshold',
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=True, is_categorical_factor=False,
        predictor=True, dv=True, mediator=True,
    ),
    'vch_pwPE': _r(
        plot_label='pwPEs',
        plot_label_verbose='VCH pwPEs',
        plot_label_lay_audience=np.nan,
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_pwPE_negative': _r(
        plot_label='pwPEs (-)',
        plot_label_verbose='VCH pwPEs (-)',
        plot_label_lay_audience='Visual Belief Updates (-)',
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_pwPE_positive': _r(
        plot_label='pwPEs (+)',
        plot_label_verbose='VCH pwPEs (+)',
        plot_label_lay_audience='Visual Belief Updates (+)',
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_pwPE_negative_0': _r(
        plot_label='pwPEs (-) No Target',
        plot_label_verbose='VCH pwPEs (-) No Target',
        plot_label_lay_audience='Visual Belief Updates (-) No Target',
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_pwPE_positive_0': _r(
        plot_label='pwPEs (+) No Target',
        plot_label_verbose='VCH pwPEs (+) No Target',
        plot_label_lay_audience='Visual Belief Updates (+) No Target',
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_pwPE_negative_75': _r(
        plot_label='pwPEs (-) 75%',
        plot_label_verbose='VCH pwPEs (-) 75% Trials',
        plot_label_lay_audience='Visual Belief Updates (-) 75% Trials',
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_pwPE_positive_75': _r(
        plot_label='pwPEs (+) 75%',
        plot_label_verbose='VCH pwPEs (+) 75% Trials',
        plot_label_lay_audience='Visual Belief Updates (+) 75% Trials',
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_pwPE_negative_present': _r(
        plot_label='pwPEs (-) Target',
        plot_label_verbose='VCH pwPEs (-) Target Present',
        plot_label_lay_audience='Visual Belief Updates (-) Target Present',
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_pwPE_positive_present': _r(
        plot_label='pwPEs (+) Target',
        plot_label_verbose='VCH pwPEs (+) Target Present',
        plot_label_lay_audience='Visual Belief Updates (+) Target Present',
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    # VCH pwPE bias variants
    'vch_pwPE_bias_0_median': _r(
        plot_label='pwPE Bias (0%)',
        plot_label_verbose='VCH pwPE Bias (0%)',
        plot_label_lay_audience=np.nan,
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_pwPE_bias_75_median': _r(
        plot_label='pwPE Bias (75%)',
        plot_label_verbose='VCH pwPE Bias (75%)',
        plot_label_lay_audience=np.nan,
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_pwPE_bias_median': _r(
        plot_label='pwPE Bias',
        plot_label_verbose='VCH pwPE Bias',
        plot_label_lay_audience=np.nan,
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_pwPE_ch_median': _r(
        plot_label='pwPE (CH)',
        plot_label_verbose='VCH pwPE (CH)',
        plot_label_lay_audience=np.nan,
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_pwPE_negative_0_median': _r(
        plot_label='pwPEs (-) 0%',
        plot_label_verbose='VCH pwPEs (-) 0% Trials',
        plot_label_lay_audience=np.nan,
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_pwPE_positive_0_median': _r(
        plot_label='pwPEs (+) 0%',
        plot_label_verbose='VCH pwPEs (+) 0% Trials',
        plot_label_lay_audience=np.nan,
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_pwPE_negative_75_median': _r(
        plot_label='pwPEs (-) 75%',
        plot_label_verbose='VCH pwPEs (-) 75% Trials',
        plot_label_lay_audience=np.nan,
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_pwPE_positive_75_median': _r(
        plot_label='pwPEs (+) 75%',
        plot_label_verbose='VCH pwPEs (+) 75% Trials',
        plot_label_lay_audience=np.nan,
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    # 3-level Julia HGF parameters
    'vch_beta_3lev': _r(
        plot_label='Decision Precision (β, 3L)',
        plot_label_verbose='VCH Decision Precision (β, Julia 3-level)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=True,
    ),
    'vch_nu_3lev': _r(
        plot_label='Prior weighting (ν, 3L)',
        plot_label_verbose='VCH Prior Weighting (ν, Julia 3-level)',
        plot_label_lay_audience=np.nan,
        distribution='gamma', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=True,
    ),
    'vch_nu_3lev_log': _r(
        plot_label='Prior weighting log (ν, 3L)',
        plot_label_verbose='VCH Prior Weighting log (ν, Julia 3-level)',
        plot_label_lay_audience=np.nan,
        distribution='gamma', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_omega_3lev': _r(
        plot_label='Contingency Belief Evolution Rate (ω₂, 3L)',
        plot_label_verbose='VCH Belief Evolution Rate (ω₂, Julia 3-level)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_omega3_3lev': _r(
        plot_label='Volatility (ω₃, 3L)',
        plot_label_verbose='VCH Volatility (ω₃, Julia 3-level)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    # 3-level MATLAB HGF parameters
    'vch_short_psychedelic_bl_beta': _r(
        plot_label='Decision Precision (β, MATLAB)',
        plot_label_verbose='VCH Decision Precision (β, MATLAB)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=True,
    ),
    'vch_short_psychedelic_bl_nu': _r(
        plot_label='Prior weighting (ν, MATLAB)',
        plot_label_verbose='VCH Prior Weighting (ν, MATLAB)',
        plot_label_lay_audience=np.nan,
        distribution='gamma', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=True,
    ),
    'vch_short_psychedelic_bl_omega2': _r(
        plot_label='Contingency Belief Evolution Rate (ω₂, MATLAB)',
        plot_label_verbose='VCH Belief Evolution Rate (ω₂, MATLAB)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_short_psychedelic_bl_omega3': _r(
        plot_label='Volatility (ω₃, MATLAB)',
        plot_label_verbose='VCH Volatility (ω₃, MATLAB)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    # Nominal HGF condition
    'vch_nu_nominal': _r(
        plot_label='Prior Weighting (ν)',
        plot_label_verbose='VCH Prior Weighting — nominal (ν)',
        plot_label_lay_audience=np.nan,
        distribution='gamma', data_type='continuous',
        need_non_normalized=True, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=True,
    ),
    'vch_beta_nominal': _r(
        plot_label='Decision Precision (β)',
        plot_label_verbose='VCH Decision Precision (β)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=True, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_omega_nominal': _r(
        plot_label='Contingency Belief Evolution Rate — nominal (ω)',
        plot_label_verbose='VCH Belief Evolution Rate — nominal (ω)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=True, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    # Average HGF condition (per-variable averages from a separate analysis)
    'vch_nu_avg': _r(
        plot_label='Prior weighting — avg (ν)',
        plot_label_verbose='VCH Prior Weighting — average (ν)',
        plot_label_lay_audience=np.nan,
        distribution='gamma', data_type='continuous',
        need_non_normalized=True, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=True,
    ),
    'vch_beta_avg': _r(
        plot_label='Decision Precision — avg (β)',
        plot_label_verbose='VCH Decision Precision — average (β)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=True, is_categorical_factor=False,
        predictor=True, dv=True, mediator=True,
    ),
    'vch_omega_avg': _r(
        plot_label='Contingency Belief Evolution Rate — avg (ω)',
        plot_label_verbose='VCH Belief Evolution Rate — average (ω)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=True, is_categorical_factor=False,
        predictor=True, dv=True, mediator=True,
    ),
    # xprob/belief nominal variants
    'vch_xprob_median_nominal': _r(
        plot_label='Contingency Belief — nominal',
        plot_label_verbose='VCH Contingency Belief — nominal (xprob)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_pwPE_median_nominal': _r(
        plot_label='pwPE — nominal',
        plot_label_verbose='VCH pwPE — nominal',
        plot_label_lay_audience=np.nan,
        distribution='lognormal', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    # xprob/belief median variants
    'vch_xprob_block_1': _r(
        plot_label='xprob (Block 1)',
        plot_label_verbose='VCH xprob (Block 1)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_xprob_block_12': _r(
        plot_label='xprob (Block 12)',
        plot_label_verbose='VCH xprob (Block 12)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=True, mediator=False,
    ),
    'vch_belief_block_1': _r(
        plot_label='Belief (Block 1)',
        plot_label_verbose='VCH Belief (Block 1)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'vch_belief_block_12': _r(
        plot_label='Belief (Block 12)',
        plot_label_verbose='VCH Belief (Block 12)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),

    # ── IVs / Predictors ──────────────────────────────────────────────────────
    'vasdose_bl': _r(
        plot_label='Last Dose (VAS)',
        plot_label_verbose='Last SP Dose (Subjective VAS)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Now based on your own subjective rating of dose strength, how strong would you consider your last psychedelic dose to be, using this slider scale on a scale of 0 to 100?',
        redcap_choices=np.nan,  # continuous VAS slider (0–100) — no discrete response options
    ),
    'motivation': _r(
        plot_label='Reason for Use',
        plot_label_verbose='Reason for SP use',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=True,
        predictor=True, dv=False, mediator=False,
        redcap_field_label='Do you mainly use serotonergic psychedelics for recreational, therapeutic, or spiritual purposes?',
        redcap_choices={
            1: 'Recreational',
            2: 'Therapeutic',
            3: 'Both (neither more than the other)',
            4: 'I do not use psychedelics',
            5: 'Spiritual or Religious',
            6: 'All Three (no one more than the others)',
        },
    ),
    'psychedelic_age': _r(
        plot_label='Age at First Use (Years)',
        plot_label_verbose='Age at First SP Use',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=False, mediator=False,
        redcap_field_label='At what age did you first use a serotonergic psychedelic?',
        redcap_choices=np.nan,  # free-entry numeric — no discrete response options
    ),
    'psycheduse_life_nomic': _r(
        plot_label='Lifetime SP Uses (Count)',
        plot_label_verbose='Lifetime SP Uses',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=False, mediator=False,
        redcap_field_label='NOT including microdoses(!) -- How many times in your life have you used serotonergic psychedelics?',
        redcap_choices=np.nan,  # free-entry numeric — no discrete response options
    ),
    'avg_life_dose': _r(
        plot_label='Avg. Dose (LSD μg eq.)',
        plot_label_verbose='Average Lifetime SP Dose (LSD μg equivalents)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=False, mediator=False,
    ),
    'life_exposure': _r(
        plot_label='Total Exposure (LSD μg eq.)',
        plot_label_verbose='Total SP Exposure (LSD μg equivalents)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'psyched_lastuse_dose': _r(
        plot_label='Last Dose (LSD μg eq.)',
        plot_label_verbose='Last SP Dose (5 pt rating)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='What subjective strength would you consider your last serotonergic psychedelic dose to be?',
        redcap_choices={
            1: 'Microdose (20-30ug; barely detectable or "threshold" drug effect)',
            2: 'Low dose (40-60ug; detectable but very mild drug effect; hallucinations barely present)',
            3: 'Medium dose (90-250ug; full psychedelic effect; many hallucinations)',
            4: 'Heavy dose (300-500ug; strong psychedelic effect)',
            5: 'Very heavy dose (700ug+; out-of-body experiences and beyond)',
        },
    ),
    'psych_dayslastuse_nomicro': _r(
        plot_label='Days Since Last Use',
        plot_label_verbose='Days Since Last SP Use',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='About how many days has it been since you last used a serotonergic psychedelic NOT including microdoses?',
        redcap_choices=np.nan,  # free-entry numeric — no discrete response options
    ),
    'psycheduse_recency_cutoff_hppd': _r(
        plot_label='Recency of last SP (ordinal)',
        plot_label_verbose='Recency of last SP (ordinal)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='ordinal',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'psychedelic_use_PC1': _r(
        plot_label='SP Use PC1',
        plot_label_verbose='SP Use Primary Component',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=False, mediator=False,
    ),
    'psychedelic_rank_use_PC1': _r(
        plot_label='SP Use PC1',
        plot_label_verbose='SP Use Primary Component',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=False, mediator=False,
    ),
    'age_v2': _r(
        plot_label='Age',
        plot_label_verbose='Current Age',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=False, mediator=False,
        redcap_field_label='What is your age? (in years)',
        redcap_choices={
            17: 'under 18',
            **{i: str(i) for i in range(18, 66)},
            66: 'over 65',
        },
    ),
    'sex_v2': _r(
        plot_label='Sex',
        plot_label_verbose='Sex (M/F)',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=True,
        predictor=True, dv=False, mediator=False,
        redcap_field_label='What sex were you assigned at birth?',
        redcap_choices={1: 'Male', 2: 'Female', 3: 'Other'},
    ),
    'mental_illness2_v2': _r(
        plot_label='Mental Illness Hx',
        plot_label_verbose='Mental Illness History',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=True,
        predictor=True, dv=False, mediator=False,
        redcap_field_label='Have you ever received (or been extremely confident that you should have received) a mental illness diagnosis (besides those listed above)?',
        redcap_choices=np.nan,  # free-text entry — no coded response options
    ),
    'psych_spectrum_v2': _r(
        plot_label='Psychosis Spectrum Dx',
        plot_label_verbose='Psychosis Spectrum Diagnosis',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=True,
        predictor=True, dv=False, mediator=False,
        redcap_field_label=(
            'Have you ever been diagnosed with a psychotic spectrum disorder such as Schizophrenia, '
            'Schizophreniform disorder, Schizoaffective disorder, Schizotypal personality disorder, '
            'Bipolar with psychotic features, Major Depressive Disorder with psychotic features, '
            'brief psychotic disorder, delusional disorder or psychosis NOS (not otherwise specified)?'
        ),
        redcap_choices=np.nan,  # free-text entry — no coded response options
    ),
    'psychedelic_primary': _r(
        plot_label='Primary SP',
        plot_label_verbose='Primary SP',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=True,
        predictor=True, dv=False, mediator=False,
        redcap_field_label='What psychedelic do you most often use?',
        redcap_choices={
            1: 'Psilocybin (magic mushrooms)',
            2: 'LSD',
            3: 'Mescaline (peyote, san pedro)',
            4: 'DMT (Ayahuasca, smoked)',
            5: '5-MeO-DMT (bufo, toads)',
            6: 'Other',
        },
    ),
    'lastdose_recency': _r(
        plot_label='Last Dose x Recency',
        plot_label_verbose='Last Dose x Recency',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'amph_lifetime': _r(
        plot_label='Amphetamine Use (Y/N)',
        plot_label_verbose='Amphetamine Use (Y/N)',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=True,
        predictor=True, dv=False, mediator=False,
        redcap_field_label='Amphetamines or stimulants (speed, meth, Adderall, Ritalin, etc.) -- USED WITHOUT A PRESCRIPTION AND/OR FOR NON-PRESCRIBED PURPOSE',
        redcap_choices={1: 'Yes', 2: 'No'},
    ),
    'atypicals_6mo': _r(
        plot_label='Past 6 Month Non-serotonergic Psychedelic Uses',
        plot_label_verbose='Past 6 Month Non-serotonergic Psychedelic Uses',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    'race_bipoc': _r(
        plot_label='BIPOC',
        plot_label_verbose='BIPOC',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=True,
        predictor=True, dv=False, mediator=False,
    ),
    'inhalants_lifetime': _r(
        plot_label='Inhalant Use (Y/N)',
        plot_label_verbose='Inhalant Use (Y/N)',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=True,
        predictor=True, dv=False, mediator=False,
        redcap_field_label='Inhalants (hydrocarbons only -- glue, gasoline, paint thinner, etc. NOT poppers or laughing gas/NO)',
        redcap_choices={1: 'Yes', 2: 'No'},
    ),
    'coke_lifetime': _r(
        plot_label='Cocaine Use (Y/N)',
        plot_label_verbose='Cocaine Use (Y/N)',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=True,
        predictor=True, dv=False, mediator=False,
        redcap_field_label='Cocaine (coke, crack, etc.)',
        redcap_choices={1: 'Yes', 2: 'No'},
    ),
    'race_asian': _r(
        plot_label='Asian',
        plot_label_verbose='Asian',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=True,
        predictor=True, dv=False, mediator=False,
    ),
    'raven_total': _r(
        plot_label='RAVEN Score (/9)',
        plot_label_verbose='RAVEN Score (/9)',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=False, mediator=False,
    ),
    'psycheduse_month_nomic': _r(
        plot_label=np.nan,
        plot_label_verbose=np.nan,
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Approximately how many times in the past month (30 days) did you use a serotonergic psychedelic NOT including microdosing?',
        redcap_choices=np.nan,  # free-entry numeric — no discrete response options
    ),
    'psycheduse_6month_nomic': _r(
        plot_label=np.nan,
        plot_label_verbose=np.nan,
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Approximately how many times in the past 6 months (180 days) did you use a serotonergic psychedelic NOT including microdosing?',
        redcap_choices=np.nan,  # free-entry numeric — no discrete response options
    ),
    'psycheduse_year_nomic': _r(
        plot_label=np.nan,
        plot_label_verbose=np.nan,
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Approximately how many times in the past year did you use a serotonergic psychedelic NOT including microdosing?',
        redcap_choices=np.nan,  # free-entry numeric — no discrete response options
    ),
    'lastdose_recency': _r(
        plot_label='Last Dose x Recency',
        plot_label_verbose='Last Dose x Recency',
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
    ),
    # ── Past-month drug use — binary (Y/N) ────────────────────────────────────

    'alc_month_yn': _r(
        plot_label='Alcohol (Past Mo.)',
        plot_label_verbose='Alcohol Use (Past Month)',
        plot_label_lay_audience='Drank Alcohol in the Past Month',
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=True,
        predictor=True, dv=False, mediator=False,
        redcap_field_label=np.nan, redcap_choices=np.nan,
    ),
    'ghb_month_yn': _r(
        plot_label='Sed.-Hypnotic (Past Mo.)',
        plot_label_verbose='Sedative-Hypnotic Use (Past Month)',
        plot_label_lay_audience='Used a Sedative or Hypnotic in the Past Month',
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=True,
        predictor=True, dv=False, mediator=False,
        redcap_field_label=np.nan, redcap_choices=np.nan,
    ),
    'opioids_month_yn': _r(
        plot_label='Opioid (Past Mo.)',
        plot_label_verbose='Opioid Use (Past Month)',
        plot_label_lay_audience='Used an Opioid in the Past Month',
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=True,
        predictor=True, dv=False, mediator=False,
        redcap_field_label=np.nan, redcap_choices=np.nan,
    ),
    'mj_month_yn': _r(
        plot_label='Cannabis (Past Mo.)',
        plot_label_verbose='Cannabis Use (Past Month)',
        plot_label_lay_audience='Used Cannabis in the Past Month',
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=True,
        predictor=True, dv=False, mediator=False,
        redcap_field_label=np.nan, redcap_choices=np.nan,
    ),
    'atypicals_month_yn': _r(
        plot_label='Atyp. Psychedelic (Past Mo.)',
        plot_label_verbose='Atypical Psychedelics Use (Past Month)',
        plot_label_lay_audience='Used a Non-Serotonergic Psychedelic in the Past Month',
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=True,
        predictor=True, dv=False, mediator=False,
        redcap_field_label=np.nan, redcap_choices=np.nan,
    ),
    'stimulants_month_yn': _r(
        plot_label='Stimulant (Past Mo.)',
        plot_label_verbose='Stimulants Use (Past Month)',
        plot_label_lay_audience='Used a Stimulant in the Past Month',
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=True,
        predictor=True, dv=False, mediator=False,
        redcap_field_label=np.nan, redcap_choices=np.nan,
    ),
    # GABAergic composite: shipped as `sedatives_month = alc_month + ghb_month`.
    # Both act at GABA receptors — alcohol is a GABA-A positive allosteric
    # modulator, and the `ghb_*` REDCap item is the whole sedative-hypnotic class
    # ("Benzodiazepines or sleeping pills — Benzos, Xanax, Valium, GHB, Ambien"),
    # i.e. GABA-A (benzodiazepines, Z-drugs) and GABA-B (GHB) agents.
    #
    # Opioids are deliberately NOT in this composite: they act at mu-opioid
    # receptors, not GABA receptors.  `depressants_month` (= alc + ghb + opioids)
    # also ships and is the WRONG column for this covariate set.
    'sedatives_month_yn': _r(
        plot_label='GABAergic (Past Mo.)',
        plot_label_verbose='GABAergic Drug Use (Past Month)',
        plot_label_lay_audience='Used Alcohol or a Sedative in the Past Month',
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=True,
        predictor=True, dv=False, mediator=False,
        redcap_field_label=np.nan, redcap_choices=np.nan,
    ),
    'atypicals_6month': _r(
        plot_label=np.nan,
        plot_label_verbose=np.nan,
        plot_label_lay_audience=np.nan,
        distribution='student_t', data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=False, mediator=False,
    ),
    'highest_education_3level_final': _r(
        plot_label=np.nan,
        plot_label_verbose=np.nan,
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='ordinal',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=False, mediator=False,
    ),
    'ses_employ_3level': _r(
        plot_label=np.nan,
        plot_label_verbose=np.nan,
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='ordinal',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=False, mediator=False,
    ),
    'race_v2': _r(
        plot_label=np.nan,
        plot_label_verbose=np.nan,
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=True,
        predictor=True, dv=False, mediator=False,
        redcap_field_label='What is your race?',
        redcap_choices={
            1: 'American Indian/Alaska Native',
            2: 'Asian',
            3: 'Native Hawaiian or Other Pacific Islander',
            4: 'Black or African American',
            5: 'White',
            6: 'More than one race',
            7: 'Unknown/I prefer not to say',
        },
    ),

    # ── Empirically-motivated covariate set variables (added 2026-06-13) ────────
    # These variables support the empirical_covariates set, selected from:
    #   (1) group-comparison differences in the combined HPPD+CAPS split tables
    #       (results/descriptive/tables/clinical_table_hppd_split_caps_split.png,
    #        results/descriptive/tables/demographics_table_hppd_split_caps_split.png)
    #   (2) Spearman correlations with caps_vision in SP users with valid caps_vision
    #       (n=130; see results/supplement/caps_vision_confounds/spearman_table.csv
    #        and 04_visualizations/supplement/caps_vision_confounds_spearman.py)
    # ⚠️  PLOT LABELS BELOW ARE PROPOSALS — confirm with Max before use in figures.
    # ────────────────────────────────────────────────────────────────────────────

    # Spearman ρ = −0.18 with caps_vision (p = .046).
    # Numeric ordinal (1–9 survey scale); Gelman-normalized as continuous in R
    # (is_categorical_factor=False, same treatment as highest_education_3level_final).
    # Column corrected in data_prep.py (record 61: NaN → 6 "College student").
    # Raw numeric column is highest_education; education_ordinal is the string-label
    # Categorical version (for display only — use highest_education in models).
    'highest_education': _r(
        plot_label='Education Level',
        plot_label_verbose='Education Level (Ordinal 1–9)',
        plot_label_lay_audience='Education Level',
        distribution=np.nan, data_type='ordinal',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=False, mediator=False,
        redcap_field_label='What is the highest level of school you have completed or the highest degree you have received?',
        redcap_choices={
            1: 'Less than 8th grade education',
            2: 'Some High School',
            3: 'Highschool diploma or GED',
            4: 'Some College',
            5: '2-year degree (e.g., Associate of Arts)',
            6: 'Currently a college student',
            7: '4-year degree (e.g., B.A. or B.S.)',
            8: "Master's degree",
            9: 'Doctoral degree or 4 year professional degree',
        },
    ),

    # 6-level balanced collapse of highest_education (see data_prep.py for map).
    # Collapses sparse tails: levels 2+3 → 1, levels 8+9 → 6.
    # All cells have n ≥ 11 in the caps_vision analysis sample (n = 130 SP users).
    # Used with mo() (monotonic effects) in brms via MONOTONIC_COVARIATES.
    # Spearman ρ = −0.167, p = .057 with caps_vision.
    'highest_education_balanced': _r(
        plot_label='Education Level',
        plot_label_verbose='Education Level (6-level balanced)',
        plot_label_lay_audience='Education Level',
        distribution=np.nan, data_type='ordinal',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=True, dv=False, mediator=False,
        redcap_field_label='What is the highest level of school you have completed or the highest degree you have received?',
        redcap_choices={
            1: 'HS diploma or GED or less (original levels 2–3)',
            2: 'Some college (original level 4)',
            3: '2-year degree (original level 5)',
            4: 'College student, currently enrolled (original level 6)',
            5: "4-year bachelor's degree (original level 7)",
            6: "Master's or higher (original levels 8–9)",
        },
    ),

    # Spearman ρ = −0.14 with caps_vision (p = .111, non-significant) but
    # included based on group differences in the clinical CAPS-split table.
    # Binary 0/1 (0 = no current psych medication, 1 = any current psych medication).
    # Treated as factor in R (is_categorical_factor=True, auto-derived from
    # data_type='categorical' and predictor=True in CATEGORICAL_FACTOR_VARS).
    'medication_current_v2': _r(
        plot_label='Psych Rx.',
        plot_label_verbose='Current Psychiatric Medication',
        plot_label_lay_audience='Currently Taking Psychiatric Medication',
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=True,
        predictor=True, dv=False, mediator=False,
        redcap_field_label='Are you currently taking any of the following psychiatric medications?',
        redcap_choices=np.nan,  # composite binary derived from multiple medication items
    ),

    # Spearman ρ = +0.28 with caps_vision (p = .001).
    # Binary 0/1 derived when the analysis dataframe is built, from
    # outside_us_v2 using the canonical European country-code list.
    # Treated as factor in R (is_categorical_factor=True).
    'location_europe': _r(
        plot_label='European',
        plot_label_verbose='European',
        plot_label_lay_audience='Located in Europe',
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=True,
        predictor=True, dv=False, mediator=False,
        redcap_field_label=np.nan,  # derived column — not a direct REDCap item
        redcap_choices=np.nan,
    ),

    # ── VCH Task — Quality Control / Engagement Metrics ──────────────────────
    # Used in: 04_visualizations/supplement/vch_beta_qc_scatter_supplement.py.
    # None of these enter statistical models (distribution=np.nan throughout).
    # effort_qc and distraction_qc are 1–5 Likert items from the post-task
    # self-report; n_timeouts_total/n_timeout_trials are derived from task JSON
    # (each "timeout" = a stimulus that was re-presented because no response was
    # detected within the response window; hence "Representations").
    # z_composite_rt/accuracy are EXPERIMENTAL multiplicative + directional
    # disengagement scores.  They are not reported and enter no figure or model;
    # they are retained only so the columns have a documented definition.
    # Formula:
    #   1. z-score delta components (z_delta_rt, z_delta_accuracy) across SP users
    #   2. rectify to enforce direction: max(0, −z_delta_rt), max(0, −z_delta_accuracy)
    #   3. multiply by raw longest_streak_length (not z-scored, avoids sign ambiguity)
    #   4. z-score the product across SP users
    # Higher score = BOTH a long same-response streak AND a meaningful directional drop in
    # RT / accuracy during that streak vs. outside it (lower RT = faster; lower acc = worse).
    'effort_qc': _r(
        plot_label='"I did the best I could on the tasks"',
        plot_label_verbose='Self-Reported Effort (1–5 Likert)',
        plot_label_lay_audience='Effort Rating',
        distribution=np.nan, data_type='ordinal',
        need_non_normalized=np.nan, inplace_normalized=np.nan, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='I did the best I could on the tasks',
        redcap_choices=np.nan,  # 1–5: Strongly Disagree → Strongly Agree
    ),
    'distraction_qc': _r(
        plot_label='"I was not distracted during the tasks"',
        plot_label_verbose='Self-Reported Non-Distraction (1–5 Likert)',
        plot_label_lay_audience='Distraction Rating',
        distribution=np.nan, data_type='ordinal',
        need_non_normalized=np.nan, inplace_normalized=np.nan, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='I was not distracted by my phone or other things during the surveys and tasks',
        redcap_choices=np.nan,  # 1–5: Strongly Disagree → Strongly Agree
    ),
    'n_timeouts_total': _r(
        plot_label='Total Trial Repeats',
        plot_label_verbose='Total Stimulus Re-presentations (sum across 4 blocks)',
        plot_label_lay_audience='Total Times Stimulus Was Repeated',
        distribution=np.nan, data_type='continuous',
        need_non_normalized=np.nan, inplace_normalized=np.nan, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,  # derived from task JSON: Σ(n_stim − n_resp) per block
        redcap_choices=np.nan,
    ),
    'n_timeout_trials': _r(
        plot_label='Trials with ≥1 Trial Repeat',
        plot_label_verbose='Distinct Trials with At Least One Stimulus Re-presentation (Trial Repeat)',
        plot_label_lay_audience='Trials Where Stimulus Was Repeated',
        distribution=np.nan, data_type='continuous',
        need_non_normalized=np.nan, inplace_normalized=np.nan, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,  # derived from task JSON: count of trials with stim_count ≥ 2
        redcap_choices=np.nan,
    ),
    'threshold_empiric_v_nominal': _r(
        plot_label='Empiric vs. QUEST-ideal Threshold',
        plot_label_verbose='Absolute Deviation of Empiric 75% Hit Rate from QUEST Ideal (|vch_bl_yes_75 − 0.75|)',
        plot_label_lay_audience='How Well the Task Was Calibrated',
        distribution=np.nan, data_type='continuous',
        need_non_normalized=np.nan, inplace_normalized=np.nan, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,  # computed: |vch_bl_yes_75 − 0.75|
        redcap_choices=np.nan,
    ),
    'z_composite_rt': _r(
        plot_label='Response Perseveration × RT Drop',
        plot_label_verbose='Z-scored composite (multiplicative): streak_length × max(0, −z_delta_rt); '
                           'higher = long streak WITH faster RT inside vs. outside [EXPERIMENTAL]',
        plot_label_lay_audience='Button-mashing Indicator (Reaction Time)',
        distribution=np.nan, data_type='continuous',
        need_non_normalized=np.nan, inplace_normalized=np.nan, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,  # computed: Z(longest_streak_length × max(0, −z_delta_rt))
        redcap_choices=np.nan,
    ),
    'z_composite_accuracy': _r(
        plot_label='Response Perseveration × Accuracy Drop',
        plot_label_verbose='Z-scored composite (multiplicative): streak_length × max(0, −z_delta_accuracy); '
                           'higher = long streak WITH lower accuracy inside vs. outside [EXPERIMENTAL]',
        plot_label_lay_audience='Button-mashing Indicator (Accuracy)',
        distribution=np.nan, data_type='continuous',
        need_non_normalized=np.nan, inplace_normalized=np.nan, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,  # computed: Z(longest_streak_length × max(0, −z_delta_accuracy))
        redcap_choices=np.nan,
    ),

    # ── Hardware control ──────────────────────────────────────────────────────
    # 3-level display-class grouping of the free-text `monitor_check` REDCap
    # field: 'Mac' | 'Windows/Other Laptop' | 'External Monitor'.
    #
    # ADDED 2026-08-31 in response to a reviewer who argued that unmeasured
    # display properties (panel type, brightness, refresh rate, viewing
    # distance) inject noise into visual contrast thresholding.  Used as a
    # covariate in the model type `nice_covariates_spusers_hardware_control`.
    #
    # DERIVATION — hand coding of the free-text `monitor_check` REDCap field.
    # It ships as a materialised column in data/final/df_public_*.csv and is read
    # from there; nothing recomputes it at run time.  Do NOT attempt to derive it
    # from `monitor_check` with keyword matching, which will not reproduce the
    # hand coding.  Nothing in this repository derives it; to change the coding,
    # change it where the dataframe is built and re-export.  The rules are
    # documented in 04_visualizations/supplement/README.md.
    #
    # NOTE: string-valued, so it MUST reach R as a factor and must never be
    # Gelman-normalized.  data_type='categorical' + predictor=True puts it in
    # CATEGORICAL_FACTOR_VARS automatically (that list is derived, not
    # hand-maintained), which is what excludes it from normalization and emits
    # it into the auto-generated categorical_factor_vars.R.
    #
    # Empirical relationship to the VCH measures (N = 193, non-responders
    # dropped; 04_visualizations/supplement/hardware_keydown_check.py):
    #   vch_threshold   Kruskal-Wallis H(2) = 0.29,  p = .864, eps^2 = 0.002
    #   d_prime_overall Kruskal-Wallis H(2) = 11.25, p = .0036, eps^2 = 0.059
    #                   (laptop 2.30 > Mac 2.11 > external monitor 1.91)
    'monitor_check_operationalized_final': _r(
        plot_label='Display Class',
        plot_label_verbose='Display Class (Mac / Windows-Other Laptop / External Monitor)',
        plot_label_lay_audience='Type of Screen Used',
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False,
        is_categorical_factor=True,
        predictor=True, dv=False, mediator=False,
        redcap_field_label='What kind of monitor are you using?',
        redcap_choices=np.nan,   # free text, hand-coded; see derivation note above
    ),

    # ══════════════════════════════════════════════════════════════════════════
    # Pipeline-support variables (QC, screening, CONSORT, recruitment, raw task
    # payloads and long-format VCH counters).
    #
    # These are read by pipeline code but are NOT analysis variables: none is a
    # predictor, outcome or mediator in any model, which is why `distribution`
    # is np.nan throughout.  They are registered so that VARIABLE_REGISTRY names
    # every column the code actually reads -- that is what lets the shipped
    # df_public_<date>.csv be cut down to the registry's key set without
    # breaking a script.  Do not set predictor=True on any of these without
    # also supplying a real brms family: CATEGORICAL_FACTOR_VARS is derived as
    # (data_type in binary/categorical AND predictor), so flipping that flag
    # silently changes R-side as.factor() handling.
    #
    # data_type follows the REDCap data-dictionary field type rather than the
    # values observed in df_public, because df_public is a 228-row
    # post-exclusion subset in which screening-stage fields are near-constant.
    # np.nan means free text, a date, or a base64 task payload.
    # ══════════════════════════════════════════════════════════════════════════
    'activecannabisuse_lastuse': _r(
        plot_label='activecannabisuse_lastuse',
        plot_label_verbose='How many days has it been since you last used cannabis? (if more than...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='How many days has it been since you last used cannabis? (if more than...',
        redcap_choices=np.nan,
    ),
    'addiction': _r(
        plot_label='addiction',
        plot_label_verbose='addiction',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'age_qc': _r(
        plot_label='What is your age? (in years)',
        plot_label_verbose='What is your age? (in years)',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='What is your age? (in years)',
        redcap_choices=np.nan,
    ),
    'alc_lifetime': _r(
        plot_label='alc_lifetime',
        plot_label_verbose='Alcoholic beverages (beer, wine, liquor, etc.)',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Alcoholic beverages (beer, wine, liquor, etc.)',
        redcap_choices='1, Yes | 2, No',
    ),
    'alc_month': _r(
        plot_label='alc_month',
        plot_label_verbose='Approximately how many times have you used alcohol in the past month (...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Approximately how many times have you used alcohol in the past month (...',
        redcap_choices=np.nan,
    ),
    'antidepressants': _r(
        plot_label='antidepressants',
        plot_label_verbose='antidepressants',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'antipsychotic': _r(
        plot_label='antipsychotic',
        plot_label_verbose='antipsychotic',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'anxiety_disorder': _r(
        plot_label='anxiety_disorder',
        plot_label_verbose='anxiety_disorder',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'attentionq_qc': _r(
        plot_label='attentionq_qc',
        plot_label_verbose='I paid close attention to the questions.',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='I paid close attention to the questions.',
        redcap_choices='1, strongly disagree | 2, disagree | 3, neither agree nor disagree | 4, agree | 5, strongl...',
    ),
    'attn_check_surveybl': _r(
        plot_label='attn_check_surveybl',
        plot_label_verbose='Returns 1 if they did the attention check above correctly',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Returns 1 if they did the attention check above correctly',
        redcap_choices='if([attn_check_animal]=\'2\',1,0)',
    ),
    'attn_check_surveybl2': _r(
        plot_label='attn_check_surveybl2',
        plot_label_verbose='Please select "Very Often or Always True" to show you\'re paying atten...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Please select "Very Often or Always True" to show you\'re paying atten...',
        redcap_choices='1, Very Often or Always True | 2, Often True | 3, Sometimes True | 4, Rarely True | 5, Nev...',
    ),
    'attn_check_surveybl3': _r(
        plot_label='attn_check_surveybl3',
        plot_label_verbose='Please select "Yes" to show you are paying attention :)',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Please select "Yes" to show you are paying attention :)',
        redcap_choices=np.nan,
    ),
    'attn_check_surveybl4': _r(
        plot_label='attn_check_surveybl4',
        plot_label_verbose='Please enter "1" to show you are still paying attention :)',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Please enter "1" to show you are still paying attention :)',
        redcap_choices=np.nan,
    ),
    'atypical_since_sp': _r(
        plot_label='atypical_since_sp',
        plot_label_verbose='Returns 1 if they\'ve used an atypical AND it was more recent than whe...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Returns 1 if they\'ve used an atypical AND it was more recent than whe...',
        redcap_choices='if( ( ([pcp_lifetime] = \'1\' and ([pcp_dayslastuse] < [sp_lastuse_days_screen])) or ([mdm...',
    ),
    'atypicals_life_yn': _r(
        plot_label='atypicals_life_yn',
        plot_label_verbose='atypicals_life_yn',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'atypicals_month': _r(
        plot_label='atypicals_month',
        plot_label_verbose='atypicals_month',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'benzos': _r(
        plot_label='benzos',
        plot_label_verbose='benzos',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'cannabis_frequency': _r(
        plot_label='cannabis_frequency',
        plot_label_verbose='How often do you usually use cannabis?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='How often do you usually use cannabis?',
        redcap_choices='1, I do not use cannabis | 2, less than once a year | 3, once a year | 4, once every 3-6 m...',
    ),
    'cognition_screener_v2': _r(
        plot_label='cognition_screener_v2',
        plot_label_verbose='Have you been diagnosed with any neurological or medical problem that...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Have you been diagnosed with any neurological or medical problem that...',
        redcap_choices=np.nan,
    ),
    'consent_baseline_complete': _r(
        plot_label='consent_baseline_complete',
        plot_label_verbose='consent_baseline_complete',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'continue_date': _r(
        plot_label='continue_date',
        plot_label_verbose='Date when they will be 3 months out from SP:',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type=np.nan,
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Date when they will be 3 months out from SP:',
        redcap_choices=np.nan,
    ),
    'correct_answer_v2': _r(
        plot_label='Correct answer?',
        plot_label_verbose='Correct answer?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Correct answer?',
        redcap_choices='if([raven_resp_1_v2]=\'3\',1,0)',
    ),
    'd_prime': _r(
        plot_label='d_prime',
        plot_label_verbose='d_prime',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'depressants_month': _r(
        plot_label='depressants_month',
        plot_label_verbose='depressants_month',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'eating_disorder': _r(
        plot_label='eating_disorder',
        plot_label_verbose='eating_disorder',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'english_fluency': _r(
        plot_label='Are you fluent in English?',
        plot_label_verbose='Are you fluent in English?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Are you fluent in English?',
        redcap_choices=np.nan,
    ),
    'geo_crit': _r(
        plot_label='geo_crit',
        plot_label_verbose='Please click this button to continue :)',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Please click this button to continue :)',
        redcap_choices='1, Click here!',
    ),
    'ghb_lifetime': _r(
        plot_label='ghb_lifetime',
        plot_label_verbose='Sedative-Hypnotics like Benzodiazepines or sleeping pills ("Benzos", X...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Sedative-Hypnotics like Benzodiazepines or sleeping pills ("Benzos", X...',
        redcap_choices='1, Yes | 2, No',
    ),
    'ghb_month': _r(
        plot_label='ghb_month',
        plot_label_verbose='Approximately how many times have you used benzodiazepines, depressant...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Approximately how many times have you used benzodiazepines, depressant...',
        redcap_choices=np.nan,
    ),
    'headphones_check': _r(
        plot_label='headphones_check',
        plot_label_verbose='What model or type of headphones are you using?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type=np.nan,
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='What model or type of headphones are you using?',
        redcap_choices=np.nan,
    ),
    'hits': _r(
        plot_label='hits',
        plot_label_verbose='hits',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'honesty_qc': _r(
        plot_label='honesty_qc',
        plot_label_verbose='I answered the questions truthfully.',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='I answered the questions truthfully.',
        redcap_choices='1, strongly disagree | 2, disagree | 3, neither agree nor disagree | 4, agree | 5, strongl...',
    ),
    'howtheyfoundus': _r(
        plot_label='How did you hear about our study?',
        plot_label_verbose='How did you hear about our study?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='How did you hear about our study?',
        redcap_choices='1, Yale Intro to Psychology | 2, Bluelight.org | 3, Shroomery.org | 4, Reddit | 5, Contact...',
    ),
    'hppd_current': _r(
        plot_label='hppd_current',
        plot_label_verbose='Are you currently living with Hallucinogen Persisting Perception Disor...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Are you currently living with Hallucinogen Persisting Perception Disor...',
        redcap_choices=np.nan,
    ),
    'hppd_ever': _r(
        plot_label='hppd_ever',
        plot_label_verbose='Have you ever been diagnosed with Hallucinogen Persisting Perception D...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Have you ever been diagnosed with Hallucinogen Persisting Perception D...',
        redcap_choices=np.nan,
    ),
    'ineligibile_reason': _r(
        plot_label='ineligibile_reason',
        plot_label_verbose='ineligibile_reason',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'intox_screen_v2': _r(
        plot_label='intox_screen_v2',
        plot_label_verbose='Have you used any psychoactive drugs (including alcohol) today? This d...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Have you used any psychoactive drugs (including alcohol) today? This d...',
        redcap_choices=np.nan,
    ),
    'location_summary': _r(
        plot_label='location_summary',
        plot_label_verbose='location_summary',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'mental_illness_dx_current_30': _r(
        plot_label='mental_illness_dx_current_30',
        plot_label_verbose='mental_illness_dx_current_30',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'mental_illness_dx_current_5': _r(
        plot_label='mental_illness_dx_current_5',
        plot_label_verbose='mental_illness_dx_current_5',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'mental_illness_dx_current_6': _r(
        plot_label='mental_illness_dx_current_6',
        plot_label_verbose='mental_illness_dx_current_6',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'mental_illness_dx_current_9': _r(
        plot_label='mental_illness_dx_current_9',
        plot_label_verbose='mental_illness_dx_current_9',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'mj_6month': _r(
        plot_label='mj_6month',
        plot_label_verbose='Approximately how many times have you used Cannabis in the past 6 mont...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Approximately how many times have you used Cannabis in the past 6 mont...',
        redcap_choices=np.nan,
    ),
    'mj_lifetime': _r(
        plot_label='MJ use repeated',
        plot_label_verbose='MJ use repeated',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='MJ use repeated',
        redcap_choices='if([activecannabisuse_yn]<>"",[activecannabisuse_yn],[mj_life_yn])',
    ),
    'mj_month': _r(
        plot_label='mj_month',
        plot_label_verbose='Approximately how many times have you used Cannabis in the past month?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Approximately how many times have you used Cannabis in the past month?',
        redcap_choices=np.nan,
    ),
    'monitor_check': _r(
        plot_label='monitor_check',
        plot_label_verbose='What type of computer monitor are you using or type of laptop?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type=np.nan,
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='What type of computer monitor are you using or type of laptop?',
        redcap_choices=np.nan,
    ),
    'mood_disorder': _r(
        plot_label='mood_disorder',
        plot_label_verbose='mood_disorder',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'no_computer': _r(
        plot_label='no_computer',
        plot_label_verbose='Are you ONLY able to use a tablet or mobile device for this study?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Are you ONLY able to use a tablet or mobile device for this study?',
        redcap_choices=np.nan,
    ),
    'nonbenzo_anxiolytics': _r(
        plot_label='nonbenzo_anxiolytics',
        plot_label_verbose='nonbenzo_anxiolytics',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'nonbenzo_sedatives': _r(
        plot_label='nonbenzo_sedatives',
        plot_label_verbose='nonbenzo_sedatives',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'num_psych_diagnoses': _r(
        plot_label='num_psych_diagnoses',
        plot_label_verbose='num_psych_diagnoses',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'onlinehonesty_qc': _r(
        plot_label='onlinehonesty_qc',
        plot_label_verbose='I find it easier to answer senstive questions honestly online, compare...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='I find it easier to answer senstive questions honestly online, compare...',
        redcap_choices='1, strongly disagree | 2, disagree | 3, neither agree nor disagree | 4, agree | 5, strongl...',
    ),
    'opioid_antagonists': _r(
        plot_label='opioid_antagonists',
        plot_label_verbose='opioid_antagonists',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'opioids_lifetime': _r(
        plot_label='opioids_lifetime',
        plot_label_verbose='Opioids (heroin, morphine, methadone, codeine, etc.) -- USED WITHOUT A...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Opioids (heroin, morphine, methadone, codeine, etc.) -- USED WITHOUT A...',
        redcap_choices='1, Yes | 2, No',
    ),
    'opioids_month': _r(
        plot_label='opioids_month',
        plot_label_verbose='opioids_month',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'outside_us_v2': _r(
        plot_label='Where do you live?',
        plot_label_verbose='Where do you live?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Where do you live?',
        redcap_choices='1, Afghanistan|2, Albania|3, Algeria|4, Andorra|5, Angola|6, Antigua and Barbuda|7, Argent...',
    ),
    'pc_screen_yale1': _r(
        plot_label='pc_screen_yale1',
        plot_label_verbose='1. I think that I have felt that there are odd or unusual things going...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='1. I think that I have felt that there are odd or unusual things going...',
        redcap_choices='0, Definitely disagree | 1, Somewhat disagree | 2, Slightly disagree | 3, Not sure | 4, Sl...',
    ),
    'pc_screen_yale2': _r(
        plot_label='pc_screen_yale2',
        plot_label_verbose='2. I think that I might be able to predict the future.',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='2. I think that I might be able to predict the future.',
        redcap_choices='0, Definitely disagree | 1, Somewhat disagree | 2, Slightly disagree | 3, Not sure | 4, Sl...',
    ),
    'perceived_benefit': _r(
        plot_label='perceived_benefit',
        plot_label_verbose='On a scale of 0 to 100, with 50 being no effect, how has your psychede...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='On a scale of 0 to 100, with 50 being no effect, how has your psychede...',
        redcap_choices='Harmful | No Effect | Beneficial',
    ),
    'persist_vis_current': _r(
        plot_label='persist_vis_current',
        plot_label_verbose='persist_vis_current',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'persistvis_distress': _r(
        plot_label='persistvis_distress',
        plot_label_verbose='Did these visual effects ever cause you any distress or make it hard t...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Did these visual effects ever cause you any distress or make it hard t...',
        redcap_choices='1, Yes | 2, No | 3, NA',
    ),
    'persistvis_most': _r(
        plot_label='persistvis_most',
        plot_label_verbose='Which effect was the most vivid/intense?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Which effect was the most vivid/intense?',
        redcap_choices='13, I have never experienced any of the above visual effects | 1, Halos or auras around th...',
    ),
    'persistvis_psychdoses': _r(
        plot_label='persistvis_psychdoses',
        plot_label_verbose='About how many times had you taken psychedelics before you first notic...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='About how many times had you taken psychedelics before you first notic...',
        redcap_choices=np.nan,
    ),
    'persistvis_txseek': _r(
        plot_label='persistvis_txseek',
        plot_label_verbose='Did you ever seek help for these visual effects?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Did you ever seek help for these visual effects?',
        redcap_choices='1, Yes | 2, No | 3, NA',
    ),
    'personality_disorder': _r(
        plot_label='personality_disorder',
        plot_label_verbose='personality_disorder',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'phone_number': _r(
        plot_label='phone_number',
        plot_label_verbose='phone_number',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'psyched_yearsofuse': _r(
        plot_label='psyched_yearsofuse',
        plot_label_verbose='How many years in total have you used serotonergic psychedelic ( EXCLU...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='How many years in total have you used serotonergic psychedelic ( EXCLU...',
        redcap_choices=np.nan,
    ),
    'psychedelicuse_lifetimetot': _r(
        plot_label='psychedelicuse_lifetimetot',
        plot_label_verbose='How many times in your life have you used serotonergic psychedelics ?...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='How many times in your life have you used serotonergic psychedelics ?...',
        redcap_choices=np.nan,
    ),
    'psycheduse_yn': _r(
        plot_label='psycheduse_yn',
        plot_label_verbose='Have you ever used a Serotonergic psychedelic (magic mushrooms, LSD, p...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Have you ever used a Serotonergic psychedelic (magic mushrooms, LSD, p...',
        redcap_choices='1, Yes -- not just a microdose | 2, No | 3, Yes but only microdoses (non-detectable effect...',
    ),
    'qc_bad_data': _r(
        plot_label='qc_bad_data',
        plot_label_verbose='Did they fail the major QC suggesting repeat survey taking and need th...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Did they fail the major QC suggesting repeat survey taking and need th...',
        redcap_choices=np.nan,
    ),
    'qc_notes': _r(
        plot_label='qc_notes',
        plot_label_verbose='Any additional notes on QC for this participant?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type=np.nan,
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Any additional notes on QC for this participant?',
        redcap_choices=np.nan,
    ),
    'qc_passed': _r(
        plot_label='qc_passed',
        plot_label_verbose='Quality Check passed?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Quality Check passed?',
        redcap_choices=np.nan,
    ),
    'race_qc': _r(
        plot_label='What is your race?',
        plot_label_verbose='What is your race?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='What is your race?',
        redcap_choices='1, American Indian/Alaska Native | 2, Asian | 3, Native Hawaiian or Other Pacific Islander...',
    ),
    'raven_total_score_v2': _r(
        plot_label='Number of correct answers',
        plot_label_verbose='Number of correct answers',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Number of correct answers',
        redcap_choices='[correct_answer4_v2]+[correct_answer5_v2]+[correct_answer6_v2]+[correct_answer7_v2]+[corre...',
    ),
    'record_id': _r(
        plot_label='Record ID',
        plot_label_verbose='Record ID',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Record ID',
        redcap_choices=np.nan,
    ),
    'salvage_yn': _r(
        plot_label='salvage_yn',
        plot_label_verbose='If they never finished but have SOME data -- do they otherwise pass QC...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='If they never finished but have SOME data -- do they otherwise pass QC...',
        redcap_choices=np.nan,
    ),
    'schizophrenia_spectrum': _r(
        plot_label='schizophrenia_spectrum',
        plot_label_verbose='schizophrenia_spectrum',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'screening_pass': _r(
        plot_label='screening_pass',
        plot_label_verbose='Did the participant pass screening ? IE [phone_number] is not a Sinch...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Did the participant pass screening ? IE [phone_number] is not a Sinch...',
        redcap_choices=np.nan,
    ),
    'screening_survey_complete': _r(
        plot_label='screening_survey_complete',
        plot_label_verbose='screening_survey_complete',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'sedatives_6month': _r(
        plot_label='sedatives_6month',
        plot_label_verbose='sedatives_6month',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'sedatives_life_yn': _r(
        plot_label='sedatives_life_yn',
        plot_label_verbose='sedatives_life_yn',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'sedatives_month': _r(
        plot_label='sedatives_month',
        plot_label_verbose='sedatives_month',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'seizure_hx_v2': _r(
        plot_label='seizure_hx_v2',
        plot_label_verbose='Do you have a known seizure disorder?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Do you have a known seizure disorder?',
        redcap_choices=np.nan,
    ),
    'si_2_v2': _r(
        plot_label='Do you need a hearing aid?',
        plot_label_verbose='Do you need a hearing aid?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Do you need a hearing aid?',
        redcap_choices=np.nan,
    ),
    'simulant_medication': _r(
        plot_label='simulant_medication',
        plot_label_verbose='simulant_medication',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'sleep_disorder': _r(
        plot_label='sleep_disorder',
        plot_label_verbose='sleep_disorder',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'sp_experience': _r(
        plot_label='sp_experience',
        plot_label_verbose='sp_experience',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'sp_fraud_5meo': _r(
        plot_label='sp_fraud_5meo',
        plot_label_verbose='What is your preferred method of use of the above psychedelic?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='What is your preferred method of use of the above psychedelic?',
        redcap_choices='1, Injected | 2, Blotter/tab | 3, Pill | 4, Suppository | 5, Other',
    ),
    'sp_fraud_dmt': _r(
        plot_label='sp_fraud_dmt',
        plot_label_verbose='What is your preferred method of use of the above psychedelic?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='What is your preferred method of use of the above psychedelic?',
        redcap_choices='1, Pill + Monoamine oxidase inhibitor (MAOI; eg. Harmine) | 2, Pill (pure DMT) | 3, Inject...',
    ),
    'sp_fraud_lsd': _r(
        plot_label='sp_fraud_lsd',
        plot_label_verbose='What is your preferred method of use of the above psychedelic?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='What is your preferred method of use of the above psychedelic?',
        redcap_choices='1, Smoked | 2, Vaporized | 3, Snorted | 4, Injected | 5, Other',
    ),
    'sp_fraud_mesc': _r(
        plot_label='sp_fraud_mesc',
        plot_label_verbose='What is your preferred method of use of the above psychedelic?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='What is your preferred method of use of the above psychedelic?',
        redcap_choices='1, Smoked | 2, Snorted | 3, Injected | 4, Lozenge | 5, Other',
    ),
    'sp_fraud_psi': _r(
        plot_label='sp_fraud_psi',
        plot_label_verbose='What is your preferred method of use of the above psychedelic?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='What is your preferred method of use of the above psychedelic?',
        redcap_choices='1, Smoked | 2, Vaporized | 3, Snorted | 4, Injected | 5, Blotter/tab | 6, Other',
    ),
    'sp_naiive': _r(
        plot_label='sp_naiive',
        plot_label_verbose='Returns a 1 if they have used nonmicrodose SP doses OR haven\'t used S...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Returns a 1 if they have used nonmicrodose SP doses OR haven\'t used S...',
        redcap_choices='if(([psycheduse_yn]=\'2\' and [interested_spstudy]<>"") or ([psycheduse_yn]=\'1\') or ([ps...',
    ),
    'sp_type_recent': _r(
        plot_label='sp_type_recent',
        plot_label_verbose='Which serotonergic psychedelic did you use most recently ?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Which serotonergic psychedelic did you use most recently ?',
        redcap_choices='1, Psilocybin (magic mushrooms) | 2, LSD | 3, Mescaline (peyote, san pedro) | 4, DMT (Ayah...',
    ),
    'sp_type_recent_qc': _r(
        plot_label='sp_type_recent_qc',
        plot_label_verbose='Which psychedelic did you use most recently ?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='categorical',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Which psychedelic did you use most recently ?',
        redcap_choices='1, Psilocybin (magic mushrooms) | 2, LSD | 3, Mescaline (peyote, san pedro) | 4, DMT (Ayah...',
    ),
    'space_junkqc': _r(
        plot_label='Passed space junk QC?',
        plot_label_verbose='Passed space junk QC?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Passed space junk QC?',
        redcap_choices=np.nan,
    ),
    'stimulants_6month': _r(
        plot_label='stimulants_6month',
        plot_label_verbose='stimulants_6month',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'stimulants_life_yn': _r(
        plot_label='stimulants_life_yn',
        plot_label_verbose='stimulants_life_yn',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'stimulants_month': _r(
        plot_label='stimulants_month',
        plot_label_verbose='stimulants_month',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'student_yn': _r(
        plot_label='student_yn',
        plot_label_verbose='Are you a Yale student taking this study as part of Introduction to Ps...',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Are you a Yale student taking this study as part of Introduction to Ps...',
        redcap_choices=np.nan,
    ),
    'task_data_ach_task_short_baseline': _r(
        plot_label='task_data_ach_task_short_baseline',
        plot_label_verbose='task_data_ach_task_short_baseline',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type=np.nan,
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'task_data_prltask': _r(
        plot_label='task_data_prltask',
        plot_label_verbose='task_data_prltask',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type=np.nan,
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'task_data_vch_short_psychedelic_bl': _r(
        plot_label='task_data_vch_short_psychedelic_bl',
        plot_label_verbose='task_data_vch_short_psychedelic_bl',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type=np.nan,
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'timestamp_survey_bl': _r(
        plot_label='timestamp_survey_bl',
        plot_label_verbose='timestamp_survey_bl',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type=np.nan,
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'total_vch_correct_rejects': _r(
        plot_label='total_vch_correct_rejects',
        plot_label_verbose='total_vch_correct_rejects',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'total_vch_trials_0': _r(
        plot_label='total_vch_trials_0',
        plot_label_verbose='total_vch_trials_0',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'total_vch_trials_25': _r(
        plot_label='total_vch_trials_25',
        plot_label_verbose='total_vch_trials_25',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'total_vch_trials_50': _r(
        plot_label='total_vch_trials_50',
        plot_label_verbose='total_vch_trials_50',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'total_vch_trials_75': _r(
        plot_label='total_vch_trials_75',
        plot_label_verbose='total_vch_trials_75',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'trauma_disorder': _r(
        plot_label='trauma_disorder',
        plot_label_verbose='trauma_disorder',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'us_loc_v2': _r(
        plot_label='Do you live in the US?',
        plot_label_verbose='Do you live in the US?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Do you live in the US?',
        redcap_choices=np.nan,
    ),
    'vch_d_prime': _r(
        plot_label='vch_d_prime',
        plot_label_verbose='vch_d_prime',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'vch_false_alarm_rate': _r(
        plot_label='vch_false_alarm_rate',
        plot_label_verbose='vch_false_alarm_rate',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'vch_hits': _r(
        plot_label='vch_hits',
        plot_label_verbose='vch_hits',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'vch_hits_25': _r(
        plot_label='vch_hits_25',
        plot_label_verbose='vch_hits_25',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'vch_hits_50': _r(
        plot_label='vch_hits_50',
        plot_label_verbose='vch_hits_50',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'vch_miss_rate': _r(
        plot_label='vch_miss_rate',
        plot_label_verbose='vch_miss_rate',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'vch_misses': _r(
        plot_label='vch_misses',
        plot_label_verbose='vch_misses',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='continuous',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label=np.nan,
        redcap_choices=np.nan,
    ),
    'waiting_emailed_yn': _r(
        plot_label='waiting_emailed_yn',
        plot_label_verbose='Has this person been emailed by us yet?',
        plot_label_lay_audience=np.nan,
        distribution=np.nan, data_type='binary',
        need_non_normalized=False, inplace_normalized=False, is_categorical_factor=False,
        predictor=False, dv=False, mediator=False,
        redcap_field_label='Has this person been emailed by us yet?',
        redcap_choices=np.nan,
    ),
}

# Remove duplicate key that was defined twice (lastdose_recency was repeated above)
# Python dicts silently overwrite; the last definition wins — already handled.


###############################################################################
# SECTION 2: PIPELINE CONFIGURATION
# (copied verbatim from hpc_config.py and variable_labels.py)
###############################################################################

# ── Categorical / factor variables ────────────────────────────────────────────
# Derived from VARIABLE_REGISTRY: all variables that are ever used as predictors
# AND have data_type 'binary' or 'categorical'. These must never be
# Gelman-normalized or dummy-coded in Python; R converts them with as.factor().
# To add a new categorical predictor: set data_type='binary'/'categorical' and
# predictor=True in its VARIABLE_REGISTRY entry — do NOT extend this list directly.
CATEGORICAL_FACTOR_VARS = sorted(
    k for k, v in VARIABLE_REGISTRY.items()
    if v['data_type'] in ('binary', 'categorical') and v['predictor']
)

# ── R-side subsetting keywords ────────────────────────────────────────────────
# Row filters a model-type name may request.  The R scripts match these by
# grepl() on the model name and apply the subset before fitting.
#
# This list is the single source of truth, and its membership is load-bearing:
# _covs_for_variant() strips these suffixes to resolve a model name to a
# covariate set, so a name containing an unlisted keyword fails to resolve and
# raises at generation time.  Every keyword here is reachable from a configured
# model type; an unused keyword is an unused way to subset the sample, which
# this pipeline should not carry.
#
# Outlier exclusion is handled by '_iqr', which fences on the raw scale.  A
# threshold in Gelman-normalized space is not definable as a row filter, because
# normalization runs in R *after* subsetting — at filter time there is no
# normalized column and no fixed reference sample.
#
# No keyword may be a substring of another; grepl() would
# otherwise match both.
R_SIDE_SUFFIXES = [
    '_spusers',        # psycheduse_yn == "Yes" — the primary analysis sample
    '_iqr',            # IQR fence on the focal predictor (outlier sensitivity)
    '_nocurrenthppd',  # persist_vis_current == 0
    '_nopsychosis',    # excludes psychotic-spectrum diagnoses
    '_nonan_caps',     # participants with non-missing caps_vision
]

# ── Canonical base covariate sets ─────────────────────────────────────────────
BASE_COVARIATE_SETS = {
    'nice_covariates': [
        'age_v2', 'sex_v2', 'raven_total', 'mental_illness2_v2'
    ],
    'age_control': ['age_v2'],
    'true_univariate': [],
    'main': [
        'race_v2', 'raven_total', 'highest_education_3level_final',
        'sex_v2', 'ses_employ_3level', 'motivation',
    ],
    'empiric_covariates': [
        'amph_lifetime', 'race_bipoc', 'sex_v2', 'inhalants_lifetime',
        'coke_lifetime', 'age_v2', 'atypicals_6month', 'raven_total',
    ],

    # Covariate set chosen to carry the strongest available confounders for
    # caps_vision, identified from the group-difference tables and the Spearman
    # confound screen (Supplementary Table S1).  Reported as a sensitivity
    # specification alongside the primary model, which is nice_covariates.
    #
    # highest_education_balanced is wrapped in mo() when building brms formula
    # strings; see MONOTONIC_COVARIATES below and get_covs() in
    # 03_hpc/generate_hpc_jobs.py.
    'empirical_covariates': [
        'age_v2',
        'highest_education_balanced',
        'medication_current_v2',
    ],


}

# ── Derived / alias covariate sets ────────────────────────────────────────────
# Defined after the dict literal so they can reference other BASE_COVARIATE_SETS
# keys.  All downstream consumers (generate_hpc_jobs.py, get_covs()) pull from
# BASE_COVARIATE_SETS — do NOT duplicate these anywhere else.

# Backwards-compat alias: old job files may reference 'univariate' instead of
# 'true_univariate'.  Kept so existing .txt files continue to resolve correctly.
BASE_COVARIATE_SETS['univariate'] = BASE_COVARIATE_SETS['true_univariate'].copy()


# ── Past-month drug-use sensitivity sets ──────────────────────────────────────
# Reviewer request: show the primary associations hold net of RECENT drug use.
# Both extend `nice_covariates` and are used with the `_spusers` R-side suffix
# (-> nice_covariates key names below are the HPC results directory names).
#
# Screening evidence for both sets — collinearity is NOT the constraint here
# (max VIF 1.86, dropping to 1.43 when the depressant classes are collapsed);
# sparsity is:
#   04_visualizations/supplement/drugs_month_sensitivity_analysis.py
#   results/supplement/drugs_month_sensitivity_analysis/
#
# `drugs_month` — every remaining class kept separate, as reported in the
# descriptive table.  Still the messy version: of 186 SP users only 17 endorse
# ghb, 17 atypicals and 18 stimulants, too sparse for their individual
# coefficients to be interpretable on their own.
#
# `opioids_month_yn` is deliberately absent.  Only 7 of 186 SP users endorse it,
# and 2 of the 113 in the caps_vision analysis sample.  The models would not
# converge with it in the covariate set, including at adapt_delta = 0.999.
# Past-month opioid use is still reported descriptively; it is unusable only as
# a model covariate.
BASE_COVARIATE_SETS['drugs_month'] = BASE_COVARIATE_SETS['nice_covariates'] + [
    'alc_month_yn', 'ghb_month_yn', 'mj_month_yn',
    'atypicals_month_yn', 'stimulants_month_yn',
]

# `drugs_trimmed_month` — restricted to the classes that are actually estimable,
# and grouped by receptor pharmacology rather than by survey item.
#   sedatives_month_yn  GABAergic: alcohol (GABA-A PAM) + the sedative-hypnotic
#                       item (benzodiazepines and Z-drugs at GABA-A, GHB at
#                       GABA-B).  107 positive of 186.
#   mj_month_yn         cannabinoid
#   stimulants_month_yn monoaminergic
# Opioids are dropped rather than folded in: they act at mu-opioid receptors, so
# adding them to the GABAergic composite would make it a grab-bag rather than a
# receptor class.  Atypical psychedelics are dropped as too sparse (n = 17) to
# estimate on their own.
#
# NOTE: use `sedatives_month_yn` (= alc + ghb).  A `depressants_month` count
# column (= alc + ghb + opioids) also ships and is the wrong grouping here; it
# has no `_yn` companion and nothing in this pipeline should derive one.
BASE_COVARIATE_SETS['drugs_trimmed_month'] = BASE_COVARIATE_SETS['nice_covariates'] + [
    'sedatives_month_yn', 'mj_month_yn'
]

# Exposure-control variants: append SP-exposure variables when modelling VCH params
# → PPA risk so that VCH effects are estimated net of psychedelic exposure history.
# dict.fromkeys deduplicates while preserving order.
_EXPOSURE_ADDON = ['avg_life_dose', 'psychedelic_age', 'psycheduse_life_nomic']


# ── Hardware-control sensitivity set ──────────────────────────────────────────
# Reviewer request 2026-08-31: show the primary associations hold net of the
# DISPLAY the participant used.  `nice_covariates` + the 3-level display-class
# grouping, used with the `_spusers` R-side suffix, i.e. model type
# `nice_covariates_spusers_hardware_control`.
#
# NAMING: the R-side suffix is stripped by _strip_r_side() / _covs_for_variant()
# before lookup, so `nice_covariates_spusers_hardware_control` resolves to the
# key `nice_covariates_hardware_control` below.  'hardware_control' shares no
# substring with any keyword in R_SIDE_SUFFIXES or with any grepl() pattern in
# nonsp_predictors.R / hpc_mediation.R -- checked 2026-08-31, and required by
# the naming rule above (no keyword may be a substring of another).
#
# The covariate is a STRING factor, not a count or flag: it is in
# CATEGORICAL_FACTOR_VARS via its VARIABLE_REGISTRY entry, so R as.factor()s it
# and the normalizer skips it.  Verify that after any job-gen run -- a silently
# normalized string column is the failure mode this guards against.
BASE_COVARIATE_SETS['nice_covariates_hardware_control'] = (
    BASE_COVARIATE_SETS['nice_covariates'] + ['monitor_check_operationalized_final']
)


# ── Labels for different covariate sets (used in results/supplement/sensitivity_analyses figures) ──
COVARIATE_SET_LABELS = {
    'nice_covariates':
       "Age + Sex + IQ + Mental Illness Dx\n (+ 42 SP-naïve)",
    'nice_covariates_spusers':
       "Age + Sex + IQ + Mental Illness Dx\n(Primary)",
    'nice_covariates_spusers_iqr':
       "Age + Sex + IQ + Mental Illness Dx\n(IQR outlier removal)",
    'age_control': "Age",
    'age_control_spusers': "Age\n(SP users)",
    'true_univariate': "Univariate",
    'true_univariate_spusers': "Univariate\n(SP users)",
    'empirical_covariates': "Age + Education + Psych Rx + European Location",
    'empirical_covariates_spusers': "Age + Education + Psych Rx\n(SP users)",
    'nice_covariates_beta_spusers':
       "β + Age + Sex + IQ + Mental Illness Dx",
    'nice_covariates_beta_spusers_iqr':
       "β + Age + Sex + IQ + Mental Illness Dx\n(IQR outlier removal)",
    'nice_covariates_spusers_nonan_caps':
       "Age + Sex + IQ + Mental Illness Dx\n(No Missing CAPS)",
    'drugs_month_spusers':
       "Age + Sex + IQ + Mental Illness Dx\n+ Past-Month Drug Use (5 classes)",
    'drugs_trimmed_month_spusers':
       "Age + Sex + IQ + Mental Illness Dx\n+ Past-Month Sedative/Cannabis Use",
    'nice_covariates_spusers_hardware_control':
       "Age + Sex + IQ + Mental Illness Dx\n+ Display Class",}

# ── Monotonic-effects covariates ──────────────────────────────────────────────
# Variables in this set are wrapped in mo() when building brms formula strings
# via get_covs() in 03_hpc/generate_hpc_jobs.py.  All variables here must be
# integer-valued ordered scales suitable for brms monotonic effects.
# Do NOT include continuous or binary variables — mo() is for ordinal predictors only.
MONOTONIC_COVARIATES = {
    'highest_education_balanced',   # 6-level balanced education ordinal; see data_prep.py
}

# ── Normalization conventions ─────────────────────────────────────────────────
# Both lists are DERIVED from VARIABLE_REGISTRY, the same way CATEGORICAL_FACTOR_VARS
# is above.  To change how a variable is normalized, set `need_non_normalized=True`
# or `inplace_normalized=True` on its registry entry — do NOT extend these lists.
#
# NEED_NON_NORMALIZED — the variable's brms family constrains its scale (gamma →
#   strictly positive; beta / zero_inflated_beta → (0,1); hurdle / negbinomial →
#   non-negative integer).  Gelman-normalizing the raw column would push values
#   outside the family's support, so the raw column is preserved and a separate
#   `{col}_normalized` column carries the normalized values for use as a predictor.
#   caps_vision is the canonical case: a count DV that also appears as a covariate.
#
# INPLACE_NORMALIZED — the family is an unbounded real (student_t / gaussian), so
#   the raw column is itself Gelman-normalized and `{col}_normalized` is an alias
#   of it.  The same column name is therefore valid in both the mediator formula
#   and the DV formula.
#
# A registry key whose name already ends in `_normalized` names a column that holds
# normalized values by construction; leave both flags False on those entries so they
# are not re-normalized into a `{col}_normalized_normalized` duplicate.
NEED_NON_NORMALIZED = sorted(
    k for k, v in VARIABLE_REGISTRY.items()
    if v['need_non_normalized'] is True
)

INPLACE_NORMALIZED = sorted(
    k for k, v in VARIABLE_REGISTRY.items()
    if v['inplace_normalized'] is True
)

# ── VCH quality-control flag columns ──────────────────────────────────────────
# ── Analysis variable groupings (from variable_labels.py) ────────────────────

iv_type_dict = {}

iv_type_dict["clinical_predictors_continuous"] = [
    'age_v2', 'absorption_tot_bl', 'asi_tot', 'pdi_total', 'lshs_total',
    'phq9_tot', 'aaq_tot_bl', 'ffmq_total_bl',
]
iv_type_dict["clinical_predictors_categorical"] = ['mental_illness2_v2', 'sex_v2']

iv_type_dict["sp_predictors_continuous"] = [
    'psychedelic_age', 'psycheduse_life_nomic', 'avg_life_dose',
    'life_exposure', 'psychedelic_use_PC1',
]
iv_type_dict["other_sp_predictors"] = [
    'psychedelic_rank_use_PC1', 'psycheduse_recency_cutoff_hppd',
    'psyched_lastuse_dose', 'vasdose_bl', 'lastdose_recency',
    'psycheduse_month_nomic', 'psycheduse_6month_nomic', 'psycheduse_year_nomic',
]
iv_type_dict["recent_sp_uses"] = [
    'psycheduse_month_nomic', 'psycheduse_6month_nomic', 'psycheduse_year_nomic',
]
iv_type_dict["sp_predictors_categorical"] = ['psychedelic_primary', 'motivation']

severity_vars = [
    "persistvis_time", "hppd_true_chronicity", "persistvis_duration",
    "baggot_total", "asi_tot", "pdi_total", "lshs_total", "absorption_tot_bl",
]

iv_type_dict["sp_predictors"] = ['psychedelic_age', 'psycheduse_life_nomic', 'avg_life_dose']

iv_type_dict["vch_predictors"] = [
    'vch_threshold', 'vch_bl_yes_0', 'vch_nu', 'vch_beta', 'vch_bl_yes_75', 'vch_omega',
]
iv_type_dict["vch_main_vars"] = [
    'vch_threshold', 'vch_bl_yes_0', 'vch_nu', 'vch_beta', 'vch_omega', 'vch_bl_yes_75',
]
iv_type_dict["vch_behavior"] = ['vch_threshold', 'vch_bl_yes_75', 'vch_bl_yes_0']
iv_type_dict["vch_computations"] = ['vch_nu', 'vch_beta', 'vch_omega']

pwpe_cols = [
    'vch_pwPE_median', 'vch_pwPE_negative_median', 'vch_pwPE_positive_median',
    'vch_pwPE_negative_0_median', 'vch_pwPE_positive_0_median', 'vch_pwPE_ch_median',
    'vch_pwPE_negative_75_median', 'vch_pwPE_positive_75_median',
    'vch_pwPE_bias_0_median', 'vch_pwPE_bias_75_median', 'vch_pwPE_bias_median',
]

xprob_cols = [
    'vch_xprob_median', 'vch_xprob_change', 'vch_xprob_block_1', 'vch_xprob_block_12',
]
belief_cols = [
    'vch_belief_median', 'vch_belief_change', 'vch_belief_block_1', 'vch_belief_block_12',
]
xprob_belief_cols = ['vch_xprob_median', 'vch_belief_change']

iv_type_dict["vch_computations_states"] = xprob_belief_cols + pwpe_cols
iv_type_dict["vch_computations_expanded"] = (
    iv_type_dict["vch_computations"] + xprob_belief_cols + pwpe_cols
)

iv_type_dict["vch_comps_3lev"] = [
    'vch_beta_3lev', 'vch_nu_3lev', 'vch_nu_3lev_log', 'vch_omega_3lev', 'vch_omega3_3lev',
]
iv_type_dict["vch_comp_mat"] = [
    'vch_short_psychedelic_bl_omega3', 'vch_short_psychedelic_bl_omega2',
    'vch_short_psychedelic_bl_nu', 'vch_short_psychedelic_bl_beta',
]
iv_type_dict["vch_comp_nominal"] = ['vch_nu_nominal', 'vch_beta_nominal', 'vch_omega_nominal']
iv_type_dict["vch_comp_extended_nominal"] = (
    iv_type_dict["vch_comp_nominal"]
    + ['vch_xprob_median_nominal', 'vch_belief_change', 'vch_pwPE_median_nominal']
)

# ── Alternative HGF-parameter estimates ──────────────────────────────────────
# vch_comp_avg: per-variable average HGF parameters from a separate analysis
#   (f"{var}_avg" for each var in vch_computations; e.g. 'vch_nu_avg').
# vch_comp_nominal: HGF parameters fitted under nominal (non-optimised) priors.
iv_type_dict["vch_comp_avg"] = [f"{v}_avg" for v in iv_type_dict["vch_computations"]]
# → ['vch_nu_avg', 'vch_beta_avg', 'vch_omega_avg']

# ── SDT and metacognition variables ─────────────────────────────────────────
# All computed when the analysis dataframe is built, and cached
# to data/final/sdt_metacog_cache.csv.  No external CSV needed.

iv_type_dict["sdt_hppd"] = [
    'criterion_overall', 'd_prime_overall', 'mean_conf_fas',
]

all_predictors = (
    iv_type_dict["clinical_predictors_continuous"]
    + iv_type_dict["sp_predictors_continuous"]
    + iv_type_dict["vch_predictors"]
    + iv_type_dict["clinical_predictors_categorical"]
    + iv_type_dict["sp_predictors_categorical"]
)

main_predictors = (
    iv_type_dict["sp_predictors"]
    + iv_type_dict["vch_behavior"]
    + iv_type_dict["vch_computations"]
)

# Backward-compat alias (same list as NEED_NON_NORMALIZED)
need_non_normalized_versions = NEED_NON_NORMALIZED

caps_types = {}
caps_types["endorsement"] = ['caps_total', 'caps_vision', 'caps_intensity', 'caps_temporal_lobe']
caps_types["frequency"] = ['caps_frequency', 'caps_vision_frequency', 'caps_intensity_frequency']
caps_types["max_frequency"] = [
    'caps_maximum_frequency', 'caps_vision_maximum_frequency', 'caps_intensity_maximum_frequency',
]
caps_types["max_weighted_frequency"] = [
    'caps_max_weighted_frequency', 'caps_vision_max_weighted_frequency',
    'caps_intensity_max_weighted_frequency',
]

hppd_variables = [
    'persist_vis_yn',
    'hppd_true_chronicity',
    'persistvis_duration',
    'baggot_total',
]

# ── Labels for variable groupings (eg. to be used in the Regression tables) ────────────────────
iv_type_group_labels = {"sp_predictors": "SP Use Patterns",
                        "vch_behavior": "VCH Task Behavior",
                        "vch_computations": "HGF Estimates"}


# ── Outcome bundles ────────────────────────────────────────────────────────────
#
# WHAT THIS IS:
#   A reference dict that records the canonical results directory, human label,
#   and intended analysis sample for every primary outcome in the manuscript.
#   One entry per outcome (or tightly grouped set of outcomes that share a
#   results subdirectory).
#
# WHERE IT IS USED:
#   1. 04_visualizations/0X_all_figures.py — DV_CONFIGS in the CONFIG section
#      reads "results_dir" and "sample" to know where to write figures and which
#      DataFrame to pass to plotting helpers.  When you add a new DV to
#      DV_CONFIGS you MUST add a matching entry here first; mismatches between
#      outcome_bundles and DV_CONFIGS will cause figures to land in the wrong
#      directory or use the wrong sample.
#   2. Anywhere new scripts route output — always derive the path from
#      outcome_bundles[key]["results_dir"] rather than hardcoding it.
#
# SAMPLE VALUES:
#   "full"           → df       (all QC-passing subjects, n ≈ 228)
#   "sp_users"       → df_sp / df_sp_plot  (psycheduse_yn == "Yes", n ≈ 186)
#   "hppd_subsample" → SP users with persist_vis_yn == 1
#
# HOW TO EXTEND:
#   1. Add an entry here (unique key, correct results_dir, correct sample).
#   2. Add a matching entry to DV_CONFIGS in 04_visualizations/0X_all_figures.py.
#   3. Update 04_visualizations/README.md if the output path is new.
#
outcome_bundles = {
    "hppd_binary": {
        "dvs": ["hppd_binary"],
        "results_dir": "results/hppd_binary",
        "label": "HPPD History",
        "sample": "sp_users",      # SP users only (psycheduse_yn == "Yes")
    },
    "hppd_severity": {
        "dvs": ["hppd_true_chronicity", "persistvis_duration", "baggot_total"],
        "results_dir": "results/hppd_severity",
        "label": "HPPD Severity",
        "sample": "hppd_subsample",
    },
    "caps_vision": {
        "dvs": ["caps_vision"],
        "results_dir": "results/caps_vision",
        "label": "CAPS Visual",
        "sample": "full",
    },
    "caps_total": {
        "dvs": ["caps_total"],
        "results_dir": "results/caps_total",
        "label": "CAPS Total",
        "sample": "full",
    },
    "lshs_total": {
        "dvs": ["lshs_total"],
        "results_dir": "results/lshs_total",
        "label": "LSHS Total",
        "sample": "full",
    },
    "baggot_total": {
        "dvs": ["baggot_total"],
        "results_dir": "results/baggot_total",
        "label": "PPA Sx Count",
        "sample": "sp_users",      # SP users only (psycheduse_yn == "Yes")
    },
}


###############################################################################
# SECTION 3: COLOR PALETTES
###############################################################################

def build_linear_palette(vibrant_rgb, dark_rgb=(0.2, 0.2, 0.2), values=None, n_levels=7, as_hex=True):
    """Build a linear palette from dark -> vibrant using 0-100 value stops."""
    if values is None:
        values = np.linspace(0, 100, n_levels)
    norm_values = np.array(values, dtype=float) / 100.0
    cmap = LinearSegmentedColormap.from_list('custom_linear', [dark_rgb, vibrant_rgb])
    colors = [cmap(v) for v in norm_values]
    if as_hex:
        colors = [to_hex(c) for c in colors]
    return colors


def make_binary_palette(palette, low_index=0, high_index=-1):
    """Pick two colors from a palette for binary plots."""
    return [palette[low_index], palette[high_index]]


ELECTRIC_BLUE_SPEC = {
    'name': 'electric_blue',
    'vibrant_rgb': (0.0, 0.749, 1.0),
    'dark_rgb': (0.2, 0.2, 0.2),
    'values': [5, 25, 40, 55, 70, 90, 100],
}

electric_blue_palette = build_linear_palette(
    ELECTRIC_BLUE_SPEC['vibrant_rgb'],
    dark_rgb=ELECTRIC_BLUE_SPEC['dark_rgb'],
    values=ELECTRIC_BLUE_SPEC['values'],
)

caps_vision_palette = electric_blue_palette
binary_palette = make_binary_palette(electric_blue_palette, low_index=0, high_index=-2)


###############################################################################
# SECTION 4: FLAT LABEL DICTS (backward compatibility)
###############################################################################

# Build flat dicts from VARIABLE_REGISTRY
def _is_nan(v):
    try:
        return isinstance(v, float) and np.isnan(v)
    except Exception:
        return False

dv_to_lab_short = {
    k: v['plot_label']
    for k, v in VARIABLE_REGISTRY.items()
    if not _is_nan(v['plot_label'])
}

dv_to_lab = {
    k: v['plot_label_verbose']
    for k, v in VARIABLE_REGISTRY.items()
    if not _is_nan(v['plot_label_verbose'])
}

dv_to_lab_gen = {
    k: v['plot_label_lay_audience']
    for k, v in VARIABLE_REGISTRY.items()
    if not _is_nan(v['plot_label_lay_audience'])
}

# Add entries from variable_labels.py that were appended after dict construction
dv_to_lab['vch_nu_nominal']       = 'VCH Prior Weighting — nominal (ν)'
dv_to_lab['vch_beta_nominal']     = 'VCH Decision Precision (β)'
dv_to_lab['vch_omega_nominal']    = 'VCH Belief Evolution Rate — nominal (ω)'
dv_to_lab_short['vch_nu_nominal']    = 'Prior Weighting (ν)'
dv_to_lab_short['vch_beta_nominal']  = 'Decision Precision (β)'
dv_to_lab_short['vch_omega_nominal'] = 'Contingency Belief Evolution Rate — nominal (ω)'

# Add _normalized variants for all keys
for _k in list(dv_to_lab_short.keys()):
    if not _k.endswith('_normalized'):
        dv_to_lab_short[f'{_k}_normalized'] = dv_to_lab_short[_k]
for _k in list(dv_to_lab.keys()):
    if not _k.endswith('_normalized'):
        dv_to_lab[f'{_k}_normalized'] = dv_to_lab[_k]

# Backward-compat aliases
dv_to_lab_supershort = dv_to_lab_short   # old name
dv_to_lab_organized  = dv_to_lab         # old name
dv_to_lab_concise    = dv_to_lab_short   # old name

# dv_to_lab_organized is exposed flat (via the alias above) rather than as the
# nested dict some older call sites used.  Nothing in this repository relies on
# the nested form.

