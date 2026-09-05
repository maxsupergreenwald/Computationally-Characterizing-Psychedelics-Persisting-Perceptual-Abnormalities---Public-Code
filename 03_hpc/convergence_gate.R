################################################################################
# convergence_gate.R — block summary output from a model that did not converge
################################################################################
#
# WHAT THIS IS FOR
#   A model whose posterior has not converged does not produce numbers worth
#   reporting.  This file decides, from the fit alone, whether a model's SUMMARY
#   RESULT TABLES may be written.  When the gate fails, the tables are withheld
#   and the reason is printed, so a run that produced nothing usable is visible
#   in the SLURM .out file rather than showing up later as an unexplained gap.
#
#   Sourced by both nonsp_predictors.R and hpc_mediation.R.  Copied into every
#   job bundle by generate_hpc_jobs.py; edit 03_hpc/convergence_gate.R, never a
#   copy under output/.
#
# WHAT IS GATED, AND WHAT IS NOT
#   Withheld on failure — the summary result tables:
#     nonsp_predictors.R   results/summary_dfs/{dv}.csv
#                          results/summary_dfs/{dv}_counterfactual.csv
#     hpc_mediation.R      path_coefficients_summary.csv
#                          hu_paths_summary.csv
#                          mc_mediation_summary.csv
#                          path_counterfactual_summary.csv
#
#   Still written on failure — everything needed to work out WHY:
#     the diagnostics CSV, the fit .RData, convergence / PP-check / DHARMa PNGs,
#     and (hpc_mediation.R) the files process_and_save_mediation_model() writes.
#
#   So a failing model leaves a complete forensic record and no reportable
#   numbers.  That is the intent: the absence of a summary table IS the signal.
#
# THRESHOLDS
#   These match 04_hpc/mcmc_validity_review.py (RHAT_MAX / ESS_MIN / DIV_MAX) and
#   the project brms skill.  They are one set of numbers in three places; if you
#   change one, change all three.
#
#     Rhat                   < 1.01     (Vehtari et al. 2021)
#     bulk ESS               >= 1000
#     tail ESS               >= 1000    usually the binding one; it is what makes
#                                       the 94% HDI bounds stable, not just the
#                                       point estimate
#     divergent transitions  == 0       post-warmup, as nuts_params() reports.
#                                       Divergences DURING warmup are normal, are
#                                       not reported by Stan, and are not gated.
#
# WHICH PARAMETERS
#   Rhat and ESS are checked on the FOCAL COEFFICIENTS ONLY — the spvar and
#   mediator paths the model exists to estimate — not on every parameter in the
#   fit.  The caller names them explicitly; see gate_parameters_* below.
#
#   This is deliberately narrower than mcmc_validity_review.py, which reviews
#   every coefficient in the results CSVs.  The two are complementary: this gate
#   stops an unusable model from producing a reportable number at all, and the
#   review still runs over whatever does get written.  A model can therefore pass
#   this gate and still be flagged downstream on a covariate — that is expected.
#
#   The divergence check is not parameter-specific and applies to every model,
#   including covariate-only runs that have no focal coefficient.
#
# A NAMED PARAMETER THAT IS ABSENT IS A FAILURE, NOT A SKIP
#   If a caller names a parameter that is not in the draws, the gate fails with
#   "not found in draws".  A missing focal coefficient means either a naming bug
#   or a model that did not fit what it was asked to fit; neither should quietly
#   produce a summary table.  Callers must therefore only name parameters that
#   genuinely apply to the model at hand — e.g. the hu-submodel paths only when
#   the DV family has a varying hu submodel.
################################################################################

CONVERGENCE_RHAT_MAX <- 1.01   # fail at >= this
CONVERGENCE_ESS_MIN  <- 1000   # fail below this, for bulk AND tail
CONVERGENCE_DIV_MAX  <- 0      # fail above this


# gate_parameters_nonsp(): focal draw-column names for a single-response fit.
#   spvar  focal predictor, verbatim (underscores preserved — brms strips them
#          from RESPONSE names only, and a single-response fit has no prefix)
#   has_hu TRUE when the family carries a varying hu submodel
# Returns character(0) for a covariate-only run (spvar == ""), which leaves the
# Rhat/ESS half of the gate vacuous while the divergence check still applies.
gate_parameters_nonsp <- function(spvar, has_hu = FALSE) {
  if (is.null(spvar) || !nzchar(spvar)) return(character(0))
  pars <- paste0("b_", spvar)
  if (isTRUE(has_hu)) pars <- c(pars, paste0("b_hu_", spvar))
  pars
}


# gate_parameters_mediation(): focal draw-column names for a bivariate fit.
#   a  path  spvar        -> mediator
#   c' path  spvar        -> DV
#   b  path  mediator     -> DV
# plus the hu-submodel counterparts of c' and b when hu varies.
#
# mediator_name_in_brms / dv_brms are the brms-INTERNAL response names (see
# brms_response_map() in hpc_mediation.R).  Predictor names keep their
# underscores; only response names are stripped.
#
# mediator_name_in_outcome_model is the mediator as it appears as a PREDICTOR in
# the outcome model, which is not always the same string as its response name --
# hence the two separate arguments.
gate_parameters_mediation <- function(mediator_name_in_brms, dv_brms, spvar,
                                      mediator_name_in_outcome_model,
                                      has_hu = FALSE) {
  pars <- c(
    paste0("b_", mediator_name_in_brms, "_", spvar),           # a
    paste0("b_", dv_brms, "_", spvar),                          # c'
    paste0("b_", dv_brms, "_", mediator_name_in_outcome_model)  # b
  )
  if (isTRUE(has_hu)) {
    pars <- c(pars,
      paste0("b_hu_", dv_brms, "_", spvar),                          # c' on the hurdle
      paste0("b_hu_", dv_brms, "_", mediator_name_in_outcome_model)) # b  on the hurdle
  }
  pars
}


# check_convergence_gate(): evaluate the gate.
#   fit         a fitted brmsfit
#   parameters  focal draw-column names, from one of the helpers above
# Returns a list:
#   pass     TRUE if the model may write its summary tables
#   reasons  character vector, one entry per failure, empty when pass
#   detail   data.frame of the per-parameter numbers (for the console report)
#   n_div    post-warmup divergent transitions
check_convergence_gate <- function(fit, parameters) {
  reasons <- character(0)

  ## Divergences first: not parameter-specific, and applies to every model.
  nuts  <- brms::nuts_params(fit)
  n_div <- sum(nuts$Parameter == "divergent__" & nuts$Value == 1)
  if (n_div > CONVERGENCE_DIV_MAX) {
    reasons <- c(reasons, sprintf(
      "%d divergent transition%s after warmup (allowed: %d) — the sampler could not explore the posterior geometry",
      n_div, if (n_div == 1L) "" else "s", CONVERGENCE_DIV_MAX))
  }

  detail <- NULL
  if (length(parameters) > 0) {
    draws   <- posterior::as_draws_df(fit)
    present <- parameters %in% names(draws)

    if (any(!present)) {
      reasons <- c(reasons, sprintf(
        "focal parameter%s not found in draws: %s",
        if (sum(!present) == 1L) "" else "s",
        paste(parameters[!present], collapse = ", ")))
    }

    if (any(present)) {
      ## summarise_draws() on just the focal columns. rhat / ess_bulk / ess_tail
      ## come straight from posterior; nothing is recomputed here.
      detail <- posterior::summarise_draws(
        posterior::subset_draws(draws, variable = parameters[present]),
        "rhat", "ess_bulk", "ess_tail")

      bad_rhat <- detail$variable[!is.na(detail$rhat) & detail$rhat >= CONVERGENCE_RHAT_MAX]
      if (length(bad_rhat) > 0) {
        reasons <- c(reasons, sprintf(
          "Rhat >= %s on %s (max %.4f) — chains have not mixed",
          format(CONVERGENCE_RHAT_MAX), paste(bad_rhat, collapse = ", "),
          max(detail$rhat[detail$variable %in% bad_rhat])))
      }

      bad_bulk <- detail$variable[!is.na(detail$ess_bulk) & detail$ess_bulk < CONVERGENCE_ESS_MIN]
      if (length(bad_bulk) > 0) {
        reasons <- c(reasons, sprintf(
          "bulk ESS < %d on %s (min %.0f) — too few effective draws for a stable point estimate",
          CONVERGENCE_ESS_MIN, paste(bad_bulk, collapse = ", "),
          min(detail$ess_bulk[detail$variable %in% bad_bulk])))
      }

      bad_tail <- detail$variable[!is.na(detail$ess_tail) & detail$ess_tail < CONVERGENCE_ESS_MIN]
      if (length(bad_tail) > 0) {
        reasons <- c(reasons, sprintf(
          "tail ESS < %d on %s (min %.0f) — the 94%% HDI bounds are not stable",
          CONVERGENCE_ESS_MIN, paste(bad_tail, collapse = ", "),
          min(detail$ess_tail[detail$variable %in% bad_tail])))
      }
    }
  }

  list(pass = length(reasons) == 0, reasons = reasons, detail = detail, n_div = n_div)
}


# report_convergence_gate(): print the verdict.  Always called, pass or fail, so
# the .out file records the numbers either way.  `emit` is print() in
# nonsp_predictors.R and message() in hpc_mediation.R, matching each script's
# existing logging.
report_convergence_gate <- function(gate, label, emit = message) {
  bar <- paste(rep("*", 78), collapse = "")

  if (!is.null(gate$detail)) {
    emit(sprintf("[convergence gate] focal parameters for %s:", label))
    for (i in seq_len(nrow(gate$detail))) {
      emit(sprintf("[convergence gate]   %-46s rhat=%.4f  bulk_ESS=%7.0f  tail_ESS=%7.0f",
                   gate$detail$variable[i], gate$detail$rhat[i],
                   gate$detail$ess_bulk[i], gate$detail$ess_tail[i]))
    }
  }
  emit(sprintf("[convergence gate] divergent transitions after warmup: %d", gate$n_div))

  if (isTRUE(gate$pass)) {
    emit(sprintf("[convergence gate] PASS — %s: summary tables will be written.", label))
    return(invisible(TRUE))
  }

  emit(bar)
  emit(sprintf("** CONVERGENCE GATE FAILED — %s", label))
  emit("** SUMMARY RESULT TABLES WILL NOT BE WRITTEN. Reasons:")
  for (r in gate$reasons) emit(paste0("**   - ", r))
  emit("** Thresholds: Rhat < 1.01, bulk and tail ESS >= 1000, 0 divergent transitions.")
  emit("** The fit, the diagnostics CSV and all diagnostic plots WERE still written.")
  emit("** Re-fit this model; do not report it. See the brms-convergence skill for what to change.")
  emit(bar)
  invisible(FALSE)
}
