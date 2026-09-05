############################################################################################################################################################################################################################
### HDI helper — 94% highest density interval via bayestestR::hdi()
############################################################################################################################################################################################################################
# Returns the 94% highest density interval of a posterior as a named pair,
# c(lower = ..., upper = ...).
#
# Thin wrapper over bayestestR::hdi(). It exists for two reasons: callers here
# want a plain named numeric pair rather than the one-row data.frame bayestestR
# returns, and callers want a non-finite draw to be an error rather than
# something quietly dropped.
#
# This replaces a hand-rolled Kruschke sliding-window implementation that was
# duplicated across several scripts. The two agree exactly on every posterior
# draw count this pipeline produces: verified bit-identical (max abs diff 0)
# across 79 fixed-effect parameters from 18 saved fits spanning the bernoulli,
# hurdle_negbinomial, and zero_inflated_negbinomial families. The old code
# sized its window with floor(cred_mass * n) where bayestestR uses ceiling(),
# so the two can differ by one draw when length(samples) is not a multiple of
# 50; brms draw counts here are always multiples of 1000, so full posteriors
# are unaffected.
#
# Fails loudly on non-finite input. Posterior draws should never be NA, NaN, or
# Inf; when they are, it signals an upstream problem with the fit (for example a
# shape parameter collapsing toward 0 in a hurdle model). Dropping those draws
# silently would report an interval computed from a smaller, different posterior
# than the one the model actually produced.
# Condition class carried by every compute_hdi_94() failure. Call sites that wrap
# model steps in tryCatch() re-raise on this class rather than logging and moving
# on, so a bad posterior fails the job instead of silently producing no output.
HDI_ERROR_CLASS <- "hdi_integrity_error"

compute_hdi_94 <- function(samples, cred_mass = 0.94) {
  if (!is.numeric(samples)) {
    stop(errorCondition(
      paste0("compute_hdi_94(): expected a numeric vector, got ",
             paste(class(samples), collapse = "/"), "."),
      class = HDI_ERROR_CLASS))
  }

  n_bad <- sum(!is.finite(samples))
  if (n_bad > 0L) {
    stop(errorCondition(sprintf(
      paste0("compute_hdi_94(): %d of %d draws are NA/NaN/Inf. Posterior draws ",
             "should always be finite, so this signals an upstream problem with ",
             "the fit. Investigate the model rather than filtering the draws."),
      n_bad, length(samples)), class = HDI_ERROR_CLASS))
  }

  if (length(samples) < 2L) {
    stop(errorCondition(sprintf(
      "compute_hdi_94(): need at least 2 draws to form an interval, got %d.",
      length(samples)), class = HDI_ERROR_CLASS))
  }

  h <- bayestestR::hdi(samples, ci = cred_mass)

  # bayestestR warns and returns NA when ceiling(ci * n) >= n. That is a real
  # failure for our purposes, so promote it to an error.
  if (!is.finite(h$CI_low) || !is.finite(h$CI_high)) {
    stop(errorCondition(sprintf(
      paste0("compute_hdi_94(): bayestestR::hdi() returned a non-finite interval for ",
             "n = %d draws at ci = %s. The sample is too short for the requested ",
             "credible mass."),
      length(samples), format(cred_mass)), class = HDI_ERROR_CLASS))
  }

  c(lower = h$CI_low, upper = h$CI_high)
}

############################################################################################################################################################################################################################
### Master Processing Function for BRMS Mediation Models
############################################################################################################################################################################################################################
process_and_save_mediation_model <- function(
  fit,
  model_name,
  predictor,
  mediator,
  mediator_non_normalized = NULL,  # For cases where mediator has underscores
  dv,
  spvar = NULL,
  spvar2 = NULL,
  dataframe,
  mediator_name_in_dv = NULL,
  results_base_dir = NULL        # Override output base dir; defaults to RESULTS_DIR

) {

  # If mediator_non_normalized not provided, assume same as mediator
  if (is.null(mediator_non_normalized)) {
    mediator_non_normalized <- mediator
  }

  # Create base directory
  if (is.null(results_base_dir)) results_base_dir <- RESULTS_DIR
  base_dir <- file.path(results_base_dir, model_name)
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  
  message(paste0("\n", paste(rep("=", 80), collapse = "")))
  message(paste0("PROCESSING MODEL: ", model_name))
  message(paste0(paste(rep("=", 80), collapse = "")))
  
  # Get response names (BRMS removes underscores)
  dv_clean <- gsub("_", "", dv)
  mediator_clean <- gsub("_", "", mediator_non_normalized)
  
  ###########################################################
  # 1. SAVE THE FIT
  ###########################################################
  message("\n[1/5] Saving model fit object...")
  tryCatch({
    save(fit, file = file.path(base_dir, paste0("fit_", model_name, ".RData")))
    message("✓ Fit saved successfully")
  }, error = function(e) {
    message(paste0("✗ Failed to save fit: ", e$message))
  })
  
  ###########################################################
  # 2. DHARMA DIAGNOSTICS PLOTS
  ###########################################################
  message("\n[2/5] Creating DHARMa diagnostic plots...")
  
  # Determine whether DV family is integer-valued so DHARMa applies the
  # correct randomized quantile residual computation. Without integer = TRUE,
  # DHARMa treats count/binary responses as continuous, producing incorrect
  # QQ plots and p-values for count and binary DVs.
  INTEGER_DV_FAMILIES <- c("hurdle_negbinomial", "bernoulli", "negbinomial",
                            "zero_inflated_negbinomial", "zero_inflated_poisson", "poisson")
  dv_fam_str <- tryCatch(
    fit$family[[dv_clean]]$family,
    error = function(e) tryCatch(family(fit)$family, error = function(e2) "unknown")
  )
  integer_dv <- isTRUE(any(sapply(INTEGER_DV_FAMILIES,
                                  function(f) grepl(f, dv_fam_str, ignore.case = TRUE))))
  message(paste0("  DV family: ", dv_fam_str, "  →  integer_dv = ", integer_dv))

  # DV diagnostics
  sim_residuals_dv <- NULL
  tryCatch({
    message(paste0("  Creating DHARMa residuals for DV: ", dv_clean))
    sim_residuals_dv <- dh_check_brms(fit, resp = dv_clean, integer = integer_dv)
    message("  ✓ DHARMa residuals created for DV")
    
    # Save comprehensive DHARMa diagnostic plot for DV
    tryCatch({
      png(file.path(base_dir, paste0("dharma_comprehensive_", dv_clean, "_", model_name, ".png")), 
          width = 1200, height = 800)
      plot(sim_residuals_dv)
      dev.off()
      message(paste0("  ✓ Comprehensive DHARMa plot saved for DV: ", dv_clean))
    }, error = function(e) {
      dev.off()
      message(paste0("  ✗ Failed to save comprehensive DHARMa plot for DV: ", e$message))
    })
    
    # Plot vs spvar
    if (!is.null(spvar)) {
      tryCatch({
        png(file.path(base_dir, paste0("dharma_diagnostics_", dv_clean, "_", spvar, "_", model_name, ".png")), 
            width = 800, height = 600)
        plotResiduals(sim_residuals_dv, fit$data[[spvar]])
        dev.off()
        message(paste0("  ✓ DHARMa plot saved: ", dv_clean, " vs ", spvar))
      }, error = function(e) {
        dev.off()
        message(paste0("  ✗ Failed to create DHARMa plot for ", dv_clean, " vs ", spvar, ": ", e$message))
      })
    }
    
    # Plot vs mediator
    tryCatch({
      png(file.path(base_dir, paste0("dharma_diagnostics_", dv_clean, "_", mediator, "_", model_name, ".png")), 
          width = 800, height = 600)
      plotResiduals(sim_residuals_dv, fit$data[[mediator]])
      dev.off()
      message(paste0("  ✓ DHARMa plot saved: ", dv_clean, " vs ", mediator))
    }, error = function(e) {
      dev.off()
      message(paste0("  ✗ Failed to create DHARMa plot for ", dv_clean, " vs ", mediator, ": ", e$message))
    })
    
    # Plot vs spvar2 if provided
    if (!is.null(spvar2)) {
      tryCatch({
        png(file.path(base_dir, paste0("dharma_diagnostics_", dv_clean, "_", spvar2, "_", model_name, ".png")), 
            width = 800, height = 600)
        plotResiduals(sim_residuals_dv, fit$data[[spvar2]])
        dev.off()
        message(paste0("  ✓ DHARMa plot saved: ", dv_clean, " vs ", spvar2))
      }, error = function(e) {
        dev.off()
        message(paste0("  ✗ Failed to create DHARMa plot for ", dv_clean, " vs ", spvar2, ": ", e$message))
      })
    }
    
  }, error = function(e) {
    message(paste0("  ✗ Failed to create DHARMa residuals for DV: ", e$message))
  })
  
  # Mediator diagnostics
  sim_residuals_mediator <- NULL
  tryCatch({
    message(paste0("  Creating DHARMa residuals for mediator: ", mediator_clean))
    sim_residuals_mediator <- dh_check_brms(fit, resp = mediator_clean)
    message("  ✓ DHARMa residuals created for mediator")
    
    # Save comprehensive DHARMa diagnostic plot for mediator
    tryCatch({
      png(file.path(base_dir, paste0("dharma_comprehensive_", mediator_clean, "_", model_name, ".png")), 
          width = 1200, height = 800)
      plot(sim_residuals_mediator)
      dev.off()
      message(paste0("  ✓ Comprehensive DHARMa plot saved for mediator: ", mediator_clean))
    }, error = function(e) {
      dev.off()
      message(paste0("  ✗ Failed to save comprehensive DHARMa plot for mediator: ", e$message))
    })
    
    # Plot vs spvar
    if (!is.null(spvar)) {
      tryCatch({
        png(file.path(base_dir, paste0("dharma_diagnostics_", mediator_clean, "_", spvar, "_", model_name, ".png")), 
            width = 800, height = 600)
        plotResiduals(sim_residuals_mediator, fit$data[[spvar]])
        dev.off()
        message(paste0("  ✓ DHARMa plot saved: ", mediator_clean, " vs ", spvar))
      }, error = function(e) {
        dev.off()
        message(paste0("  ✗ Failed to create DHARMa plot for ", mediator_clean, " vs ", spvar, ": ", e$message))
      })
    }
    
    # Plot vs spvar2 if provided
    if (!is.null(spvar2)) {
      tryCatch({
        png(file.path(base_dir, paste0("dharma_diagnostics_", mediator_clean, "_", spvar2, "_", model_name, ".png")), 
            width = 800, height = 600)
        plotResiduals(sim_residuals_mediator, fit$data[[spvar2]])
        dev.off()
        message(paste0("  ✓ DHARMa plot saved: ", mediator_clean, " vs ", spvar2))
      }, error = function(e) {
        dev.off()
        message(paste0("  ✗ Failed to create DHARMa plot for ", mediator_clean, " vs ", spvar2, ": ", e$message))
      })
    }
    
  }, error = function(e) {
    message(paste0("  ✗ Failed to create DHARMa residuals for mediator: ", e$message))
  })
  
  ###########################################################
  # 3. SUMMARY DATAFRAMES WITH DIAGNOSTICS
  ###########################################################
  message("\n[3/5] Creating summary dataframes with diagnostics...")
  
  # Get basic summary
  summarized <- summary(fit)
  
  # Initialize diagnostics for DV
  diagnostics_list_dv <- list()
  diagnostics_list_dv$response <- dv
  diagnostics_list_dv$mediator <- mediator
  diagnostics_list_dv$N <- nobs(fit)
  diagnostics_list_dv$spvar <- ifelse(is.null(spvar), NA, spvar)
  # Programmatic divergent-transition count via brms::nuts_params() — more reliable than
  # rstan::get_num_divergent(fit$fit), which is deprecated in newer rstan.
  .n_div <- sum(
    brms::nuts_params(fit)$Parameter == "divergent__" &
    brms::nuts_params(fit)$Value      == 1
  )
  diagnostics_list_dv$num_divergents <- .n_div

  # Initialize diagnostics for mediator
  diagnostics_list_mediator <- list()
  diagnostics_list_mediator$response <- mediator
  diagnostics_list_mediator$mediator <- NA
  diagnostics_list_mediator$N <- nobs(fit)
  diagnostics_list_mediator$spvar <- ifelse(is.null(spvar), NA, spvar)
  diagnostics_list_mediator$num_divergents <- .n_div  # same fit object — same divergent count
  
  # DHARMa tests for DV
  if (!is.null(sim_residuals_dv)) {
    tryCatch({
      ks_test <- testUniformity(sim_residuals_dv, plot = FALSE)
      diagnostics_list_dv$dharma_ks_pval <- ks_test$p.value
    }, error = function(e) {
      diagnostics_list_dv$dharma_ks_pval <- NA
    })
    
    tryCatch({
      outlier_test <- testOutliers(sim_residuals_dv, plot = FALSE)
      diagnostics_list_dv$dharma_outlier_pval <- outlier_test$p.value
    }, error = function(e) {
      diagnostics_list_dv$dharma_outlier_pval <- NA
    })
    
    tryCatch({
      dispersion_test <- testDispersion(sim_residuals_dv, plot = FALSE)
      diagnostics_list_dv$dharma_dispersion_pval <- dispersion_test$p.value
    }, error = function(e) {
      diagnostics_list_dv$dharma_dispersion_pval <- NA
    })
    
    tryCatch({
      zeroinflation_test <- testZeroInflation(sim_residuals_dv, plot = FALSE)
      diagnostics_list_dv$dharma_zeroinflation_pval <- zeroinflation_test$p.value
    }, error = function(e) {
      diagnostics_list_dv$dharma_zeroinflation_pval <- NA
    })
    
    tryCatch({
      # DHARMa 0.4.7: testQuantiles(plot=F) does NOT rank-transform the predictor,
      # but plotResiduals() always does before calling testQuantiles internally.
      # A direct call with the raw fittedPredictedResponse gives a p-value that
      # diverges from what any plot would show. Fix: rank-transform fitted values
      # explicitly before passing, replicating what plotResiduals does internally.
      fitted_rk <- rank(sim_residuals_dv$fittedPredictedResponse,
                        ties.method = "average") /
                   length(sim_residuals_dv$fittedPredictedResponse)
      heteroscedasticity_test <- testQuantiles(sim_residuals_dv,
                                               predictor = fitted_rk, plot = FALSE)
      # Combined BH-adjusted p + per-quantile (q25/q50/q75) individual pvals.
      # DHARMa colors a quantile line RED when its individual p < 0.05, even if
      # the combined BH-adjusted p is > 0.05. Saving all three makes the CSV
      # reflect exactly what the image flags visually.
      diagnostics_list_dv$dharma_heteroscedasticity_pval      <- heteroscedasticity_test$p.value
      diagnostics_list_dv$dharma_heteroscedasticity_q25_pval  <- if (!is.null(heteroscedasticity_test$pvals)) heteroscedasticity_test$pvals[1] else NA_real_
      diagnostics_list_dv$dharma_heteroscedasticity_q50_pval  <- if (!is.null(heteroscedasticity_test$pvals)) heteroscedasticity_test$pvals[2] else NA_real_
      diagnostics_list_dv$dharma_heteroscedasticity_q75_pval  <- if (!is.null(heteroscedasticity_test$pvals)) heteroscedasticity_test$pvals[3] else NA_real_
    }, error = function(e) {
      diagnostics_list_dv$dharma_heteroscedasticity_pval     <- NA
      diagnostics_list_dv$dharma_heteroscedasticity_q25_pval <- NA
      diagnostics_list_dv$dharma_heteroscedasticity_q50_pval <- NA
      diagnostics_list_dv$dharma_heteroscedasticity_q75_pval <- NA
    })

    # testQuantiles vs spvar — tests whether residual quantile lines deviate
    # when plotted against the predictor (matches compiled_resid_dv_sp.png).
    # This is a DIFFERENT test from dharma_heteroscedasticity_pval (which tests
    # vs fitted values).  Both are needed because a model can show no issues vs
    # fitted values but significant quantile deviation vs a specific predictor.
    if (!is.null(spvar) && spvar %in% colnames(fit$data)) {
      tryCatch({
        spvar_rk <- rank(fit$data[[spvar]], ties.method = "average") /
                    length(fit$data[[spvar]])
        qt_spvar <- testQuantiles(sim_residuals_dv, predictor = spvar_rk, plot = FALSE)
        diagnostics_list_dv$dharma_quantiles_vs_spvar_pval     <- qt_spvar$p.value
        diagnostics_list_dv$dharma_quantiles_vs_spvar_q25_pval <- if (!is.null(qt_spvar$pvals)) qt_spvar$pvals[1] else NA_real_
        diagnostics_list_dv$dharma_quantiles_vs_spvar_q50_pval <- if (!is.null(qt_spvar$pvals)) qt_spvar$pvals[2] else NA_real_
        diagnostics_list_dv$dharma_quantiles_vs_spvar_q75_pval <- if (!is.null(qt_spvar$pvals)) qt_spvar$pvals[3] else NA_real_
      }, error = function(e) {
        diagnostics_list_dv$dharma_quantiles_vs_spvar_pval     <- NA
        diagnostics_list_dv$dharma_quantiles_vs_spvar_q25_pval <- NA
        diagnostics_list_dv$dharma_quantiles_vs_spvar_q50_pval <- NA
        diagnostics_list_dv$dharma_quantiles_vs_spvar_q75_pval <- NA
      })
    } else {
      diagnostics_list_dv$dharma_quantiles_vs_spvar_pval     <- NA
      diagnostics_list_dv$dharma_quantiles_vs_spvar_q25_pval <- NA
      diagnostics_list_dv$dharma_quantiles_vs_spvar_q50_pval <- NA
      diagnostics_list_dv$dharma_quantiles_vs_spvar_q75_pval <- NA
    }

    # testQuantiles vs mediator — tests whether DV residual quantile lines deviate
    # when plotted against the mediator predictor (matches compiled_resid_dv_med.png).
    {
      med_col_for_qt <- if (mediator %in% colnames(fit$data)) mediator else NULL
      if (!is.null(med_col_for_qt)) {
        tryCatch({
          med_rk <- rank(fit$data[[med_col_for_qt]], ties.method = "average") /
                    length(fit$data[[med_col_for_qt]])
          qt_med <- testQuantiles(sim_residuals_dv, predictor = med_rk, plot = FALSE)
          diagnostics_list_dv$dharma_quantiles_vs_mediator_pval     <- qt_med$p.value
          diagnostics_list_dv$dharma_quantiles_vs_mediator_q25_pval <- if (!is.null(qt_med$pvals)) qt_med$pvals[1] else NA_real_
          diagnostics_list_dv$dharma_quantiles_vs_mediator_q50_pval <- if (!is.null(qt_med$pvals)) qt_med$pvals[2] else NA_real_
          diagnostics_list_dv$dharma_quantiles_vs_mediator_q75_pval <- if (!is.null(qt_med$pvals)) qt_med$pvals[3] else NA_real_
        }, error = function(e) {
          diagnostics_list_dv$dharma_quantiles_vs_mediator_pval     <- NA
          diagnostics_list_dv$dharma_quantiles_vs_mediator_q25_pval <- NA
          diagnostics_list_dv$dharma_quantiles_vs_mediator_q50_pval <- NA
          diagnostics_list_dv$dharma_quantiles_vs_mediator_q75_pval <- NA
        })
      } else {
        diagnostics_list_dv$dharma_quantiles_vs_mediator_pval     <- NA
        diagnostics_list_dv$dharma_quantiles_vs_mediator_q25_pval <- NA
        diagnostics_list_dv$dharma_quantiles_vs_mediator_q50_pval <- NA
        diagnostics_list_dv$dharma_quantiles_vs_mediator_q75_pval <- NA
      }
    }
  } else {
    diagnostics_list_dv$dharma_ks_pval <- NA
    diagnostics_list_dv$dharma_outlier_pval <- NA
    diagnostics_list_dv$dharma_dispersion_pval <- NA
    diagnostics_list_dv$dharma_zeroinflation_pval <- NA
    diagnostics_list_dv$dharma_heteroscedasticity_pval <- NA
    diagnostics_list_dv$dharma_quantiles_vs_spvar_pval <- NA
    diagnostics_list_dv$dharma_quantiles_vs_mediator_pval <- NA
  }

  # DHARMa tests for mediator
  if (!is.null(sim_residuals_mediator)) {
    tryCatch({
      ks_test <- testUniformity(sim_residuals_mediator, plot = FALSE)
      diagnostics_list_mediator$dharma_ks_pval <- ks_test$p.value
    }, error = function(e) {
      diagnostics_list_mediator$dharma_ks_pval <- NA
    })
    
    tryCatch({
      outlier_test <- testOutliers(sim_residuals_mediator, plot = FALSE)
      diagnostics_list_mediator$dharma_outlier_pval <- outlier_test$p.value
    }, error = function(e) {
      diagnostics_list_mediator$dharma_outlier_pval <- NA
    })
    
    tryCatch({
      dispersion_test <- testDispersion(sim_residuals_mediator, plot = FALSE)
      diagnostics_list_mediator$dharma_dispersion_pval <- dispersion_test$p.value
    }, error = function(e) {
      diagnostics_list_mediator$dharma_dispersion_pval <- NA
    })
    
    tryCatch({
      zeroinflation_test <- testZeroInflation(sim_residuals_mediator, plot = FALSE)
      diagnostics_list_mediator$dharma_zeroinflation_pval <- zeroinflation_test$p.value
    }, error = function(e) {
      diagnostics_list_mediator$dharma_zeroinflation_pval <- NA
    })
    
    tryCatch({
      # DHARMa 0.4.7: testQuantiles(plot=F) does NOT rank-transform the predictor,
      # but plotResiduals() always does before calling testQuantiles internally.
      # Fix: rank-transform fitted values explicitly so the CSV p-value matches the plot.
      fitted_rk_med <- rank(sim_residuals_mediator$fittedPredictedResponse,
                            ties.method = "average") /
                       length(sim_residuals_mediator$fittedPredictedResponse)
      heteroscedasticity_test <- testQuantiles(sim_residuals_mediator,
                                               predictor = fitted_rk_med, plot = FALSE)
      # Combined BH-adjusted p + per-quantile (q25/q50/q75) individual pvals.
      diagnostics_list_mediator$dharma_heteroscedasticity_pval      <- heteroscedasticity_test$p.value
      diagnostics_list_mediator$dharma_heteroscedasticity_q25_pval  <- if (!is.null(heteroscedasticity_test$pvals)) heteroscedasticity_test$pvals[1] else NA_real_
      diagnostics_list_mediator$dharma_heteroscedasticity_q50_pval  <- if (!is.null(heteroscedasticity_test$pvals)) heteroscedasticity_test$pvals[2] else NA_real_
      diagnostics_list_mediator$dharma_heteroscedasticity_q75_pval  <- if (!is.null(heteroscedasticity_test$pvals)) heteroscedasticity_test$pvals[3] else NA_real_
    }, error = function(e) {
      diagnostics_list_mediator$dharma_heteroscedasticity_pval     <- NA
      diagnostics_list_mediator$dharma_heteroscedasticity_q25_pval <- NA
      diagnostics_list_mediator$dharma_heteroscedasticity_q50_pval <- NA
      diagnostics_list_mediator$dharma_heteroscedasticity_q75_pval <- NA
    })

    # testQuantiles vs spvar — tests whether mediator residual quantile lines
    # deviate when plotted against the predictor (matches compiled_resid_med_sp.png).
    if (!is.null(spvar) && spvar %in% colnames(fit$data)) {
      tryCatch({
        spvar_rk_med <- rank(fit$data[[spvar]], ties.method = "average") /
                        length(fit$data[[spvar]])
        qt_spvar_med <- testQuantiles(sim_residuals_mediator, predictor = spvar_rk_med, plot = FALSE)
        diagnostics_list_mediator$dharma_quantiles_vs_spvar_pval     <- qt_spvar_med$p.value
        diagnostics_list_mediator$dharma_quantiles_vs_spvar_q25_pval <- if (!is.null(qt_spvar_med$pvals)) qt_spvar_med$pvals[1] else NA_real_
        diagnostics_list_mediator$dharma_quantiles_vs_spvar_q50_pval <- if (!is.null(qt_spvar_med$pvals)) qt_spvar_med$pvals[2] else NA_real_
        diagnostics_list_mediator$dharma_quantiles_vs_spvar_q75_pval <- if (!is.null(qt_spvar_med$pvals)) qt_spvar_med$pvals[3] else NA_real_
      }, error = function(e) {
        diagnostics_list_mediator$dharma_quantiles_vs_spvar_pval     <- NA
        diagnostics_list_mediator$dharma_quantiles_vs_spvar_q25_pval <- NA
        diagnostics_list_mediator$dharma_quantiles_vs_spvar_q50_pval <- NA
        diagnostics_list_mediator$dharma_quantiles_vs_spvar_q75_pval <- NA
      })
    } else {
      diagnostics_list_mediator$dharma_quantiles_vs_spvar_pval     <- NA
      diagnostics_list_mediator$dharma_quantiles_vs_spvar_q25_pval <- NA
      diagnostics_list_mediator$dharma_quantiles_vs_spvar_q50_pval <- NA
      diagnostics_list_mediator$dharma_quantiles_vs_spvar_q75_pval <- NA
    }
  } else {
    diagnostics_list_mediator$dharma_ks_pval <- NA
    diagnostics_list_mediator$dharma_outlier_pval <- NA
    diagnostics_list_mediator$dharma_dispersion_pval <- NA
    diagnostics_list_mediator$dharma_zeroinflation_pval <- NA
    diagnostics_list_mediator$dharma_heteroscedasticity_pval <- NA
    diagnostics_list_mediator$dharma_quantiles_vs_spvar_pval <- NA
  }

    # Create results dataframes
    results_df_dv <- summarized$fixed
    results_df_dv$var <- rownames(results_df_dv)
    results_df_dv$response <- dv
    results_df_dv$mediator <- mediator
    
    # Filter to only DV rows (those starting with dv_clean)
    results_df_dv <- results_df_dv[grepl(paste0("^", dv_clean, "_"), results_df_dv$var), ]
    
    results_df_mediator <- summarized$fixed
    results_df_mediator$var <- rownames(results_df_mediator)
    results_df_mediator$response <- mediator
    results_df_mediator$mediator <- NA
    
    # Filter to only mediator rows (those starting with mediator_clean)
    results_df_mediator <- results_df_mediator[grepl(paste0("^", mediator_clean, "_"), results_df_mediator$var), ]
    
    # Add probabilities from posterior draws
    tryCatch({
      # Fixed-effect draw names.  BOTH prefixes are needed:
      #   b_{resp}_{term}     ordinary population-level coefficients
      #   bsp_{resp}_mo{term} "special" coefficients — monotonic mo() terms
      # Until 2026-09-03 this grepped "^b_" only, which silently excluded every
      # mo() covariate: its row still reached the summary CSV (it comes from
      # summary(fit)$fixed) but arrived with NA for prob_below_0, prob_above_0
      # and BOTH 94% HDI bounds, because the left_join below found no match.
      # The 94% HDI is this project's canonical interval, so that was a real gap
      # rather than a cosmetic one.  nonsp_predictors.R:398,592 already used
      # "^b_|^bsp_"; this brings the mediation helper into line.
      #
      # Stripping the prefix yields the compound name summary(fit)$fixed uses as
      # its rowname ("hppdbinary_mohighest_education_balanced"), which is the
      # join key.  The draw name itself is kept for indexing so a bsp_ parameter
      # is not looked up under a b_ name that does not exist.
      fixed_effect_params <- grep("^b_|^bsp_", variables(fit), value = TRUE)

      # Extract posterior draws
      draws <- as_draws_df(fit)

      # Calculate probabilities and 94% HDI for each fixed effect
      results <- map_df(fixed_effect_params, ~ {
        samples  <- draws[[.x]]
        hdi_vals <- compute_hdi_94(samples)
        data.frame(
          var          = str_remove(.x, "^(bsp_|b_)"),  # compound name; join key
          prob_below_0 = mean(samples < 0),
          prob_above_0 = mean(samples > 0),
          hdi_lower_94 = hdi_vals[["lower"]],
          hdi_upper_94 = hdi_vals[["upper"]]
        )
      })
      
      # Merge with result dataframes by var (compound names will match)
      results_df_dv <- left_join(results_df_dv, results, by = "var")
      results_df_mediator <- left_join(results_df_mediator, results, by = "var")
      
      message("  ✓ Probabilities added to summary dataframes")
    }, error = function(e) {
      # A compute_hdi_94() failure means the posterior itself is bad. Re-raise it so
      # the job fails, instead of logging and exiting 0 with no output written.
      if (inherits(e, HDI_ERROR_CLASS)) stop(e)
      message(paste0("  ✗ Failed to add probabilities: ", e$message))
    })
  
  # Merge with diagnostics
  diagnostics_df_dv <- as.data.frame(diagnostics_list_dv, stringsAsFactors = FALSE)
  diagnostics_df_mediator <- as.data.frame(diagnostics_list_mediator, stringsAsFactors = FALSE)
  
  # Add diagnostics columns to results
  for (col in names(diagnostics_df_dv)) {
    if (!(col %in% names(results_df_dv))) {
      results_df_dv[[col]] <- diagnostics_df_dv[[col]][1]
    }
  }
  
  for (col in names(diagnostics_df_mediator)) {
    if (!(col %in% names(results_df_mediator))) {
      results_df_mediator[[col]] <- diagnostics_df_mediator[[col]][1]
    }
  }
  
  # Save summary dataframes with actual variable names
  tryCatch({
    write.csv(results_df_dv, 
              file.path(base_dir, paste0("summary_", dv, "_", model_name, ".csv")), 
              row.names = FALSE)
    message(paste0("  ✓ Summary dataframe saved: ", dv))
  }, error = function(e) {
    message(paste0("  ✗ Failed to save summary for ", dv, ": ", e$message))
  })
  
  tryCatch({
    write.csv(results_df_mediator, 
              file.path(base_dir, paste0("summary_", mediator, "_", model_name, ".csv")), 
              row.names = FALSE)
    message(paste0("  ✓ Summary dataframe saved: ", mediator))
  }, error = function(e) {
    message(paste0("  ✗ Failed to save summary for ", mediator, ": ", e$message))
  })
  
  ###########################################################
  # 4. MEDIATION EFFECTS
  ###########################################################
  message("\n[4/5] Calculating mediation effects...")
  
  mediation_results_list <- list()
  
  # Calculate mediation for main predictor
  tryCatch({
    # Use mediator_name_in_dv if provided, otherwise use mediator_clean
    med_name_for_dv <- if (!is.null(mediator_name_in_dv)) {
      mediator_name_in_dv  # Don't strip underscores - use as-is
    } else {
      mediator_clean
    }
    message(paste0("  Calculating mediation: ", predictor, " → ", mediator, " → ", dv))
    med_out <- calculate_mediation_effect(
      fit = fit,
      predictor = predictor,
      mediator = mediator_clean,
      dv = dv_clean,
      mediator_name_in_dv = med_name_for_dv
    )
    
    # Save results
    write.csv(med_out$results, 
              file.path(base_dir, paste0("mediation_results_", predictor, "_", model_name, ".csv")), 
              row.names = FALSE)
    
    # Save plot
    ggsave(file.path(base_dir, paste0("mediation_plot_", predictor, "_", model_name, ".png")),
           plot = med_out$plot, width = 10, height = 8, dpi = 300)
    
    mediation_results_list[[predictor]] <- med_out
    message(paste0("  ✓ Mediation results saved for ", predictor))
    
  }, error = function(e) {
    message(paste0("  ✗ Failed to calculate mediation for ", predictor, ": ", e$message))
  })
  
  # Calculate for spvar2 if provided
  if (!is.null(spvar2) && spvar2 != predictor) {
    tryCatch({
      message(paste0("  Calculating mediation: ", spvar2, " → ", mediator, " → ", dv))
      med_out <- calculate_mediation_effect(
        fit = fit,
        predictor = spvar2,
        mediator = mediator_clean,
        dv = dv_clean,
        mediator_name_in_dv = med_name_for_dv
      )
      
      # Save results
      write.csv(med_out$results, 
                file.path(base_dir, paste0("mediation_results_", spvar2, "_", model_name, ".csv")), 
                row.names = FALSE)
      
      # Save plot
      ggsave(file.path(base_dir, paste0("mediation_plot_", spvar2, "_", model_name, ".png")),
             plot = med_out$plot, width = 10, height = 8, dpi = 300)
      
      mediation_results_list[[spvar2]] <- med_out
      message(paste0("  ✓ Mediation results saved for ", spvar2))
      
    }, error = function(e) {
      message(paste0("  ✗ Failed to calculate mediation for ", spvar2, ": ", e$message))
    })
  }
  
  ###########################################################
  # 6. POSTERIOR PREDICTIVE CHECKS
  ###########################################################
  message("\n[5/5] Creating posterior predictive check plots...")
  
  # PP check for DV
  tryCatch({
    pp_dv <- pp_check(fit, resp = dv_clean)
    ggsave(file.path(base_dir, paste0("pp_check_", dv, "_", model_name, ".png")),
           plot = pp_dv, width = 10, height = 8, dpi = 300)
    message(paste0("  ✓ PP check plot saved for ", dv))
  }, error = function(e) {
    message(paste0("  ✗ Failed to create PP check for ", dv, ": ", e$message))
  })
  
  # PP check for mediator
  tryCatch({
    pp_mediator <- pp_check(fit, resp = mediator_clean)
    ggsave(file.path(base_dir, paste0("pp_check_", mediator, "_", model_name, ".png")),
           plot = pp_mediator, width = 10, height = 8, dpi = 300)
    message(paste0("  ✓ PP check plot saved for ", mediator))
  }, error = function(e) {
    message(paste0("  ✗ Failed to create PP check for ", mediator, ": ", e$message))
  })
  
  ###########################################################
  # DONE!
  ###########################################################
  message(paste0("\n", paste(rep("=", 80), collapse = "")))
  message(paste0("✓ PROCESSING COMPLETE FOR: ", model_name))
  message(paste0("All outputs saved to: ", base_dir))
  message(paste0(paste(rep("=", 80), collapse = ""), "\n"))
  
  return(list(
    mediation_results = mediation_results_list,
    diagnostics_dv = diagnostics_df_dv,
    diagnostics_mediator = diagnostics_df_mediator
  ))
}
