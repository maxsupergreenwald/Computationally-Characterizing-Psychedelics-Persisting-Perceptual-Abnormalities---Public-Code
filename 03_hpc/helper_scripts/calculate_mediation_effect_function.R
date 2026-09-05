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
### Mediation Effect Calculator
############################################################################################################################################################################################################################
calculate_mediation_effect <- function(fit, predictor, mediator, dv, mediator_name_in_dv = NULL) {
  
  # If mediator_name_in_dv is not provided, use the same name as mediator
  if (is.null(mediator_name_in_dv)) {
    mediator_name_in_dv <- mediator
  }
  
  message("===== CALCULATING MEDIATION EFFECT =====")
  message(paste0("Predictor: ", predictor))
  message(paste0("Mediator: ", mediator))
  message(paste0("Mediator name in DV formula: ", mediator_name_in_dv))
  message(paste0("DV: ", dv))
  
  # Extract posterior draws
  draws <- as_draws_df(fit)
  
  # Get all coefficient names
  all_coefs <- grep("^b_", names(draws), value = TRUE)
  
  # Helper function to find coefficient by fuzzy matching
  # BRMS converts variable names: e.g., "vch_nu_normalized" stays as "vch_nu_normalized"
  # but user might provide "vchnunormalized"
  find_coef <- function(response, predictor_var, all_names) {
    # Try exact match first
    exact_name <- paste0("b_", response, "_", predictor_var)
    if (exact_name %in% all_names) {
      return(exact_name)
    }
    
    # Try pattern matching - look for coefficients that:
    # 1. Start with "b_response_"
    # 2. Contain all the characters from predictor_var (ignoring underscores)
    pattern_start <- paste0("^b_", response, "_")
    candidates <- grep(pattern_start, all_names, value = TRUE)
    
    # Remove underscores and compare
    predictor_clean <- gsub("_", "", tolower(predictor_var))
    
    for (candidate in candidates) {
      # Extract just the predictor part after "b_response_"
      pred_part <- sub(pattern_start, "", candidate)
      pred_clean <- gsub("_", "", tolower(pred_part))
      
      if (pred_clean == predictor_clean) {
        return(candidate)
      }
    }
    
    return(NULL)
  }
  
  # Find the coefficients
  a_path_name <- find_coef(mediator, predictor, all_coefs)
  c_prime_name <- find_coef(dv, predictor, all_coefs)

  # Construct b_path using dv (already cleaned, no underscores) 
  # and mediator_name_in_dv (try WITH underscores first)
  b_path_name_with_underscores <- paste0("b_", dv, "_", mediator_name_in_dv)

  # Try to find the b_path coefficient with fallback logic
  b_path_name <- NULL
  if (b_path_name_with_underscores %in% names(draws)) {
    # Success with underscores preserved (this is the expected case)
    b_path_name <- b_path_name_with_underscores
    message(paste0("  b path found (with underscores): ", b_path_name))
  } else {
    # Fallback: try with underscores removed from mediator name
    mediator_name_in_dv_clean <- gsub("_", "", mediator_name_in_dv)
    b_path_name_clean <- paste0("b_", dv, "_", mediator_name_in_dv_clean)
    
    if (b_path_name_clean %in% names(draws)) {
      # Success with cleaned version
      b_path_name <- b_path_name_clean
      message(paste0("  b path found (no underscores): ", b_path_name))
    } else {
      # Neither version exists - print debug info and fail
      message("\n===== B PATH COEFFICIENT NOT FOUND =====")
      message(paste0("Tried with underscores: ", b_path_name_with_underscores))
      message(paste0("Tried without underscores: ", b_path_name_clean))
      message(paste0("\nAvailable coefficients starting with 'b_", dv, "_':"))
      dv_coefs <- grep(paste0("^b_", dv, "_"), names(draws), value = TRUE)
      message(paste0("  ", paste(dv_coefs, collapse = "\n  ")))
      message(paste0("\nAll available coefficients:"))
      message(paste0("  ", paste(all_coefs, collapse = "\n  ")))
      
      stop(paste0("Could not find 'b' path coefficient for: ", dv, " ~ ", mediator_name_in_dv, 
                  "\nTried: ", b_path_name_with_underscores, " and ", b_path_name_clean))
    }
  }
  
  message(paste0("Found coefficients:"))
  message(paste0("  a path: ", ifelse(is.null(a_path_name), "NOT FOUND", a_path_name)))
  message(paste0("  b path: ", ifelse(is.null(b_path_name), "NOT FOUND", b_path_name)))
  message(paste0("  c' path: ", ifelse(is.null(c_prime_name), "NOT FOUND", c_prime_name)))
  
  # Check if coefficients exist
  if (is.null(a_path_name)) {
    stop(paste0("Could not find 'a' path coefficient for: ", mediator, " ~ ", predictor, "\n",
                "Available coefficients: ", paste(all_coefs, collapse = ", ")))
  }
  if (is.null(b_path_name)) {
    stop(paste0("Could not find 'b' path coefficient for: ", dv, " ~ ", mediator_name_in_dv, "\n",
                "Available coefficients: ", paste(all_coefs, collapse = ", ")))
  }
  if (is.null(c_prime_name)) {
    stop(paste0("Could not find 'c' path coefficient for: ", dv, " ~ ", predictor, "\n",
                "Available coefficients: ", paste(all_coefs, collapse = ", ")))
  }
  
  # Try to find the b_path coefficient with fallback logic
  b_path_name <- NULL
  if (b_path_name_with_underscores %in% names(draws)) {
    # Success with underscores preserved (this is the expected case)
    b_path_name <- b_path_name_with_underscores
    message(paste0("  b path found (with underscores): ", b_path_name))
  } else {
    # Fallback: try with underscores removed from mediator name
    mediator_name_in_dv_clean <- gsub("_", "", mediator_name_in_dv)
    b_path_name_clean <- paste0("b_", dv, "_", mediator_name_in_dv_clean)
    
    if (b_path_name_clean %in% names(draws)) {
      # Success with cleaned version
      b_path_name <- b_path_name_clean
      message(paste0("  b path found (no underscores): ", b_path_name))
    } else {
      # Neither version exists - print debug info and fail
      message("\n===== B PATH COEFFICIENT NOT FOUND =====")
      message(paste0("Tried with underscores: ", b_path_name_with_underscores))
      message(paste0("Tried without underscores: ", b_path_name_clean))
      message(paste0("\nAvailable coefficients starting with 'b_", dv, "_':"))
      dv_coefs <- grep(paste0("^b_", dv, "_"), names(draws), value = TRUE)
      message(paste0("  ", paste(dv_coefs, collapse = "\n  ")))
      message(paste0("\nAll available coefficients:"))
      message(paste0("  ", paste(all_coefs, collapse = "\n  ")))
      
      stop(paste0("Could not find 'b' path coefficient for: ", dv, " ~ ", mediator_name_in_dv, 
                  "\nTried: ", b_path_name_with_underscores, " and ", b_path_name_clean))
    }
  }
    # Try to find the b_path coefficient with fallback logic
  b_path_name <- NULL
  if (b_path_name_with_underscores %in% names(draws)) {
    # Success with underscores preserved (this is the expected case)
    b_path_name <- b_path_name_with_underscores
    message(paste0("  b path found (with underscores): ", b_path_name))
  } else {
    # Fallback: try with underscores removed from mediator name
    mediator_name_in_dv_clean <- gsub("_", "", mediator_name_in_dv)
    b_path_name_clean <- paste0("b_", dv, "_", mediator_name_in_dv_clean)
    
    if (b_path_name_clean %in% names(draws)) {
      # Success with cleaned version
      b_path_name <- b_path_name_clean
      message(paste0("  b path found (no underscores): ", b_path_name))
    } else {
      # Neither version exists - print debug info and fail
      message("\n===== B PATH COEFFICIENT NOT FOUND =====")
      message(paste0("Tried with underscores: ", b_path_name_with_underscores))
      message(paste0("Tried without underscores: ", b_path_name_clean))
      message(paste0("\nAvailable coefficients starting with 'b_", dv, "_':"))
      dv_coefs <- grep(paste0("^b_", dv, "_"), names(draws), value = TRUE)
      message(paste0("  ", paste(dv_coefs, collapse = "\n  ")))
      message(paste0("\nAll available coefficients:"))
      message(paste0("  ", paste(all_coefs, collapse = "\n  ")))
      
      stop(paste0("Could not find 'b' path coefficient for: ", dv, " ~ ", mediator_name_in_dv, 
                  "\nTried: ", b_path_name_with_underscores, " and ", b_path_name_clean))
    }
  }
  # Extract the paths
  a_path <- draws[[a_path_name]]
  b_path <- draws[[b_path_name]]
  c_prime <- draws[[c_prime_name]]
  # Try to find the b_path coefficient with fallback logic
  b_path_name <- NULL
  if (b_path_name_with_underscores %in% names(draws)) {
    # Success with underscores preserved (this is the expected case)
    b_path_name <- b_path_name_with_underscores
    message(paste0("  b path found (with underscores): ", b_path_name))
  } else {
    # Fallback: try with underscores removed from mediator name
    mediator_name_in_dv_clean <- gsub("_", "", mediator_name_in_dv)
    b_path_name_clean <- paste0("b_", dv, "_", mediator_name_in_dv_clean)
    
    if (b_path_name_clean %in% names(draws)) {
      # Success with cleaned version
      b_path_name <- b_path_name_clean
      message(paste0("  b path found (no underscores): ", b_path_name))
    } else {
      # Neither version exists - print debug info and fail
      message("\n===== B PATH COEFFICIENT NOT FOUND =====")
      message(paste0("Tried with underscores: ", b_path_name_with_underscores))
      message(paste0("Tried without underscores: ", b_path_name_clean))
      message(paste0("\nAvailable coefficients starting with 'b_", dv, "_':"))
      dv_coefs <- grep(paste0("^b_", dv, "_"), names(draws), value = TRUE)
      message(paste0("  ", paste(dv_coefs, collapse = "\n  ")))
      message(paste0("\nAll available coefficients:"))
      message(paste0("  ", paste(all_coefs, collapse = "\n  ")))
      
      stop(paste0("Could not find 'b' path coefficient for: ", dv, " ~ ", mediator_name_in_dv, 
                  "\nTried: ", b_path_name_with_underscores, " and ", b_path_name_clean))
    }
  }
  # Calculate indirect effect (a * b)
  indirect_effect <- a_path * b_path
  
  # Calculate total effect (indirect + direct)
  total_effect <- indirect_effect + c_prime
  
  # Calculate proportion mediated
  prop_mediated <- indirect_effect / total_effect
  
  # Summarize results
  results <- data.frame(
      predictor = predictor,
      mediator = mediator,
      dv = dv,
      effect = c("a_path", "b_path", "c_prime_direct", "indirect_ab", "total_effect", "prop_mediated"),
      estimate = c(
        mean(a_path),
        mean(b_path),
        mean(c_prime),
        mean(indirect_effect),
        mean(total_effect),
        mean(prop_mediated, na.rm = TRUE)
      ),
      lower_95 = c(
        quantile(a_path, 0.025),
        quantile(b_path, 0.025),
        quantile(c_prime, 0.025),
        quantile(indirect_effect, 0.025),
        quantile(total_effect, 0.025),
        quantile(prop_mediated, 0.025, na.rm = TRUE)
      ),
      upper_95 = c(
        quantile(a_path, 0.975),
        quantile(b_path, 0.975),
        quantile(c_prime, 0.975),
        quantile(indirect_effect, 0.975),
        quantile(total_effect, 0.975),
        quantile(prop_mediated, 0.975, na.rm = TRUE)
      ),
      lower_94_hdi = c(
        compute_hdi_94(a_path)[["lower"]],
        compute_hdi_94(b_path)[["lower"]],
        compute_hdi_94(c_prime)[["lower"]],
        compute_hdi_94(indirect_effect)[["lower"]],
        compute_hdi_94(total_effect)[["lower"]],
        compute_hdi_94(prop_mediated)[["lower"]]
      ),
      upper_94_hdi = c(
        compute_hdi_94(a_path)[["upper"]],
        compute_hdi_94(b_path)[["upper"]],
        compute_hdi_94(c_prime)[["upper"]],
        compute_hdi_94(indirect_effect)[["upper"]],
        compute_hdi_94(total_effect)[["upper"]],
        compute_hdi_94(prop_mediated)[["upper"]]
      ),
      prob_above_0 = c(
        mean(a_path > 0),
        mean(b_path > 0),
        mean(c_prime > 0),
        mean(indirect_effect > 0),
        mean(total_effect > 0),
        mean(prop_mediated > 0, na.rm = TRUE)
      ),
      prob_below_0 = c(
        mean(a_path < 0),
        mean(b_path < 0),
        mean(c_prime < 0),
        mean(indirect_effect < 0),
        mean(total_effect < 0),
        mean(prop_mediated < 0, na.rm = TRUE)
      )
    )
  
  message("\n===== MEDIATION RESULTS =====")
  print(results)
  
  # Plot the indirect effect distribution
  # p <- ggplot(data.frame(indirect = indirect_effect), aes(x = indirect)) +
  #   geom_density(fill = "skyblue", alpha = 0.5) +
  #   geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  #   geom_vline(xintercept = mean(indirect_effect), color = "blue", linewidth = 1) +
  #   labs(
  #     title = paste0("Indirect Effect Distribution: ", predictor, " → ", mediator, " → ", dv),
  #     x = "Indirect Effect (a × b)",
  #     y = "Density",
  #     subtitle = paste0("Mean = ", round(mean(indirect_effect), 4), 
  #                       " | 95% CI [", round(quantile(indirect_effect, 0.025), 4), 
  #                       ", ", round(quantile(indirect_effect, 0.975), 4), "]")
  #   ) +
  #   theme_minimal()
  # 
  # print(p)
  # Plot the indirect effect distribution with enhanced styling
  p <- ggplot(data.frame(x = indirect_effect), aes(x = x)) +
    {
      # Calculate density
      dens <- density(indirect_effect)
      dens_df <- data.frame(x = dens$x, y = dens$y)
      ci_lower <- quantile(indirect_effect, 0.025)
      ci_upper <- quantile(indirect_effect, 0.975)
      
      # Separate by sign
      dens_neg <- dens_df[dens_df$x < 0, ]
      dens_pos <- dens_df[dens_df$x >= 0, ]
      
      # Further separate into CI vs tails
      dens_neg_tail <- dens_neg[dens_neg$x < ci_lower, ]
      dens_neg_ci <- dens_neg[dens_neg$x >= ci_lower, ]
      dens_pos_tail <- dens_pos[dens_pos$x > ci_upper, ]
      dens_pos_ci <- dens_pos[dens_pos$x <= ci_upper, ]
      
      list(
        # Shade negative side - TAIL (light red)
        if(nrow(dens_neg_tail) > 0) 
          geom_area(data = dens_neg_tail, aes(x = x, y = y), 
                    fill = "red", alpha = 0.2),
        
        # Shade negative side - CI (darker red)
        if(nrow(dens_neg_ci) > 0) 
          geom_area(data = dens_neg_ci, aes(x = x, y = y), 
                    fill = "red", alpha = 0.5),
        
        # Shade positive side - TAIL (light blue)
        if(nrow(dens_pos_tail) > 0)
          geom_area(data = dens_pos_tail, aes(x = x, y = y), 
                    fill = "blue", alpha = 0.2),
        
        # Shade positive side - CI (darker blue)
        if(nrow(dens_pos_ci) > 0)
          geom_area(data = dens_pos_ci, aes(x = x, y = y), 
                    fill = "blue", alpha = 0.5)
      )
    } +
    
    # Add density line on top
    geom_density(fill = NA, color = "black", linewidth = 0.8) +
    
    # Add vertical lines
    geom_vline(xintercept = mean(indirect_effect), linetype = "dashed", color = "darkblue", linewidth = 1) +
    geom_vline(xintercept = 0, color = "black", linewidth = 1.2) +
    geom_vline(xintercept = quantile(indirect_effect, 0.025), linetype = "dotted", color = "darkgray", linewidth = 0.8) +
    geom_vline(xintercept = quantile(indirect_effect, 0.975), linetype = "dotted", color = "darkgray", linewidth = 0.8) +
    
    # Add probability annotations
    annotate("text", 
             x = min(density(indirect_effect)$x), 
             y = max(density(indirect_effect)$y) * 0.9,
             label = paste0("P < 0: ", sprintf("%.4f", mean(indirect_effect < 0))),
             hjust = 0, size = 5, fontface = "bold", color = "darkred") +
    
    annotate("text", 
             x = max(density(indirect_effect)$x), 
             y = max(density(indirect_effect)$y) * 0.9,
             label = paste0("P > 0: ", sprintf("%.4f", mean(indirect_effect > 0))),
             hjust = 1, size = 5, fontface = "bold", color = "darkblue") +
    
    # Add CI annotation
    annotate("text",
             x = max(density(indirect_effect)$x) * 0.65,
             y = max(density(indirect_effect)$y) * 0.5,
             label = paste0("95% CI: [", sprintf("%.4f", quantile(indirect_effect, 0.025)), 
                            ", ", sprintf("%.4f", quantile(indirect_effect, 0.975)), "]"),
             hjust = 0.5, size = 4, color = "black") +
    
    labs(
      title = paste0("Indirect Effect Distribution: ", predictor, " → ", mediator, " → ", dv),
      x = "Indirect Effect (a × b)",
      y = "Density"
    ) +
    theme_minimal() +
    theme(text = element_text(size = 12))
  
  print(p)
  
  return(list(
    results = results,
    draws = data.frame(
      a_path = a_path,
      b_path = b_path,
      c_prime = c_prime,
      indirect_effect = indirect_effect,
      total_effect = total_effect,
      prop_mediated = prop_mediated
    ),
    plot = p
  ))
}

