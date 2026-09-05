############################################################################################################################################################################################################################
### Enhanced posterior distribution plot with 95% CI shading
############################################################################################################################################################################################################################
# Enhanced posterior distribution plot with 95% CI shading
plot_posterior_with_ci <- function(posterior_samples, 
                                   title = "Posterior Distribution",
                                   subtitle = NULL,
                                   xlabel = "Effect Size") {
  
  # Calculate statistics
  mean_val <- mean(posterior_samples)
  ci_lower <- quantile(posterior_samples, 0.025)
  ci_upper <- quantile(posterior_samples, 0.975)
  prob_above <- mean(posterior_samples > 0)
  prob_below <- mean(posterior_samples < 0)
  
  # Get density for plotting
  dens <- density(posterior_samples)
  dens_df <- data.frame(x = dens$x, y = dens$y)
  
  # Separate by sign (for coloring)
  dens_neg <- dens_df[dens_df$x < 0, ]
  dens_pos <- dens_df[dens_df$x >= 0, ]
  
  # Further separate into CI vs tails for each side
  # Negative side
  dens_neg_tail <- dens_neg[dens_neg$x < ci_lower, ]  # Outside CI
  dens_neg_ci <- dens_neg[dens_neg$x >= ci_lower, ]   # Inside CI
  
  # Positive side
  dens_pos_tail <- dens_pos[dens_pos$x > ci_upper, ]  # Outside CI
  dens_pos_ci <- dens_pos[dens_pos$x <= ci_upper, ]   # Inside CI
  
  # Create plot
  p <- ggplot(data.frame(x = posterior_samples), aes(x = x)) +
    
    # Shade negative side - TAIL (light red)
    {if(nrow(dens_neg_tail) > 0) 
      geom_area(data = dens_neg_tail, aes(x = x, y = y), 
                fill = "red", alpha = 0.2)} +
    
    # Shade negative side - CI (darker red)
    {if(nrow(dens_neg_ci) > 0) 
      geom_area(data = dens_neg_ci, aes(x = x, y = y), 
                fill = "red", alpha = 0.5)} +
    
    # Shade positive side - TAIL (light blue)
    {if(nrow(dens_pos_tail) > 0)
      geom_area(data = dens_pos_tail, aes(x = x, y = y), 
                fill = "blue", alpha = 0.2)} +
    
    # Shade positive side - CI (darker blue)
    {if(nrow(dens_pos_ci) > 0)
      geom_area(data = dens_pos_ci, aes(x = x, y = y), 
                fill = "blue", alpha = 0.5)} +
    
    # Add density line on top
    geom_density(fill = NA, color = "black", linewidth = 0.8) +
    
    # Add vertical lines
    geom_vline(xintercept = mean_val, linetype = "dashed", color = "darkblue", linewidth = 1) +
    geom_vline(xintercept = 0, color = "black", linewidth = 1.2) +
    geom_vline(xintercept = ci_lower, linetype = "dotted", color = "darkgray", linewidth = 0.8) +
    geom_vline(xintercept = ci_upper, linetype = "dotted", color = "darkgray", linewidth = 0.8) +
    
    # Add probability annotations
    annotate("text", 
             x = min(dens_df$x), 
             y = max(dens$y) * 0.9,
             label = paste0("P < 0: ", sprintf("%.4f", prob_below)),
             hjust = 0, size = 5, fontface = "bold", color = "darkred") +
    
    annotate("text", 
             x = max(dens_df$x), 
             y = max(dens$y) * 0.9,
             label = paste0("P > 0: ", sprintf("%.4f", prob_above)),
             hjust = 1, size = 5, fontface = "bold", color = "darkblue") +
    
    # Add CI annotation
    annotate("text",
             #x = mean_val, # og placement over line
             x = max(dens_df$x)*0.65,
             y = max(dens$y) * 0.5,
             label = paste0("95% CI: [", sprintf("%.4f", ci_lower), 
                            ", ", sprintf("%.4f", ci_upper), "]"),
             hjust = 0.5, size = 4, color = "black") +
    
    labs(
      title = title,
      subtitle = subtitle,
      x = xlabel,
      y = "Density"
    ) +
    theme_minimal() +
    theme(text = element_text(size = 12))
  
  return(p)
}
