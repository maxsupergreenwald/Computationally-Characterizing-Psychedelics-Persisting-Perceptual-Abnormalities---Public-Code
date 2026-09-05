############################################################################################################################################################################################################################
### function for quickly running relevant dharma diagnostics
############################################################################################################################################################################################################################
dharma_check <- function(fit,resp,mediator,spvar,spvar2=NULL){
  sim_residuals <- dh_check_brms(fit, resp) #persistvisyn #vchthreshold
  
  plotResiduals(sim_residuals, fit$data[[spvar]])
  plotResiduals(sim_residuals, fit$data[[mediator]])
  plotResiduals(sim_residuals, fit$data[[spvar2]])
}
