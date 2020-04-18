# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ----
# fitBest Function                                                          ----
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# ************************************************************************* ----
# Main Function                                                             ----
# >                                                                         ----

fitBest <- function(DATA.Fit.Predict, forecast.lead, eq.number, sample.number, 
                    Y.name, equations.rhs, ...) {
  
  # TODO: CHECK ARGUMENTS OF GLM (...)
  
  # |__ Estimate the Best Equations ============================================
  
    tmp_lhs      <- paste0(Y.name, "_plus_", forecast.lead)
    tmp_equation <- paste0(tmp_lhs, equations.rhs[eq.number[forecast.lead]])
    
    tmp_fit <- glm(
      formula = as.formula(tmp_equation),
      family  = glm.family,
      data    = DATA.Fit.Predict[[forecast.lead]]$data_fit[[sample.number]],
      ...
    )
    
    return(tmp_fit)
    
}
