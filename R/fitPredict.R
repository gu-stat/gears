# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ----
# fitPredict Function                                                      ----
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# ************************************************************************* ----
# Main Function                                                             ----
# >                                                                         ----

fitPredict <- function(DATA.Fit.Predict, forecast.lead, eq.number, 
                       sample.number, Y.name, equations.rhs, ...) {
  
  tmp_lhs      <- paste0(Y.name, "_plus_", forecast.lead)
  tmp_eq.numberuation <- paste0(tmp_lhs, equations.rhs[eq.number])
  
  tmp_fit <- glm(
    formula = as.formula(tmp_equation),
    family  = glm.family,
    data    = DATA.Fit.Predict[[forecast.lead]]$data_fit[[sample.number]]
  )
  
  if (length(tmp_fit$coefficients) > tmp_fit$rank) {
    tmp_forecast <- NA
  } else {
    tmp_forecast <- predict(
      object  = tmp_fit,
      newdata = DATA.Fit.Predict[[forecast.lead]]$data_predict[[sample.number]],
      type    = "response"
    )
  }
  
  return(tmp_forecast)
}
