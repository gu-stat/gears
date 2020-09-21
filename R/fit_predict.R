#' Fit and Predict Function
#'
#' For each rolling sample, \code{fit_predict} uses the training data set to
#' fit a GLM. Then, it uses the testing data set to get the predictions.
#'
#' @param forecast.lead A numeric value with the number of the forecast lead.
#' @param eq.number A numeric value with the number of the equation being
#'     estimated.
#' @param sample.number A numeric value with the number of the rolling sample
#'     being used in the estimation/prediction.
#' @param glm.family A description of the error distribution to be used in the
#'     model. Inherits from \code{gears}.
#' @param Y.name Name of the y variable (variable to be forecasted) created by
#'    \code{gears}.
#' @param all.equations.rhs Vector with all possible combinations of right-hand
#'    side variables.
#' @param DF.Fit.Predict a data frame with training data and test data for the
#'    estimation estage.
#' @param ... Other parameters passed on the \link[stats]{glm}.
#'
#' @return Returns a numerical value giving the prediction for the selected
#'     equation, rolling sample, and forecast lead.
#'
#' @keywords internal
#'
fit_predict <- function(forecast.lead, eq.number, sample.number, glm.family,
                        Y.name, all.equations.rhs, DF.Fit.Predict, ...) {

  tmp_lhs      <- paste0(Y.name, "_plus_", forecast.lead)
  tmp_equation <- paste0(tmp_lhs, all.equations.rhs[eq.number])

  tmp_fit <- stats::glm(
    formula = stats::as.formula(tmp_equation),
    family  = glm.family,
    data    = DF.Fit.Predict[[forecast.lead]]$data_fit[[sample.number]],
    ...
  )

  if (length(tmp_fit$coefficients) > tmp_fit$rank) {
    tmp_forecast <- NA
  } else {
    tmp_forecast <- stats::predict(
      object  = tmp_fit,
      newdata = DF.Fit.Predict[[forecast.lead]]$data_predict[[sample.number]],
      type    = "response"
    )
  }

  return(tmp_forecast)
}
