#' Estimates the Best Equation =================================================
#'
#' Given the number of the equation that minimizes the selected error measure,
#' \code{fit_best} estimates the equation for each rolling sample and forecast
#' lead.
#'
#' @param forecast.lead Number of the forecast lead.
#' @param eq.number Number (position)of the equation in the vector of all
#'     possible equations.
#' @param sample.number Number of the specific random sample that will be
#'    evaluated.
#' @param glm.family A description of the error distribution to be used in the
#'     model. Inherits from \code{gears}.
#' @param Y.name Name of the y variable (variable to be forecasted) created by
#'    \code{gears}.
#' @param all.equations.rhs Vector with all possible combinations of right-hand
#'    side variables.
#' @param DF.Fit.Predict a data frame with training data and test data for the
#'    estimation estage.
#' @param ... Additional parameters passed to \link[stats]{glm}.
#'
#' @return Returns an object of class inheriting from "glm" (see
#'     \link[stats]{glm}).
#'
#' @keywords internal
#'
fit_best <- function(forecast.lead, eq.number, sample.number, glm.family,
                     Y.name, all.equations.rhs, DF.Fit.Predict, ...) {

  tmp_lhs      <- paste0(Y.name, "_plus_", forecast.lead)
  tmp_equation <- paste0(tmp_lhs, all.equations.rhs[eq.number])

  tmp_fit <- stats::glm(
    formula = stats::as.formula(tmp_equation),
    family  = glm.family,
    data    = DF.Fit.Predict[[forecast.lead]]$data_fit[[sample.number]],
    ...
  )

  return(tmp_fit)

}
