#' Estimates the Best Equation =================================================
#'
#' Given the number of the equation that minimizes the selected error measure,
#' \code{fit_best} estimates the equation for each rolling sample and forecast
#' lead.
#'
#' @param forecast.lead
#' @param eq.number
#' @param sample.number
#' @param ...
#'
#' @return Returns an object of class inheriting from "glm" (see
#'     \link[stats]{glm}).
#'
#' @keywords internal
#'
fit_best <- function(forecast.lead, eq.number, sample.number, ...) {

  tmp_lhs      <- paste0(Y_name, "_plus_", forecast.lead)
  tmp_equation <- paste0(tmp_lhs, all.equations.rhs[eq.number])

  tmp_fit <- stats::glm(
    formula = stats::as.formula(tmp_equation),
    family  = glm.family,
    data    = DF.Fit.Predict[[forecast.lead]]$data_fit[[sample.number]]#,
    #...
  )

  return(tmp_fit)

}
