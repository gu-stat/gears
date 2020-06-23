#' Returns the Forecast Errors for the GEARS Method.
#'
#' Calculates the forecast errors using all methods in
#' \link[gears]{error_measures}: "MSE", "MAD", "sMAPE", "MASE", and "OWA".
#'
#' @param DATA A data frame or a univariate time series.
#' @param forecast.horizon Numeric value giving the number of periods for
#'     forecasting.
#' @param Y_name Character string with the name of the Y (left-hand side)
#'     variable.
#' @param total.equations A numeric value giving the number of equations to be
#'     estimated (i.e., the number of equations from the output of
#'     \link[gears]{all_models_rhs})
#' @param number.rs A numeric value giving the number of rolling samples.
#' @param DF.Fit.Predict A data frame with the training data and testing data.
#'     It is the output of \link[gears]{create_DF_Fit_Predict}.
#' @param forecasts.gears A vector containing the forecasted values using the
#'     GEARS method. It is the output of \code{prediction.gears} inside
#'     \code{gears}.
#' @param names_measures Character vector with the names of the error measures
#'     used in \link[gears]{error_measures}.
#'
#' @return A list. First level of the list represents the different forecast
#'     leads. The table inside each level returns the equation/model number
#'     (row) and the error measure (column), and the values inside the cells
#'     are the forecast errors.
#' @keywords internal
#'
forecast_errors <- function(DATA,
                            forecast.horizon,
                            Y_name,
                            total.equations,
                            number.rs,
                            DF.Fit.Predict,
                            forecasts.gears,
                            names_measures = c("mse", "mad", "smape",
                                               "mase", "owa")
                            ) {
  ## \__ Error for each forecast lead ----------------------------------------
  ## First level of the list represents the different forecast leads.
  ## The table inside each level returns the equation/model number (row) and
  ## the error measure (column), and the values inside the cells are the
  ## forecast errors.
  lapply(
    X = 1:forecast.horizon,
    function(H = X) {

      tmp_lhs <- paste0(Y_name, "_plus_", H)

      tmp_ts_freq <- stats::frequency(DATA)

      #### \____ Error for each equation/model ----
      ### The level of the list represents the different equations.
      ### The table inside each level returns the number of the rolling sample
      ### (column) and the error measure (row), and the value inside the cell
      ### is the forecast error.
      tmp <- lapply(
        X = 1:total.equations,
        function(EQ = X){
          sapply(
            X = 1:number.rs,
            function(RS = X) {

              tmp_in <- DF.Fit.Predict[[H]]$data_fit[[RS]][,tmp_lhs]

              tmp_out <- DF.Fit.Predict[[H]]$data_predict[[RS]][,tmp_lhs]

              tmp_forecasts <- forecasts.gears[[H]][RS, EQ]

              names(tmp_forecasts) <- NULL

              error_measures(
                forecasts        = tmp_forecasts,
                outsample        = tmp_out,
                insample         = tmp_in,
                ts.frequency     = tmp_ts_freq,
                forecast.horizon = H
              )
            }
          )
        }
      )

      ## Create table by equation/model number (row) and error measure (column)
      tmp_errors_by_eq <- sapply(
        X = names_measures,
        function(E = X){
          sapply(X = 1:total.equations, function(X) mean(unlist(tmp[[X]][E, ])))
        }
      )

      ## RETURN ----
      tmp_errors_by_eq
    }
  )
}
