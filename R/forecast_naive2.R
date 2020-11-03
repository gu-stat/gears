#' NAIVE2 Forecasts
#'
#' Used to return the point forecasts using the Naive2 method. Code modified
#' from the M4 Competition's original code.
#'
#' @param ts.data A numeric vector of time series observations (a ts object).
#' @param ts.frequency The frequency of a ts object.
#' @param forecast.horizon A numeric value with the length of the forecast lead.
#' @param alpha.level A numeric value with the alpha level to be used in the
#'     test to detect seasonality. Default is 0.05.
#'
#' @return A ts object of forecasted values obtained using the Naive2 method.
#' @export
#'
#' @examples
#' forecast_naive2(
#' ts.data          = datasets::WWWusage,
#' ts.frequency     = stats::frequency(datasets::WWWusage),
#' forecast.horizon = 10,
#' alpha.level      = 0.05
#' )
#'
#' forecast_naive2(
#'   ts.data          = datasets::AirPassengers,
#'   ts.frequency     = stats::frequency(datasets::AirPassengers),
#'   forecast.horizon = 15,
#'   alpha.level      = 0.05
#' )
#'
#' forecast_naive2(
#'   ts.data          = datasets::EuStockMarkets[, "DAX"],
#'   ts.frequency     = stats::frequency(datasets::EuStockMarkets[, "DAX"]),
#'   forecast.horizon = 40,
#'   alpha.level      = 0.05
#' )
forecast_naive2 <- function(ts.data, ts.frequency, forecast.horizon,
                            alpha.level = 0.05){

  tmp.seasonality.test <- FALSE

  if (stats::is.ts(ts.data) == FALSE){
    tmp.ts.data <- stats::ts(ts.data, frequency = ts.frequency)
  } else {
    tmp.ts.data <- ts.data
  }


  # |__ Seasonality Test =======================================================

  if (ts.frequency > 1){

    tmp.seasonality.test <- seasonality_test(
      ts.data      = tmp.ts.data,
      ts.frequency = ts.frequency,
      alpha.level  = alpha.level
    )

  }

  # |__ Seasonal Adjustment ====================================================

  if (isTRUE(tmp.seasonality.test)){

    tmp.decomp <- stats::decompose(tmp.ts.data, type = "multiplicative")

    deseason.ts.data <- tmp.ts.data / tmp.decomp$seasonal

    tmp.start <- length(tmp.decomp$seasonal) - ts.frequency + 1
    tmp.end   <- length(tmp.decomp$seasonal)

    SIout <- utils::head(
      rep(tmp.decomp$seasonal[tmp.start:tmp.end], forecast.horizon),
      forecast.horizon
    )

  } else {

    deseason.ts.data <- tmp.ts.data

    SIout <- rep(1, forecast.horizon)

  }

  # |__ Forecasts: Naive =======================================================

  last.deseason.ts.data <- utils::tail(deseason.ts.data, 1)

  tmp.forecast.horizon <- forecast.horizon - 1

  tmp.end.original   <- stats::tsp(deseason.ts.data)[2]
  tmp.start.forecast <- tmp.end.original + (1/ts.frequency)

  tmp.end.forecast  <-
    tmp.start.forecast + (tmp.forecast.horizon * 1 / ts.frequency)

  forecast.naive <- stats::ts(
    data      = last.deseason.ts.data,
    start     = tmp.start.forecast,
    end       = tmp.end.forecast,
    frequency = ts.frequency
  )

  # |__ Forecasts: Naive2 ======================================================

  forecast.naive2 <- forecast.naive * SIout

  # |__ RETURN =================================================================
  return(forecast.naive2)
}
