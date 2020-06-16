# > NAIVE2 Forecasts ###################################################### ----

# Used to return the point forecasts using the Naive2 method.
# Code modified from the M4 Competition's original code.

#' Title
#'
#' @param ts.data A numeric vector of time series observations (a ts object).
#' @param forecast.horizon A numeric value with the length of the forecast lead.
#' @param alpha.level A numeric value with the alpha level to be used in the
#'     test to detect seasonality. Default is 0.05.
#'
#' @return A ts object of forecasted values obtained using the Naive2 method.
#'
#' @examples
#' forecast_naive2(
#' ts.data          = datasets::WWWusage,
#' forecast.horizon = 10,
#' alpha.level      = 0.05
#' )
#'
#' forecast_naive2(
#'   ts.data          = datasets::AirPassengers,
#'   forecast.horizon = 15,
#'   alpha.level      = 0.05
#' )
#'
#' forecast_naive2(
#'   ts.data          = datasets::EuStockMarkets[, "DAX"],
#'   forecast.horizon = 40,
#'   alpha.level      = 0.05
#' )
forecast_naive2 <- function(ts.data, forecast.horizon, alpha.level = 0.05){

  tmp.ts.freq <- frequency(ts.data)

  tmp.seasonality.test <- FALSE

  # |__ Seasonality Test =======================================================

  if (tmp.ts.freq > 1){

    tmp.seasonality.test <- seasonality_test(
      ts.data     = ts.data,
      ts.freq     = tmp.ts.freq,
      alpha.level = alpha.level
    )

  }

  # |__ Seasonal Adjustment ====================================================

  if (isTRUE(tmp.seasonality.test)){

    tmp.decomp <- decompose(ts.data, type = "multiplicative")

    deseason.ts.data <- ts.data/tmp.decomp$seasonal

    tmp.start <- length(tmp.decomp$seasonal) - tmp.ts.freq + 1
    tmp.end   <- length(tmp.decomp$seasonal)

    SIout <- head(
      rep(tmp.decomp$seasonal[tmp.start:tmp.end], forecast.horizon),
      forecast.horizon
    )

  } else {

    deseason.ts.data <- ts.data

    SIout <- rep(1, forecast.horizon)

  }

  # |__ Forecasts: Naive =======================================================

  last.deseason.ts.data <- tail(deseason.ts.data, 1)

  tmp.forecast.horizon <- forecast.horizon - 1

  tmp.time  <- tsp(deseason.ts.data)[2]
  tmp.start <- tmp.time + (1/tmp.ts.freq)
  tmp.end   <- tmp.start + (tmp.forecast.horizon * 1 / tmp.ts.freq)

  forecast.naive <- ts(
    data      = last.deseason.ts.data,
    start     = tmp.start,
    end       = tmp.end,
    frequency = tmp.ts.freq
  )

  # |__ Forecasts: Naive2 ======================================================

  forecast.naive2 <- forecast.naive * SIout

  # |__ RETURN =================================================================
  return(forecast.naive2)
}
