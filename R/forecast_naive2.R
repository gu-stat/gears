# > NAIVE2 Forecasts ###################################################### ----

# Used to return the point forecasts using the Naive2 method.
# Code modified from the M4 Competition's original code.

forecast_naive2 <- function(ts.data, forecast.horizon, alpha.level){

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

  tmp.time  <- tsp(last.deseason.ts.data)[1]
  tmp.start <- tmp.time + (1/tmp.ts.freq)

  tmp.forecast.horizon <- forecast.horizon - 1

  forecast.naive <- ts(
    data      = last.deseason.ts.data,
    start     =  tmp.start,
    end       = tmp.start + (tmp.forecast.horizon * 1 / tmp.ts.freq),
    frequency = tmp.ts.freq
  )

  # |__ Forecasts: Naive2 ======================================================

  forecast.naive2 <- forecast.naive * SIout

  # |__ RETURN =================================================================
  return(forecast.naive2)
}

# Examples

# forecast_naive2(
#   ts.data          = insample,
#   forecast.horizon = 6,
#   alpha.level      = 0.05
# )

forecast_naive2(
  ts.data          = datasets::WWWusage,
  forecast.horizon = 10,
  alpha.level      = 0.05
)

forecast_naive2(
  ts.data          = datasets::AirPassengers,
  forecast.horizon = 10,
  alpha.level      = 0.05
)

forecast_naive2(
  ts.data          = datasets::EuStockMarkets[, "DAX"],
  forecast.horizon = 10,
  alpha.level      = 0.05
)


