#' Functions to calculate the forecast errors.
#'
#' \code{error_functions} returns the forecast errors for the selected error
#' measure.
#'
#' @param error.measure A character string with the name of the error measure.
#'     These can be 'mse', 'mad', 'smape', 'mase', and 'owa' (see Details).
#' @param forecasts A numeric vector with the point forecasts.
#' @param outsample A numeric vector with the test data set.
#' @param insample A numeric vector with the train data set.
#' @param forecast.horizon A numeric value with the length of the forecast lead.
#' @param alpha.level A numeric value with the alpha level to be used in the
#'     test to detect seasonality. Default is 0.05.
#'
#' @return A numeric value with the forecast error.
#' @export
#'
# @import forecast_naive2
#'
#' @examples
#' # Using NAIVE2:
#' # Observations until the 100th will be in the insample (train) data set.
#' tmp.cut.at <- 100
#'
#' tmp.forecast.horizon <- length(datasets::AirPassengers) - tmp.cut.at
#'
#' tmp.orig.start <- tsp(datasets::AirPassengers)[1]
#' tmp.orig.end   <- tsp(datasets::AirPassengers)[2]
#' tmp.orig.freq  <- tsp(datasets::AirPassengers)[3]
#'
#'
#' # Get train data (insample)
#'
#' tmp.train.start <- tmp.orig.start
#' tmp.train.end   <- tmp.orig.start + ((tmp.cut.at - 1) * 1 / tmp.orig.freq)
#'
#' tmp.train.data  <- window(
#'   x         = datasets::AirPassengers,
#'   start     = tmp.train.start,
#'   end       = tmp.train.end,
#'   frequency = tmp.orig.freq
#' )
#'
#' # Get test data (outsample)
#'
#' tmp.test.start <- tmp.orig.start + (tmp.cut.at * 1 / tmp.orig.freq)
#' tmp.test.end   <- tmp.orig.end
#'
#' tmp.test.data  <- window(
#'   x         = datasets::AirPassengers,
#'   start     = tmp.test.start,
#'   end       = tmp.test.end,
#'   frequency = tmp.orig.freq
#' )
#'
#' # Get Forecasts
#' tmp.forecasts <- forecast_naive2(
#'   ts.data          = tmp.train.data,
#'   forecast.horizon = tmp.forecast.horizon,
#'   alpha.level      = 0.05
#' )
#'
#' # Get MSE
#'
#' error_functions(
#'   error.measure      = "mse",
#'   forecasts          = tmp.forecasts,
#'   outsample          = tmp.test.data,
#'   #'insample           = tmp.test.data,
#'   forecast.horizon   = tmp.forecast.horizon,
#'   alpha.level        = 0.05
#' )
#'
#' # Get OWA
#'
#' error_functions(
#'   error.measure      = "owa",
#'   forecasts          = tmp.forecasts,
#'   outsample          = tmp.test.data,
#'   insample           = tmp.test.data,
#'   forecast.horizon   = tmp.forecast.horizon,
#'   alpha.level        = 0.05
#' )
error_functions <- function(error.measure,
                            forecasts,
                            outsample,
                            insample = NULL,
                            forecast.horizon,
                            alpha.level = 0.05) {

  if (error.measure == "mse") {
    tmp.error <- fcn_MSE(forecasts = forecasts, outsample = outsample)

  } else if (error.measure == "mad") {
    tmp.error <- fcn_MAD(forecasts = forecasts, outsample = outsample)

  } else if (error.measure == "smape") {
    tmp.error <- fcn_SMAPE(forecasts = forecasts, outsample = outsample)

  } else if (error.measure == "mase") {
    tmp.error <- fcn_MASE(
      forecasts = forecasts,
      outsample        = outsample,
      insample         = insample
    )

  } else if (error.measure == "owa") {
    tmp.error <- fcn_OWA(
      forecasts = forecasts,
      outsample        = outsample,
      insample         = insample,
      forecast.horizon = forecast.horizon,
      alpha.level      = alpha.level
    )
  }

  # |__ RETURN =================================================================
  return(tmp.error)
}

#' @describeIn error_functions MSE - Mean Absolute Deviation
fcn_MSE <- function(forecasts, outsample){

  tmp.forecasts <- as.numeric(forecasts)
  tmp.outsample <- as.numeric(outsample)

  mean((tmp.forecasts - tmp.outsample)**2)
}

#' @describeIn error_functions MAD - Mean Absolute Deviation
fcn_MAD <-  function(forecasts, outsample){

  tmp.forecasts <- as.numeric(forecasts)
  tmp.outsample <- as.numeric(outsample)

  mean(abs(tmp.forecasts - tmp.outsample))
}

#' @describeIn error_functions SMAPE - Symmetric Mean Absolute Percentage Error
fcn_SMAPE <-  function(forecasts, outsample){

  tmp.forecasts <- as.numeric(forecasts)
  tmp.outsample <- as.numeric(outsample)

  tmp.numerator   <- abs(tmp.outsample - tmp.forecasts) * 200
  tmp.denominator <- abs(tmp.outsample) + abs(tmp.forecasts)

  mean(tmp.numerator/tmp.denominator)
}

#' @describeIn error_functions MASE - Mean Absolute Scaled Error
fcn_MASE <- function(forecasts, outsample, insample){

  tmp.forecasts <- as.numeric(forecasts)
  tmp.outsample <- as.numeric(outsample)
  tmp.insample  <- as.numeric(insample)

  ts.freq <- stats::frequency(insample)

  forecastsNaiveSD <- rep(NA, ts.freq)

  for (j in (ts.freq + 1):length(tmp.insample)){
    forecastsNaiveSD <- c(forecastsNaiveSD, tmp.insample[j - ts.freq])
  }

  masep<-mean(abs(tmp.insample - forecastsNaiveSD), na.rm = TRUE)

  mase_t <- (abs(tmp.outsample - tmp.forecasts))/masep

  return(mean(mase_t))
}

#' @describeIn error_functions OWA - Overall Weighted Average
fcn_OWA <- function(forecasts, outsample, insample, forecast.horizon,
                    alpha.level) {

  # \__ Naive2 -----------------------------------------------------------------

  ## \____ Forecasts ----
  tmp.forecasts.naive2 <- forecast_naive2(
    ts.data          = insample,
    forecast.horizon = forecast.horizon,
    alpha.level      = alpha.level
  )

  ## \____ SMAPE ----
  tmp.SMAPE.naive2 <- fcn_SMAPE(
    forecasts = tmp.forecasts.naive2,
    outsample        = outsample
  )

  ## \____ MASE ----
  tmp.MASE.naive2 <- fcn_MASE(
    forecasts = tmp.forecasts.naive2,
    outsample        = outsample,
    insample         = insample
  )

  # \__ User's Method ----------------------------------------------------------

  ## \____ SMAPE ----
  tmp.SMAPE.method <- fcn_SMAPE(
    forecasts = forecasts,
    outsample        = outsample
  )

  ## \____ MASE ----
  tmp.MASE.method <- fcn_MASE(
    forecasts = forecasts,
    outsample        = outsample,
    insample         = insample
  )

  # \__ Final OWA --------------------------------------------------------------

  tmp.ratio.MASE  <- tmp.MASE.method / tmp.MASE.naive2
  tmp.ratio.SMAPE <- tmp.SMAPE.method / tmp.SMAPE.naive2

  final.owa <- (tmp.ratio.MASE + tmp.ratio.SMAPE) / 2

  # |__ RETURN =================================================================
  return(final.owa)
}
