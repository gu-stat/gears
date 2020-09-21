#' Functions to calculate the forecast errors.
#'
#' \code{error_measures} returns the forecast errors for the following error
#' measures: MSE, MAD, sMAPE, MASE, and OWA.
#'
# @param error.measure A character string with the name of the error measure.
#'     These can be 'mse', 'mad', 'smape', 'mase', and 'owa' (see Details).
#' @param forecasts A numeric vector with the point forecasts.
#' @param outsample A numeric vector with the test data set.
#' @param insample A numeric vector with the train data set.
#' @param ts.frequency The frequency of a ts object.
#' @param forecast.horizon A numeric value with the length of the forecast lead.
#' @param alpha.level A numeric value with the alpha level to be used in the
#'     test to detect seasonality. Default is 0.05.
#'
#' @return A 5-object list with the numeric values of the forecast error for
#'     each one of the error measures.
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
#' tmp.orig.start <- stats::tsp(datasets::AirPassengers)[1]
#' tmp.orig.end   <- stats::tsp(datasets::AirPassengers)[2]
#' tmp.orig.freq  <- stats::tsp(datasets::AirPassengers)[3]
#'
#'
#' # Get train data (insample)
#'
#' tmp.train.start <- tmp.orig.start
#' tmp.train.end   <- tmp.orig.start + ((tmp.cut.at - 1) * 1 / tmp.orig.freq)
#'
#' tmp.train.data  <- stats::window(
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
#' tmp.test.data  <- stats::window(
#'   x         = datasets::AirPassengers,
#'   start     = tmp.test.start,
#'   end       = tmp.test.end,
#'   frequency = tmp.orig.freq
#' )
#'
#' # Get Forecasts
#' tmp.forecasts <- gears::forecast_naive2(
#'   ts.data          = tmp.train.data,
#'   ts.frequency     = tmp.orig.freq,
#'   forecast.horizon = tmp.forecast.horizon,
#'   alpha.level      = 0.05
#' )
#'
#' # Get error measures for the forecast
#' error_measures(
#'   forecasts          = tmp.forecasts,
#'   outsample          = tmp.test.data,
#'   insample           = tmp.train.data,
#'   ts.frequency       = tmp.orig.freq,
#'   forecast.horizon   = tmp.forecast.horizon,
#'   alpha.level        = 0.05
#' )
error_measures <- function(forecasts,
                           outsample,
                           insample = NULL,
                           ts.frequency,
                           forecast.horizon,
                           alpha.level = 0.05) {

  tmp.MSE   <- measure_MSE(forecasts = forecasts, outsample = outsample)
  tmp.MAD   <- measure_MAD(forecasts = forecasts, outsample = outsample)
  tmp.SMAPE <- measure_SMAPE(forecasts = forecasts, outsample = outsample)

  tmp.MASE <- measure_MASE(
    forecasts = forecasts,
    outsample = outsample,
    insample  = insample
  )

  tmp.OWA <- measure_OWA(
    forecasts        = forecasts,
    outsample        = outsample,
    insample         = insample,
    ts.frequency     = ts.frequency,
    forecast.horizon = forecast.horizon,
    alpha.level      = alpha.level
  )

 tmp.error <- list(
   "mse"   = tmp.MSE,
   "mad"   = tmp.MAD,
   "smape" = tmp.SMAPE,
   "mase"  = tmp.MASE,
   "owa"   = tmp.OWA
 )


  # |__ RETURN =================================================================
  return(tmp.error)
}

#' @describeIn error_measures MSE - Mean Absolute Deviation
measure_MSE <- function(forecasts, outsample){

  tmp.forecasts <- as.numeric(forecasts)
  tmp.outsample <- as.numeric(outsample)

  mean((tmp.forecasts - tmp.outsample)**2)
}

#' @describeIn error_measures MAD - Mean Absolute Deviation
measure_MAD <-  function(forecasts, outsample){

  tmp.forecasts <- as.numeric(forecasts)
  tmp.outsample <- as.numeric(outsample)

  mean(abs(tmp.forecasts - tmp.outsample))
}

#' @describeIn error_measures SMAPE - Symmetric Mean Absolute Percentage Error
measure_SMAPE <-  function(forecasts, outsample){

  tmp.forecasts <- as.numeric(forecasts)
  tmp.outsample <- as.numeric(outsample)

  tmp.numerator   <- abs(tmp.outsample - tmp.forecasts) * 200
  tmp.denominator <- abs(tmp.outsample) + abs(tmp.forecasts)

  mean(tmp.numerator/tmp.denominator)
}

#' @describeIn error_measures MASE - Mean Absolute Scaled Error
measure_MASE <- function(forecasts, outsample, insample){

  tmp.forecasts <- as.numeric(forecasts)
  tmp.outsample <- as.numeric(outsample)
  tmp.insample  <- as.numeric(insample)

  ts.frequency <- stats::frequency(insample)

  forecastsNaiveSD <- rep(NA, ts.frequency)

  for (j in (ts.frequency + 1):length(tmp.insample)){
    forecastsNaiveSD <- c(forecastsNaiveSD, tmp.insample[j - ts.frequency])
  }

  masep <- mean(abs(tmp.insample - forecastsNaiveSD), na.rm = TRUE)

  mase_t <- (abs(tmp.outsample - tmp.forecasts))/masep

  return(mean(mase_t))
}

#' @describeIn error_measures OWA - Overall Weighted Average
measure_OWA <- function(forecasts, outsample, insample, forecast.horizon,
                        alpha.level, ts.frequency) {

  # \__ Naive2 -----------------------------------------------------------------

  ## \____ Get Forecasts ----
  tmp.forecasts.naive2 <- forecast_naive2(
    ts.data          = insample,
    ts.frequency     = ts.frequency,
    forecast.horizon = forecast.horizon,
    alpha.level      = alpha.level
  )

  ## \____ SMAPE ----
  tmp.SMAPE.naive2 <- measure_SMAPE(
    forecasts = tmp.forecasts.naive2,
    outsample = outsample
  )

  ## \____ MASE ----
  tmp.MASE.naive2 <- measure_MASE(
    forecasts = tmp.forecasts.naive2,
    outsample = outsample,
    insample  = insample
  )

  # \__ User's Method ----------------------------------------------------------

  ## \____ SMAPE ----
  tmp.SMAPE.method <- measure_SMAPE(
    forecasts = forecasts,
    outsample = outsample
  )

  ## \____ MASE ----
  tmp.MASE.method <- measure_MASE(
    forecasts  = forecasts,
    outsample  = outsample,
    insample   = insample
  )

  # \__ Final OWA --------------------------------------------------------------

  tmp.ratio.MASE  <- tmp.MASE.method / tmp.MASE.naive2
  tmp.ratio.SMAPE <- tmp.SMAPE.method / tmp.SMAPE.naive2

  final.owa <- (tmp.ratio.MASE + tmp.ratio.SMAPE) / 2

  # |__ RETURN =================================================================
  return(final.owa)
}
