# > HELPERS ############################################################### ----

# |__ MSE - Mean Absolute Deviation ============================================

fcn_MSE <- function(method.forecasts, outsample){
  mean((method.forecasts - outsample)**2)
}

# |__ MAD - Mean Absolute Deviation ============================================

fcn_MAD <-  function(method.forecasts, outsample){
  mean(abs(method.forecasts-outsample))
}

# |__ SMAPE - Symmetric Mean Absolute Percentage Error =========================

fcn_SMAPE <-  function(method.forecasts, outsample){

  tmp.numerator   <- abs(outsample - method.forecasts) * 200
  tmp.denominator <- abs(outsample) + abs(method.forecasts)

  mean(tmp.numerator/tmp.denominator)
}

# |__ MASE - Mean Absolute Scaled Error ========================================

fcn_MASE <- function(method.forecasts, outsample, insample){

  ts.freq <- frequency(insample)

  forecastsNaiveSD <- rep(NA, ts.freq)

  for (j in (ts.freq + 1):length(insample)){
    forecastsNaiveSD <- c(forecastsNaiveSD, insample[j - ts.freq])
  }

  masep<-mean(abs(insample - forecastsNaiveSD), na.rm = TRUE)

  mase <- (abs(outsample - method.forecasts))/masep

  return(mean(mase))
}

# |__ OWA - Overall Weighted Average ===========================================

fcn_OWA <- function(method.forecasts, outsample, insample, forecast.horizon,
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
    method.forecasts = tmp.forecasts.naive2,
    outsample        = outsample
  )

  ## \____ MASE ----
  tmp.MASE.naive2 <- fcn_MASE(
    method.forecasts = tmp.forecasts.naive2,
    outsample        = outsample,
    insample         = insample
  )

  # \__ User's Method ----------------------------------------------------------

  ## \____ SMAPE ----
  tmp.SMAPE.method <- fcn_SMAPE(
    method.forecasts = method.forecasts,
    outsample        = outsample
  )

  ## \____ MASE ----
  tmp.MASE.method <- fcn_MASE(
    method.forecasts = method.forecasts,
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

# > Main Function ######################################################### ----

error_functions <- function(error.measure,
                            forecasts,
                            outsample,
                            insample = NULL,
                            forecast.horizon,
                            alpha.level) {

  if (error.measure == "mse") {
    tmp.error <- fcn_MSE(method.forecasts = forecasts, outsample = outsample)

  } else if (error.measure == "mad") {
    tmp.error <- fcn_MAD(method.forecasts = forecasts, outsample = outsample)

  } else if (error.measure == "smape") {
    tmp.error <- fcn_SMAPE(method.forecasts = forecasts, outsample = outsample)

  } else if (error.measure == "mase") {
    tmp.error <- fcn_MASE(
      method.forecasts = forecasts,
      outsample        = outsample,
      insample         = insample
    )

  } else if (error.measure == "owa") {
    tmp.error <- fcn_OWA(
      method.forecasts = forecasts,
      outsample        = outsample,
      insample         = insample,
      forecast.horizon = forecast.horizon,
      alpha.level      = alpha.level
    )
  }

  # |__ RETURN =================================================================
  return(tmp.error)
}

# TODO; TEST THIS FUNCTION
