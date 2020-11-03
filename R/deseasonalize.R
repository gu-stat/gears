#' Function to deseasonalize a time series.
#'
#' This function deseasonalizes a time series using the "multiplicative"
#' decomposition only.
#'
#' @param ts.data A ts object (i.e., a vector of time series observations).
#' @param ts.frequency The frequency of a ts object.
#' @param forecast.horizon Number of periods for forecasting.
#' @param alpha.level The alpha level to be used in the test to detect
#'     seasonality. Default is 0.05.
#'
#' @return A list with the deseasonalized data (\code{"deseasonTS"}) and the
#'     seasonal component (\code{"seasonalComp"}).
#'
#' @examples
#' # Example using a TS with frequency > 1 (in this case, frequency = 12)
#' tmpAir <- deseason(
#'   ts.data          = datasets::AirPassengers,
#'   ts.frequency     = stats::frequency(datasets::AirPassengers),
#'   forecast.horizon = 20
#' )
#'
#' head(datasets::AirPassengers) # original data
#' head(tmpAir$deseasonTS)       # deseasonalized series w/o multiplicative
#'                               # component, not equal to the original
#' head(tmpAir$deseasonTS * tmpAir$seasonalComp) # this should be equal to the
#'                                               # original data
#' # Example using a TS with frequency = 1
#' tmpUsage <- deseason(
#'   ts.data          = datasets::WWWusage,
#'   ts.frequency     = stats::frequency(datasets::WWWusage),
#'   forecast.horizon = 10
#' )
#'
#' head(datasets::WWWusage)  # original data
#' head(tmpUsage$deseasonTS) # since frequency = 1, results should be the same
#'
#' @export
#'
deseason <- function(ts.data, ts.frequency, forecast.horizon,
                     alpha.level = 0.05){
  # Deseasonalize the data
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

  # Return
  list("deseasonTS" = deseason.ts.data, "seasonalComp" = SIout)
}
