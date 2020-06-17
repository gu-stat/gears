# >  ###################################################### ----

# Used to determine whether a time series is seasonal.
# Code modified from the M4 Competition's original code.

#' Seasonality Test
#'
#' Used to determine whether a time series is seasonal. Code modified from the
#' M4 Competition's original code.
#'
#' @param ts.data A ts object (i.e., a vector of time series observations).
#' @param ts.freq The frequency of a ts object.
#' @param alpha.level The alpha level to be used in the test to detect
#'     seasonality. Default is 0.05.
#'
#' @return A boolean indicating whether or not the ts object has seasonality.
#'
#' @keywords internal
# @examples
# seasonality_test(
#   ts.data = datasets::WWWusage,
#   ts.freq = frequency(datasets::WWWusage),
#   alpha.level = 0.05
# )
#
# seasonality_test(
#   ts.data = datasets::AirPassengers,
#   ts.freq = frequency(datasets::AirPassengers),
#   alpha.level = 0.05
# )
#
# seasonality_test(
#   ts.data = datasets::EuStockMarkets[, "DAX"],
#   ts.freq = frequency(datasets::EuStockMarkets[, "DAX"]),
#   alpha.level = 0.05
# )
seasonality_test <- function(ts.data, ts.freq, alpha.level = 0.05){

  tmp.tcrit <- stats::qnorm(1 - alpha.level)

  if (length(ts.data) < 3*ts.freq){

    test_seasonal <- FALSE

  } else {

    xacf <- stats::acf(ts.data, plot = FALSE)$acf[-1, 1, 1]

    clim <- tmp.tcrit/sqrt(length(ts.data)) * sqrt(cumsum(c(1, 2 * xacf^2)))

    test_seasonal <- (abs(xacf[ts.freq]) > clim[ts.freq])

    if (is.na(test_seasonal) ==  TRUE){
      test_seasonal <- FALSE
    }

  }

  # |__ RETURN =================================================================
  return(test_seasonal)
}
