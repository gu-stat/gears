# >  ###################################################### ----

# Used to determine whether a time series is seasonal.
# Code modified from the M4 Competition's original code.

#' Seasonality Test
#'
#' Used to determine whether a time series is seasonal. Code modified from the
#' M4 Competition's original code.
#'
#' @param ts.data A ts object (i.e., a vector of time series observations).
#' @param ts.frequency The frequency of a ts object.
#' @param alpha.level The alpha level to be used in the test to detect
#'     seasonality. Default is 0.05.
#'
#' @return A boolean indicating whether or not the ts object has seasonality.
#'
#' @keywords internal
# @examples
# seasonality_test(
#   ts.data = datasets::WWWusage,
#   ts.frequency = frequency(datasets::WWWusage),
#   alpha.level = 0.05
# )
#
# seasonality_test(
#  ts.data = datasets::AirPassengers,
#   ts.frequency = frequency(datasets::AirPassengers),
#   alpha.level = 0.05
# )
#
# seasonality_test(
#   ts.data = datasets::EuStockMarkets[, "DAX"],
#   ts.frequency = frequency(datasets::EuStockMarkets[, "DAX"]),
#   alpha.level = 0.05
# )
seasonality_test <- function(ts.data, ts.frequency, alpha.level = 0.05){

  tmp.tcrit <- stats::qnorm(1 - alpha.level)

  if (length(ts.data) < 3 * ts.frequency){

    test_seasonal <- FALSE

  } else {

    xacf <- stats::acf(ts.data, plot = FALSE)$acf[-1, 1, 1]

    clim <- tmp.tcrit/sqrt(length(ts.data)) * sqrt(cumsum(c(1, 2 * xacf^2)))

    test_seasonal <- (abs(xacf[ts.frequency]) > clim[ts.frequency])

    if (is.na(test_seasonal) ==  TRUE){
      test_seasonal <- FALSE
    }

  }

  # |__ RETURN =================================================================
  return(test_seasonal)
}
