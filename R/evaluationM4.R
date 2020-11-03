#' Function to evaluate the forecasts following the M4 Competition guidelines
#'
#' @param DATA A ts object (i.e., a vector of time series observations).
#' @param forecast.list A list with ts objects.
#' @param alpha Significance level
#'
#' @return A list with three values: "smape", "mase", and "owa".
#' @export
#'
evaluationM4 <- function(DATA, forecast.list,  alpha){

  # \____ Forecasts Naive ------------------------------------------------------

  tmpForecastsNaive2 <- lapply(
    X = 1:length(DATA),
    function(X) {
       gears::forecast_naive2(
        ts.data          = DATA[[X]]$x,
        ts.frequency     = stats::frequency(DATA[[X]]$x),
        forecast.horizon = DATA[[X]]$h,
        alpha.level      = alpha
      )
    }
  )
  # \____ sMAPE ----------------------------------------------------------------

  ## sMAPE for all hourly series

  smapeM4All <- lapply(
    X = 1:length(DATA),
    function(X) {

      tmp <- apply(
        forecast.list[[X]],
        MARGIN            = 1,
        FUN               = error_measures,
        outsample         = DATA[[X]]$xx,
        insample          = DATA[[X]]$x,
        ts.frequency      = stats::frequency(DATA[[X]]$x),
        forecast.horizon  = DATA[[X]]$h,
        alpha.level       = alpha,
        error.measure     = "smape"
      )

      unlist(tmp)
    }
  )

  smapeM4All           <- do.call(rbind, smapeM4All)
  #colnames(smapeM4All) <- methods.names

  smapeM4All           <- colMeans(smapeM4All)

  ## SMAPE for naive2

  smapeM4Naive2 <- lapply(
    X = 1:length(DATA),
    function(X) {

      tmp <- error_measures(
        forecasts         = tmpForecastsNaive2[[X]],
        outsample         = DATA[[X]]$xx,
        insample          = DATA[[X]]$x,
        ts.frequency      = stats::frequency(DATA[[X]]$x),
        forecast.horizon  = DATA[[X]]$h,
        alpha.level       = alpha,
        error.measure     = "smape"
      )

      unlist(tmp)
    }
  )

  smapeM4Naive2 <- mean(unlist(smapeM4Naive2))

  # \____ MASE -------------------------------------------------------------------

  ## MASE for all hourly series

  maseM4All <- lapply(
    X = 1:length(DATA),
    function(X) {

      tmp <- apply(
        forecast.list[[X]],
        MARGIN            = 1,
        FUN               = error_measures,
        outsample         = DATA[[X]]$xx,
        insample          = DATA[[X]]$x,
        ts.frequency      = stats::frequency(DATA[[X]]$x),
        forecast.horizon  = DATA[[X]]$h,
        alpha.level       = alpha,
        error.measure     = "mase"
      )

      unlist(tmp)
    }
  )

  maseM4All           <- do.call(rbind, maseM4All)
 #colnames(maseM4All) <- M4comp2018::submission_info$ID[1:25]

  maseM4All           <- colMeans(maseM4All)

  ## MASE for naive2

  maseM4Naive2 <- lapply(
    X = 1:length(DATA),
    function(X) {

      tmp <- error_measures(
        forecasts         = tmpForecastsNaive2[[X]],
        outsample         = DATA[[X]]$xx,
        insample          = DATA[[X]]$x,
        ts.frequency      = stats::frequency(DATA[[X]]$x),
        forecast.horizon  = DATA[[X]]$h,
        alpha.level       = alpha,
        error.measure     = "mase"
      )

      unlist(tmp)
    }
  )

  maseM4Naive2 <- mean(unlist(maseM4Naive2))

  # \____ OWA --------------------------------------------------------------------

  owaM4All <- rowMeans(
    cbind(smapeM4All / smapeM4Naive2, maseM4All / maseM4Naive2)
  )

  # \____ Return -----------------------------------------------------------------

  cbind("smape" = smapeM4All, "mase" = maseM4All, "owa" = owaM4All)
}
