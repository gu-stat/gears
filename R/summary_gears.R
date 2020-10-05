#' Summary method for "gears" class
#'
#' @param object An object of class "gears"
#' @param ... Other arguments passed to or from other methods
#'
#' @export summary.gears
#' @export
summary.gears <- function(object, ...){

  cat("Summary of ForecastS with GEARS\n")

  cat(paste("\nNumber of Estimated Models:", object$total_equations_estimated))

  tmpForecastHorizon <- length(object$details)

  cat(paste("\n\nBest Equations:"))
  cat(paste("\nForecast Lead", 1:tmpForecastHorizon,":",
    lapply(
      1:tmpForecastHorizon,
      function(X) object$details[[X]]$best_equation
    )
  ))

  cat("\n\nAverage In-Sample Prediction Errors:\n")
  print(object$min_prediction_errors)

  if (is.null(object$out_sample_forecasts)) {
    cat("\nNo forecasts\n")
  } else {
    cat("\nForecasts:\n")
    tmpBind <- cbind(object$out_sample_forecasts, object$lower, object$upper)

    colnames(tmpBind) <- c(
      "Point Forecasts",
      paste0("Lo ", object$level),
      paste0("Hi ", object$level)
    )

    print(tmpBind)
  }
}
