#' Summary method for "gears" class
#'
#' @param object An object of class "gears"
#' @param ... Other arguments passed to or from other methods
#'
#' @export summary.gears
#' @export
summary.gears <- function(object, ...){

  if (object$betas == "both"){
    tmpBind <- lapply(
      X = 1:2,
      function(X) {
        tmp <- cbind(object$out_sample_forecasts[[X]], object$lower[[X]], object$upper[[X]])
        colnames(tmp) <- c(
          "Point Forecasts", paste0("Lo ", object$level), paste0("Hi ", object$level)
        )
        tmp
      }
    )

    names(tmpBind) <- c("beta.selection = LAST", "beta.selection = AVERAGE")

  } else {

    tmpBind <- cbind(object$out_sample_forecasts, object$lower, object$upper)

    colnames(tmpBind) <- c(
      "Point Forecasts",
      paste0("Lo ", object$level),
      paste0("Hi ", object$level)
    )
  }

  # -------------------------------------------------------------------------- #

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
    print(tmpBind)
  }
}
