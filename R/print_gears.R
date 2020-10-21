#' Print method for "gears" class
#'
#' @param x An object of class "gears"
#' @param ... Other arguments passed to or from other methods
#'
#' @export print.gears
#' @export
print.gears <- function(x, ...){

  if (x$betas == "both"){
    tmpBind <- lapply(
      X = 1:2,
      function(X) {
         tmp <- cbind(x$out_sample_forecasts[[X]], x$lower[[X]], x$upper[[X]])
         colnames(tmp) <- c(
           "Point Forecasts", paste0("Lo ", x$level), paste0("Hi ", x$level)
         )
         tmp
      }
    )

    names(tmpBind) <- c("beta.selection = LAST", "beta.selection = AVERAGE")

  } else {

    tmpBind <- cbind(x$out_sample_forecasts, x$lower, x$upper)

    colnames(tmpBind) <- c(
      "Point Forecasts", paste0("Lo ", x$level), paste0("Hi ", x$level)
    )
  }

  print(tmpBind)
}
