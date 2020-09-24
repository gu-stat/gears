#' Print method for "gears" class
#'
#' @param x An object of class "gears"
#' @param ... Other arguments passed to or from other methods
#'
#' @export print.gears
#' @export
print.gears <- function(x, ...){

  tmpBind <- cbind(x$out_sample_forecasts, x$lower, x$upper)

  colnames(tmpBind) <- c(
    "Point Forecasts", paste0("Lo ", x$level), paste0("Hi ", x$level)
  )

  print(tmpBind)
}
