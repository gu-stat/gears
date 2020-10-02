#' Plot method for "gears" class
#'
#' @param x An object of class "gears".
#' @param col Color for the line with observed values. Default is 1 ("black").
#' @param fcol Color for the line with forecasts. Default is 4 ("blue").
#' @param ylim Limits for the Y-axis.
#' @param main Main title for the plot.
#' @param xlab Label for the X-axis.
#' @param ylab Label for the Y-axis.
#' @param type Line type for the line with observed values. See
#'     \link[graphics]{plot} for more options.
#' @param flty Line type for the line with forecasts.
#' @param flwd Width of the line with forecasts.
#' @param ... Other arguments passed to or from other methods.
#'
#' @export plot.gears
#' @export
plot.gears <- function(x, col = 1, fcol = 4, ylim = NULL, main = NULL,
                       xlab = "", ylab = "", type = "l", flty = 1, flwd = 2,
                       ...){

  # |__ Original Data ==========================================================

  tmpNForecasts <- length(x$out_sample_forecasts)
  tmpBindTS     <- cbind(x$x, x$out_sample_forecasts)

  # |__ Data - Observed and Forecasts ==========================================

  tmpStart <- stats::start(tmpBindTS)
  tmpFreq  <- stats::frequency(tmpBindTS)

  if (tmpNForecasts > 1){
    tmpEnd   <- stats::tsp(tmpBindTS)[2] + (round(tmpNForecasts / 2) / tmpFreq)

    dataTS <- stats::ts(
      data      = c(
        x$x, x$out_sample_forecasts, rep(NA, round(tmpNForecasts / 2))
      ),
      start     = tmpStart,
      end       = tmpEnd,
      frequency = tmpFreq
    )

  } else {
    tmpEnd   <- stats::tsp(tmpBindTS)[2] + (tmpNForecasts / tmpFreq)

    dataTS <- stats::ts(
      data      = c(x$x, x$out_sample_forecasts, rep(NA, tmpNForecasts)),
      start     = tmpStart,
      end       = tmpEnd,
      frequency = tmpFreq
    )
  }

  # |__ Data - Gaps ============================================================

  tmpStartGAP <- stats::tsp(x$x)[2]
  tmpEndGAP   <- stats::tsp(x$out_sample_forecasts)[1]

  dataGAP <- stats::ts(
    data      = c(utils::tail(x$x, 1), utils::head(x$out_sample_forecasts, 1)),
    start     = tmpStartGAP,
    end       = tmpEndGAP,
    frequency = tmpFreq
  )

  # |__ Dates - Intervals ======================================================

  tmpDatesPolygon <- stats::time(x$out_sample_forecasts)

  # |__ Y-Axis Limits ==========================================================

  if (is.null(ylim)) {
    ylim <- c(
      min(x$x, x$out_sample_forecasts, na.rm = TRUE),
      max(x$x, x$out_sample_forecasts, na.rm = TRUE)
    )

    ylim[1] <- min(ylim[1], x$lower)
    ylim[2] <- max(ylim[2], x$upper)
  }

  # |__ Title ==================================================================

  if (is.null(main)) {
    main <- paste0(
      "Forecasts from GEARS with ", x$level, "% prediction interval"
    )
  }

  # |__ Plot ===================================================================

  graphics::plot(
    x    = dataTS,
    ylim = ylim,
    main = main,
    col  = col,
    type = type,
    xlab = xlab,
    ylab = ylab
  )
  graphics::lines(dataGAP, col="white")

  if (tmpNForecasts > 1){
    graphics::polygon(
      x      = c(tmpDatesPolygon, rev(tmpDatesPolygon)),
      y      = c(x$lower, rev(x$upper)),
      col    = 'grey80',
      border = NA
    )
    graphics::lines(x$out_sample_forecasts, col = fcol, lty = flty, lwd = flwd)
  } else {
    tmpLowerDate <- c(tmpDatesPolygon) - (0.5 / tmpFreq)
    tmpUpperDate <- c(tmpDatesPolygon) + (0.5 / tmpFreq)
    graphics::polygon(
      x      = c(tmpLowerDate, tmpUpperDate, tmpUpperDate, tmpLowerDate),
      y      = c(rep(x$lower, 2), rep(rev(x$upper), 2)),
      col    = 'grey80',
      border = NA
    )
    graphics::points(x$out_sample_forecasts, col = fcol, lty = flty, lwd = flwd)
  }

}
