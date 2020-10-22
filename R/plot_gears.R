#' Plot method for "gears" class
#'
#' @param x An object of class "gears".
#' @param col Color for the line with observed values. Default is 1 ("black").
#' @param fcol Color for the line with forecasts. Default is 4 ("blue").
#' @param ylim Limits for the Y-axis.
#' @param main Main title for the plot.
#' @param sub Sub-title for the plot.
#' @param xlab Label for the X-axis.
#' @param ylab Label for the Y-axis.
#' @param type Line type for the line with observed values. See
#'     \link[graphics]{plot} for more options.
#' @param flty Line type for the line with forecasts.
#' @param flwd Width of the line with forecasts.
#' @param cex  Size of the main title when \code{betas.selection = "both"}.
#' @param ... Other arguments passed to or from other methods.
#'
#' @export plot.gears
#' @export
plot.gears <- function(x, col = 1, fcol = 4, ylim = NULL, main = NULL,
                       sub = NULL, xlab = "", ylab = "", type = "l", flty = 1,
                       flwd = 2, cex = 1.75, ...){

  # > PLOT FUNCTION ####################################################### ----

  fnct_plot <- function(observed, forecasts, lower, upper, level, betas,
                        main, sub, ...){

    # |__ Original Data ========================================================

    tmpForeLen <- length(forecasts)
    tmpBindTS  <- cbind(observed, forecasts)

    # |__ Data - Observed and Forecasts ========================================

    tmpStart <- stats::start(tmpBindTS)
    tmpFreq  <- stats::frequency(tmpBindTS)

    if (tmpForeLen > 1){
      tmpEnd   <- stats::tsp(tmpBindTS)[2] + (round(tmpForeLen / 2) / tmpFreq)

      dataTS <- stats::ts(
        data      = c(observed, forecasts, rep(NA, round(tmpForeLen / 2))),
        start     = tmpStart,
        end       = tmpEnd,
        frequency = tmpFreq
      )

    } else {
      tmpEnd   <- stats::tsp(tmpBindTS)[2] + (tmpForeLen / tmpFreq)

      dataTS <- stats::ts(
        data      = c(observed, forecasts, rep(NA, tmpForeLen)),
        start     = tmpStart,
        end       = tmpEnd,
        frequency = tmpFreq
      )
    }

    # |__ Data - Gaps ==========================================================

    tmpStartGAP <- stats::tsp(observed)[2]
    tmpEndGAP   <- stats::tsp(forecasts)[1]

    dataGAP <- stats::ts(
      data      = c(utils::tail(observed, 1), utils::head(forecasts, 1)),
      start     = tmpStartGAP,
      end       = tmpEndGAP,
      frequency = tmpFreq
    )

    # |__ Dates - Intervals ====================================================

    tmpDatesPolygon <- stats::time(forecasts)

    # |__ Y-Axis Limits ========================================================

    if (is.null(ylim)) {
      ylim <- c(
        min(observed, forecasts, na.rm = TRUE),
        max(observed, forecasts, na.rm = TRUE)
      )
      ylim[1] <- min(ylim[1], lower)
      ylim[2] <- max(ylim[2], upper)
    }

    # |__ Title ================================================================

    if (is.null(main)) {
      main <- paste0(
        "Forecasts from GEARS with ", level, "% prediction interval"
      )
    }

    # |__ Subtitle =============================================================

    if (is.null(sub)) {
      if (betas == "last") {
        sub <- paste0(
          "Obtained using the coefficients for the last random sample."
        )
      } else {
        sub <- paste0(
          "Obtained using the average of coefficients for all random samples."
        )
      }
    }

    # |__ Plot =================================================================

    graphics::plot(
      x    = dataTS,
      ylim = ylim,
      main = main,
      sub  = sub,
      col  = col,
      type = type,
      xlab = xlab,
      ylab = ylab
    )
    graphics::lines(dataGAP, col="white")

    if (tmpForeLen > 1){

      tmpPaddingX <- rep(utils::tail(tmpDatesPolygon, 1) + (0.1 / tmpFreq), 2)
      tmpPaddingY <- c(utils::tail(lower, 1), utils::head(rev(upper), 1))

      graphics::polygon(
        x      = c(tmpDatesPolygon, tmpPaddingX, rev(tmpDatesPolygon)),
        y      = c(lower, tmpPaddingY, rev(upper)),
        col    = 'grey80',
        border = NA
      )
      graphics::lines(forecasts, col = fcol, lty = flty, lwd = flwd)
    } else {
      lowerDate <- c(tmpDatesPolygon) - (0.5 / tmpFreq)
      upperDate <- c(tmpDatesPolygon) + (0.5 / tmpFreq)
      graphics::polygon(
        x      = c(lowerDate, upperDate, upperDate, lowerDate),
        y      = c(rep(lower, 2), rep(rev(upper), 2)),
        col    = 'grey80',
        border = NA
      )
      graphics::points(forecasts, col = fcol, lty = flty, lwd = flwd)
    }

  }

  # > CHECK BETAS ######################################################### ----

  # |__ Each Selection =========================================================

  if (x$betas != "both") {
    fnct_plot(
      observed  = x$x,
      forecasts = x$out_sample_forecasts,
      lower     = x$lower,
      upper     = x$upper,
      level     = x$level,
      betas     = x$betas,
      main      = main,
      sub       = sub,
      ...
    )
  }
  # |__ Both ===================================================================
  else {

    # \____ Current par settings -----------------------------------------------

    ## Save current par settings and return after finished

    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op))

    # \____ Set new layout -----------------------------------------------------

    graphics::layout(
      matrix(1:length(x$out_sample_forecasts), ncol = 1, byrow = TRUE)
    )

    # \____ Set new par --------------------------------------------------------

    ## Margins

    graphics::par(mar = c(3, 3, 2, 3), oma = c(3, 3, 3, 3))

    # \____ Title --------------------------------------------------------------

    ## Main ("Outer")

    if (is.null(main)) {
      tmpMain <-
        paste0("Forecasts from GEARS with ", x$level, "% prediction interval")
    } else {
      tmpMain <- main
    }

    ## "Sub-titles"

    tmpMainEach <- list(
      "Obtained using the coefficients for the last random sample.",
      "Obtained using the average of coefficients for all random samples."
    )

    # \____ Plots --------------------------------------------------------------

    tmpBetas <- list("last", "average")

    tmpPlot <- lapply(
      X = 1:length(x$out_sample_forecasts),
      function(X) {

        fnct_plot(
          observed  = x$x,
          forecasts = x$out_sample_forecasts[[X]],
          lower     = x$lower[[X]],
          upper     = x$upper[[X]],
          level     = x$level,
          betas     = tmpBetas[[X]],
          main      = tmpMainEach[[X]],
          sub       = "",
          ...
        )

        if (X == 1) {
          graphics::mtext(tmpMain, outer = T, line = 1, side = 3, cex = cex)
        }

      }
    )

  }
}
