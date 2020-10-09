#' Optimization Method for the GEARS Method for Time Series Forecasting
#'
#' This algorithm will search for the "best" size of the rolling sample and the
#' "best" number of rolling samples. "Best" in the sense of minimizing the
#' selected error measure.
#'
#' @param DATA A data frame or a univariate time series.
#' @param forecast.horizon Number of periods for forecasting.
#' @param search.size.rs Vector with the number of observations in the rolling
#'     sample to be optimized.
#' @param search.number.rs Vector with the number of rolling samples to be
#'     optimized.
#' @param level Confidence level for prediction intervals. Numeric value
#'     between 0 and 100.
#' @param details If \code{TRUE}, the function outputs the entire glm object for
#'     all random samples for each forecast lead. To save memory, the default
#'     is \code{FALSE}.
#' @param glm.family A description of the error distribution to be used in the
#'     model. See \link[stats]{glm} for details.
#' @param y.name The name of the Y (left-hand side) variable. If NULL (default),
#'     the function creates a temporary name.
#' @param y.max.lags A numeric value that gives the maximum number of lags of
#'     the Y (left-hand side) variable (see Details). Can be NULL (default) if
#'     the past values of the Y variable are not included in the model.
#' @param x.names List with names of the X (right-hand side) variables that
#'     have a maximum number of lags (see Details). Can be NULL (default) if
#'     univariate model or if your model does not have variables of this type.
#' @param x.max.lags List of numeric values that give the maximum number of lags
#'     of the X (right-hand side) variables. Can be NULL (default).
#' @param x.fixed.names List with names of the X (right-hand side) variables
#'     that have a fixed number of lags (see Details).
#'     Can be NULL (default) if univariate model or if your model does not have
#'     variables of this type.
#' @param x.fixed.lags List of numeric values that give the fixed number of lags
#'     of the variables in x.fixed.names. Can be NULL (default).
#' @param x.interaction.names List of character vectors with names of the
#'     variables to be included as interaction terms (see Details).
#'     Can be NULL (default) if your model does not have interaction terms.
#' @param x.interaction.lags List of numeric vectors with lags of the
#'     variables to be included as interaction terms. List and numeric vectors
#'     should have the same length as the ones in x.interaction.names.
#'     Can be NULL (default) if your model does not have interaction terms.
#' @param last.obs Index number of the last observation to be considered.
#' @param use.intercept If \code{use.intercept == "both"} (default), the
#'     function returns all possible model equations with and without intercept.
#'     If \code{use.intercept == "without"}, then the function returns all
#'     possible right-hand side equations without intercept. If
#'     \code{use.intercept == "with"}, then only the right-hand side equations
#'     with intercept are returned.
#' @param error.measure Error measure to be used when calculating the in-sample
#'     prediction errors.
#' @param betas.selection If \code{betas.selection == "last"}, the
#'     estimated coefficients from the last rolling sample will be used to
#'     obtain the out-of-sample forecasts. If
#'     \code{betas.selection == "average"}, then the average of the coefficients
#'     from all rolling samples will be used.
#'     If \code{betas.selection == "both"}, then two out-of-sample forecasts
#'     will be estimated (one with \code{betas.selection == "last"} and
#'     another with \code{betas.selection == "average"}).
#' @param ... Further arguments passed to \link[stats]{glm}.
#'
#' @details If \code{y.max.lags} equals to a number, then all lags of
#'     \code{y.name} up to this number (and starting at 0) will be included
#'     in the list of variables to enter the right-hand side of the model
#'     equations. For example, if \code{y.max.lags = 2}, then
#'     \ifelse{html}{\out{Y<sub>t</sub>, Y<sub>t-1</sub>, Y<sub>t-2</sub>
#'     }}{\eqn{Y_t, Y_{t-1},Y_{t-2}}} will be included in the list of variables
#'     to enter the right-hand side of the equation.
#'
#' @return An object of class "\code{gears}".
#'     An object of class "\code{gears}" is a list containing the following
#'     elements:
#'     \describe{
#'      \item{fitted_best}{A list containing information about
#'          the best fitted model for each error measure (outer level of the
#'          list), forecast horizon (inner level 1), and each rolling sample
#'          (inner level 2).}
#'      \item{out_sample_point_forecasts}{A data frame containing point
#'          forecasts for each error measure and beta-selection method.}
#'      \item{insample_point_estimates}{List containing the out-of-sample point
#'          estimates for each rolling sample (row) and each forecast horizon
#'          (column).}
#'      \item{insample_mean_error_measure}{A data frame containing the mean
#'          prediction errors for each error measure and forecast horizon.}
#'      \item{best_equations}{A data frame containing information about the best
#'          equation fitted for each error measure and forecast horizon.}
#'     }
#' @author Gustavo Varela-Alvarenga
#'
#' @keywords ts glm optimization
#'
#'#' @examples
#' # Univariate Time Series Forecasting - Data of class 'ts'.
#' gears(
#'   DATA                = datasets::WWWusage,
#'   forecast.horizon    = 12,
#'   search.size.rs      = c(18:20),
#'   search.number.rs    = c(10:12),
#'   level               = 95,
#'   details             = FALSE,
#'   glm.family          = "quasi",
#'   y.name              = NULL,
#'   y.max.lags          = 2,
#'   x.names             = NULL,
#'   x.max.lags          = NULL,
#'   x.fixed.names       = NULL,
#'   x.fixed.lags        = NULL,
#'   x.interaction.names = NULL,
#'   x.interaction.lags  = NULL,
#'   last.obs            = length(datasets::WWWusage),
#'   use.intercept       = "both",
#'   error.measure       = "mse",
#'   betas.selection     = "last"
#' )
gears_optim <- function(DATA,
                  forecast.horizon,
                  search.size.rs,
                  search.number.rs,
                  glm.family = c(
                    "gaussian", "binomial", "poisson", "Gamma", "quasi"
                  ),
                  level = 95,
                  details = FALSE,
                  y.name = NULL,
                  y.max.lags = NULL,
                  x.names = NULL,
                  x.max.lags = NULL,
                  x.fixed.names = NULL,
                  x.fixed.lags = NULL,
                  x.interaction.names = NULL,
                  x.interaction.lags = NULL,
                  last.obs = NULL,
                  use.intercept = c("both", "without", "with"),
                  error.measure = c("mse", "mae", "mase", "smape", "owa"),
                  betas.selection = c("last", "average", "both"),
                  ...) {

  # > CHECKS ############################################################## ----
  #
  # TODO: INCLUDE CHECKS FROM OLD FILE "new_gears.R"
  # TODO: FIX PROBLEMS WHEN NOT-SPECIFYING last.obs, use.intercept AND NULL
  # VARIABLES
  # TODO: insample_point_estimates RETURNS THE SAME INFO TWICE IF
  # error.measure.list IS SPECIFIED.
  #
  #use.intercept <- match.arg(use.intercept)

  # > Last Observation #################################################### ----

  if (is.null(last.obs)) {

    if (stats::is.ts(DATA)) {
      last.obs <- length(DATA)
    } else {
      last.obs <- dim(DATA)[1]
    }

  }

  # > Training and Test Data Sets ######################################### ----

  if (stats::is.ts(DATA)) {

    tmpStart <- stats::tsp(DATA)[1]
    tmpEnd   <- stats::tsp(DATA)[2]
    tmpFreq  <- stats::tsp(DATA)[3]

    if (tmpFreq > 1) {

      tmpTrainEnd    <- tmpEnd - ((forecast.horizon + tmpFreq) / tmpFreq)

      tmpTestStart   <- tmpTrainEnd + (1 / tmpFreq)
      tmpTestEnd     <- tmpEnd # or use tmpEnd - (tmpFreq / tmpFreq)

    } else {

      tmpTrainEnd    <- tmpEnd - (forecast.horizon / tmpFreq)

      tmpTestStart   <- tmpTrainEnd + (1 / tmpFreq)
      tmpTestEnd     <- tmpEnd # or use tmpEnd - (tmpFreq / tmpFreq)

    }

    trainingData <- stats::window(
      x     = DATA,
      start = tmpStart,
      end   = tmpTrainEnd
    )

    testData <- stats::window(
      x     = DATA,
      start = tmpTestStart,
      end   = tmpTestEnd
    )


  } else {

    tmpFreq <- stats::frequency(DATA)

    tmpTrainEnd <- last.obs - forecast.horizon

    tmpTestStart <- last.obs - forecast.horizon + 1
    tmpTestEnd   <- last.obs

    trainingData <- DATA[1:tmpTrainEnd, ]
    testData     <- DATA[tmpTestStart:tmpTestEnd, ]

  }

  # > All Combinations Size and Number #################################### ----

  allCombinations <- expand.grid(
    'size.rs'   = search.size.rs,
    'number.rs' = search.number.rs
  )

  nCombinations <- dim(allCombinations)[1]

  # > Run GEARS ########################################################### ----

  forecastErrorsL <- lapply(
    X = seq(1:nCombinations),
    function(X) {

      tmpGears <- gears(
        DATA      = trainingData,
        size.rs   = allCombinations[X, "size.rs"],
        number.rs = allCombinations[X, "number.rs"],
        last.obs  = tmpTrainEnd,
        forecast.horizon,
        glm.family,
        level,
        details,
        y.name,
        y.max.lags,
        x.names,
        x.max.lags,
        x.fixed.names,
        x.fixed.lags,
        x.interaction.names,
        x.interaction.lags,
        use.intercept,
        error.measure,
        betas.selection
      )

      # Forecast Error
      error_measures(
        forecasts         = tmpGears$out_sample_forecasts,
        outsample         = testData,
        insample          = trainingData,
        ts.frequency      = tmpFreq,
        forecast.horizon  = forecast.horizon,
        alpha.level       = 1 - (level / 100),
        error.measure     = error.measure
      )

      # TODO: add betas.selection to see which is best, also intercept

    }
  )

  # > Get Minimum ######################################################### ----

  tmpMinimumError <- min(unlist(forecastErrorsL))

  tmpMinimumObs   <- which(do.call(rbind, forecastErrorsL) == tmpMinimumError)

  # > Add it to the Selected allCombinations ############################## ----

  allCombinationsForecasts <- cbind(
    allCombinations[tmpMinimumObs, ],
    rep(tmpMinimumError, length(tmpMinimumObs))
  )

  names(allCombinationsForecasts)[3] <- error.measure

  # > Return ############################################################## ----

  allCombinationsForecasts

}
