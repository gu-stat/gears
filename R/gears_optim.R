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
#' @param use.parallel If parallel computing should be used. Default is
#'     \code{FALSE}.
#' @param num.cores Number of cores to use if parallel computing is used.
#'     If \code{use.parallel == TRUE}, then the default is
#'     \code{detectCores() - 1}. For more details, see
#'     \link[parallel]{detectCores}.
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
#' @return A data frame with five columns: "size.rs", "number.rs", "intercept",
#'    "betas" and the last one named after the selected error measure. The last
#'    column returns the minimum value for the selected error measure, and the
#'    columns "size.rs" and "number.rs" return the sample sizes and the number
#'    of rolling samples used to get this minimum. The column "intercept" gives
#'    information on whether the intercept was used ("with"), or not
#'    ("without") on the best model. The "betas" column tells which
#'    \code{betas.selection} returns the best results.
#'
#' @author Gustavo Varela-Alvarenga
#'
#' @keywords ts glm optimization
#'
#' @examples
#' # Univariate Time Series Forecasting - Data of class 'ts'.
#' # Using betas.selection = "last"
#' gears_optim(
#'   DATA                = datasets::WWWusage,
#'   forecast.horizon    = 12,
#'   search.size.rs      = c(20, 30),
#'   search.number.rs    = c(10, 12),
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
#'   betas.selection     = "last",
#'   use.parallel        = FALSE,
#'    num.cores          = NULL
#' )
#'
#' # Univariate Time Series Forecasting - Data of class 'ts'.
#' # Using betas.selection = "both"
#' gears_optim(
#'   DATA                = datasets::WWWusage,
#'   forecast.horizon    = 12,
#'   search.size.rs      = c(20, 30),
#'   search.number.rs    = c(10, 12),
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
#'   betas.selection     = "both",
#'   use.parallel        = FALSE,
#'   num.cores           = NULL
#' )
#'
#' gears_optim(
#'   DATA                = datasets::WWWusage,
#'   forecast.horizon    = 12,
#'   search.size.rs      = c(20, 30),
#'   search.number.rs    = c(10, 12),
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
#'   use.intercept       = "with",
#'   error.measure       = "mse",
#'   betas.selection     = "both",
#'   use.parallel        = FALSE,
#'   num.cores           = NULL
#' )
#'
#' # With Parallel Computing
#' # Univariate Time Series Forecasting - Data of class 'ts'.
#' ## If num.cores = NULL, the function uses
#' ## num.cores = parallel::detectCores() - 1
#' gears_optim(
#'   DATA                = datasets::WWWusage,
#'   forecast.horizon    = 12,
#'   search.size.rs      = c(20, 30),
#'   search.number.rs    = c(10, 12),
#'   level               = 95,
#'   glm.family          = "quasi",
#'   y.max.lags          = 2,
#'   use.intercept       = "both",
#'   error.measure       = "mse",
#'   betas.selection     = "last",
#'   use.parallel        = TRUE,
#'   num.cores           = NULL
#' )
#'
#' @export
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
                  use.parallel = FALSE,
                  num.cores = NULL,
                  ...) {

  # > CHECKS ############################################################## ----
  #
  # TODO: FIX PROBLEMS WHEN NOT-SPECIFYING NULL VARIABLES
  # TODO: Check problems with checks under Univariate TS Inputs
  # > CHECK ARGUMENTS ##################################################### ----

  # |__ GLM Family =============================================================

  if (missing(glm.family)){
    glm.family <- "quasi"
  } else {
    glm.family <- match.arg(glm.family)
  }

  # |__ Last Observation =======================================================

  if (is.null(last.obs)) {
    if (stats::is.ts(DATA)) {
      last.obs <- length(DATA)
    } else {
      last.obs <- dim(DATA)[1]
    }
  }

  # |__ Intercept ==============================================================

  if (missing(use.intercept)){
    use.intercept <- "both"
  } else {
    use.intercept <- match.arg(use.intercept)
  }

  # |__ Error Measure ==========================================================

  if (missing(error.measure)){
    error.measure <- "mse"
  } else {
    error.measure <- match.arg(error.measure)
  }

  # |__ Betas Selection ========================================================

  if (missing(betas.selection)){
    betas.selection <- "both"
  } else {
    betas.selection <- match.arg(betas.selection)
  }

  # |__ Checks that STOP the function ==========================================

  checks(
    DATA,
    forecast.horizon,
    size.rs = max(search.size.rs),
    number.rs = max(search.number.rs),
    glm.family,
    level,
    y.name,
    y.max.lags,
    x.names,
    x.max.lags,
    x.fixed.names,
    x.fixed.lags,
    x.interaction.names,
    x.interaction.lags,
    last.obs,
    use.intercept,
    error.measure,
    betas.selection
  )

  # # |__ Univariate TS Inputs ===================================================
  #
  # if (class(DATA)[1] == "ts") {
  #
  #   if (is.null(y.max.lags) & !is.null(x.max.lags)) {
  #     warning(paste0(
  #       "DATA input is an univariate time series. Function will use the highest ",
  #       "value in 'x.max.lags' as the maximum number of lags of Y."
  #     ))
  #
  #     y.max.lags <- max(unlist(x.max.lags))
  #     x.max.lags <- NULL
  #
  #   } else if (!is.null(y.max.lags) & !is.null(x.max.lags)) {
  #
  #     warning(paste0(
  #       "DATA input is an univariate time series. Function will use ",
  #       "'y.max.lags' as the maximum number of lags of Y and will disregard ",
  #       "'x.max.lags'."
  #     ))
  #
  #     x.max.lags <- NULL
  #
  #     if (is.null(y.name)) {
  #       if (is.null(x.names)) {
  #         y.name <- "Y"
  #       } else {
  #         if (length(x.names) > 1) {
  #           warning(paste0(
  #             "DATA input is an univariate time series. Function will ",
  #             "use 'y.max.lags' as the maximum number of lags of Y and will ",
  #             "disregard 'x.max.lags'."
  #           ))
  #         }
  #       }
  #     }
  #
  #     if (length(y.max.lags) > 1) {
  #       warning(paste0(
  #         "Variable y.max.lags must have only one value. ",
  #         "Function will use the highest value in y.max.lags as the ",
  #         "maximum number of lags of Y."
  #       ))
  #       y.max.lags <- max(y.max.lags)
  #     }
  #
  #   } else if (is.null(y.max.lags) & is.null(x.max.lags)) {
  #
  #     warning(paste0(
  #       "DATA input is an univariate time series. Function will use ",
  #       "'x.fixed.lags' as the number of lags of Y."
  #     ))
  #
  #     if (is.null(x.fixed.names)) {
  #       if (is.null(y.name)) y.name <- "Y"
  #       x.fixed.names <- lapply(1:length(x.fixed.lags), function(X) y.name)
  #     }
  #   }
  # }
  #
  # # |__ X variable name ========================================================
  #
  # if (class(DATA)[1] == "data.frame") {
  #   if (is.null(x.names) & !is.null(y.max.lags)) {
  #     x.names <- y.name
  #   }
  # }

  # > DATA MANIPULATION ################################################### ----

  # > Training and Test Data Sets ######################################### ----

  if (stats::is.ts(DATA)) {

    tmpStart <- stats::tsp(DATA)[1]
    tmpEnd   <- stats::tsp(DATA)[2]
    tmpFreq  <- stats::tsp(DATA)[3]

    if (tmpFreq > 1) {

      tmpTrainEnd    <- tmpEnd - ((forecast.horizon + tmpFreq) / tmpFreq)

      tmpTestStart   <- tmpTrainEnd + (1 / tmpFreq)
      tmpTestEnd     <- tmpEnd # or use tmpEnd - (tmpFreq / tmpFreq)

      tmpLastObs     <- which(stats::time(DATA) == tmpTrainEnd)

    } else {

      tmpTrainEnd    <- tmpEnd - (forecast.horizon / tmpFreq)

      tmpTestStart   <- tmpTrainEnd + (1 / tmpFreq)
      tmpTestEnd     <- tmpEnd # or use tmpEnd - (tmpFreq / tmpFreq)

      tmpLastObs     <- which(stats::time(DATA) == tmpTrainEnd)
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

    tmpLastObs   <- dim(trainingData)[1]

  }

  # > All Combinations Size and Number #################################### ----

  allCombinations <- expand.grid(
    'size.rs'   = search.size.rs,
    'number.rs' = search.number.rs
  )

  nCombinations <- dim(allCombinations)[1]

  # > Run GEARS ########################################################### ----

  # |_ Check for Parallel ======================================================

  if (isTRUE(use.parallel)){

    # Number of cores

    ## Do this to overcome CRAN's problem with more than 2 cores:
    if (is.null(num.cores)) {
      tmpCheck <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
      if (nzchar(tmpCheck) && tmpCheck == "TRUE") {
        # use 2 cores in CRAN/Travis/AppVeyor
        no_cores <- 2L
      } else {
        # use all cores in devtools::test()
        no_cores <- parallel::detectCores() - 1
      }
    } else {
      no_cores <- num.cores
    }

    # Initiate cluster

    if (Sys.info()[['sysname']] != "Windows") {
      tmpCluster <- parallel::makeCluster(no_cores, type = "FORK")
      on.exit(parallel::stopCluster(tmpCluster), add = TRUE)
    } else {

      tmpCluster <- parallel::makeCluster(no_cores)
      on.exit(parallel::stopCluster(tmpCluster), add = TRUE)

      parallel::clusterExport(
        cl      = tmpCluster,
        envir   = environment(),
        varlist = list(
          "nCombinations",
          "allCombinations",
          "trainingData",
          "testData",
          "tmpTrainEnd",
          "tmpFreq",
          "tmpLastObs"
        )
      )
      parallel::clusterEvalQ(tmpCluster, library("gears"))

    }

    # Temporary lapply function
    tmFcnPar <- function(...) parallel::parLapply(cl = tmpCluster, ...)

  } else {
    tmFcnPar <- function(...) lapply(...)
  }

  forecastErrorsL <- tmFcnPar(
    X = seq(1:nCombinations),
    function(X) {

      tmpSize   <- allCombinations[X, "size.rs"]
      tmpNumber <- allCombinations[X, "number.rs"]

      if (use.intercept == "both") {

        tmpList <- list("with", "without")

        tmpForecastErrorBoth <- lapply(
          X = seq_along(tmpList),
          function(X) {
            tmpGearsBoth <- gears(
              DATA      = trainingData,
              size.rs   = tmpSize,
              number.rs = tmpNumber,
              last.obs  = tmpLastObs,
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
              use.intercept = tmpList[[X]],
              error.measure,
              betas.selection
            )

            # Forecast Error

            if (betas.selection == "both") {
              tmpErrorBetas <- sapply(
                X = 1:2,
                function(X){
                  error_measures(
                    forecasts         = tmpGearsBoth$out_sample_forecasts[[X]],
                    outsample         = testData,
                    insample          = trainingData,
                    ts.frequency      = tmpFreq,
                    forecast.horizon  = forecast.horizon,
                    alpha.level       = 1 - (level / 100),
                    error.measure     = error.measure
                  )
                }
              )

              names(tmpErrorBetas) <- c("last", "average")

              tmpErrorBetas

            } else {
              tmpErrorBetas <- error_measures(
                forecasts         = tmpGearsBoth$out_sample_forecasts,
                outsample         = testData,
                insample          = trainingData,
                ts.frequency      = tmpFreq,
                forecast.horizon  = forecast.horizon,
                alpha.level       = 1 - (level / 100),
                error.measure     = error.measure
              )

              names(tmpErrorBetas) <- betas.selection

              tmpErrorBetas
            }

          }
        )

        names(tmpForecastErrorBoth) <- c(unlist(tmpList))

        tmpUnlistedErrors <- unlist(tmpForecastErrorBoth)

        tmpMinError       <- min(tmpUnlistedErrors)

        tmpForecastError  <- tmpUnlistedErrors[tmpUnlistedErrors==tmpMinError]

        tmpIntercept <- gsub("^(.*?)\\..*$", "\\1", names(tmpForecastError))
        tmpFinalBeta <- gsub(".*\\.","",names(tmpForecastError))

      } else {

        tmpGearsSingle <- gears(
          DATA      = trainingData,
          size.rs   = allCombinations[X, "size.rs"],
          number.rs = allCombinations[X, "number.rs"],
          last.obs  = tmpLastObs,
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

        if (betas.selection == "both") {
          tmpErrorBetas <- sapply(
            X = 1:2,
            function(X){
              error_measures(
                forecasts         = tmpGearsSingle$out_sample_forecasts[[X]],
                outsample         = testData,
                insample          = trainingData,
                ts.frequency      = tmpFreq,
                forecast.horizon  = forecast.horizon,
                alpha.level       = 1 - (level / 100),
                error.measure     = error.measure
              )
            }
          )

          names(tmpErrorBetas) <- c("last", "average")

        } else {
          tmpErrorBetas <- error_measures(
            forecasts         = tmpGearsSingle$out_sample_forecasts,
            outsample         = testData,
            insample          = trainingData,
            ts.frequency      = tmpFreq,
            forecast.horizon  = forecast.horizon,
            alpha.level       = 1 - (level / 100),
            error.measure     = error.measure
          )

          names(tmpErrorBetas) <- betas.selection

        }

        tmpUnlistedErrors <- unlist(tmpErrorBetas)

        tmpMinError       <- min(tmpUnlistedErrors)

        tmpForecastError  <- tmpUnlistedErrors[tmpUnlistedErrors==tmpMinError]

        tmpFinalBeta <- names(tmpForecastError)
        tmpIntercept <- rep(use.intercept, length(tmpFinalBeta))

      }

      # Return

      tmpReturn <- do.call(
        cbind, list(tmpIntercept, tmpFinalBeta, tmpForecastError)
      )

      colnames(tmpReturn) <- c("intercept", "betas", error.measure)
      rownames(tmpReturn) <- NULL

      tmpReturn

    }
  )

  # > Get Minimum ######################################################### ----

  tmpMinimumError <- min(unlist(forecastErrorsL))

  tmpMinimumObs   <- which(
    do.call(rbind, forecastErrorsL)[, error.measure] == tmpMinimumError
  )

  # > Add it to the Selected allCombinations ############################## ----

  allCombinationsForecasts <- cbind(
    allCombinations[tmpMinimumObs, ],
    do.call(rbind, lapply(X = tmpMinimumObs, function(X) forecastErrorsL[[X]]))
  )

  # > Return ############################################################## ----

  allCombinationsForecasts

}
