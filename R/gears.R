#' GEARS Method for Time Series Forecasting
#'
#' Use the GEARS Method for forecasting using data from univariate time series
#' objects or from data frame objects.
#'
#' @param DATA A data frame or a univariate time series.
#' @param forecast.horizon Number of periods for forecasting.
#' @param size.rs Number of observations in the rolling sample.
#' @param number.rs Number of rolling samples.
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
#' @return An object of class "\code{gears}". The function \code{summary} is
#'     used to obtain and print a summary of the results,
#'     while the function \code{plot} produces a plot of the forecasts and
#'     prediction intervals.
#'     An object of class "\code{gears}" is a list containing the following
#'     elements:
#'     \describe{
#'      \item{out_sample_point_forecasts}{If \code{betas.selection != "both"},
#'          the output is a \code{ts} object containing the out-of-sample
#'          point forecasts obtained from the best model.
#'          Otherwise, it is a list with two \code{ts} objects,
#'          each containing the out-of-sample point forecasts for the two cases.
#'          The first one is for when \code{betas.selection == "last"}, and the
#'          second one is for when \code{betas.selection == "average"}.}
#'      \item{lower}{If \code{betas.selection != "both"},
#'          the output is a \code{ts} object containing the lower limits for
#'          prediction intervals. Otherwise, it is a list with two
#'          \code{ts} objects, each containing the the lower limits for the two
#'          cases. The first one is for when \code{betas.selection == "last"},
#'          and the second one is for when \code{betas.selection == "average"}.}
#'      \item{upper}{If \code{betas.selection != "both"},
#'          the output is a \code{ts} object containing the upper limits for
#'          prediction intervals. Otherwise, it is a list with two
#'          \code{ts} objects, each containing the the upper limits for the two
#'          cases. The first one is for when \code{betas.selection == "last"},
#'          and the second one is for when \code{betas.selection == "average"}.}
#'      \item{x}{A \code{ts} object with the original time series.
#'          If \code{DATA} is a data.frame, it returns the dependent variable as
#'          a \code{ts} object.}
#'      \item{x.fitted}{A \code{ts} object of \code{length = number.rs} with
#'          the in-sample predictions.}
#'      \item{total_models}{The total number of models evaluated by
#'          the function.}
#'      \item{total_equations_estimated}{The total number of equations estimated
#'          by the function. It is the value of \code{total_models} times the
#'          number of rolling samples (\code{number.rs}) and the forecast
#'          horizon (\code{forecast.horizon}).}
#'      \item{betas}{The \code{beta.selection} used.}
#'      \item{error_measure}{The \code{error.measure} used.}
#'      \item{min_prediction_errors}{A matrix object with two columns - one for
#'          the forecast lead number, and the other for the values of the
#'          minimum prediction error for the given measure -, and number of rows
#'          equal to \code{number.rs}.}
#'      \item{mse_out_sample_prediction}{A matrix object with two columns - one
#'          for the forecast lead number, and the other for the values of the
#'          MSE for the best model -, and number of rows
#'          equal to \code{number.rs}.}
#'      \item{details}{A list containing information about
#'          the best fitted model for each forecast lead. If
#'          \code{details == FALSE}, the list contains information on the best
#'          equation estimated for that particular forecast lead, and its
#'          associated in-sample predictions (for all rolling samples) and their
#'          associated prediction errors. If \code{details == TRUE}, it also
#'          includes all the information from the \link[stats]{glm} fit for each
#'          forecast lead.}
#'     \item{model}{A list with the model call used in this function.}
#'     }
#' @author Gustavo Varela-Alvarenga
#'
#' @keywords ts glm
#'
#' @examples
#' # Univariate Time Series Forecasting - Data of class 'ts'.
#' gears(
#'   DATA                = datasets::WWWusage,
#'   forecast.horizon    = 12,
#'   size.rs             = 20,
#'   number.rs           = 10,
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
#'
#' # Using Parallel Computing - Univariate Time Series Forecasting - Data of
#' # class 'ts'.
#' gears:::gears(
#'   DATA                = datasets::WWWusage,
#'   forecast.horizon    = 10,
#'   size.rs             = 20,
#'   number.rs           = 12,
#'   glm.family          = "quasi",
#'   y.max.lags          = 2,
#'   error.measure       = "mse",
#'   betas.selection     = "last",
#'   use.parallel        = TRUE,
#'   num.cores           = 2
#' )
#'
#' # Univariate Time Series Forecasting - Data from a data frame.
#' gears(
#'   DATA                = gears::commodities_prices,
#'   forecast.horizon    = 5,
#'   size.rs             = 20,
#'   number.rs           = 12,
#'   level               = 95,
#'   details             = FALSE,
#'   glm.family          = "quasi",
#'   y.name              = "PORK_PRICE",
#'   y.max.lags          = 2,
#'   x.names             = NULL,
#'   x.max.lags          = NULL,
#'   x.fixed.names       = NULL,
#'   x.fixed.lags        = NULL,
#'   x.interaction.names = NULL,
#'   x.interaction.lags  = NULL,
#'   last.obs            = 100,
#'   use.intercept       = "both",
#'   error.measure       = "mse",
#'   betas.selection     = "last"
#' )
#'
#' # Multivariate Time Series Forecasting
#' gears(
#'   DATA                = gears::commodities_prices,
#'   forecast.horizon    = 1,
#'   size.rs             = 12,
#'   number.rs           = 5,
#'   level               = 95,
#'   details             = FALSE,
#'   glm.family          = "quasi",
#'   y.name              = "PORK_PRICE" ,
#'   y.max.lags          = 1,
#'   x.names             = list("BEEF_PRICE", "WHEAT_PRICE"),
#'   x.max.lags          = list(2, 1),
#'   x.fixed.names       = NULL,
#'   x.fixed.lags        = NULL,
#'   x.interaction.names = NULL,
#'   x.interaction.lags  = NULL,
#'   last.obs            = 150,
#'   use.intercept       = "both",
#'   error.measure       = "mse" ,
#'   betas.selection     = "last"
#' )
#'
#' # Multivariate Time Series Forecasting - With Interactions
#' gears(
#'   DATA                = gears::commodities_prices,
#'   forecast.horizon    = 1,
#'   size.rs             = 25,
#'   number.rs           = 5,
#'   level               = 95,
#'   details             = FALSE,
#'   glm.family          = "quasi",
#'   y.name              = "PORK_PRICE",
#'   y.max.lags          = 1,
#'   x.names             = list("BEEF_PRICE"),
#'   x.max.lags          = list(2),
#'   x.fixed.names       = list("CORN_PRICE", "POULTRY_PRICE"),
#'   x.fixed.lags        = list(4, 5),
#'   x.interaction.names = list("CORN_PRICE*CORN_PRICE", "WHEAT_PRICE*BEEF_PRICE"),
#'   x.interaction.lags  = list(c(1, 1), c(2,3)),
#'   last.obs            = 150,
#'   use.intercept       = "without",
#'   error.measure       = "mse",
#'   betas.selection     = "last"
#' )
#'
#' @export
gears <- function(DATA,
                  forecast.horizon,
                  size.rs,
                  number.rs,
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
                  ...) { # ... TO ACCOUNT FOR OTHER OPTIONS PASSED TO glm


  # > TODO ################################################################ ----
  #
  # TODO: FIX PROBLEMS WHEN NOT-SPECIFYING NULL VARIABLES
  # TODO: Check problems with checks under Univariate TS Inputs
  #
  # > CHECK ARGUMENTS ##################################################### ----

  # |__ Y variable name ========================================================

  if (is.null(y.name)) {
    yName  <- "Y_t"
  } else {
    yName  <- paste0(y.name, "_t")
  }

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
    size.rs,
    number.rs,
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

  # |__ DF Fit Predict =========================================================
  ## Used as argument for the "fit_predict", "fit_best", and "prediction_erros"
  ## functions.

  dfFitPredict <- create_DF_Fit_Predict(
    DATA,
    forecast.horizon,
    size.rs,
    number.rs,
    y.name,
    y.max.lags,
    x.names,
    x.max.lags,
    x.fixed.names,
    x.fixed.lags,
    x.interaction.names,
    x.interaction.lags,
    last.obs
  )

  # |__ DF Forecast ============================================================

  dfForecast <- create_DF_Forecast(
    DATA,
    forecast.horizon,
    y.name,
    y.max.lags,
    x.names,
    x.max.lags,
    x.fixed.names,
    x.fixed.lags,
    x.interaction.names,
    x.interaction.lags
  )

  # > ALL POSSIBLE MODELS/EQUATIONS ####################################### ----

  # |__ All equations rhs ======================================================

  allEquationsRhs <- all_models_rhs(
    y.name,
    y.max.lags,
    x.names,
    x.max.lags,
    x.fixed.names,
    x.fixed.lags,
    x.interaction.names,
    x.interaction.lags,
    use.intercept
  )

  # |__ Total Number of Equations ==============================================

  totalEquations <- length(allEquationsRhs)

  # |__ Total Number of Equations Estimated ====================================

  totalEquationsEstimated <- totalEquations * number.rs * forecast.horizon

  # > PARALLEL ############################################################ ----
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

    if (Sys.info()[['sysname']] == "Darwin") {
      tmpCluster <- parallel::makeCluster(no_cores, type = "FORK")
      on.exit(parallel::stopCluster(tmpCluster), add = TRUE)
    } else if (Sys.info()[['sysname']] == "Windows"){

      tmpCluster <- parallel::makeCluster(no_cores)
      on.exit(parallel::stopCluster(tmpCluster), add = TRUE)

      parallel::clusterExport(
        cl      = tmpCluster,
        envir   = environment(),
        varlist = list(
          "dfFitPredict",
          "dfForecast",
          "yName",
          "allEquationsRhs",
          "totalEquations"
        )
      )

      parallel::clusterCall(
        cl = tmpCluster,
        assign,
        "fit_predict",
        fit_predict,
        envir = .GlobalEnv
      )
    } else {
      tmpCluster <- parallel::makeForkCluster(nnodes = no_cores)
      on.exit(parallel::stopCluster(tmpCluster), add = TRUE)
    }

    # Temporary lapply function
    tmFcnPar <- function(...) parallel::parLapply(cl = tmpCluster, ...)

  } else {
    tmFcnPar <- function(...) lapply(...)
  }

  # > ESTIMATION ########################################################## ----

  # |__ Fit and Predict ========================================================
  ## Returns a list where the first level is the forecast lead, and inside each
  ## level there is a table that returns the forecast by equation/model number
  ## (column) and rolling sample (row).

  predictionGEARS <- tmFcnPar(
    X   = 1:forecast.horizon,
    function(X) {
      h <- X
      sapply(
        X = 1:totalEquations,
        #function(Eq = X) {
        function(X) {   # Using this instead of what's on the previous line to
          Eq <- X       # address the issue: "no visible binding for global
                        # variable" generated by "R CMD check"
          sapply(
            X                 = 1:number.rs,
            FUN               = fit_predict,
            forecast.lead     = h,
            eq.number         = Eq,
            glm.family        = glm.family,
            Y.name            = yName,
            all.equations.rhs = allEquationsRhs,
            DF.Fit.Predict    = dfFitPredict,
            ... # further arguments passed on to glm
          )
        }
      )
    }
  )

  # |__ Get Prediction Errors ==================================================

  ## \____ Selected Measure ----------------------------------------------------
  ## Returns a list where the first level is the forecast lead, and inside each
  ## level there is a table that returns the prediction errors by equation/
  ## model number (row) and error measure (column)

  predictionErrorsUser <- prediction_errors(
    DATA             = DATA,
    forecast.horizon = forecast.horizon,
    Y.name           = yName,
    total.equations  = totalEquations,
    number.rs        = number.rs,
    DF.Fit.Predict   = dfFitPredict,
    forecasts.gears  = predictionGEARS,
    names.measures   = error.measure
  )

  predictionErrors <- sapply(
    X = 1:forecast.horizon,
    function(X) predictionErrorsUser[[X]][, error.measure]
  )

  ## \____ MSE -----------------------------------------------------------------

  ## This will be used for the confidence intervals

  predictionMSE <- prediction_errors(
    DATA             = DATA,
    forecast.horizon = forecast.horizon,
    Y.name           = yName,
    total.equations  = totalEquations,
    number.rs        = number.rs,
    DF.Fit.Predict   = dfFitPredict,
    forecasts.gears  = predictionGEARS,
    names.measures   = "mse"
  )

  # > BEST EQUATION ####################################################### ----

  # |__ Get Equations that Min. Selected Error Function ========================

  equationsMinError <- apply(predictionErrors, 2, which.min)

  # |__ Estimate Best Model ====================================================

  ## Returns a list of size equal to forecast.horizon. The first level is the
  ## forecast lead. The second level is given by the number of the particular
  ## rolling sample. The output is the glm object associated with the best
  ## equation (that minimizes the prediction errors), for each rolling sample
  ## and each forecast lead.
  ## E.g., best_model[[2]][[1]]  returns the glm object for the second
  ##                       ...   forecast lead, and for the first equation
  ##                       ...   number.
  ##                       [[number.rs]
  ##

  bestModel <- lapply(
    X   = 1:forecast.horizon,
    #function(h = X) {
    function(X) {   # Using this instead of what's on the previous line to
      h <- X        # address the issue: "no visible binding for global
      # variable" generated by "R CMD check"
      lapply(
        X                 = 1:number.rs,
        forecast.lead     = h,
        eq.number         = equationsMinError[h],
        FUN               = fit_best,
        glm.family        = glm.family,
        Y.name            = yName,
        all.equations.rhs = allEquationsRhs,
        DF.Fit.Predict    = dfFitPredict,
        ...
      )
    }
  )

  # > Out-of-Sample Forecasts ############################################# ----

  # |__ Get Forecasts for all random samples ===================================

  outSampleForecastsAll <- sapply(
    X   = 1:forecast.horizon,
    #function(h = X) {
    function(X) {   # Using this instead of what's on the previous line to
      h <- X        # address the issue: "no visible binding for global
                    # variable" generated by "R CMD check"
      sapply(
        X = 1:number.rs,
        function(X) {
          stats::predict(
            object  = bestModel[[h]][[X]],
            newdata = dfForecast[last.obs, ],
            type    = "response"
          )
        }
      )
    }
  )

  ## \____ Get the Forecasts based on the user's Betas selection ---------------

  tmpOutSampleForecasts <- list(
    "betas.last"    = outSampleForecastsAll[number.rs, ],
    "betas.average" = colMeans(outSampleForecastsAll)
  )

  if (betas.selection == "last") {
    outSampleForecasts <- tmpOutSampleForecasts$betas.last
  } else if (betas.selection == "average") {
    outSampleForecasts <- tmpOutSampleForecasts$betas.average
  } else {
    outSampleForecasts <- tmpOutSampleForecasts
  }

  # > Results GEARS ####################################################### ----

  # |_ Best Fitted Models ======================================================

  #> bestModel

  # |_ Best Equation for Each Forecast Lead ====================================

  bestEquations <- as.character(lapply(
    X = 1: forecast.horizon,
    function(X) stats::formula(bestModel[[X]][[1]])
  ))

  # |_ Out-of-sample Predictions - Best Model ==================================
  ## \hat{y}_{89}, ..., \hat{y}_100
  ## Calling it "in-sample predictions to avoid confusion"

  inSamplePredictions <- sapply(
    X = 1:forecast.horizon,
    #function(h = X) {
    function(X) {   # Using this instead of what's on the previous line to
      h <- X        # address the issue: "no visible binding for global
                    # variable" generated by "R CMD check"
      eq <- equationsMinError[h]
      predictionGEARS[[h]][, eq]
    }
  )

  # |_ Mean Out-of-sample Prediction Errors - Best Model =======================

  meanPredictionErrors <- sapply(
    X = 1:forecast.horizon,
    #function(h = X) {
    function(X) {   # Using this instead of what's on the previous line to
      h <- X        # address the issue: "no visible binding for global
      # variable" generated by "R CMD check"
      eq <- equationsMinError[h]
      predictionErrorsUser[[h]][eq, ]
    }
  )

  # meanPredictionErrorsUser <- cbind(
  #   seq(1:forecast.horizon),
  #   meanPredictionErrors[error.measure, ]
  # )

  meanPredictionErrorsUser <- cbind(
    seq(1:forecast.horizon),
    as.numeric(meanPredictionErrors)
  )

  colnames(meanPredictionErrorsUser) <- c("forecast.lead", error.measure)

  # |_ MSE Out-of-sample Prediction - Best Model ===============================

  bestMSE <- sapply(
    X = 1:forecast.horizon,
    #function(h = X) {
    function(X) {   # Using this instead of what's on the previous line to
      h <- X        # address the issue: "no visible binding for global
      # variable" generated by "R CMD check"
      eq <- equationsMinError[h]
      predictionMSE[[h]][eq, ]
    }
  )

  bestOutMSE <- cbind(
    seq(1:forecast.horizon),
    as.numeric(bestMSE)
  )

  colnames(bestOutMSE) <- c("forecast.lead", "mse")


  ## |_ Get Prediction Interval ================================================

  ## TODO: Use df = best_model[[1]][[1]]$df.residual to get df based on degrees
  ## of freedom from the model (depends on the number of parameters)
  ## OR
  ## Use = size.rs - 2 or Use number.rs - 2 since numbers.rs is what was used to
  ## calculate the MSE?
  ## Going with number.rs - 2 for now.

  tCrit <- stats::qt(
    p = 0.5 - level / 200,
    df = number.rs - 2,
    lower.tail = FALSE
  )

  fct_lower <- function(X, forecasts.gears) {
    forecasts.gears[X] - tCrit * sqrt(as.numeric(bestMSE[X]))
  }

  fct_upper <- function(X, forecasts.gears) {
    forecasts.gears[X] + tCrit * sqrt(as.numeric(bestMSE[X]))
  }

  if (betas.selection == "both"){
    lower <- lapply(
      X = 1:length(outSampleForecasts),
      function(X) {
        tmp <- sapply(1:forecast.horizon, fct_lower, outSampleForecasts[[X]])
        names(tmp) <- NULL
        tmp
      }
    )

    upper <- lapply(
      X = 1:length(outSampleForecasts),
      function(X) {
        tmp <- sapply(1:forecast.horizon, fct_upper, outSampleForecasts[[X]])
        names(tmp) <- NULL
        tmp
      }
    )

  } else {
    lower <- sapply(1:forecast.horizon, fct_lower, outSampleForecasts)
    upper <- sapply(1:forecast.horizon, fct_upper, outSampleForecasts)

    names(lower) <- names(upper) <- NULL
  }



  ## |_ Transform Results to ts objects  =======================================

  if (stats::is.ts(DATA)) {
    tmpFreq     <- stats::frequency(DATA)
    tmpStartOut <- stats::tsp(DATA)[2] + (1 / tmpFreq)
    tmpStartX   <- stats::tsp(DATA)[2] + ((1 - number.rs) / tmpFreq)
    tmpEndX     <- stats::tsp(DATA)[2]

    tmpX <- stats::window(
      x         = DATA,
      start     = tmpStartX,
      end       = tmpEndX,
      frequency = tmpFreq
    )

    tmpXAll <- DATA

  } else {
    tmpFreq     <- 1
    tmpStartOut <- last.obs + (1 / tmpFreq)
    tmpStartX   <- last.obs + ((1 - number.rs) / tmpFreq)
    tmpEndX     <- last.obs

    tmpX <- stats::ts(
      data      = DATA[, y.name],
      start     = tmpStartX,
      end       = tmpEndX,
      frequency = tmpFreq
    )

    tmpXAll <- stats::ts(
      data      = DATA[, y.name],
      start     = 1,
      end       = tmpEndX,
      frequency = tmpFreq
    )
  }

  if (betas.selection == "both"){
    outSampleForecastsTS <- lapply(
      X = 1:length(outSampleForecasts),
      function(X) {
        stats::ts(
          data      = outSampleForecasts[[X]],
          start     = tmpStartOut,
          frequency = tmpFreq
        )
      }
    )

    lowerTS <- lapply(
      X = 1:length(lower),
      function(X) {
        stats::ts(data = lower[[X]], start = tmpStartOut, frequency = tmpFreq)
      }
    )

    upperTS <- lapply(
      X = 1:length(upper),
      function(X) {
        stats::ts(data = upper[[X]], start = tmpStartOut, frequency = tmpFreq)
      }
    )

    names(outSampleForecastsTS) <- c("betas.last", "betas.average")
    names(lowerTS) <- names(upperTS) <- c("betas.last", "betas.average")

  } else {

    outSampleForecastsTS <- stats::ts(
      data      = outSampleForecasts,
      start     = tmpStartOut,
      frequency = tmpFreq
    )

    lowerTS <- stats::ts(data = lower, start = tmpStartOut, frequency = tmpFreq)
    upperTS <- stats::ts(data = upper, start = tmpStartOut, frequency = tmpFreq)

  }

  inSamplePredictionsTS <- lapply(
    X = 1:forecast.horizon,
    function(X) {
      tmpFit <- stats::ts(as.numeric(inSamplePredictions[ , X]))
      stats::tsp(tmpFit) <- stats::tsp(tmpX)
      tmpFit
    }
  )

  predictionErrorTS <- lapply(
    X = 1:forecast.horizon,
    function(X) inSamplePredictionsTS[[X]] - tmpX
  )

  ## |_ Create Details Object ==================================================

  if (isTRUE(details)) {
    detailsGEARS <- lapply(
      X = 1:forecast.horizon,
      function(X) {
        list(
          best_equation          = bestEquations[X],
          best_models            = c(bestModel[[X]]),
          in_sample_predictions  = inSamplePredictionsTS[[X]],
          prediction_errors      = predictionErrorTS[[X]]
        )
      }
    )
  } else {
    detailsGEARS <- lapply(
      X = 1:forecast.horizon,
      function(X) {
        list(
          best_equation         = bestEquations[X],
          in_sample_predictions = inSamplePredictionsTS[[X]],
          prediction_errors     = predictionErrorTS[[X]]
        )
      }
    )
  }

  tmpDetailsNames <- paste0("forecast_lead_", 1:forecast.horizon)

  names(detailsGEARS) <- tmpDetailsNames

  # > Return Results ###################################################### ----
  outGEARS <- list(
    out_sample_forecasts      = outSampleForecastsTS,
    lower                     = lowerTS,
    upper                     = upperTS,
    x                         = tmpXAll,
    x.fitted                  = tmpX,
    total_models              = totalEquations,
    total_equations_estimated = totalEquationsEstimated,
    level                     = level,
    betas                     = betas.selection,
    error_measure             = error.measure,
    min_prediction_errors     = meanPredictionErrorsUser,
    mse_out_sample_prediction = bestOutMSE,
    details                   = detailsGEARS
  )

  outGEARS$model$call <- match.call()

  return(structure(outGEARS, class = "gears"))
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ----
}
# END ----
