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
#'     forecast errors.
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
#' @return An object of class "forecast.gears".
#'     An object of class "forecast.gears" is a list containing the following
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
#'
#' @export
#'
#' @examples
#' # Univariate Time Series Forecasting - Data of class 'ts'.
#' gears(
#'   DATA                = datasets::WWWusage,
#'   forecast.horizon    = 10,
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
#'   x.names             = list("BEEF_PRICE", "WHEAT_PRICE") ,
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
                  ...) { # ... TO ACCOUNT FOR OTHER OPTIONS PASSED TO glm


  # > CHECKS ############################################################## ----
  #
  # TODO: INCLUDE CHECKS FROM OLD FILE "new_gears.R"
  # TODO: FIX PROBLEMS WHEN NOT-SPECIFYING last.obs, use.intercept AND NULL
  # VARIABLES
  # TODO: insample_point_estimates RETURNS THE SAME INFO TWICE IF
  # error.measure.list IS SPECIFIED.
  #
  use.intercept <- match.arg(use.intercept)

  # > Helpers ############################################################# ----

  # > DATA ################################################################ ----

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

  # |__ Y variable name ========================================================

  if (is.null(y.name)) {
    yName <- "Y_t"
  } else {
    yName <- paste0(y.name, "_t")
  }

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

  # > ESTIMATION ########################################################## ----

  # |__ Fit and Predict ========================================================
  ## Returns a list where the first level is the forecast lead, and inside each
  ## level there is a table that returns the forecast by equation/model number
  ## (column) and rolling sample (row).

  # TODO: https://www.r-bloggers.com/no-visible-binding-for-global-variable/
  predictionGEARS <- lapply(
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
            DF.Fit.Predict    = dfFitPredict#,
            #... # further arguments passed on to glm
          )
        }
      )
    }
  )

  # |__ Get Prediction Errors ==================================================

  ## \____ All measures --------------------------------------------------------
  ## Returns a list where the first level is the forecast lead, and inside each
  ## level there is a table that returns the prediction errors by equation/
  ## model number (row) and error measure (column)

  predictionErrorsAll <- prediction_errors(
    DATA,
    forecast.horizon,
    Y.name          = yName,
    total.equations = totalEquations,
    number.rs,
    DF.Fit.Predict  = dfFitPredict,
    forecasts.gears = predictionGEARS
  )

  ## \____ User's measure-------------------------------------------------------
  ## Returns a table with the prediction errors by equation/model number (row)
  ## and forecast lead (column), only for the selected error measures.

  predictionErrors <- sapply(
    X = 1:forecast.horizon,
    function(X) predictionErrorsAll[[X]][, error.measure]
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

  if (betas.selection == "last") {
    outSampleForecasts <- outSampleForecastsAll[number.rs, ]
  } else if (betas.selection == "average") {
    outSampleForecasts <- colMeans(outSampleForecastsAll)
  } else {
    outSampleForecasts <- rbind(
      outSampleForecastsAll[number.rs, ],
      colMeans(outSampleForecastsAll)
    )
  }

  # > Results GEAR ######################################################## ----

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
      predictionErrorsAll[[h]][eq, ]
    }
  )

  meanPredictionErrorsUser <- cbind(
    seq(1:forecast.horizon),
    meanPredictionErrors[error.measure, ]
  )

  colnames(meanPredictionErrorsUser) <- c("forecast.lead", error.measure)

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

  tmp_lower_fct <- function(X, forecasts.gears) {
    forecasts.gears[X] - tCrit * meanPredictionErrors["mse", X]
  }

  tmp_upper_fct <- function(X, forecasts.gears) {
    forecasts.gears[X] + tCrit * meanPredictionErrors["mse", X]
  }

  if (betas.selection != "both"){
    lower <- sapply(1:forecast.horizon, tmp_lower_fct, outSampleForecasts)
    upper <- sapply(1:forecast.horizon, tmp_upper_fct, outSampleForecasts)
  } else {
    lower <- rbind(
      sapply(1:forecast.horizon, tmp_lower_fct, outSampleForecasts[1, ]),
      sapply(1:forecast.horizon, tmp_lower_fct, outSampleForecasts[2, ])
    )

    upper <- rbind(
      sapply(1:forecast.horizon, tmp_upper_fct, outSampleForecasts[1, ]),
      sapply(1:forecast.horizon, tmp_upper_fct, outSampleForecasts[2, ])
    )
  }

  names(lower) <- names(upper) <- NULL

  ## |_ Create Details Object ==================================================

  if (isTRUE(details)) {
    detailsGEARS <- lapply(
      X = 1:forecast.horizon,
      function(X) {
        list(
          best_equation = bestEquations[X],
          best_models   = c(bestModel[[X]]),
          fitted        = inSamplePredictions[ , X],
          residuals     = inSamplePredictions[ , X] -
            dfForecast[(last.obs - number.rs + 1):last.obs, yName]
        )
      }
    )
  } else {
    detailsGEARS <- lapply(
      X = 1:forecast.horizon,
      function(X) {
        list(
          best_equation = bestEquations[X],
          fitted        = inSamplePredictions[ , X],
          residuals     = inSamplePredictions[ , X] -
            dfForecast[(last.obs - number.rs + 1):last.obs, yName]
        )
      }
    )
  }

  tmpDetailsNames <- paste0("forecast_lead_", 1:forecast.horizon)

  names(detailsGEARS) <- tmpDetailsNames

  # > Return Results ###################################################### ----

  outGEARS <- list(
    out_sample_forecasts  = outSampleForecasts,
    lower                 = lower,
    upper                 = upper,
    level                 = level,
    betas                 = betas.selection,
    error_measure         = error.measure,
    min_prediction_errors = meanPredictionErrorsUser,
    details               = detailsGEARS
  )

  outGEARS$model$call <- match.call()

  return(structure(outGEARS, class = "gears"))
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ----
} # UNCOMMENT THIS WHEN DONE WITH FUNCTION, TO CLOSE MAIN BRACKET
print.gears <- function(x, ...){
  #print(paste0("Best GEARS Model:  ",x$best.title))
  print(cbind(x$out_sample_forecasts, x$lower, x$upper))
}
# END ----
