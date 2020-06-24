#' GEARS Method for Time Series Forecasting
#'
#' Use the GEARS Method for forecasting using data from univariate time series
#' objects or from data frame objects.
#'
#' @param DATA A data frame or a univariate time series.
#' @param forecast.horizon Number of periods for forecasting.
#' @param size.rs Number of observations in the rolling sample.
#' @param number.rs Number of rolling samples.
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
#'   size.rs             = 12,
#'   number.rs           = 10,
#'   glm.family          = "quasi",
#'   y.name              = "PORK_PRICE" ,
#'   y.max.lags          = 2,
#'   x.names             = NULL,
#'   x.max.lags          = NULL,
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
#' # Multivariate Time Series Forecasting
#' gears(
#'   DATA                = gears::commodities_prices,
#'   forecast.horizon    = 1,
#'   size.rs             = 12,
#'   number.rs           = 5,
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

  # |__ DF.Fit.Predict =========================================================

  DF.Fit.Predict <- gears:::create_DF_Fit_Predict(
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

  # |__ DF.Forecast ============================================================

  DF.Forecast <- gears:::create_DF_Forecast(
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
    Y_name <- "Y_t"
  } else {
    Y_name <- paste0(y.name, "_t")
  }

  # |__ All equations rhs ======================================================

  all.equations.rhs <- gears:::all_models_rhs(
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

  total.equations <- length(all.equations.rhs)

  # |__ Total Number of Equations Estimated ====================================

  total.equations.estimated <- total.equations * number.rs * forecast.horizon

  # > ESTIMATION ########################################################## ----

  # |__ Fit and Predict ========================================================
  ## Creates list where the first level is the forecast lead, and inside each
  ## level there is a table returning the forecast by equation/model's number
  ## (column) and rolling sample (row).
  prediction.gears <- lapply(
    X   = 1:forecast.horizon,
    function(h = X) {
      sapply(
        X = 1:total.equations,
        function(Eq = X) {
          sapply(
            X                = 1:number.rs,
            FUN              = gears:::fit_predict,
            forecast.lead    = h,
            eq.number        = Eq,
            ... # further arguments passed on to glm
          )
        }
      )
    }
  )

  # |__ Get Prediction Errors ==================================================

  ## \____ All measures --------------------------------------------------------
  ## Creates a list where the first level is the forecast lead, and inside each
  ## level there is a table returning the prediction errors by equation/model's
  ## number (row) and error measure (column)
  prediction.errors.all <- gears:::forecast_errors(
    DATA,
    forecast.horizon,
    Y_name,
    total.equations,
    number.rs,
    DF.Fit.Predict,
    forecasts.gears = prediction.gears
  )

  ## \____ User's measure-------------------------------------------------------
  ## Returns a table with the prediction errors by equation/model's number (row)
  ## and forecast lead (column).
  prediction.errors <- sapply(
    X = 1:forecast.horizon,
    function(X) prediction.errors.all[[X]][, error.measure]
  )
  # > BEST EQUATION ####################################################### ----

  # |__ Get Equations that Min. Selected Error Function ========================

  equations.min <- apply(prediction.errors, 2, which.min)

  # |__ Estimate Best Model ====================================================
  ## Returns a list of size equal to forecast.horizon. Each level returns the
  ## glm object associated with the forecast lead and its best equation.
  estimates.best <- lapply(
    X   = 1:forecast.horizon,
    function(h = X) {
      lapply(
        X                = 1:number.rs,
        forecast.lead    = h,
        eq.number        = equations.min[h],
        FUN              = gears:::fit_best,
        ...
      )
    }
  )

  # > Out-of-Sample Forecasts ############################################# ----

  if (betas.selection == "both") {

    if (error.measure == "list") {

      out.forecasts1 <- cbind.data.frame(
        "error.measure" = do.call(rbind, error.measure.list),
        "betas.selection" = "average",
        do.call(rbind, out_forecasts(betas = "average", estimates.best))
      )

      out.forecasts2 <- cbind.data.frame(
        "error.measure" = do.call(rbind, error.measure.list),
        "betas.selection" = "last",
        do.call(rbind, out_forecasts(betas = "last", estimates.best))
      )

      out.forecasts <- rbind(out.forecasts1, out.forecasts2)

    } else {

      out.forecasts1 <- cbind.data.frame(
        error.measure,
        "betas.selection" = "average",
        do.call(rbind, out_forecasts(betas = "average", estimates.best))
      )

      out.forecasts2 <- cbind.data.frame(
        error.measure,
        "betas.selection" = "last",
        do.call(rbind, out_forecasts(betas = "last", estimates.best))
      )

      out.forecasts <- rbind(out.forecasts1, out.forecasts2)

    }

  } else {
    if (error.measure == "list") {

      out.forecasts <- cbind.data.frame(
        "error.measure" = do.call(rbind, error.measure.list),
        "betas.selection" = betas.selection,
        do.call(rbind, out_forecasts(betas = betas.selection, estimates.best))
      )

    } else {

      out.forecasts <- cbind.data.frame(
        error.measure,
        "betas.selection" = betas.selection,
        do.call(rbind, out_forecasts(betas = betas.selection, estimates.best))
      )

    }
  }

  # > Results GEAR ######################################################## ----

  # |_ Best Fitted Models ======================================================

  fitted_best <- estimates.best

  # |_ Out-of-Sample Point Forecasts ===========================================

  out_sample_point_forecasts <- out.forecasts

  # |_ In-sample Predictions ===================================================

  ## \__ Point Estimates -------------------------------------------------------

  insample_point_estimates <- lapply(
    X = 1:length(prediction.errors),
    function(Er = X) {
      sapply(
        X = 1:forecast.horizon,
        function(H = X) {
          EQ <- equations.min[[Er]][H]
          prediction.gears[[H]][, EQ]
        }
      )
    }
  )

  ## \__ Errors ----------------------------------------------------------------
  # TODO: GET INSAMPLE ERRORS

  ## \__ Mean Error Measure ----------------------------------------------------

  tmp_insample_mean_error_measure <- lapply(
    X = 1:length(prediction.errors),
    function(Er = X) {
      sapply(
        X = 1:forecast.horizon,
        function(X) prediction.errors[[Er]][equations.min[[Er]][X], X]
      )
    }
  )

  if (error.measure == "list") {
    tmp_insample_mean_error_measure <- cbind.data.frame(
      do.call(rbind, error.measure.list),
      do.call(rbind, tmp_insample_mean_error_measure)
    )

  } else {
    tmp_insample_mean_error_measure <- cbind.data.frame(
      error.measure,
      do.call(rbind, tmp_insample_mean_error_measure)
    )
  }

  names(tmp_insample_mean_error_measure) <- c(
    "error.measure",
    paste0("mean.value.for.h.", 1:forecast.horizon)
  )

  best_equations <- lapply(
    X = 1:length(prediction.errors),
    function(Er = X) {
      paste0(
        paste0(Y_name, "_plus_", 1:forecast.horizon),
        all.equations.rhs[equations.min[[Er]]]
      )
    }
  )

  if (error.measure == "list") {
    best_equations <- cbind.data.frame(
      do.call(rbind, error.measure.list),
      do.call(rbind, best_equations)
    )

  } else {
    best_equations <- cbind.data.frame(
      error.measure,
      do.call(rbind, best_equations)
    )
  }

  names(best_equations) <- c(
    "error.measure",
    paste0("best.eq.for.h.", 1:forecast.horizon)
  )

  # > Return Results ###################################################### ----

  return(list(
    "fitted_best"                 = fitted_best,
    "out_sample_point_forecasts"  = out_sample_point_forecasts,
    "insample_point_estimates"    = insample_point_estimates,
    "insample_mean_error_measure" = tmp_insample_mean_error_measure,
    "best_equations"              = best_equations
  ))

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ----
} # UNCOMMENT THIS WHEN DONE WITH FUNCTION, TO CLOSE MAIN BRACKET
# END ----
