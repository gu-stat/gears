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
#'     prediction errors. Use \code{error.measure == "list"} if more than one
#'     measure should be used. Then, add the list with measures to the
#'     \code{error.measure.list} argument.
#' @param error.measure.list List with the names of the error measures to be
#'     used when calculating the in-sample prediction errors. Only use this if
#'     \code{error.measure == "list"}.
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
#'   error.measure       = "list",
#'   error.measure.list  = list("mse", "mad"),
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
                  error.measure = c("mse", "mae", "mase", "smape", "list"),
                  error.measure.list = NULL,
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

  # Fit-Predict ================================================================

  fitPredict <- function(forecast.lead, eq.number, sample.number, ...) {

    tmp_lhs      <- paste0(Y_name, "_plus_", forecast.lead)
    tmp_equation <- paste0(tmp_lhs, all.equations.rhs[eq.number])

    tmp_fit <- glm(
      formula = as.formula(tmp_equation),
      family  = glm.family,
      data    = DF.Fit.Predict[[forecast.lead]]$data_fit[[sample.number]],
      ...
    )

    if (length(tmp_fit$coefficients) > tmp_fit$rank) {
      tmp_forecast <- NA
    } else {
      tmp_forecast <- predict(
        object  = tmp_fit,
        newdata = DF.Fit.Predict[[forecast.lead]]$data_predict[[sample.number]],
        type    = "response"
      )
    }

    return(tmp_forecast)
  }

  # |__ Estimate the Best Equations ============================================

  fitBest <- function(forecast.lead, eq.number, sample.number, ...) {

    tmp_lhs      <- paste0(Y_name, "_plus_", forecast.lead)
    tmp_equation <- paste0(tmp_lhs, all.equations.rhs[eq.number[forecast.lead]])

    tmp_fit <- glm(
      formula = as.formula(tmp_equation),
      family  = glm.family,
      data    = DF.Fit.Predict[[forecast.lead]]$data_fit[[sample.number]],
      ...
    )

    return(tmp_fit)

  }

  # |__ Get Prediction Errors ==================================================

  predictionError <- function(error.measure, ...) {

    if (error.measure == "mase") {

      prediction.errors <- sapply(
        X = 1:forecast.horizon,
        function(H = X) {

          tmp_lhs <- paste0(Y_name, "_plus_", H)

          tmp <- sapply(
            X = 1:total.equations,
            function(EQ = X, ERROR.MEASURE = error.measure){
              sapply(
                X = 1:number.rs,
                function(RS = X) {

                  tmp_in <- DF.Fit.Predict[[H]]$data_fit[[RS]][,tmp_lhs]

                  tmp_out <- DF.Fit.Predict[[H]]$data_predict[[RS]][,tmp_lhs]

                  tmp_forecasts <- prediction.gears[[H]][RS, EQ]

                  error_functions(
                    error.measure = ERROR.MEASURE,
                    forecasts     = tmp_forecasts,
                    outsample     = tmp_out,
                    insample      = tmp_in
                  )
                }
              )
            }
          )

          tmp_avg_by_eq <- apply(tmp, 2, mean)

          tmp_avg_by_eq
        }
      )

    } else if (error.measure %in% c("mse", "mad", "smape")) {

      prediction.errors <- sapply(
        X = 1:forecast.horizon,
        function(H = X) {
          tmp <- lapply(
            X = 1:total.equations,
            function(EQ = X, ERROR.MEASURE = error.measure){
              tmp_lhs <- paste0(Y_name, "_plus_", H)

              tmp_forecasts <- prediction.gears[[H]][, EQ]

              tmp_out <- do.call(rbind, DF.Fit.Predict[[H]]$data_predict)[,tmp_lhs]

              tmp_in <- NULL

              error_functions(
                error.measure = ERROR.MEASURE,
                forecasts     = tmp_forecasts,
                outsample     = tmp_out,
                insample      = tmp_in
              )

            }
          )

          tmp

        }
      )

    } else if (error.measure == "owa") {

      prediction.errors <- sapply(
        X = 1:forecast.horizon,
        function(H = X) {

          tmp_lhs <- paste0(Y_name, "_plus_", H)

          sapply(
            X = 1:total.equations,
            function(EQ = X){
              TMP <- lapply(
                X = 1:number.rs,
                function(RS = X) {

                  tmp_in <- DF.Fit.Predict[[H]]$data_fit[[RS]][,tmp_lhs]

                  tmp_out <- DF.Fit.Predict[[H]]$data_predict[[RS]][,tmp_lhs]

                  tmp_forecasts <- prediction.gears[[H]][RS, EQ]

                  fcn_owa(
                    forecasts.values = tmp_forecasts,
                    insample         = tmp_in,
                    outsample        = tmp_out,
                    forecast.horizon = H
                  )
                }
              )

              TMP2 <- do.call(rbind, TMP)[, 1:4]

              TMP3 <- apply(TMP2, 2, mean)

              mean_mase_gears  <- TMP3["MASE"]
              mean_mase_naive2 <- TMP3["MASE_NAIVE2"]

              mean_smape_gears  <- TMP3["SMAPE"]
              mean_smape_naive2 <- TMP3["sMAPE_NAIVE2"]

              # FINAL OWA
              tmp_owa <- ((mean_mase_gears/mean_mase_naive2) +
                            (mean_smape_gears/mean_smape_naive2))/2
              names(tmp_owa) <- NULL

              tmp_owa
            }
          )
        }
      )
    }

    return(prediction.errors)

  }

  # |__ Get Prediction Errors ==================================================

  outForecasts <- function(betas, ...) {

    # |__ Betas Selection ------------------------------------------------------

    if (betas == "last") {

      fcn.betas <- function(X, errorMeasure, ...) {

        tmp_fit <- estimates.best[[errorMeasure]][[X]][[number.rs]]

        tmp_forecast <- predict(
          object  = tmp_fit,
          newdata = DF.Forecast[last.obs, ],
          type    = "response"
        )

        names(tmp_forecast) <- sum(last.obs + X)

        return(tmp_forecast)

      }

    } else if (betas == "average") {

      fcn.betas <- function(X, errorMeasure, ...) {

        tmp_h <- X

        tmp_coeff <- lapply(
          X = 1:number.rs,
          function(X) estimates.best[[errorMeasure]][[tmp_h]][[X]]$coefficients
        )

        tmp_coeff    <- do.call(rbind, tmp_coeff)
        avg_coeff    <- apply(tmp_coeff, 2, mean)

        if ("(Intercept)" %in% names(avg_coeff)) {

          tmp_names <- names(avg_coeff)[which(names(avg_coeff)!="(Intercept)")]

          tmp_forecast <- sum(
            avg_coeff[tmp_names] * DF.Forecast[last.obs, tmp_names]
          )

          tmp_forecast <- sum(tmp_forecast, avg_coeff["(Intercept)"])

        } else {
          tmp_forecast <- sum(
            avg_coeff * DF.Forecast[last.obs, names(avg_coeff)]
          )
        }

        names(tmp_forecast) <- sum(last.obs + tmp_h)

        return(tmp_forecast)

      }

    }

    out.forecasts <- lapply(
      X = 1:length(prediction.errors),
      function(Er = X) {
        sapply(
          X = 1:forecast.horizon,
          errorMeasure = Er,
          FUN = fcn.betas)
      }
    )

    return(out.forecasts)

  }

  # > DATA ################################################################ ----

  # |__ DF.Fit.Predict =========================================================

  DF.Fit.Predict <- create_DF_Fit_Predict(
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

  DF.Forecast <- create_DF_Forecast(
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

  all.equations.rhs <- all_models_rhs(
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

  prediction.gears <- lapply(
    X   = 1:forecast.horizon,
    function(h = X) {
      sapply(
        X = 1:total.equations,
        function(Eq = X) {
          sapply(
            X                = 1:number.rs,
            FUN              = fitPredict,
            forecast.lead    = h,
            eq.number        = Eq,
            ...
          )
        }
      )
    }
  )

  # |__ Get Prediction Errors ==================================================

  if (error.measure == "list") {

    if (is.null(error.measure.list)) {
      warning(paste0(
        "If error.measure = 'list', you must provide a list of measures with",
        " error.measure.list = list()"
      ))
    }

    prediction.errors <- lapply(
      error.measure.list, predictionError, forecast.horizon, ...
    )

  } else {
    prediction.errors <- list(predictionError(
      error.measure, forecast.horizon, ...
    ))
  }

  # > BEST EQUATION ####################################################### ----

  # |__ Get Equations that Min. Selected Error Function ========================

  equations.min <- lapply(
    X = 1:length(prediction.errors),
    function(X) apply(prediction.errors[[X]], 2, which.min)
  )

  # |__ Estimate Best Model ====================================================

  estimates.best <- lapply(
    X = 1:length(prediction.errors),
    function(Er = X) {
      lapply(
        X   = 1:forecast.horizon,
        function(h = X) {
          lapply(
            X                = 1:number.rs,
            forecast.lead    = h,
            eq.number        = equations.min[[Er]],
            FUN              = fitBest,
            ...
          )
        }
      )
    }
  )

  # > Out-of-Sample Forecasts ############################################# ----

  if (betas.selection == "both") {

    if (error.measure == "list") {

      out.forecasts1 <- cbind.data.frame(
        "error.measure" = do.call(rbind, error.measure.list),
        "betas.selection" = "average",
        do.call(rbind, outForecasts(betas = "average", estimates.best))
      )

      out.forecasts2 <- cbind.data.frame(
        "error.measure" = do.call(rbind, error.measure.list),
        "betas.selection" = "last",
        do.call(rbind, outForecasts(betas = "last", estimates.best))
      )

      out.forecasts <- rbind(out.forecasts1, out.forecasts2)

    } else {

      out.forecasts1 <- cbind.data.frame(
        error.measure,
        "betas.selection" = "average",
        do.call(rbind, outForecasts(betas = "average", estimates.best))
      )

      out.forecasts2 <- cbind.data.frame(
        error.measure,
        "betas.selection" = "last",
        do.call(rbind, outForecasts(betas = "last", estimates.best))
      )

      out.forecasts <- rbind(out.forecasts1, out.forecasts2)

    }

  } else {
    if (error.measure == "list") {

      out.forecasts <- cbind.data.frame(
        "error.measure" = do.call(rbind, error.measure.list),
        "betas.selection" = betas.selection,
        do.call(rbind, outForecasts(betas = betas.selection, estimates.best))
      )

    } else {

      out.forecasts <- cbind.data.frame(
        error.measure,
        "betas.selection" = betas.selection,
        do.call(rbind, outForecasts(betas = betas.selection, estimates.best))
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
