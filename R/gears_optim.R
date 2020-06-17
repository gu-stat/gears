#' Optimization Algorithm for the GEARS Method for Time Series Forecasting
#'
#' Finds the sample size and number of samples that minimizes the in-sample
#' (ex-ante) forecasting error.
#'
#' @param DATA A data frame or a univariate time series.
#' @param forecast.horizon Number of periods for forecasting.
#' @param search.size.rs A vector or a numeric value that indicates the
#'     possible rolling sample sizes that the algorithm will search through.
#' @param search.number.rs A vector or a numeric value that indicates the
#'     possible number of rolling sample that the algorithm will search through.
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
#' @param parallel.use Boolean. Whether or not parallel computing should be
#'     used. \code{gears_optim} uses the
#'     \url{https://github.com/HenrikBengtsson/future.apply}{future.apply}
#'     package.
#' @param DATA.insample Training data set. Either as a data frame object, or
#'     a time series (ts) object.
#' @param DATA.outsample Testing data set. Either as a data frame object, or
#'     a time series (ts) object.
#' @param ... Further arguments passed to \link[stats]{glm}.
#'
#' @return A list containing the following
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
#' @export
#'
gears_optim <- function(DATA,
                        forecast.horizon,
                        search.size.rs,
                        search.number.rs,
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
                        parallel.use = FALSE,
                        DATA.insample,
                        DATA.outsample,
                        ...) { # ... TO ACCOUNT FOR OTHER OPTIONS PASSED TO glm

  # > Helpers ############################################################# ---

  # |__ Check future.apply  ----------------------------------------------------

  if (isTRUE(parallel.use)) {

    if (isFALSE("package:future.apply" %in% search())) {
      stop(paste0(
        "Package 'future.apply' not loaded. This package is used for parallel ",
        "computing. Load it or use 'parallel.use = FALSE'."
      ))
    }
  }

  # |__ Evaluation -------------------------------------------------------------

  fcn_eval <- function(X, out_forecasts) {

    tmp_out_forecasts <- as.numeric(out_forecasts[X, -c(1,2)])

    tmp_mad <- error_functions(
      error.measure = "mad",
      forecasts     = tmp_out_forecasts,
      outsample     = DATA.outsample
    )

    tmp_mse <- error_functions(
      error.measure = "mse",
      forecasts     = tmp_out_forecasts,
      outsample     = DATA.outsample
    )
    # TODO: VER PROBLEMAS AQUI APOS MUDANCAS EM fcn_OWA
    tmp_owa <- error_functions(
      error.measure      = "owa",
      forecasts          = tmp_out_forecasts,
      outsample          = DATA.outsample,
      insample           = DATA.insample,
      forecast.horizon   = forecast.horizon,
      alpha.level        = 0.05
    )
    # tmp_owa <- fcn_OWA(
    #   forecasts.values = tmp_out_forecasts,
    #   insample         = DATA.insample,
    #   outsample        = DATA.outsample,
    #   forecast.horizon = forecast.horizon
    # )

    tmp_results <- c(
      "Out_MAD"          = tmp_mad,
      "Out_MSE"          = tmp_mse,
      "OUT_SMAPE_N2"     = as.numeric(tmp_owa["sMAPE_NAIVE2"]),
      "OUT_MASE_N2"      = as.numeric(tmp_owa["MASE_NAIVE2"]),
      "Out_SMAPE"        = as.numeric(tmp_owa["SMAPE"]),
      "Out_MASE"         = as.numeric(tmp_owa["MASE"]),
      "Out_OWA"          = as.numeric(tmp_owa["OWA"])
    )

    return(tmp_results)

  }

  # |__ Grid Results -----------------------------------------------------------

  fcn_grid_results <- function(TMP.SIZE.RS, TMP.NUMBER.RS,...) {

    tmp_gears <- gears(
      DATA,
      forecast.horizon,
      size.rs             = TMP.SIZE.RS,
      number.rs           = TMP.NUMBER.RS,
      glm.family,
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
      error.measure.list,
      betas.selection
    )

    tmp_out_gears <- tmp_gears$out_sample_point_forecast

    tmp_eval_gears <- cbind(
      "size.rs"         = TMP.SIZE.RS,
      "number.rs"       = TMP.NUMBER.RS,
      tmp_out_gears,
      do.call(rbind, lapply(
        X = 1:dim(tmp_out_gears)[1],
        fcn_eval,
        out_forecasts = tmp_out_gears
      ))
    )

    return(tmp_eval_gears)
  }

  # |__ Get Minimum ------------------------------------------------------------

  # > GET COMBINATIONS #################################################### ----
  all_combs <- expand.grid(
    'size.rs'   = search.size.rs,
    'number.rs' = search.number.rs
  )

  # > RUN GEARS ########################################################### ----

  if (parallel.use == TRUE) {
    # TODO: VER EFEITOS DE USAR 'multisession' E NAO multisession (SEM ASPAS)
    future.apply::plan('multisession')

    results_gears <- future.apply::future_lapply(
      X = seq_along(search.size.rs),
      function(SIZE.RS = X) {
        do.call(rbind, lapply(
          X = search.number.rs,
          function(X) fcn_grid_results(TMP.NUMBER.RS = X, TMP.SIZE.RS = SIZE.RS)
        ))
      }
    )

  } else {

    results_gears <- lapply(
      X = search.number.rs,
      function(NUMBER.RS = X) {
        do.call(rbind, lapply(
          X = search.size.rs,
          function(X) fcn_grid_results(
            TMP.NUMBER.RS = NUMBER.RS,
            TMP.SIZE.RS   = X
          )
        ))
      }
    )

  }

  results_gears_all <- do.call(rbind, results_gears)

  min_results_gears <- apply(results_gears_all[,-c(1:4)], 2, which.min)

  results_N2_smape    <-   results_gears_all[min_results_gears['OUT_SMAPE_N2'], ]
  results_N2_mase     <-   results_gears_all[min_results_gears['OUT_MASE_N2'], ]
  results_gears_mad   <-   results_gears_all[min_results_gears['Out_MAD'], ]
  results_gears_mse   <-   results_gears_all[min_results_gears['Out_MSE'], ]
  results_gears_smape <-   results_gears_all[min_results_gears['Out_SMAPE'], ]
  results_gears_mase  <-   results_gears_all[min_results_gears['Out_MASE'], ]
  results_gears_owa   <-   results_gears_all[min_results_gears['Out_OWA'], ]

  return(list(
    'results_N2_smape'    = results_N2_smape,
    'results_N2_mase'     = results_N2_mase,
    'results_gears_mad'   = results_gears_mad,
    'results_gears_mse'   = results_gears_mse,
    'results_gears_smape' = results_gears_smape,
    'results_gears_mase'  = results_gears_mase,
    'results_gears_owa'   = results_gears_owa
  ))

  # return(list(results_gears))

}
