#' Create the Fit & Predict Data Frame
#'
#' This function creates a data frame with training data and test data for the
#' estimation stage.
#'
#' @param DATA A data frame or a univariate time series.
#' @param forecast.horizon Number of periods for forecasting.
#' @param size.rs Number of observations in the rolling sample.
#' @param number.rs Number of rolling samples.
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
#'
#' @return A data frame.
#' @keywords internal
#'
# @examples
# create_DF_Fit_Predict(
#   DATA                = gears::commodities_prices,
#   forecast.horizon    = 2,
#   size.rs             = 5,
#   number.rs           = 3,
#   y.name              = "PORK_PRICE" ,
#   y.max.lags          = 2,
#   x.names             = list("BEEF_PRICE") ,
#   x.max.lags          = list(3),
#   x.fixed.names       = list("BEEF_PRICE") ,
#   x.fixed.lags        = list(4),
#   x.interaction.names = list("CORN_PRICE*CORN_PRICE", "BEEF_PRICE*WHEAT_PRICE"),
#   x.interaction.lags  = list(c(1, 1), c(3,3)),
#   last.obs            = 100
# )
create_DF_Fit_Predict <- function(DATA,
                                  forecast.horizon,
                                  size.rs,
                                  number.rs,
                                  y.name = NULL,
                                  y.max.lags = NULL,
                                  x.names = NULL,
                                  x.max.lags = NULL,
                                  x.fixed.names = NULL,
                                  x.fixed.lags = NULL,
                                  x.interaction.names = NULL,
                                  x.interaction.lags = NULL,
                                  last.obs = NULL) {

  # > Helpers ############################################################# ----

  fcn.lag <- function(df, n.lags) c(rep(NA, n.lags), df)[1:length(df)]

  fcn.cut.DF <- function(h) {

    tmp.DF.Fit.h <- cbind(
      DF.Y.fit[[h]],
      DF.X_Lagged.fit[[h]],
      DF.X_Fixed.fit[[h]],
      DF.Inter.fit[[h]]
    )

    tmp.unique.vars <- unique(colnames(tmp.DF.Fit.h))

    tmp.DF.Fit.h <- tmp.DF.Fit.h[, tmp.unique.vars]

    tmp.DF.Fit.h <- subset(
      tmp.DF.Fit.h,
      select = !colnames(tmp.DF.Fit.h) %in% c("NA.")
    )

    tmp.DF.Fit.h <- cbind.data.frame(
      tmp_n_t_plus = 1:dim(tmp.DF.Fit.h)[1],
      tmp_n_t = c(rep(NA, h),1:(dim(tmp.DF.Fit.h)[1] - h)),
      tmp.DF.Fit.h
    )

    # Cut fit

    tmp.fcn.cut.fit <- function(rs) {

      max.number.rs <- last.obs - (h + size.rs)

      rs.start <- max.number.rs - (rs - 1)
      rs.end <- rs.start + (size.rs - 1)

      subset(tmp.DF.Fit.h, tmp_n_t %in% c(rs.start:rs.end))

      #tmp.DF.Fit.h[tmp.DF.Fit.h$tmp_n_t %in% c(rs.start:rs.end), ]
    }

    # Cut predict

    tmp.fcn.cut.predict <- function(rs) {

      max.number.rs <- last.obs - (2 * h) - (size.rs - 1)

      rs.start <- max.number.rs - (rs - 1) + h
      rs.end <- rs.start + (size.rs - 1)

      subset(tmp.DF.Fit.h, tmp_n_t == rs.end)
      #tmp.DF.Fit.h[tmp.DF.Fit.h$tmp_n_t == rs.end, ]
    }

    # tmp.DF.cut

    tmp.DF.cut <- data.frame(
        forecast.h = h,

        rolling.sample = number.rs:1,

        data_fit = I(lapply(X = number.rs:1, FUN = tmp.fcn.cut.fit)),

        data_predict = I(lapply(X = number.rs:1, FUN = tmp.fcn.cut.predict)),

        stringsAsFactors = FALSE
    )

    return(tmp.DF.cut)
  }

  # > DATA MANIPULATION ################################################### ----

  # |__ Y name =================================================================

  if (is.null(y.name)) {
    Y_name <- "Y_t"
  } else {
    Y_name <- paste0(y.name, "_t")
  }

  # |__ DF.mod =================================================================

  # \____ Class: ts ------------------------------------------------------------

  if (class(DATA) == "ts") {

    # Define Y variable name

    DF.mod.Y <- data.frame("Y_t" = as.numeric(DATA), stringsAsFactors = FALSE)
    colnames(DF.mod.Y) <- Y_name

    DF.mod.X.lag   <- NULL
    DF.mod.X.fixed <- NULL
    DF.mod.X.inter <- NULL

  } else {

    # \____ Class: data.frame---------------------------------------------------

    # Define Y variable name

    DF.mod.Y <- as.data.frame(DATA[, y.name])

    ## Add _t to variable name (e.g., Y becomes Y_t)

    Y_name <- paste0(y.name, "_t")

    ## Change Y variable in original dataset to Y_t

    colnames(DF.mod.Y) <- Y_name

    ## Add _minus_lag# to variable name (e.g., Y_t becomes Y_t_minus_1),
    ## and create a vector with all Y variables.

    if (is.null(y.max.lags)) {
      y.rhs.var.vec <- c(Y_name)
    } else {
      y.rhs.var.vec <- c(Y_name, paste0(Y_name, "_minus_", 1:y.max.lags))
    }

    # ........................................................................ #

    # Define X-lagged names

    if (is.null(x.max.lags)) {
      x.var.vec <- c(x.names)
      DF.mod.X.lag <- NULL
    } else {

      DF.mod.X.lag <- as.data.frame(DATA[, unlist(x.names)])

      ## Add _t to variable name (e.g., X becomes X_t)

      tmp_x_lag_names <- unlist(x.names)

      X_lagged_names_list <- lapply(
        X   = 1:length(x.names),
        FUN = function(X) paste0(tmp_x_lag_names[[X]], "_t")
      )

      ## Change X variables in original dataset to X_t

      X_lagged_names <- unlist(X_lagged_names_list)

      colnames(DF.mod.X.lag) <- X_lagged_names

      ## Add _minus_lag# to variable name (e.g., X_t becomes X_t_minus_1),
      ## and create a vector with all X variables.

      x.var.vec <- unlist(sapply(
        X = 1:length(x.names),
        #function(i = X) {
        function(X) {   # Using this instead of what's on the previous line to
          i <- X        # address the issue: "no visible binding for global
                        # variable" generated by "R CMD check"
          sapply(
            X = 1:x.max.lags[[i]],
            function(X) {
              paste0(x.names[[i]], "_t_minus_", X)
            }
          )
        }
      ))

    }

    # ........................................................................ #

    # Define X-fixed names
    if (is.null(x.fixed.lags)) {

      x.fixed.var.vec <- c(x.fixed.names)
      DF.mod.X.fixed <- NULL

    } else {

      DF.mod.X.fixed <- as.data.frame(DATA[, unlist(x.fixed.names)])

      x.fixed.var.vec <- paste0(
        unlist(x.fixed.names), "_t_minus_", unlist(x.fixed.lags)
      )

      colnames(DF.mod.X.fixed) <- x.fixed.var.vec

    }

    # ........................................................................ #

    # Define interactions names

    if (is.null(x.interaction.names)) {
      x.interaction.var.vec <- c(x.interaction.names)
      DF.mod.X.inter <- NULL
    } else {

      ## Add _t to variable name (e.g., X becomes X_t)

      tmp_Inter_names <- strsplit(x = unlist(x.interaction.names), split ="\\*")

      tmp_unique_inter_names <- unique(unlist(tmp_Inter_names))

      DF.mod.X.inter <- as.data.frame(DATA[, unlist(tmp_unique_inter_names)])

      Inter_names_list <- lapply(
        X   = 1:length(x.interaction.names),
        FUN = function(X) paste0(tmp_Inter_names[[X]], "_t")
      )

      ## Change X variables in original dataset to X_t

      Inter_names <- unlist(Inter_names_list)

      unique_Inter_names <- unique(Inter_names)

      colnames(DF.mod.X.inter) <- unique_Inter_names

      ## Add _minus_lag# to variable name (e.g., X_t becomes X_t_minus_1),
      ## and create a vector with all X1*X2 variables.

      x.interaction.vars <- sapply(
        X = 1:length(x.interaction.names),
        function(X) {
          paste(
            paste0(tmp_Inter_names[[X]], "_t_minus_", x.interaction.lags[[X]]),
            collapse = "*"
          )
        }
      )

      ## Create a vector with all X1, X2 variables.

      x.interaction.var.vec <- unlist(strsplit(
        x = x.interaction.vars,
        split = "\\*"
      ))

    }

  }

  # |__ DF.Y.fit ===============================================================

  DF.Y.fit <- lapply(
    X = 1:forecast.horizon,
    #function(H = X) {
    function(X) {   # Using this instead of what's on the previous line to
      H <- X        # address the issue: "no visible binding for global
                    # variable" generated by "R CMD check"

      # Add Lags Due to Forecast Horizon

      tmp.df.lagged <- fcn.lag(DF.mod.Y[, Y_name], n.lags = H)

      tmp.h <- cbind(DF.mod.Y[, Y_name], tmp.df.lagged)

      colnames(tmp.h) <- c(paste0(Y_name, "_plus_", H), Y_name)

      # Add Lags Due to y.max.lags

      if (!is.null(y.max.lags)) {

        tmp.l <- sapply(X = 1:y.max.lags, FUN = fcn.lag, df = tmp.df.lagged)

        colnames(tmp.l) <- paste0(Y_name, "_minus_", 1:y.max.lags)

        return(cbind(tmp.h, tmp.l))

      } else {
        return(tmp.h)
      }

    }
  )

  # |__ DF.X_Lagged.fit ========================================================

  if (is.null(DF.mod.X.lag)) {
    DF.X_Lagged.fit <- lapply(1:forecast.horizon, function(X) data.frame(NA))
  } else {

    DF.X_Lagged.fit <- lapply(
      X = 1:forecast.horizon,
      #function(H = X) {
      function(X) {   # Using this instead of what's on the previous line to
        H <- X        # address the issue: "no visible binding for global
                      # variable" generated by "R CMD check"

        # Add Lags Due to Forecast Horizon

        if (!is.null(dim(DF.mod.X.lag[, X_lagged_names])[2])) {

          tmp.h <- apply(
            X      = DF.mod.X.lag[, X_lagged_names],
            MARGIN = 2,
            FUN    = fcn.lag,
            n.lags = H
          )

        } else {

          tmp.h <- as.data.frame(
            fcn.lag(DF.mod.X.lag[, X_lagged_names], n.lags = H)
          )

          colnames(tmp.h) <- X_lagged_names

        }

        # Add Lags Due to x.max.lags

        tmp.l <- do.call(
          cbind,
          lapply(
            X = 1:length(X_lagged_names),
            function(j = X) {
              do.call(
                cbind,
                sapply(
                  X = 1:unlist(x.max.lags)[j],
                  #function(L = X) {
                  function(X) {

                    L <- X

                    tmp.lags <- as.data.frame(fcn.lag(tmp.h[, j], n.lags = L))

                    colnames(tmp.lags) <- paste0(
                      colnames(tmp.h)[j], "_minus_", L
                    )

                    return(tmp.lags)

                  }
                )
              )
            }
          )
        )

        return(cbind(tmp.h, tmp.l))

      }
    )

  }

  # |__ X-Fixed ================================================================
  if (is.null(DF.mod.X.fixed)) {
    DF.X_Fixed.fit <- lapply(1:forecast.horizon, function(X) data.frame(NA))
  } else {
    DF.X_Fixed.fit <- lapply(
      X = 1:forecast.horizon,
      #function(H = X) {
      function(X) {   # Using this instead of what's on the previous line to
        H <- X        # address the issue: "no visible binding for global
                      # variable" generated by "R CMD check"

        # Add Lags Due to Forecast Horizon

        if (!is.null(dim(DF.mod.X.fixed)[2])) {

          tmp.h <- apply(
            X      = DF.mod.X.fixed,
            MARGIN = 2,
            FUN    = fcn.lag,
            n.lags = H
          )

        } else {

          tmp.h <- as.data.frame(
            fcn.lag(DF.mod.X.lag, n.lags = H)
          )

          colnames(tmp.h) <- x.fixed.var.vec

        }

        # Add Lags Due to x.fixed.lags

        tmp.l <- do.call(
          cbind,
          lapply(
            X = 1:length(x.fixed.var.vec),
            function(X) {
              j <- X
              do.call(
                cbind,
                sapply(
                  X = unlist(x.fixed.lags)[j],
                  function(X) {
                    L <- X

                    tmp.lags <- as.data.frame(fcn.lag(tmp.h[, j], n.lags = L))

                    colnames(tmp.lags) <- paste0(colnames(tmp.h)[j])

                    return(tmp.lags)

                  }
                )
              )
            }
          )
        )

        return(tmp.l)

      }
    )
  }

  # |__ DF.Inter.fit ===========================================================

  if (is.null(DF.mod.X.inter)) {
    DF.Inter.fit <- lapply(1:forecast.horizon, function(X) data.frame(NA))
  } else {

    DF.Inter.fit <- lapply(
      X = 1:forecast.horizon,
      function(X) {
        H <- X

        # Add Lags Due to Forecast Horizon

        tmp.h <- apply(
          X      = DF.mod.X.inter[, Inter_names],
          MARGIN = 2,
          FUN    = fcn.lag,
          n.lags = H
        )

        # Add Lags Due to x.interaction.lags

        tmp.l <- do.call(
          cbind,
          lapply(
            X = 1:length(Inter_names),
            function(X) {
              j <- X
              do.call(
                cbind,
                sapply(
                  X = 1:unlist(x.interaction.lags)[j],
                  function(X) {

                    L <- X

                    tmp.lags <- as.data.frame(fcn.lag(tmp.h[, j], n.lags = L))

                    colnames(tmp.lags) <- paste0(
                      colnames(tmp.h)[j], "_minus_", L
                    )

                    return(tmp.lags)

                  }
                )
              )
            }
          )
        )

        # Remove duplicate columns (from cases like X1*X1)

        tmp.l.unique <- as.data.frame(tmp.l[, unique(x.interaction.var.vec)])
        colnames(tmp.l.unique) <- unique(x.interaction.var.vec)

        return(tmp.l.unique)

      }
    )

  }

  # |__ DF.Fit.Predict =========================================================

  DF.Fit.Predict <- lapply(X = 1:forecast.horizon, FUN = fcn.cut.DF)

  # > Return Results ###################################################### ----

  return(DF.Fit.Predict)

}
