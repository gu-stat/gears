#' Checks the arguments for the GEARS function.
#'
#' This function checks the arguments and stops the program if encounters an
#' error. An error message is returned to the user.
#'
#' @param DATA A data frame or a univariate time series.
#' @param forecast.horizon Number of periods for forecasting.
#' @param size.rs Number of observations in the rolling sample.
#' @param number.rs Number of rolling samples.
#' @param glm.family A description of the error distribution to be used in the
#'     model. See \link[stats]{glm} for details.
#' @param level Confidence level for prediction intervals. Numeric value
#'     between 0 and 100.
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
#' @param use.intercept Possible values \code{c("both", "with", "without")}.
#' @param error.measure Possible values
#'     \code{c("mse", "mae", "mase", "smape", "owa")}.
#' @param betas.selection Possible values \code{c("last", "average", "both")}.
#'
#' @return Stops the program and returns an error message.
#' @keywords internal
#'
checks <- function(DATA, forecast.horizon, size.rs, number.rs, glm.family,
                   level,
                   y.name = NULL,
                   y.max.lags = NULL,
                   x.names = NULL,
                   x.max.lags = NULL,
                   x.fixed.names = NULL,
                   x.fixed.lags = NULL,
                   x.interaction.names = NULL,
                   x.interaction.lags = NULL,
                   last.obs, use.intercept, error.measure, betas.selection){

  # |__ GLM Family =============================================================

  if(!missing(glm.family) & length(glm.family)>1) {
    stop("Only one 'glm.family' allowed.")
  }

  # |__ Last Obs ===============================================================

  if (!is.null(last.obs)) {

    if (stats::is.ts(DATA)) {
      checkLastObs <- length(DATA)
    } else {
      checkLastObs <- dim(DATA)[1]
    }

    if (last.obs > checkLastObs) {
      stop(paste0(
        "The value for 'last.obs' is greater than the size/length of your data."
      ))
    }

  }

  # |__ Intercept ==============================================================

  if(!missing(use.intercept) & length(use.intercept)>1) {
    stop("Only one 'use.intercept' allowed.")
  }

  # |__ Error Measure ==========================================================

  if(!missing(error.measure) & length(error.measure)>1) {
    stop("Only one 'error.measure' allowed.")
  }

  # |__ Betas Selection ========================================================

  if(!missing(betas.selection) & length(betas.selection)>1) {
    stop("Only one 'error.measure' allowed.")
  }

  # |__ DATA ===================================================================

  if (missing(DATA)) {
    stop(paste0("Missing DATA input. Please, check if you defined DATA = "))
  }

  if (!(class(DATA)[1] %in% c("data.frame", "ts"))) {
    stop(paste0(
      "DATA not of valid type. Please, transform your data to one of",
      "these formats: data.frame, ts"
    ))
  }

  # |__ Forecast Horizon =======================================================

  if (missing(forecast.horizon)) stop(paste0("Missing forecast.horizon value."))

  # |__ Size Random Sample =====================================================

  if (missing(size.rs)) {
    stop(paste0(
      "Missing size.rs input. Please, indicate the size of your rolling sample."
    ))
  }

  # |__ Number Random Samples ==================================================

  if (missing(number.rs)) {
    stop(paste0(
      "Missing number.rs input. Please, indicate the number of rolling samples."
    ))
  }

  # |__ Confidence Level =======================================================

  if (level > 100 | level <= 0) {
    stop(paste0(
      "Incorrect value for the confidence level. Please, provide a value for",
      "'level' that is between 0 and 100."
    ))
  }

  # |__ Number of Lags =========================================================

  if (is.null(y.max.lags) & is.null(x.max.lags) & is.null(x.fixed.lags)) {
    stop(paste0(
      "You need to define the number of lags to be included in your model. ",
      "Add a numeric value to 'y.max.lags' or add a list of values to either ",
      "'x.max.lags' or 'x.fixed.lags'.\n",
      "For example, if you want your model to include up to 2 past values of ",
      "y but not of x, set y.max.lags = 2 and do not include 'x.max.lags' nor ",
      "'x.fixed.lags'. Or if you want your model to consider only lags 3 and ",
      "4 of Y, set x.fixed.lags = list(3,4) and do not include y.max.lags nor ",
      "'x.max.lags'."
    ))
  }

  # |__ X vars types ===========================================================

  if (!is.null(x.names) & is.list(x.names) == FALSE) {
    stop("The argument 'x.names' should be a list.")
  }

  if (!is.null(x.max.lags) & is.list(x.max.lags) == FALSE) {
    stop("The argument 'x.max.lags' should be a list.")
  }

  if (!is.null(x.fixed.names) & is.list(x.fixed.names) == FALSE) {
    stop("The argument 'x.fixed.names' should be a list.")
  }

  if (!is.null(x.fixed.lags) & is.list(x.fixed.lags) == FALSE) {
    stop("The argument 'x.fixed.lags' should be a list.")
  }

  if (!is.null(x.interaction.names) & is.list(x.interaction.names) == FALSE) {
    stop("The argument 'x.interaction.names' should be a list.")
  }

  if (!is.null(x.interaction.lags) & is.list(x.interaction.lags) == FALSE) {
    stop("The argument 'x.interaction.lags' should be a list.")
  }

  # |__ Data Frame Inputs ======================================================

  if (class(DATA)[1] == "data.frame") {

    ## \______ Y variable ----
    if (is.null(y.name)) {
      stop(paste0(
        "Missing name of your y (dependent) variable. Please, ",
        "supply y.name and make sure it matches the one in DATA."
      ))
    }

    if (!(y.name %in% names(DATA))) {
      stop(paste0(
        "Name of Y variable (y.name) not found among the columns ",
        "of DATA. Please, check the spelling of y.name variable."
      ))
    }

    if (is.null(y.max.lags) & !is.null(x.max.lags)) {
      stop(paste0(
        "Missing name of your y (dependent) variable. Please, ",
        "supply 'y.name' and make sure it matches the one in DATA."
      ))
    }

    if (!(y.name %in% names(DATA))) {
      stop(paste0(
        "Name of Y variable (y.name) not found among the columns ",
        "of DATA. Please, check the spelling of y.name variable."
      ))
    }

    ## \______ X variables ----

    if (!is.null(x.names)) {
      if (sum(!(x.names %in% names(DATA))) >0) {
        stop(paste0(
          "Name of at least one X variable ('x.names') not found ",
          "among the columns of DATA. Please, check the spelling of ",
          "'x.names' variable."
        ))
      }
    }

    if (!is.null(x.names) & is.null(x.max.lags)) {
      stop(paste0(
        "You have defined a value for 'x.names' but 'x.max.lags' is ",
        "missing. Please, set 'x.max.lags' to be list with a number greater ",
        "than zero or remove the value for 'x.names'."
      ))
    }

    if (is.null(x.names) & !is.null(x.max.lags)) {
      stop(paste0(
        "You have defined a value for x.max.lags but 'x.names' is ",
        "missing. Please, indicate the name(s) of your x ",
        "variable(s) or remove the value for 'x.max.lags'."))
    }

    ## \____ X Variables - Fixed ----

    if (!is.null(x.fixed.names)) {
      if (sum(!(x.fixed.names %in% names(DATA))) >0) {
        stop(paste0(
          "Name of at least one X variable ('x.fixed.names') not found ",
          "among the columns of DATA. Please, check the spelling of ",
          "'x.fixed.names' variable."
        ))
      }
    }

    if (!is.null(x.fixed.names) & is.null(x.fixed.lags)) {
      stop(paste0(
        "You have defined a value for 'x.fixed.names' but ",
        "with 'x.fixed.lags' is missing. Please, set 'x.fixed.lags' to be a ",
        "list with a number greater than zero or remove the value for ",
        "'x.fixed.names'."))
    }

    if (is.null(x.fixed.names) & !is.null(x.fixed.lags)) {
      stop(paste0(
        "You have defined a value for 'x.fixed.lags' but 'x.fixed.names' ",
        "is missing. Please, indicate the name(s) of your x ",
        "variable(s) with fixed lags or remove the value for ",
        "'x.fixed.lags'."))
    }

    ## \____ X Variables - Interactions ----

    if (!is.null(x.interaction.names)) {

      tmp_Inter_names <- strsplit(x = unlist(x.interaction.names), split ="\\*")

      tmp_unique_inter_names <- unique(unlist(tmp_Inter_names))

      if (sum(!(tmp_unique_inter_names %in% names(DATA))) >0) {
        stop(paste0(
          "Name of at least one X variable ('x.interaction.names') not found ",
          "among the columns of DATA. Please, check the spelling of ",
          "'x.interaction.names' variable."
        ))
      }
    }

    if (!is.null(x.interaction.names) & is.null(x.interaction.lags)) {
      stop(paste0(
        "You have defined a value for x.interaction.names but ",
        "x.interaction.lags is missing. ",
        "Please, set 'x.interaction.lags' to be a list with numbers greater ",
        "than zero or remove the value for x.interaction.names."
      ))
    }

    if (!is.null(x.interaction.lags) & is.null(x.interaction.names)) {
      stop(paste0(
        "You have defined a value for x.interaction.lags but ",
        "x.interaction.names is missing. ",
        "Please, indicate the name(s) of your interaction variables ",
        "using x.interaction.names or remove the value for ",
        "x.interaction.lags"
      ))
    }
  }

  # |__ Minimum Size of Rolling Sample =========================================

  ## Number of Y parameters (Y_{t-0},...,Y_{t-y.max.lags})

  # if (!is.null(y.max.lags) | y.max.lags > 0) {
  #   parameters.y <- sum(y.max.lags, 1)
  # } else {
  #   parameters.y <- 0
  # }
  if (!is.null(y.max.lags)) {
    parameters.y <- sum(y.max.lags, 1)
  } else {
    parameters.y <- 0
  }

  ## Number of X parameters (X_{t-0},...,X_{t-x.max.lags})
  # if (!is.null(x.max.lags) | sum(unlist(x.max.lags)) > 0) {
  #   parameters.x <- sum(unlist(x.max.lags), 1)
  # } else {
  #   parameters.x <- 0
  # }
  if (!is.null(x.max.lags)) {
    parameters.x <- sum(unlist(x.max.lags), 1)
  } else {
    parameters.x <- 0
  }

  ## Number of FIXED X parameters (X_{t-x.fixed.lags})
  # if (!is.null(x.fixed.lags) | sum(unlist(x.fixed.lags)) > 0) {
  #   parameters.x.fixed <- sum(unlist(x.fixed.lags))
  # } else {
  #   parameters.x.fixed <- 0
  # }
  if (!is.null(x.fixed.lags)) {
    parameters.x.fixed <- sum(unlist(x.fixed.lags))
  } else {
    parameters.x.fixed <- 0
  }

  ## Number of INTERACTED X parameters (X_{t-x.interaction.lags))
  # if (!is.null(x.interaction.lags) | sum(unlist(x.interaction.lags)) > 0) {
  #   parameters.x.interaction <- sum(unlist(x.interaction.lags))
  # } else {
  #   parameters.x.interaction <- 0
  # }
  if (!is.null(x.interaction.lags)) {
    parameters.x.interaction <- sum(unlist(x.interaction.lags))
  } else {
    parameters.x.interaction <- 0
  }

  ## Total Number of parameters with no Intercept
  total.parameters <- sum(
    parameters.y, parameters.x, parameters.x.fixed, parameters.x.interaction
  )

  ## Add Intercept
  if (use.intercept != "without") {
    total.parameters <- sum(total.parameters, 1)
  }

  ## Minimum Size --------------------------------------------------------------

  min.size.rs <- sum(1, total.parameters)

  if (min.size.rs > size.rs) {
    stop(paste0(
      "The size of your random sample is too small given the number ",
      "of variables and their lags. We suggest increasing the size",
      " of random samples or decreasing the maximum number of lags.",
      "\n Minimum sample size required given the current number of ",
      "variables and their lags: ", min.size.rs
    ))
  }

  # |__ Max Number of Random Samples Given Min Size of Rolling Sample ==========

  if (!is.null(y.max.lags)) tmpYLags <- y.max.lags else tmpYLags <- 0

  if (!is.null(x.max.lags)) {
    tmpXVLags <- sum(unlist(x.max.lags))
  } else tmpXVLags <- 0

  if (!is.null(x.fixed.lags)) {
    tmpXFLags <- sum(unlist(x.fixed.lags))
  } else tmpXFLags <- 0

  if (!is.null(x.interaction.lags)) {
    tmpXILags <- sum(unlist(x.interaction.lags))
  } else tmpXILags <- 0

  y.x.max.lags <- max(tmpYLags, tmpXVLags)
  y.x.max.lags <- max(y.x.max.lags, tmpXFLags)
  y.x.max.lags <- max(y.x.max.lags, tmpXILags)

  max.number.rs.given.min.size.rs <-
    last.obs - (2*forecast.horizon) - min.size.rs - y.x.max.lags + 1

  if (max.number.rs.given.min.size.rs < number.rs) {
    stop(paste0(
      "The number of random samples is too large given the size of ",
      "your random sample, the number of lags of your variables, ",
      "and the size of the forecast length. We suggest decreasing ",
      "the number of random samples, or the size of your forecast ",
      "horizon, or the maximum number of lags of your variables.",
      "\n Maximum number of random samples giving current inputs: ",
      max.number.rs.given.min.size.rs
    ))
  }

  # |__ Max Size of Random Samples Given Min Size of Rolling Sample ============

  max.size.rs <- last.obs - (2*forecast.horizon) - number.rs - y.x.max.lags + 1

  if (max.size.rs < size.rs) {
    stop(paste0(
      "The size of your random sample is too big given the total ",
      "size of your entire sample, the number of rolling samples ",
      "you selected, the size of your forecast horizon, and the ",
      "number of lags of your variables. We suggest decreasing ",
      "the number of random samples, or the size of your forecast ",
      "horizon, or the maximum number of lags of your variables.",
      "\n Maximum size of random samples giving current inputs: ",
      max.size.rs
    ))
  }

  # |__ Number of Lags Matches Number of X vars ================================

  ## X (variable lag)
  if (!is.null(x.names) & !is.null(x.max.lags)) {
    len.x.names <- length(x.names)
    len.x.max.lags <- length(x.max.lags)
    if (len.x.names != len.x.max.lags) {
      stop("Length of 'x.names' does not match the length of 'x.max.lags'.")
    }
  }

  ## X (fixed lag)

  if (!is.null(x.fixed.names) & !is.null(x.fixed.lags)) {
    len.x.fixed.names <- length(x.fixed.names)
    len.x.fixed.lags <- length(x.fixed.lags)
    if (len.x.fixed.names != len.x.fixed.lags) {
      stop(
        "Length of 'x.fixed.names' does not match the length of 'x.fixed.lags'."
      )
    }
  }

  ## X (interactions)

  if (!is.null(x.interaction.names) & !is.null(x.interaction.lags)) {
    len.x.interaction.names <- length(x.interaction.names)
    len.x.interaction.lags <- length(x.interaction.lags)
    if (len.x.interaction.names != len.x.interaction.lags) {
      stop(paste0(
        "Length of 'x.interaction.names' does not match the length of ",
        "'x.interaction.lags'."
      ))
    }
  }

}



