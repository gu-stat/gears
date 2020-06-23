#' Generates all possible combinations of right-hand side variables.
#'
#' Takes in the names and lags of the variables, and generates all possible
#' combinations of these to form the right-hand side of the model equations.
#'
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
#' @param use.intercept If \code{use.intercept == "both"} (default), the
#'     function returns all possible model equations with and without intercept.
#'     If \code{use.intercept == "without"}, then the function returns all
#'     possible right-hand side equations without intercept. If
#'     \code{use.intercept == "with"}, then only the right-hand side equations
#'     with intercept are returned.
#'
#' @return A character vector containing all possible combinations of the input
#'     variables and lags that form the right-hand side of the model equations.
#'
#' @details If \code{y.max.lags} equals to a number, then all lags of
#'     \code{y.name} up to this number (and starting at 0) will be included
#'     in the list of variables to enter the right-hand side of the model
#'     equations. For example, if \code{y.max.lags = 2}, then
#'     \ifelse{html}{\out{Y<sub>t</sub>, Y<sub>t-1</sub>, Y<sub>t-2</sub>
#'     }}{\eqn{Y_t, Y_{t-1},Y_{t-2}}} will be included in the list of variables
#'     to enter the right-hand side of the equation.
#'
#' @keywords internal
#' @examples
#' all_models_rhs(y.name = "PORK_PRICE", y.max.lags = 2, use.intercept = "both")
#'
#' all_models_rhs(
#'   y.name              = "PORK_PRICE" ,
#'   y.max.lags          = NULL,
#'   x.names             = list("BEEF_PRICE") ,
#'   x.max.lags          = list(3),
#'   x.fixed.names       = list("CORN_PRICE", "CORN_PRICE") ,
#'   x.fixed.lags        = list(4, 5),
#'   x.interaction.names = list("CORN_PRICE*CORN_PRICE", "BEEF_PRICE*WHEAT_PRICE"),
#'   x.interaction.lags  = list(c(1, 1), c(3,3)),
#'   use.intercept       = "without"
#' )
all_models_rhs <- function(y.name = NULL,
                           y.max.lags = NULL,
                           x.names = NULL,
                           x.max.lags = NULL,
                           x.fixed.names = NULL,
                           x.fixed.lags = NULL,
                           x.interaction.names = NULL,
                           x.interaction.lags = NULL,
                           use.intercept = c("both", "without", "with")
                           ) {

  # > Helpers ############################################################# ----

  # Create Equations' RHS

  fcn.eqs.with.int <- function(x) {
    paste0(
      "~",
      apply(utils::combn(unique.vars.vec, x), 2, paste, collapse = "+")
    )
  }


  fcn.eqs.without.int <- function(x) {
    paste0(
      "~-1 +",
      apply(utils::combn(unique.vars.vec, x), 2, paste, collapse = "+")
    )
  }

  # > Variables ########################################################### ----

  # Define Y variable name

  if (is.null(y.name)) {
    Y_name <- "Y_t"
  } else {
    Y_name <- paste0(y.name, "_t")
  }

  if (is.null(y.max.lags)) {
    y.rhs.var.vec <- NULL
  } else {
    y.rhs.var.vec <- c(Y_name, paste0(Y_name, "_minus_", 1:y.max.lags))
  }

  # Define X-lagged names

  if (is.null(x.max.lags)) {
    x.var.vec <- c(x.names)
  } else {

    tmp_x_lag_names <- unlist(x.names)

    X_lagged_names_list <- lapply(
      X   = 1:length(x.names),
      FUN = function(X) paste0(tmp_x_lag_names[[X]], "_t")
    )

    ## Change X variables in original dataset to X_t

    X_lagged_names <- unlist(X_lagged_names_list)

    x.var.vec <- c(
      X_lagged_names,
      unlist(sapply(
        X = 1:length(x.names),
        function(i = X) {
          sapply(
            X = 1:x.max.lags[[i]],
            function(X) {
              paste0(x.names[[i]], "_t_minus_", X)
            }
          )
        }
      ))
    )

  }

  # Define X-fixed names
  if (is.null(x.fixed.lags)) {
    x.fixed.var.vec <- c(x.fixed.names)
  } else {
    x.fixed.var.vec <- paste0(
      unlist(x.fixed.names), "_t_minus_", unlist(x.fixed.lags)
    )
  }

  if (is.null(x.interaction.names)) {
    x.interaction.var.vec <- c(x.interaction.names)
    x.interaction.vars <- c(x.interaction.names)
  } else {

    ## Add _t to variable name (e.g., X becomes X_t)

    tmp_Inter_names <- strsplit(x = unlist(x.interaction.names), split ="\\*")

    x.interaction.vars <- sapply(
      X = 1:length(x.interaction.names),
      function(X) {
        paste(
          paste0(tmp_Inter_names[[X]], "_t_minus_", x.interaction.lags[[X]]),
          collapse = "*"
        )
      }
    )
  }

  # |__ All Equations ==========================================================

  var.vec.no.int <- c(y.rhs.var.vec, x.var.vec, x.fixed.var.vec)

  unique.vars.vec <- c(unique(var.vec.no.int), x.interaction.vars)

  # \____ RHS  -----------------------------------------------------------------

  if (use.intercept == "with") {
    all.equations.rhs <- unlist(
      sapply(X = 1:length(unique.vars.vec), FUN = fcn.eqs.with.int)
    )
  } else if (use.intercept == "without") {
    all.equations.rhs <- unlist(
      sapply(X = 1:length(unique.vars.vec), FUN = fcn.eqs.without.int)
    )
  } else {
    all.equations.rhs <- c(
      unlist(sapply(X = 1:length(unique.vars.vec), FUN = fcn.eqs.with.int)),
      unlist(sapply(X = 1:length(unique.vars.vec), FUN = fcn.eqs.without.int))
    )
  }

  return(all.equations.rhs)

}
