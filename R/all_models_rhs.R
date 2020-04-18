# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ----
# all_models_rhs Function                                                   ----
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

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
    paste0("~", apply(combn(unique.vars.vec, x), 2, paste, collapse = "+"))
  }
  
  
  fcn.eqs.without.int <- function(x) {
    paste0("~-1 +", apply(combn(unique.vars.vec, x), 2, paste, collapse = "+"))
  }
  
  # > Variables ########################################################### ----
  
  # Define Y variable name
  
  if (is.null(y.name)) {
    Y_name <- "Y_t"
  } else {
    Y_name <- paste0(y.name, "_t")
  }
  
  if (is.null(y.max.lags)) {
    y.rhs.var.vec <- c(Y_name)
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

# ## TEST
# 
# teste <- all_models_rhs(
  # y.name              = "PORK_PRICE" ,
  # y.max.lags          = 2,
  # x.names             = NULL,
  # x.max.lags          = NULL,
  # x.fixed.names       = NULL,
  # x.fixed.lags        = NULL,
  # x.interaction.names = NULL,
  # x.interaction.lags  = NULL,
  # use.intercept = "both"
# )
# 
# head(teste)
# tail(teste)
# 
# teste <- all_models_rhs(
#   y.name              = "PORK_PRICE" ,
#   y.max.lags          = 2,
#   x.names             = list("BEEF_PRICE") ,
#   x.max.lags          = list(3),
#   x.fixed.names       = list("CORN_PRICE", "CORN_PRICE") ,
#   x.fixed.lags        = list(4, 5),
#   x.interaction.names = list("CORN_PRICE*CORN_PRICE", "BEEF_PRICE*WHEAT_PRICE"),
#   x.interaction.lags  = list(c(1, 1), c(3,3)),
#   use.intercept = "both"
# )
# 
# head(teste)
# tail(teste)