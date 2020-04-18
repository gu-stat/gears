# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ----
# create_DF_Forecast Function                                               ----
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# DATA.PORK <- read.table("pork-360.txt", header = TRUE)
# 
# DATA.PORK <- DATA.PORK[1:10, ]
# 
# DATA                = DATA.PORK
# forecast.horizon    = 2
# y.name              = "PORK_PRICE"
# y.max.lags          = 2
# x.names             = list("BEEF_PRICE")
# x.max.lags          = list(3)
# x.fixed.names       = NULL
# x.fixed.lags        = NULL
# x.interaction.names = list("CORN_PRICE*CORN_PRICE", "BEEF_PRICE*WHEAT_PRICE")
# x.interaction.lags  = list(c(1, 1), c(3,3))


# ************************************************************************* ----
# Main Function                                                             ----
# >                                                                         ----

create_DF_Forecast <- function(DATA, 
                               forecast.horizon,
                               y.name = NULL, 
                               y.max.lags = NULL,
                               x.names = NULL, 
                               x.max.lags = NULL,
                               x.fixed.names = NULL, 
                               x.fixed.lags = NULL,
                               x.interaction.names = NULL, 
                               x.interaction.lags = NULL) {
  
  # > Helpers ############################################################# ----
  
  fcn.lag <- function(df, n.lags) c(rep(NA, n.lags), df)[1:length(df)]
  
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
    
    # Get Frequency
    
    input.frequency <- frequency(DATA)
    
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
        function(i = X) {
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
  
  if (!is.null(y.max.lags)) {
    
    # |__ DF.Y.forecast ========================================================
    
    DF.Y.forecast <- cbind(
      DF.mod.Y[, Y_name],
      sapply(X = 1:y.max.lags, FUN = fcn.lag, df = DF.mod.Y[, Y_name])
    )
    
    colnames(DF.Y.forecast) <- c(
      Y_name, 
      paste0(Y_name, "_minus_", 1:y.max.lags)
    )
    
  } else {
    DF.Y.forecast <- data.frame(NA)
  }
  
  # |__ DF.X_Lagged.forecast ===================================================
  
  if (is.null(DF.mod.X.lag)) {
    DF.X_Lagged.forecast <- data.frame(NA)
  } else {
    
    DF.X_Lagged.forecast <- cbind(
      DF.mod.X.lag,
      do.call(
        cbind,
        lapply(
          X = 1:length(X_lagged_names),
          function(j = X) {
            do.call(
              cbind,
              sapply(
                X = 1:unlist(x.max.lags)[j], 
                function(L = X) {
                  
                  tmp.lags <- as.data.frame(
                    fcn.lag(DF.mod.X.lag[, X_lagged_names[j]], n.lags = L)
                  )
                  
                  colnames(tmp.lags) <- paste0(X_lagged_names[j], "_minus_", L)
                  
                  return(tmp.lags)
                  
                }
              )
            )
          }
        )
      )
    )
    
  }
  
  # |__ DF.X_Fixed.forecsat ====================================================
  if (is.null(DF.mod.X.fixed)) {
    DF.X_Fixed.forecast <- data.frame(NA)
  } else {
    DF.X_Fixed.forecast <- do.call(
      cbind,
      lapply(
        X = 1:length(x.fixed.var.vec),
        function(j = X) {
          do.call(
            cbind,
            sapply(
              X = unlist(x.fixed.lags)[j], 
              function(L = X) {
                
                tmp.lags <- as.data.frame(
                  fcn.lag(DF.mod.X.fixed[, x.fixed.var.vec[j]], n.lags = L)
                )
                
                colnames(tmp.lags) <- paste0(x.fixed.var.vec[j])
                
                return(tmp.lags)
                
              }
            )
          )
        }
      )
    )
  }
  
  # |__ DF.Inter.forecast ======================================================
  
  if (is.null(DF.mod.X.inter)) {
    DF.Inter.forecast <- data.frame(NA)
  } else {
    
    tmp.DF.Inter.forecast <- do.call(
      cbind,
      lapply(
        X = 1:length(Inter_names),
        function(j = X) {
          do.call(
            cbind,
            sapply(
              X = 1:unlist(x.interaction.lags)[j], 
              function(L = X) {
                
                tmp.lags <- as.data.frame(
                  fcn.lag(DF.mod.X.inter[, Inter_names][, j], n.lags = L)
                )
                
                colnames(tmp.lags) <- paste0(
                  colnames(DF.mod.X.inter[, Inter_names])[j], "_minus_", L
                )
                
                return(tmp.lags)
                
              }
            )
          )
        }
      )
    )
    
    DF.Inter.forecast <- as.data.frame(
      tmp.DF.Inter.forecast[, unique(x.interaction.var.vec)]
    )
    
    colnames(DF.Inter.forecast) <- unique(x.interaction.var.vec)
    
  }
  
  # |__ DF.Forecast ============================================================
  
  DF.Forecast <- cbind(
    DF.Y.forecast,
    DF.X_Lagged.forecast,
    DF.X_Fixed.forecast,
    DF.Inter.forecast
  )
  
  # \____ Remove Duplicates ====================================================
  
  tmp.unique.vars <- unique(colnames(DF.Forecast))
  
  DF.Forecast <- DF.Forecast[, tmp.unique.vars]
  
  DF.Forecast <- subset(
    DF.Forecast, 
    select = !colnames(DF.Forecast) %in% c("NA.")
  )
  
  # > Return Results ###################################################### ----
  
  return(DF.Forecast)

}

## TESTS ----
# 
# DATA.PORK <- read.table("pork-360.txt", header = TRUE)
# 
# DATA.PORK <- DATA.PORK[1:10, ]
# 
# create_DF_Forecast(
#   DATA                = DATA.PORK,
#   forecast.horizon    = 2,
#   y.name              = "PORK_PRICE" ,
#   y.max.lags          = 2,
#   x.names             = list("BEEF_PRICE") ,
#   x.max.lags          = list(3),
#   x.fixed.names       = list("CORN_PRICE", "CORN_PRICE") ,
#   x.fixed.lags        = list(4, 5),
#   x.interaction.names = list("CORN_PRICE*CORN_PRICE", "BEEF_PRICE*WHEAT_PRICE"),
#   x.interaction.lags  = list(c(1, 1), c(3,3))
# )
# 
# create_DF_Forecast(
#   DATA                = DATA.PORK,
#   forecast.horizon    = 2,
#   y.name              = "PORK_PRICE" ,
#   y.max.lags          = 2,
#   x.names             = list("BEEF_PRICE") ,
#   x.max.lags          = list(3),
#   x.fixed.names       = NULL,
#   x.fixed.lags        = NULL,
#   x.interaction.names = list("CORN_PRICE*CORN_PRICE"),
#   x.interaction.lags  = list(c(1, 1))
# )
# 
# create_DF_Forecast(
#   DATA                = DATA.PORK,
#   forecast.horizon    = 2,
#   y.name              = "PORK_PRICE" ,
#   y.max.lags          = 2,
#   x.names             = list("BEEF_PRICE", "WHEAT_PRICE") ,
#   x.max.lags          = list(3, 2),
#   x.fixed.names       = NULL,
#   x.fixed.lags        = NULL,
#   x.interaction.names = list("CORN_PRICE*CORN_PRICE"),
#   x.interaction.lags  = list(c(1, 1))
# )
# 
# create_DF_Forecast(
#   DATA                = DATA.PORK,
#   forecast.horizon    = 2,
#   y.name              = "PORK_PRICE" ,
#   y.max.lags          = 2,
#   x.names             = list("BEEF_PRICE", "WHEAT_PRICE") ,
#   x.max.lags          = list(3, 2),
#   x.fixed.names       = NULL,
#   x.fixed.lags        = NULL,
#   x.interaction.names = list("CORN_PRICE*CORN_PRICE", "BEEF_PRICE*WHEAT_PRICE"),
#   x.interaction.lags  = list(c(1, 1), c(1, 1))
# )
# 
# create_DF_Forecast(
#   DATA                = DATA.PORK,
#   forecast.horizon    = 2,
#   y.name              = "PORK_PRICE" ,
#   y.max.lags          = 2,
#   x.names             = list("BEEF_PRICE", "WHEAT_PRICE") ,
#   x.max.lags          = list(3, 2),
#   x.fixed.names       = NULL,
#   x.fixed.lags        = NULL,
#   x.interaction.names = list("CORN_PRICE*CORN_PRICE", "BEEF_PRICE*WHEAT_PRICE"),
#   x.interaction.lags  = list(c(1, 1), c(3, 3))
# )
# 
# create_DF_Forecast(
#   DATA                = DATA.PORK,
#   forecast.horizon    = 2,
#   y.name              = "PORK_PRICE" ,
#   y.max.lags          = 2,
#   x.names             = NULL,
#   x.max.lags          = NULL,
#   x.fixed.names       = NULL,
#   x.fixed.lags        = NULL,
#   x.interaction.names = list("CORN_PRICE*CORN_PRICE", "BEEF_PRICE*WHEAT_PRICE"),
#   x.interaction.lags  = list(c(1, 1), c(3, 3))
# )
# 
# create_DF_Forecast(
#   DATA                = DATA.PORK,
#   forecast.horizon    = 2,
#   y.name              = "PORK_PRICE" ,
#   y.max.lags          = 2,
#   x.names             = NULL,
#   x.max.lags          = NULL,
#   x.fixed.names       = NULL,
#   x.fixed.lags        = NULL,
#   x.interaction.names = list("CORN_PRICE*CORN_PRICE"),
#   x.interaction.lags  = list(c(1, 1))
# )
# 
# #/***********************/#
#   
# create_DF_Forecast(
#   DATA                = DATA.PORK,
#   forecast.horizon    = 2,
#   y.name              = "PORK_PRICE" ,
#   y.max.lags          = NULL,
#   x.names             = list("BEEF_PRICE") ,
#   x.max.lags          = list(3),
#   x.fixed.names       = NULL,
#   x.fixed.lags        = NULL,
#   x.interaction.names = list("CORN_PRICE*CORN_PRICE", "BEEF_PRICE*WHEAT_PRICE"),
#   x.interaction.lags  = list(c(1, 1), c(3,3))
# )
# 
# create_DF_Forecast(
#   DATA                = DATA.PORK,
#   forecast.horizon    = 2,
#   y.name              = "PORK_PRICE" ,
#   y.max.lags          = NULL,
#   x.names             = list("BEEF_PRICE") ,
#   x.max.lags          = list(3),
#   x.fixed.names       = NULL,
#   x.fixed.lags        = NULL,
#   x.interaction.names = list("CORN_PRICE*CORN_PRICE"),
#   x.interaction.lags  = list(c(1, 1))
# )
# 
# create_DF_Forecast(
#   DATA                = DATA.PORK,
#   forecast.horizon    = 2,
#   y.name              = "PORK_PRICE" ,
#   y.max.lags          = NULL,
#   x.names             = list("BEEF_PRICE", "WHEAT_PRICE") ,
#   x.max.lags          = list(3, 2),
#   x.fixed.names       = NULL,
#   x.fixed.lags        = NULL,
#   x.interaction.names = list("CORN_PRICE*CORN_PRICE"),
#   x.interaction.lags  = list(c(1, 1))
# )
# 
# create_DF_Forecast(
#   DATA                = DATA.PORK,
#   forecast.horizon    = 2,
#   y.name              = "PORK_PRICE" ,
#   y.max.lags          = NULL,
#   x.names             = list("BEEF_PRICE", "WHEAT_PRICE") ,
#   x.max.lags          = list(3, 2),
#   x.fixed.names       = NULL,
#   x.fixed.lags        = NULL,
#   x.interaction.names = list("CORN_PRICE*CORN_PRICE", "BEEF_PRICE*WHEAT_PRICE"),
#   x.interaction.lags  = list(c(1, 1), c(1, 1))
# )
# 
# create_DF_Forecast(
#   DATA                = DATA.PORK,
#   forecast.horizon    = 2,
#   y.name              = "PORK_PRICE" ,
#   y.max.lags          = NULL,
#   x.names             = list("BEEF_PRICE", "WHEAT_PRICE") ,
#   x.max.lags          = list(3, 2),
#   x.fixed.names       = NULL,
#   x.fixed.lags        = NULL,
#   x.interaction.names = list("CORN_PRICE*CORN_PRICE", "BEEF_PRICE*WHEAT_PRICE"),
#   x.interaction.lags  = list(c(1, 1), c(3, 3))
# )
# 
# create_DF_Forecast(
#   DATA                = DATA.PORK,
#   forecast.horizon    = 2,
#   y.name              = "PORK_PRICE" ,
#   y.max.lags          = NULL,
#   x.names             = NULL,
#   x.max.lags          = NULL,
#   x.fixed.names       = NULL,
#   x.fixed.lags        = NULL,
#   x.interaction.names = list("CORN_PRICE*CORN_PRICE", "BEEF_PRICE*WHEAT_PRICE"),
#   x.interaction.lags  = list(c(1, 1), c(3, 3))
# )
# 
# create_DF_Forecast(
#   DATA                = DATA.PORK,
#   forecast.horizon    = 2,
#   y.name              = "PORK_PRICE" ,
#   y.max.lags          = NULL,
#   x.names             = NULL,
#   x.max.lags          = NULL,
#   x.fixed.names       = NULL,
#   x.fixed.lags        = NULL,
#   x.interaction.names = list("CORN_PRICE*CORN_PRICE"),
#   x.interaction.lags  = list(c(1, 1))
# )
# 
# create_DF_Forecast(
#   DATA                = DATA.PORK,
#   forecast.horizon    = 2,
#   y.name              = "PORK_PRICE" ,
#   y.max.lags          = 2,
#   x.names             = list("BEEF_PRICE", "WHEAT_PRICE") ,
#   x.max.lags          = list(3, 2),
#   x.fixed.names       = NULL,
#   x.fixed.lags        = NULL,
#   x.interaction.names = NULL,
#   x.interaction.lags  = NULL
# )
# 
# create_DF_Forecast(
#   DATA                = DATA.PORK,
#   forecast.horizon    = 2,
#   y.name              = "PORK_PRICE" ,
#   y.max.lags          = 2,
#   x.names             = NULL,
#   x.max.lags          = NULL,
#   x.fixed.names       = NULL,
#   x.fixed.lags        = NULL,
#   x.interaction.names = NULL,
#   x.interaction.lags  = NULL
# )
# 
# create_DF_Forecast(
#   DATA                = DATA.PORK,
#   forecast.horizon    = 2,
#   y.name              = "PORK_PRICE" ,
#   y.max.lags          = NULL,
#   x.names             = list("BEEF_PRICE", "WHEAT_PRICE") ,
#   x.max.lags          = list(3, 2),
#   x.fixed.names       = NULL,
#   x.fixed.lags        = NULL,
#   x.interaction.names = NULL,
#   x.interaction.lags  = NULL
# )