# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ----
# predictionError Function                                                  ----
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# ************************************************************************* ----
# Main Function                                                             ----
# >                                                                         ----

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
  
  # tmp_final_result <- lapply(
  #   X = 1:forecast.horizon,
  #   function(X) unlist(prediction.errors[,X])
  # )
  # 
  # #tmp_final_result
  
  return(prediction.errors)
  
}

