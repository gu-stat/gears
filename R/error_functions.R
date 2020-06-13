# ERROR MEASURES

error_functions <- function(error.measure, forecasts, outsample, insample = NULL) {
  
  # HELPERS
  
  if (error.measure == "mse") {
    
    tmp.error <- mean((forecasts-outsample)**2)
    
  } else if (error.measure == "mad") {
    
    tmp.error <- mean(abs(forecasts-outsample))
    
  } else if (error.measure == "smape") {
    
    smape_cal <- function(outsample, forecasts){
      #Used to estimate sMAPE
      outsample <- as.numeric(outsample) ; forecasts<-as.numeric(forecasts)
      smape <- (abs(outsample-forecasts)*200)/(abs(outsample)+abs(forecasts))
      return(smape)
    }
    
    tmp.error <- mean(smape_cal(outsample = outsample, forecasts = forecasts))
    
  } else if (error.measure == "mase") {
    
    mase_cal <- function(insample, outsample, forecasts){
      #Used to estimate MASE
      frq <- frequency(insample)
      forecastsNaiveSD <- rep(NA,frq)
      for (j in (frq+1):length(insample)){
        forecastsNaiveSD <- c(forecastsNaiveSD, insample[j-frq])
      }
      masep<-mean(abs(insample-forecastsNaiveSD),na.rm = TRUE)
      
      outsample <- as.numeric(outsample) ; forecasts <- as.numeric(forecasts)
      mase <- (abs(outsample-forecasts))/masep
      return(mase)
    }
    
    tmp.error <- mean(mase_cal(
      outsample = outsample, 
      forecasts = forecasts,
      insample  = insample 
    ))
    
  }
  
  return(tmp.error)
  
}