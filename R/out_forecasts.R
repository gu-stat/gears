# |__ Get Out-of-sample (ex-post) Forecasts ==================================

out_forecasts <- function(betas, ...) {

  # |__ Betas Selection ------------------------------------------------------

  if (betas == "last") {

    fcn.betas <- function(X, errorMeasure, ...) {

      tmp_fit <- estimates.best[[errorMeasure]][[X]][[number.rs]]

      tmp_forecast <- stats::predict(
        object  = tmp_fit,
        newdata = DF.Forecast[last.obs, ],
        type    = "response"
      )

      #names(tmp_forecast) <- sum(last.obs + X)

      names(tmp_forecast) <- NULL

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
