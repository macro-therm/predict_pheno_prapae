# GitHub repo: github.com/dario-ssm/phenodev_predict

# Aim: functions script

# 1. TPC Model Fitting ----

# functions developed based on devRate formulas
wang <- function(temp, k, r, topt, tmin, tmax, a){
  est <- (k/(1 + exp(-r * (temp - topt)))) * (1 - exp(-(temp - tmin)/a)) *
    (1 - exp(-(tmax - temp)/a))
  return(est)
}
mod_polynomial <- function(temp, a_0, a_1, a_2, a_3, a_4){
  est <- a_0 + a_1 * temp + a_2 * temp^2 + a_3 * temp^3 + a_4 * temp^4
  return(est)
}

briere1 <- function(temp, tmin, tmax, a) {
  est <- a * temp * (temp - tmin) * (tmax - temp)^(1/2)
  return(est)
}

lactin1 <- function(temp, a, tmax, delta_t) {
  est <- exp(a * temp) - exp(a * tmax - (tmax - temp)/delta_t)
  return(est)
}

regniere <- function(temp,tmin, tmax, phi, delta_b, delta_m, b) {
  est <- phi* (exp(b * (temp - tmin)) - ((tmax - temp)/(tmax - tmin)) * exp(-b *
                                                                              (temp - tmin)/delta_b) - ((temp - tmin)/(tmax - tmin)) * exp(b * (tmax - tmin) - (tmax - temp)/delta_m))
  return(est)
}


# 2. Predictions ----------------------------------------------


## a) linear degree-days ---------------------------------------------------

lm_pred <- function(lm_fit, tavg) {
  alpha <- coef(lm_fit)[1]
  beta <- coef(lm_fit)[2]
  preds_lm <-  alpha + beta*tavg
  return(preds_lm)
}

logic_dd <- function(heat_units, rate_cumsum){
  logic_dd = if_else(rate_cumsum < heat_units,
                     1,
                     0)
  return(logic_dd)
}

calc_dd <- function(lm_fit, tavg) {
  ldt <- -coef(lm_fit)[1]/coef(lm_fit)[2]
  time_id_dds <- case_when(tavg >= ldt ~ tavg-ldt,
                           .default = 0)
  return(time_id_dds)
}


## b) nonlinear rate summation ---------------------------------------------------

rate_pred <- function(model2ratesum, temperature, fitted_params, res) {
  fitted_param_tbl_i <- fitted_params |> 
    filter(model_name == model2ratesum)
  params_i <- coef(fitted_param_tbl_i[1, ]$model_fit[[1]])
  available_models_i <- available_models |> 
    filter(model_name == model2ratesum)
  temperature_i <- temperature
  working_formula_i <- available_models_i$working_formula
  if(res == "daily") {
    rate_preds <- temperature_i |> 
      mutate(rate_pred = map_dbl(.x = daily_tavg,
                                 .f = reformulate(working_formula_i)),
             model_name = model2ratesum)
  } else if (res == "hourly") {
    rate_preds <- temperature_i |> 
      mutate(rate_pred = map_dbl(.x = temperature,
                                 .f = reformulate(working_formula_i)),
             model_name = model2ratesum)
  } else {
    stop("`res` must be either `daily` or `hourly`")
  }
  
  return(rate_preds)
}

logic_ratesum <- function(rate_cumsum){
  logic_rate_summation = if_else(rate_cumsum < 100,
                     1,
                     0)
  return(logic_rate_summation)
}



