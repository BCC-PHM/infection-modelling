library(dplyr)
library(ggplot2)

source("gamma_params.R")


SEIR_model <- function(
  params,
  inits,
  timesteps = 100,
  t0 = 1
  ) {
  
  # Check and fill missing elements with zero
  for (elem in c("s", "e", "i", "r")) {
    if (!(elem %in% names(inits))) {
      inits[[elem]] <- 0
    }
  }
  
  # Calculate reciprocal values
  params$gamma = 1/params$`1/gamma`
  params$sigma = 1/params$`1/sigma`
  params$omega = 1/params$`1/omega`
  
  # Calculate infection 
  if (is.null(params$beta)) {
    params$beta = params$R0 * params$gamma
  }
  
  if (is.null(params$hosp)) {
    params$hosp = 0
  }
  
  
  s = rep(0, timesteps)
  e = rep(0, timesteps)
  i = rep(0, timesteps)
  r = rep(0, timesteps)
  d = rep(0, timesteps)
  admissions = rep(0, timesteps)
  
  infected_sum = rep(0, timesteps)
  
  # initial infected people
  i[1] = inits$i
  s[1] = inits$s
  e[1] = inits$e
  r[1] = inits$r
  
  for (t in 1:(timesteps-1)){
    
    s[t + 1] = s[t] +
     params$mu_b - # births
     params$beta * s[t] * i[t] + # infection
     params$omega * r[t] - # lost immunity
     params$vac * s[t] - # vaccination
     params$mu_d * s[t] # death
    
    e[t + 1] = e[t] +
     params$beta * s[t] * i[t] - # infection
     params$sigma * e[t] - # latency
     params$mu_d * e[t] # death
    
    i[t + 1] = i[t] +
     params$sigma * e[t] - # latency
     params$gamma * i[t] - # recovery
     (params$mu_d + params$alpha) * i[t] # death
    
    r[t + 1] = r[t] +
     params$gamma * i[t] + # recovery
     params$vac * s[t] - # vaccination
     params$omega * r[t] - # lost immunity
     params$mu_d * r[t] # death
    
    d[t+1] = d[t] +
     params$mu_d * s[t] +# death of susceptible
     params$mu_d * e[t] +# death of exposed
     (params$mu_d + params$alpha) * i[t] +# death of infected
     params$mu_d * r[t] # death of recovered
    
    infected_sum[t+1] = infected_sum[t] +
       params$beta * s[t] * i[t]
    
    admissions[t+1] <- params$hosp * i[t]
  }
  
  df <- data.frame(
    "Time" = t0:(t0 + timesteps + -1),
    "Susceptible" = s,
    "Exposed" = e,
    "Infectious" = i,
    "Recovered" = r,
    "Died" = d,
    "Infected_sum" = infected_sum,
    "Admissions" = admissions,
    "Admissions_sum" = cumsum(admissions)
  ) 
  
  return(df)
  }
  
SEIHR_model <- function(
    params,
    inits,
    timesteps = 100,
    t0 = 1
) {
  
  # Check and fill missing elements with zero
  for (elem in c("s", "e", "i", "r")) {
    if (!(elem %in% names(inits))) {
      inits[[elem]] <- 0
    }
  }
  
  # Calculate reciprocal values
  params$gamma = 1/params$`1/gamma`
  params$sigma = 1/params$`1/sigma`
  params$omega = 1/params$`1/omega`
  
  # Calculate infection 
  params$beta = params$R0 * params$gamma
  
  s = rep(0, timesteps)
  e = rep(0, timesteps)
  i = rep(0, timesteps)
  h = rep(0, timesteps)
  r = rep(0, timesteps)
  d = rep(0, timesteps)
  addmissions = rep(0, timesteps)
  
  infected_sum = rep(0, timesteps)
  
  # initial infected people
  i[1] = inits$i
  s[1] = inits$s
  e[1] = inits$e
  r[1] = inits$r
  h[1] = inits$h
  
  for (t in 1:(timesteps-1)){
    
    s[t + 1] = s[t] +
      params$mu_b - # births
      params$beta * s[t] * i[t] + # infection
      params$omega * r[t] - # lost immunity
      params$vac * s[t] - # vaccination
      params$mu_d * s[t] # death
    
    e[t + 1] = e[t] +
      params$beta * s[t] * i[t] - # infection
      params$sigma * e[t] - # latency
      params$mu_d * e[t] # death
    
    i[t + 1] = i[t] +
      params$sigma * e[t] - # latency
      params$gamma * i[t] - # recovery
      params$hosp * i[t] - # hospitalization
      (params$mu_d + params$alpha) * i[t] # death
    
    admissions[t+1] <- params$hosp * i[t]
    
    h[t+1] = h[t] + 
      admissions[t+1] - # hospitalization
      params$gamma * h[t] - # recovery 
      params$mu_d * h[t] # death
    
    r[t + 1] = r[t] +
      params$gamma * i[t] + # infected recovery
      params$gamma * h[t] + # hospital recovery
      params$vac * s[t] - # vaccination
      params$omega * r[t] - # lost immunity
      params$mu_d * r[t] # death
    
    d[t+1] = d[t] +
      params$mu_d * s[t] +# death of susceptible
      params$mu_d * e[t] +# death of exposed
      (params$mu_d + params$alpha) * i[t] +# death of infected
      params$mu_d * r[t] # death of recovered
    
    infected_sum[t+1] = infected_sum[t] +
      params$beta * s[t] * i[t]
  }
  
  df <- data.frame(
    "Time" = t0:(t0 + timesteps + -1),
    "Susceptible" = s,
    "Exposed" = e,
    "Infectious" = i,
    "Recovered" = r,
    "Hospital" = h,
    "Died" = d,
    "Infected_sum" = infected_sum,
    "Admissions" = admissions
  ) 
  
  return(df)
}

param_sample <- function(
    params_df
) {
  # calculate a random sample for each parameter in the data frame.
  
  # Get Poisson rate and shape for each parameter
  shape <- c()
  rate  <- c()
  for (i in 1:nrow(params_df)) {
    params <- gamma_params(
      params_df$Value[i], 
      params_df$lowerCI[i], 
      params_df$upperCI[i],
      params_df$Distribution[i])
    shape[i] <- params$shape
    rate[i] <- params$rate
  }

  params_df$shape = shape
  params_df$rate = rate
  
  # I feel like theres a better way to get rid of these warnings
  suppressWarnings(
  params_df <- params_df %>%
    rowwise() %>%
    mutate(
      Value = case_when (
        Distribution == "Gamma" ~ 
          rgamma(1, rate = rate, shape = shape),
        Distribution == "Uniform" ~
          runif(1, lowerCI, upperCI),
        Distribution == "Fixed" ~ Value,
        TRUE ~ -1
      )
    ) %>%
    ungroup() %>%
    select(c(Parameter, Value))
  )
  
  # Convert to list
  output_list <- setNames(as.list(params_df$Value),    # Convert vectors to named list
                          params_df$Parameter)
  
  return(output_list)
}







SEIR_surge_model <- function(
    params,
    inits,
    timesteps = 100,
    t0 = 1,
    t_surge=0,    #specific time step to surge
    surge_factor = 1 #the surge factor
) {
  
  # Check and fill missing elements with zero
  for (elem in c("s", "e", "i", "r")) {
    if (!(elem %in% names(inits))) {
      inits[[elem]] <- 0
    }
  }
  
  # Calculate reciprocal values
  params$gamma = 1/params$`1/gamma`
  params$sigma = 1/params$`1/sigma`
  params$omega = 1/params$`1/omega`
  
  # Calculate infection 
  if (is.null(params$beta)) {
    params$beta = params$R0 * params$gamma
  }
  
  if (is.null(params$hosp)) {
    params$hosp = 0
  }
  

  
  s = rep(0, timesteps)
  e = rep(0, timesteps)
  i = rep(0, timesteps)
  r = rep(0, timesteps)
  d = rep(0, timesteps)
  admissions = rep(0, timesteps)
  
  infected_sum = rep(0, timesteps)
  
  # initial infected people
  i[1] = inits$i
  s[1] = inits$s
  e[1] = inits$e
  r[1] = inits$r
  
  for (t in 1:(timesteps-1)){
    
    #surge modelling
    
    if (t != t_surge) {
      beta = params$beta
    } else {
      beta = surge_factor * params$beta
    }
    
    
    s[t + 1] = s[t] +
      params$mu_b - # births
      beta * s[t] * i[t] + # infection
      params$omega * r[t] - # lost immunity
      params$vac * s[t] - # vaccination
      params$mu_d * s[t] # death
    
    e[t + 1] = e[t] +
      beta * s[t] * i[t] - # infection
      params$sigma * e[t] - # latency
      params$mu_d * e[t] # death
    
    i[t + 1] = i[t] +
      params$sigma * e[t] - # latency
      params$gamma * i[t] - # recovery
      (params$mu_d + params$alpha) * i[t] # death
    
    r[t + 1] = r[t] +
      params$gamma * i[t] + # recovery
      params$vac * s[t] - # vaccination
      params$omega * r[t] - # lost immunity
      params$mu_d * r[t] # death
    
    d[t+1] = d[t] +
      params$mu_d * s[t] +# death of susceptible
      params$mu_d * e[t] +# death of exposed
      (params$mu_d + params$alpha) * i[t] +# death of infected
      params$mu_d * r[t] # death of recovered
    
    infected_sum[t+1] = infected_sum[t] +
      beta * s[t] * i[t]
    
    admissions[t+1] <- params$hosp * i[t]
  }
  
  df <- data.frame(
    "Time" = t0:(t0 + timesteps + -1),
    "Susceptible" = s,
    "Exposed" = e,
    "Infectious" = i,
    "Recovered" = r,
    "Died" = d,
    "Infected_sum" = infected_sum,
    "Admissions" = admissions,
    "Admissions_sum" = cumsum(admissions)
  ) 
  
  return(df)
}
































