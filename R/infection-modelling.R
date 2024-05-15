library(dplyr)
library(ggplot2)

source("gamma_params.R")


SEIR_model <- function(
  params,
  inits,
  timesteps = 100
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
  r = rep(0, timesteps)
  d = rep(0, timesteps)
  
  infected_sum = rep(0, timesteps)
  
  # initial infected people
  i[1] = inits$i
  s[1] = inits$s
  e[1] = inits$e
  r[1] = inits$r
  
  print(params)
  
  for (t in 1:(timesteps-1)){
    
    s[t + 1] = s[t] +
     params$mu - # births
     params$beta * s[t] * i[t] + # infection
     params$omega * r[t] - # lost immunity
     params$vac * s[t] - # vaccination
     params$mu * s[t] # death
    
    e[t + 1] = e[t] +
     params$beta * s[t] * i[t] - # infection
     params$sigma * e[t] - # latency
     params$mu * e[t] # death
    
    i[t + 1] = i[t] +
     params$sigma * e[t] - # latency
     params$gamma * i[t] - # recovery
     (params$mu + params$alpha) * i[t] # death
    
    r[t + 1] = r[t] +
     params$gamma * i[t] + # recovery
     params$vac * s[t] - # vaccination
     params$omega * r[t] - # lost immunity
     params$mu * r[t] # death
    
    d[t+1] = d[t] +
     params$mu * s[t] +# death of susceptible
     params$mu * e[t] +# death of exposed
     (params$mu + params$alpha) * i[t] +# death of infected
     params$mu * r[t] # death of recovered
    
    infected_sum[t+1] = infected_sum[t] +
       params$beta * s[t] * i[t]
  }
  
  df <- data.frame(
    "time" = 1:timesteps,
    "susceptible" = s,
    "exposed" = e,
    "infectious" = i,
    "recovered" = r,
    "died" = d,
    "infected sum" = infected_sum
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
      params_df$UpperCI[i])
    shape[i] <- params$shape
    rate[i] <- params$rate
  }

  params_df$shape = shape
  params_df$rate = rate
  
  params_df <- params_df %>%
    rowwise() %>%
    mutate(
      Value = case_when (
        Distribution == "Gamma" ~ 
          rgamma(1, rate = rate, shape = shape),
        Distribution == "Uniform" ~
          runif(1, lowerCI, UpperCI),
        Distribution == "Fixed" ~ Value,
        TRUE ~ -1
      )
    ) %>%
    ungroup() %>%
    select(c(Parameter, Value))
  
  # Convert to list
  output_list <- setNames(as.list(params_df$Value),    # Convert vectors to named list
                          params_df$Parameter)
  
  return(output_list)
}

