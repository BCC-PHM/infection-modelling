library(dplyr)
library(ggplot2)

source("gamma_params.R")


SEIR_model_prob <- function(
    params,
    inits,
    timesteps = 100
) {
  
  #create an empty list
  df <- list()
  
  # Iterate over each set of parameters
  for (n in seq_along(params)) {
    current_params <- params[[n]]
    
    # Calculate reciprocal values
    current_params$gamma <- 1 / current_params$`1/gamma`
    current_params$sigma <- 1 / current_params$`1/sigma`
    current_params$omega <- 1 / current_params$`1/omega`
    
    # Calculate infection 
    current_params$beta <- current_params$R0 * current_params$gamma
    
    s <- rep(0, timesteps)
    e <- rep(0, timesteps)
    i <- rep(0, timesteps)
    r <- rep(0, timesteps)
    d <- rep(0, timesteps)
    infected_sum <- rep(0, timesteps)
    
    # initial infected people
    i[1] <- inits$i
    s[1] <- inits$s
    e[1] <- inits$e
    r[1] <- inits$r
    
    for (t in 1:(timesteps - 1)) {
      
      s[t + 1] <- s[t] +
        current_params$mu_b - # births
        current_params$beta * s[t] * i[t] + # infection
        current_params$omega * r[t] - # lost immunity
        current_params$vac * s[t] - # vaccination
        current_params$mu_d * s[t] # death
      
      e[t + 1] <- e[t] +
        current_params$beta * s[t] * i[t] - # infection
        current_params$sigma * e[t] - # latency
        current_params$mu_d * e[t] # death
      
      i[t + 1] <- i[t] +
        current_params$sigma * e[t] - # latency
        current_params$gamma * i[t] - # recovery
        (current_params$mu_d + current_params$alpha) * i[t] # death
      
      r[t + 1] <- r[t] +
        current_params$gamma * i[t] + # recovery
        current_params$vac * s[t] - # vaccination
        current_params$omega * r[t] - # lost immunity
        current_params$mu_d * r[t] # death
      
      d[t + 1] <- d[t] +
        current_params$mu_d * s[t] + # death of susceptible
        current_params$mu_d * e[t] + # death of exposed
        (current_params$mu_d + current_params$alpha) * i[t] + # death of infected
        current_params$mu_d * r[t] # death of recovered
      
      infected_sum[t + 1] <- infected_sum[t] +
        current_params$beta * s[t] * i[t]
    }
    
    #return a list containing dataframe with respect to each [n] of the params sample 
    df[[n]] <- list(data = data.frame(
      "time" = 1:timesteps,
      "susceptible" = s,
      "exposed" = e,
      "infectious" = i,
      "recovered" = r,
      "died" = d,
      "infected sum" = infected_sum
    )) 
  }
  
  return(df)
}





################################################################################


  
param_sample_prob <- function(
    params_df,
    iter = 100 #numbers of param sample to be produced
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
  
  #create an empty list
  output_list <- list() 
  
  # loop over the defined iteration to generate multiple list with randomly generated parameters' values 
 for (i in 1:iter) {
    sampled_params_df <- params_df %>%
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
    
    output_list[[i]] <- setNames(as.list(sampled_params_df$Value),    # Convert vectors to named list
                                 params_df$Parameter)
}

  
  return(output_list)

}

################################################################################


plot_SEIR_prob =function(
    df=SEIR_df){
  #create a long format table 
  combined_df = bind_rows(df)
  
 
  # Add a group column
  combined_df =  combined_df %>% mutate(group = rep(1:length(df), each = nrow(df[[1]][["data"]])))
  
  # Compute median and CI for each time point
  summary_stats <- combined_df %>%
    group_by(data$time) %>%
    summarise(
      median_susceptible = median(data$susceptible),
      CI_lower = quantile(data$susceptible, 0.025),
      CI_upper = quantile(data$susceptible, 0.975)
    )
  
  
  s = ggplot() +
    geom_line(data = combined_df, aes(x = data$time, y = data$susceptible, group = as.factor(group)), color = "grey", alpha = 0.9) +
    # Add median line
    geom_line(data = summary_stats, aes(x = `data$time`, y = median_susceptible), color = "blue", size =0.9) +
    # Add CI lines
    geom_line(data = summary_stats, aes(x = `data$time`, y = CI_lower), linetype = "dotted", color = "red", size = 0.9) +
    geom_line(data = summary_stats, aes(x = `data$time`, y = CI_upper), linetype = "dotted", color = "red", size = 0.9)+
    geom_ribbon(data = summary_stats, aes(x = `data$time`, ymax = CI_upper, ymin = CI_lower), fill = "red", alpha=0.05)+
    labs(title = "Median and 95%CI of simulated susceptible",
         x= "day",
         y="Fraction"
         )+
    theme_minimal()+
    theme(plot.title = element_text(hjust=0.5))
    
    
  
  return(s)

  
}



plot_SEIR_prob(SEIR_df)



