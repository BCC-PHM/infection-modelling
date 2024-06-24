library(dplyr)
library(ggplot2)
library(egg)
library(gridExtra)

source("gamma_params.R")


SEIR_model_prob <- function(
    params,
    inits,
    timesteps = 100,
    print_progress = TRUE 
    
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
    vaccinated_sum = rep(0, timesteps)   #calculate the sum of vaccinated people
    
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
      
      vaccinated_sum[t+1] = vaccinated_sum[t]+ current_params$vac*s[t] # Update vaccinated count
    }
    
    #return a list containing dataframe with respect to each [n] of the params sample 
    df[[n]] <- list(data = data.frame(
      "time" = 1:timesteps,
      "susceptible" = s,
      "exposed" = e,
      "infectious" = i,
      "recovered" = r,
      "died" = d,
      "infected sum" = infected_sum,
      "vaccinated_sum" = vaccinated_sum
    )) 
    
    # Define progress intervals
    progress_intervals <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
    
    # Calculate progress percentage
    progress <- n / length(params)
    
    # Check if progress reaches certain intervals
    if (print_progress && progress %in% progress_intervals) {
      cat("Iteration", n, "completed.\n")
    }
    
  }
  
  return(df)
}





################################################################################


  
param_sample_prob <- function(
    params_df,
    vac=0,
    iter = 100,
    print_progress = TRUE 
    #numbers of param sample to be produced
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
    
    output_list[[i]]$vac = vac
    
    # Define progress intervals
    progress_intervals <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
    
    # Calculate progress percentage
    progress <- i / iter
    
    # Check if progress reaches certain intervals
    if (print_progress && progress %in% progress_intervals) {
      cat("Iteration", i, "completed.\n")
    }
 }

  return(output_list)

}

################################################################################


plot_SEIR_prob =function(
    df=SEIR_df,
    plot = "all"){
  #create a long format table 
  combined_df = bind_rows(df)
  
 
  # Add a group column
  combined_df =  combined_df %>% mutate(group = rep(1:length(df), each = nrow(df[[1]][["data"]])))
  
  ########################################################################
  # Compute median and CI for each time point for s 
  summary_stats_s <- combined_df %>%
    group_by(data$time) %>%
    summarise(
      median = median(data$susceptible),
      CI_lower = quantile(data$susceptible, 0.025),
      CI_upper = quantile(data$susceptible, 0.975)
    )
  
  
  s = ggplot() +
    geom_line(data = combined_df, aes(x = data$time, y = data$susceptible, group = as.factor(group)), color = "grey", alpha = 0.9) +
    # Add median line
    geom_line(data = summary_stats_s, aes(x = `data$time`, y = median), color = "red", size =0.9) +
    # Add CI lines
    geom_line(data = summary_stats_s, aes(x = `data$time`, y = CI_lower), linetype = "dotted", color = "red", size = 0.9) +
    geom_line(data = summary_stats_s, aes(x = `data$time`, y = CI_upper), linetype = "dotted", color = "red", size = 0.9)+
    geom_ribbon(data = summary_stats_s, aes(x = `data$time`, ymax = CI_upper, ymin = CI_lower), fill = "red", alpha=0.05)+
    labs(title = "Median and 95%CI of simulated susceptible",
         x= "day",
         y="Fraction"
         )+
    theme_minimal()+
    theme(plot.title = element_text(hjust=0.5))
  
  ########################################################################    

  # Compute median and CI for each time point for e
  summary_stats_e <- combined_df %>%
    group_by(data$time) %>%
    summarise(
      median = median(data$exposed),
      CI_lower = quantile(data$exposed, 0.025),
      CI_upper = quantile(data$exposed, 0.975)
    )
  
  
  e = ggplot() +
    geom_line(data = combined_df, aes(x = data$time, y = data$exposed, group = as.factor(group)), color = "grey", alpha = 0.9) +
    # Add median line
    geom_line(data = summary_stats_e, aes(x = `data$time`, y = median), color = "green", size =0.9) +
    # Add CI lines
    geom_line(data = summary_stats_e, aes(x = `data$time`, y = CI_lower), linetype = "dotted", color = "green", size = 0.9) +
    geom_line(data = summary_stats_e, aes(x = `data$time`, y = CI_upper), linetype = "dotted", color = "green", size = 0.9)+
    geom_ribbon(data = summary_stats_e, aes(x = `data$time`, ymax = CI_upper, ymin = CI_lower), fill = "green", alpha=0.05)+
    labs(title = "Median and 95%CI of simulated exposed",
         x= "day",
         y="Fraction"
    )+
    theme_minimal()+
    theme(plot.title = element_text(hjust=0.5))
  

  ########################################################################
  
  # Compute median and CI for each time point for i
  summary_stats_i <- combined_df %>%
    group_by(data$time) %>%
    summarise(
      median = median(data$infectious),
      CI_lower = quantile(data$infectious, 0.025),
      CI_upper = quantile(data$infectious, 0.975)
    )
  
  
  i = ggplot() +
    geom_line(data = combined_df, aes(x = data$time, y = data$infectious, group = as.factor(group)), color = "grey", alpha = 0.9) +
    # Add median line
    geom_line(data = summary_stats_i, aes(x = `data$time`, y = median), color = "blue", size =0.9) +
    # Add CI lines
    geom_line(data = summary_stats_i, aes(x = `data$time`, y = CI_lower), linetype = "dotted", color = "blue", size = 0.9) +
    geom_line(data = summary_stats_i, aes(x = `data$time`, y = CI_upper), linetype = "dotted", color = "blue", size = 0.9)+
    geom_ribbon(data = summary_stats_i, aes(x = `data$time`, ymax = CI_upper, ymin = CI_lower), fill = "blue", alpha=0.05)+
    labs(title = "Median and 95%CI of simulated infectious",
         x= "day",
         y="Fraction"
    )+
    theme_minimal()+
    theme(plot.title = element_text(hjust=0.5))
  
 
  
  ########################################################################
 
  
  # Compute median and CI for each time point for r
  summary_stats_r <- combined_df %>%
    group_by(data$time) %>%
    summarise(
      median = median(data$recovered),
      CI_lower = quantile(data$recovered, 0.025),
      CI_upper = quantile(data$recovered, 0.975)
    )
  
  
  r = ggplot() +
    geom_line(data = combined_df, aes(x = data$time, y = data$recovered, group = as.factor(group)), color = "grey", alpha = 0.9) +
    # Add median line
    geom_line(data =summary_stats_r, aes(x = `data$time`, y = median), color = "purple", size =0.9) +
    # Add CI lines
    geom_line(data =summary_stats_r, aes(x = `data$time`, y = CI_lower), linetype = "dotted", color = "purple", size = 0.9) +
    geom_line(data =summary_stats_r, aes(x = `data$time`, y = CI_upper), linetype = "dotted", color = "purple", size = 0.9)+
    geom_ribbon(data =summary_stats_r, aes(x = `data$time`, ymax = CI_upper, ymin = CI_lower), fill = "purple", alpha=0.05)+
    labs(title = "Median and 95%CI of simulated recovered",
         x= "day",
         y="Fraction"
    )+
    theme_minimal()+
    theme(plot.title = element_text(hjust=0.5))
  
  
  ########################################################################
  #combine all the line into a single graph 
  # Combine all summary statistics data frames into a single data frame
  combined_summary_stats <- bind_rows(
    mutate(summary_stats_s, group = "Susceptible"),
    mutate(summary_stats_e, group = "Exposed"),
    mutate(summary_stats_i, group = "Infectious"),
    mutate(summary_stats_r, group = "Recovered")
  )
  
  combined_summary_stats$group = factor(combined_summary_stats$group, levels = c("Susceptible",
                                                                                 "Exposed",
                                                                                 "Infectious",
                                                                                 "Recovered"))
  
  # Plot using color inside aes() and creating lines and ribbons for each group
 all=  ggplot(combined_summary_stats, aes(x = `data$time`, y = median, color = group, fill = group, ymin = CI_lower, ymax = CI_upper)) +
    geom_line(size = 0.9) +
    geom_ribbon(alpha = 0.05, linetype = 0) +
    labs(
      title = "SEIR model",
      x = "Day",
      y = "Fraction"
    ) +
   scale_color_manual(values = c("red", "green", "blue", "purple")) +
   scale_fill_manual(values = c("red", "green", "blue", "purple")) +
   theme_minimal()+
   theme(plot.title = element_text(hjust = 0.5))
    
  
 ##########################################################################
 # create all plot with empty  legend
 all_empty_legend=  ggplot(combined_summary_stats, aes(x = `data$time`, y = median, color = group, fill = group, ymin = CI_lower, ymax = CI_upper)) +
   geom_line(size = 0.9) +
   geom_ribbon(alpha = 0.05, linetype = 0) +
   labs(
     title = "SEIR model",
     x = "Day",
     y = "Fraction"
   ) +
   scale_color_manual(values = c("red", "green", "blue", "purple")) +
   scale_fill_manual(values = c("red", "green", "blue", "purple")) +
   theme_minimal()+
   theme(plot.title = element_text(hjust = 0.5),
         legend.position = "none")
 
 # create all plot with  legend only
 legend= ggplot(combined_summary_stats, aes(x = `data$time`, y = median, color = group, fill = group, ymin = CI_lower, ymax = CI_upper)) +
   geom_line(size = 0.9) +
   geom_ribbon(alpha = 0.05, linetype = 0) +
   lims(x = c(0,0), y = c(0,0))+
   scale_color_manual(values = c("red", "green", "blue", "purple")) +
   scale_fill_manual(values = c("red", "green", "blue", "purple")) +
   theme_void()+
   theme(legend.position = c(0.5,0.5),
         legend.key.size = unit(1, "cm"),
         legend.text = element_text(size =  12),
         legend.title = element_text(size = 15, face = "bold"))+
   guides(colour = guide_legend(override.aes = list(size=8)))


 
 ##stick all the graph together
 
 # #arrange plot
 combined_all <- ggarrange(all_empty_legend, s, e, i, r,legend,
                           ncol =3, nrow=2
                           )


 if (plot == "SEIR") {
   return(all)
 } else if (plot == "s") {
   return(s)
 } else if (plot == "e") {
   return(e)
 } else if (plot == "i") {
   return(i)
 } else if (plot == "r") {
   return(r)
 } else {
   return(combined_all)
 }
 
 
 
}









