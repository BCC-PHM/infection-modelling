library(dplyr)
library(ggplot2)

#https://github.com/XuelongSun/Dynamic-Model-of-Infectious-Diseases/blob/master/SIR.ipynb

SEIR_model <- function(
  params,
  inits,
  timesteps = 100
  ) {
  
  # Check and fill missing elements with zero
  for (elem in c("beta", "gamma", "sigma", "vac", "omega", "mu", "alpha")) {
    if (!(elem %in% names(params))) {
      params[[elem]] <- 0
    }
  }
  
  # Check and fill missing elements with zero
  for (elem in c("s", "e", "i", "r")) {
    if (!(elem %in% names(inits))) {
      inits[[elem]] <- 0
    }
  }
  
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
  
  for (t in 1:(T-1)){
    
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
    "time" = 1:T,
    "susceptible" = s,
    "exposed" = e,
    "infectious" = i,
    "recovered" = r,
    "died" = d,
    "infected sum" = infected_sum
  ) 
  
  return(df)
}

params <- list(
  # contact rate
  beta = 0.5,
  # recover rate
  gamma = 0.0821,
  # exposed period
  sigma = 1 / 4,
  # vaccination rate
  vac = 0.006,
  # immunity loss rate
  omega = 0.0001,
  # daily birth/death rate
  mu = 0.001,
  # infected loss of life
  alpha = 0.00000001
)

inits <- list(
  # susceptiable ratio
  s = 0.3,
  # exposed ratio
  e = 0,
  # infective ratio
  i = 0.01,
  # recovered ratio
  r = 0.7
)

# simuation Time / Day
tsteps = 170

df <- SEIR_model(params, inits, tsteps) %>%
  reshape2::melt(id.vars = "time", 
       variable.name = 'Group',
       value.name = 'Fraction')

plt <- ggplot(df, aes(x = time, y = Fraction, col = Group)) +
  geom_line(linewidth = 1.1) +
  theme_bw() +
  ggtitle("Vaccination rate = 0.002") + theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave("../output/figures/SEIR-test-new.png", plt, 
       width = 5, height = 3)