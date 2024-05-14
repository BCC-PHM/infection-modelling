library(dplyr)
library(ggplot2)

#https://github.com/XuelongSun/Dynamic-Model-of-Infectious-Diseases/blob/master/SIR.ipynb

SEIR_model <- function(
  params,
  inits,
  timesteps = 100
  ) {
  
  # Check and fill missing elements with zero
  for (elem in c("lamda", "gamma", "sigma", "vac")) {
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
  
  # initial infective people
  i[1] = inits$i
  s[1] = inits$s
  e[1] = inits$e
  
  for (t in 1:(T-1)){
    s[t + 1] = s[t] - lamda * s[t] * i[t] - vac * s[t]
    e[t + 1] = e[t] + lamda * s[t] * i[t] - sigma * e[t]
    i[t + 1] = i[t] + sigma * e[t] - gamma * i[t]
    r[t + 1] = r[t] + gamma * i[t] + vac * s[t]
  }
  
  df <- data.frame(
    "time" = 1:T,
    "susceptible" = s,
    "exposed" = e,
    "infectious" = i,
    "recovered" = r
  ) 
  
  return(df)
}

params <- list(
  # contact rate
  lamda = 0.5,
  # recover rate
  gamma = 0.0821,
  # exposed period
  sigma = 1 / 4,
  # vaccination rate
  vac = 0.006
)

inits <- list(
  # susceptiable ratio
  s = 0.99,
  # exposed ratio
  e = 0,
  # infective ratio
  i = 0.01,
  # recovered ratio
  r = 0
)

# simuation Time / Day
tsteps = 170

df <- SEIR_model(params, inits, tsteps) %>%
  reshape2::melt(id.vars = "time", 
       variable.name = 'Group',
       value.name = 'Number')

plt <- ggplot(df, aes(x = time, y = Number, col = Group)) +
  geom_line(linewidth = 1.2) +
  theme_bw() +
  ggtitle("Vaccination rate = 0.002") + theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave("../output/figures/SEIR-test-new.png", plt, 
       width = 5, height = 3)