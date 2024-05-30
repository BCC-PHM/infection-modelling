# Assume fixed values for 1/gamma, 1/sigma. Hence calculate R0 
# based on hospitalisation data assuming total number of infections is less than
# 100 times the number of hospitalisations

library(dplyr)
library(ggplot2)

source("infection-modelling.R")

set.seed(0)

post_chain <- read.csv("../data/Polio/polio_posterior_chain.csv") %>% 
  pivot_longer(
    cols = c("R_0","gamma","eta"),
    names_to='Parameter',
    values_to='Value') %>%
  group_by(Parameter) %>%
  summarize(
    Value = median(Value)
  )
post_chain

params <- list(
  "R0" = 6, # From Bayesian SEIHR model
  "1/gamma" = 10.75,
  "1/sigma" = 4,
  "1/omega" = 1e150,
  "vac" = 3.374e-5,
  "alpha" = 0,
  "mu_b" = 3.464e-5,
  "mu_d" = 3.449e-5,
  "hosp" = 0.0004
)

inits <- list(
  r = 0.89,
  e = 5.0e-5,
  s = 1 - 0.89 - 5.07e-5 - 2.0e-4,
  i = 7.0e-9,
  h = 0
)



#sum(unlist(inits))

result <- SEIHR_model(
    params,
    inits,
    timesteps = 365,
    t0 = 1
) 
result_long <- result %>%
  select(-c("Died", "Infected_sum")) %>% 
  pivot_longer(
  cols = c("Susceptible","Exposed","Infectious","Hospital","Recovered", "Admissions"),
  names_to='Group',
  values_to='Fraction')

ggplot(result_long, aes(x = Time, y = Fraction, color = Group)) +
  geom_line(linewidth = 1.2) +
  theme_bw()

result_long %>% 
  filter(Time == 365)

average_hosp_rate <- result %>% 
  summarise(avs_hosp_rate = sum(Hospital * 1144900))
average_hosp_rate
