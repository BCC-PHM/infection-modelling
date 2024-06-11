# Assume fixed values for 1/gamma, 1/sigma. Hence calculate R0 
# based on hospitalisation data assuming total number of infections is less than
# 100 times the number of hospitalisations

library(dplyr)
library(ggplot2)

source("infection-modelling.R")

set.seed(0)

post_chain <- read.csv("../data/Polio/polio_posterior_BSOL 1.csv") %>% 
  tidyr::pivot_longer(
    cols = c("R_0","eta"),
    names_to='Parameter',
    values_to='Value') %>%
  group_by(Parameter) %>%
  summarize(
    Median = quantile(Value, 0.5),
    lowerCI = quantile(Value, 0.025),
    upperCI = quantile(Value, 0.965)
  )
post_chain

params <- list(
  "R0" = 4.98, # 5, # From Bayesian SEIHR model
  "1/gamma" = 1/0.09, 
  "1/sigma" = 4,
  "1/omega" = 1e99,
  "vac" = 3.374e-5,
  "alpha" = 0,
  "mu_b" = 3.464e-5,
  "mu_d" = 3.449e-5,
  "hosp" = 0.000107# 3e-5
)

inits <- list(
  r = 0.8,
  e = 0.0000838,
  s = 1 - 0.8 - 0.0000838 - 0.000233,
  i = 0.000233
)
# 
# inits <- list(
#   r = 0.89,
#   e = 5.0e-5,
#   s = 1 - 0.89 - 5.07e-5 - 2.0e-4,
#   i = 7.0e-9,
#   a = 0
# )

result <- SEIR_model(
    params,
    inits,
    timesteps = 100000,
    t0 = 1
) 
result_long <- result %>%
  select(-c("Died", "Infected_sum")) %>% 
  tidyr::pivot_longer(
  cols = c("Susceptible","Exposed","Infectious","Recovered"),
  names_to='Group',
  values_to='Fraction')

ggplot(result_long, aes(x = Time, y = Fraction, color = Group)) +
  geom_line(linewidth = 1.2) +
  theme_bw() +
  scale_y_log10()

result_long %>% 
  filter(Time == max(result_long$Time))

average_hosp_rate <- result %>% 
  summarise(avs_hosp_rate = sum(Admissions * 1144900)/max(result_long$Time)*365)
average_hosp_rate




hosp <- 10^(seq(-6,-2,0.4))
hosp_rate <- list()
for (i in 1:length(hosp)) {
  params$hosp <- hosp[[i]]
  
  result <- SEIR_model(
    params,
    inits,
    timesteps = 100000,
    t0 = 1
  ) 
  
  average_hosp_rate <- result %>% 
    summarise(avs_hosp_rate = sum(Admissions * 1144900)/max(result_long$Time)*365)
  
  hosp_rate[[i]] <- average_hosp_rate$avs_hosp_rate
}

hosp_rate_sim <- data.frame(
  hosp_param = hosp,
  hosp_rate = unlist(hosp_rate)
)

ggplot(hosp_rate_sim, aes(hosp_param, hosp_rate)) +
  geom_point() +
  theme_bw() +
  scale_y_log10()+
  scale_x_log10() +
  geom_hline(yintercept = 3) +
  labs(
    x = "Hospital parameter \\eta",
    y = "Average yearly hospital admissions"
  )
