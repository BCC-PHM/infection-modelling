
library(readxl)
library(tidyverse)
library(rjags)


hospital_data =read_excel(
  paste(
    "//svwvap1126.addm.ads.brm.pri/PHSensitive$/Intelligence/2. Requests/REQ3280 - HES Inpatients data on polio, flu and rubella for Justin",
    "/REQ3280_Inpatients_data_10years_on_polio_flu_and_rubella_BirminghamResidents.xlsx",
    sep = ""
  )
)

#filter Polio only and only data in 2022/23
Poli_data  = hospital_data %>% 
  filter(Condition == "Polio" & Year == "2022/2023")

#Convert the dates
Poli_data$Date_of_Admission = as.Date(Poli_data$Date_of_Admission, format = "%d%m%Y")


Poli_data = Poli_data %>% 
  group_by(Date_of_Admission) %>% 
  summarise(case = n())

# Create an empty data frame with the date column
Polio_empty_FY = data.frame(Date_of_Admission = seq.Date(from = as.Date("2022-04-01"), 
                                                         to = as.Date("2023-03-31"), by="day"), 
                            stringsAsFactors = FALSE)


#left join the data
Poli_combined = Polio_empty_FY %>% left_join(Poli_data, by = "Date_of_Admission")


#############################################################################################################
#Prepare the data for Bayes model 
Poli_combined$case = ifelse(is.na(Poli_combined$case), 0, Poli_combined$case)

# Prepare data
time=nrow(Poli_combined)            # Time points
N <- 1144900        #population size 
case = Poli_combined$case  # Observed cases data
# Initial number of infected individuals

data <- list(
  time = time,
  N = N,
  case = case
  
)

# Write model specification
model_code <- "
model {
  # Prior distributions for parameters
  sigma =0.25     # Rate of transition from exposed to infectious dgamma(7.617034, 53)
  gamma =0.09  # Rate of recovery/removal dgamma(14, 201) 
  R_0 ~ dunif(0,10)         #dgamma(36,6)
  
  #contact rate
   beta = R_0*gamma # Transmission rate
  
  #death from infection
  alpha = 0
  
  #fixed birth and death
  birth <- 0.000035
  death <- 0.000034
  
  # Immunity loss rate
  omega = 0     # Rate of loss of immunity
  
  #fixed vacinnation rate
  vacc_rate = 0.00003373568
  
  #hospitalisation rate 
  hosp =0.00465
  hosprecover =0.09

  I_init = %s

  # Initial conditions
  S[1] <- N-I_init-858675
  E[1] <- 0
  I[1] <- I_init
  H[1] = 0
  R[1] <- 858675 # 75 percent of people are vaccinated in the population
  total_infected[1] = I_init

  # SEIR dynamics
  for (t in 2:time) {
    dS[t] <- birth * N - beta * S[t-1] * I[t-1] / N - death * S[t-1]+ omega * R[t-1] - vacc_rate * S[t-1]
    dE[t] <- beta * S[t-1] * I[t-1] / N - sigma * E[t-1] - death * E[t-1]
    dI[t] <- sigma * E[t-1] - gamma * I[t-1] - (death+alpha) * I[t-1] - hosp*I[t-1]
    dH[t] = hosp*I[t-1]- hosprecover*H[t-1]-death*H[t-1]
    dR[t] <- gamma * I[t-1] - death * R[t-1]- omega * R[t-1] + vacc_rate * S[t-1]+hosprecover*H[t-1]

    S[t] <- S[t-1] + dS[t]
    E[t] <- E[t-1] + dE[t]
    I[t] <- I[t-1] + dI[t]
    H[t] = H[t-1]+ dH[t]
    R[t] <- R[t-1] + dR[t]
    total_infected[t] = total_infected[t-1]+beta * S[t-1] * I[t-1] / N
    
    }

    for (t in 2:time) {
    case[t] ~ dpois(H[t])
  }
}"


I_inits <- 10^(seq(0, 3, 0.5))

results_list <- c()
for (i in 1:length(I_inits)) {
  I_init_i <- I_inits[[i]]
  
  # Add I(t=1) value into model
  model_code_i <- sprintf(model_code, as.character(I_init_i))

  # Compile model
  model <- jags.model(textConnection(model_code_i), data = data,  n.chains = 3, n.adapt = 1000)
  
  # sample from the posterior
  mcmc = coda.samples(model = model, variable.names = c("R_0"), n.iter = 20000, thin = 2)
  
  # Extract posterior chain
  polio_posterior_chain1 = as.data.frame(mcmc[[1]])
  
  # Calculate quantiles
  quantiles <- quantile(polio_posterior_chain1$R_0, c(0.025,0.5, 0.975))
  
  # Save R0 quantiles for this I(t=1) value
  result_i <- data.frame(
    I_init =  I_init_i,
    R0_median = quantiles[[2]],
    R0_lower = quantiles[[1]],
    R0_upper = quantiles[[3]]
  )
  
  # Store to list
  results_list[[i]] <- result_i
}

# Combine results into single table
results <- data.table::rbindlist(results_list)

# Plot R0 estimate vs I(t=1) with error bars
plt <- ggplot(results, aes(x = I_init, y = R0_median)) +
  geom_point() +
  geom_errorbar(aes(ymin=R0_lower, ymax=R0_upper,), width=.2,
                position=position_dodge(.9)) +
  theme_bw() +
  labs(
    y = "R0 Estimate",
    x = TeX("Initial Number of People Infected $I(t=1)$")
  ) +
  scale_x_log10()
plt
# Save figure
ggsave("../output/figures/polio_init_vs_r0_pred.png", plt, 
       width = 7, height = 5, dpi = 300)