
library(readxl)
library(tidyverse)
library(rjags)
library(egg)
library(bayesplot)
library(stringr)

# hospital_data =read_excel(
#   "//svwvap1126.addm.ads.brm.pri/PHSensitive$/Intelligence/PHM/Polio/Admissions with Polio A80.xlsx",
#   sheet = "Top 6 causes"
#   )

#filter Polio only and only data in 2022/23
Poli_data  = hospital_data %>% 
  filter(FYEAR %in% c(2021, 2122, 2223)) %>%
  mutate(
    date_str = as.character(epistart),
    date_str = case_when(
      str_length(date_str) == 8 ~ date_str,
      str_length(date_str) == 7 ~ paste0("0", date_str)
    ),
    Date_of_Admission = as.Date(date_str, format = "%d%m%Y")
  )

#Convert the dates
Poli_data$Date_of_Admission = as.Date(
  as.character(Poli_data$epistart), format = "%d%m%Y")


Poli_data <- Poli_data %>% 
  group_by(Date_of_Admission) %>% 
  summarise(case = sum(COUNT))

# Create an empty data frame with the date column
Polio_empty_FY <- 
  data.frame(
    Date_of_Admission = seq.Date(from = as.Date("2020-04-01"),
                                 to = as.Date("2023-03-31"), by="day"), 
    stringsAsFactors = FALSE)


#left join the data
Poli_combined <- Polio_empty_FY %>% 
  left_join(Poli_data, by = "Date_of_Admission")


#############################################################################################################
#Prepare the data for Bayes model 
Poli_combined$case = ifelse(is.na(Poli_combined$case), 0, Poli_combined$case)

# Prepare data
time=nrow(Poli_combined)            # Time points
N <- 57106000        #population size
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
  sigma = 0.25      # Rate of transition from exposed to infectious ~dgamma(7.617034, 53) 0.25
  gamma = 0.093  #  0.09
  epsilon_0 ~ dunif(0, 1e-5)         #dgamma(36,6)
  
  #fixed birth and death
  birth <- 0.000035
  death <- 0.000034
  
  #fixed vacinnation rate
  vacc_rate = 0.00003373568
  
  #hospitalisation rate 
  eta ~ dunif(0, 1e-3)   # 0.00465

  # Initial conditions
  I[1] <- 0 # No one needs to be initially infected since it's imported
  R[1] <- N*0.75 #7 5% of people are vaccinated in the population
  S[1] <- N - R[1]
  E[1] <- 0
  H[1] <- 0
  total_infected[1] = 0
  admission[1] = 0

  # SEIR dynamics
  for (t in 2:time) {
    dS[t] <- birth * N - epsilon_0 * S[t-1] - death * S[t-1] - vacc_rate * S[t-1]
    dE[t] <- epsilon_0 * S[t-1] - sigma * E[t-1] - death * E[t-1]
    dI[t] <- sigma * E[t-1] - gamma * I[t-1] - death * I[t-1] - eta*I[t-1]
    
    dH[t] = eta * sigma * E[t-1]
    
    dR[t] <- gamma * I[t-1] - death * R[t-1] + vacc_rate * S[t-1] 

    S[t] <- S[t-1] + dS[t]
    E[t] <- E[t-1] + dE[t]
    I[t] <- I[t-1] + dI[t]
    H[t] = dH[t]
    R[t] <- R[t-1] + dR[t]
    total_infected[t] = total_infected[t-1] + epsilon_0 * S[t-1]
    
    admission[t] = eta*I[t-1]
    }

    for (t in 2:time) {
    case[t] ~ dpois(admission[t])
  }
}"

# Compile model
model <- jags.model(textConnection(model_code), data = data,  n.chains = 3, n.adapt = 1000)

##sample from the posterior
mcmc = coda.samples(model = model, variable.names = c("epsilon_0", "eta"), n.iter = 20000, thin = 2)



plot(mcmc, density = TRUE)

polio_posterior_chain = as.data.frame(mcmc[[1]])


quantile(polio_posterior_chain$epsilon_0, c(0.025,0.5, 0.975))
quantile(polio_posterior_chain$eta, c(0.025,0.5, 0.975))

##using bayesplot

#set colour
color_scheme_set("mix-blue-red")

#trace plot
traceplot = mcmc_trace(mcmc, pars = c("epsilon_0", "eta"),
           facet_args = list(nrow=3))

#density plot

densityplot = mcmc_dens_overlay(mcmc, pars = c("epsilon_0", "eta"),
                  facet_args = list(nrow=3))

ggarrange(traceplot,densityplot,
          nrow = 1, ncol = 2)

write.csv(polio_posterior_chain, "../data/polio/polio_posterior_chain_import_model_2021to2223.csv")
