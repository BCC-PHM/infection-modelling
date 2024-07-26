
library(readxl)
library(tidyverse)
library(rjags)
library(egg)
library(bayesplot)


hospital_data =read_excel("~/R projects/infection-modelling/data/REQ3280_Inpatients_data_10years_on_polio_flu_and_rubella_BirminghamResidents.xlsx")

#filter Polio only and only data in 2022/23
Poli_data  = hospital_data %>% filter(Condition == "Polio" & Year == "2022/2023")

#Convert the dates
Poli_data$Date_of_Admission = as.Date(Poli_data$Date_of_Admission, format = "%d%m%Y")


Poli_data = Poli_data %>% group_by(Date_of_Admission) %>% summarise(case = n())

# Create an empty data frame with the date column
Polio_empty_FY = data.frame(Date_of_Admission =seq.Date(from = as.Date("2022-04-01"), to = as.Date("2023-03-31"), by="day"), stringsAsFactors = FALSE)


#left join the data
Poli_combined = Polio_empty_FY %>% left_join(Poli_data, by = "Date_of_Admission")



#############################################################################################################
#Prepare the data for Bayes model 
Poli_combined$case = ifelse(is.na(Poli_combined$case), 0, Poli_combined$case)
# Set seed for reproducibility
set.seed(123)

# Create a vector of 365 zeros
case <- rep(0, 365)

# Randomly select 3 positions to assign the value 1
positions <- sample(1:365, 3)

# Assign the value 1 to the selected positions
case[positions] <- 1


# Prepare data
time=nrow(Poli_combined)            # Time points
N <- 1361159       #population size 
# case = Poli_combined$case  # Observed cases data
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
  sigma =0.25      # Rate of transition from exposed to infectious ~dgamma(7.617034, 53) 0.25
  gamma =0.09  # Rate of recovery/removal ~dgamma(14, 201) 0.09
  R_0 ~ dunif(0,10)         #dgamma(36,6)
  
  #contact rate
   beta = R_0*gamma # Transmission rate
  
  #death from infection
  alpha = 0
  
  #fixed birth and death
  birth <- 0.000035
  death <- 0.000034
  
  # Immunity loss rate
  omega =0     # Rate of loss of immunity
  
  #fixed vacinnation rate
  vacc_rate = 0.00003373568
  
  #hospitalisation rate 
  eta ~ dunif(0,0.1)               #0.00465
  hosprecover =0.09

  # Initial conditions
  S[1] <- 1361159 -1088927-108.89272-27.22318
  E[1] <- 108.89272
  I[1] <- 27.22318
  A[1] = 0
  R[1] <- 1088927 #80% of people are vaccinated in the population
  total_infected[1] = I[1]


  # SEIR dynamics
  for (t in 2:time) {
    dS[t] <- birth * N - beta * S[t-1] * I[t-1] / N - death * S[t-1]+ omega * R[t-1] - vacc_rate * S[t-1]
    dE[t] <- beta * S[t-1] * I[t-1] / N - sigma * E[t-1] - death * E[t-1]
    dI[t] <- sigma * E[t-1] - gamma * I[t-1] - (death+alpha) * I[t-1] - eta*I[t-1]
    dA[t] = eta*I[t-1]-A[t-1]
    dR[t] <- gamma * I[t-1] - death * R[t-1]- omega * R[t-1] + vacc_rate * S[t-1]+A[t-1]

    S[t] <- S[t-1] + dS[t]
    E[t] <- E[t-1] + dE[t]
    I[t] <- I[t-1] + dI[t]
    A[t] = A[t-1]+ dA[t]
    R[t] <- R[t-1] + dR[t]
    total_infected[t] = total_infected[t-1]+beta * S[t-1] * I[t-1] / N
   
    }

    for (t in 2:time) {
    case[t] ~ dpois(A[t])
  }
}"


# Compile model
model <- jags.model(textConnection(model_code), data = data,  n.chains = 3, n.adapt = 1000)


##sample from the posterior
mcmc = coda.samples(model = model, variable.names = c("R_0", "gamma", "eta"), n.iter = 20000, thin = 2)



plot(mcmc, density = TRUE)

autocorr.plot(mcmc)


##using bayesplot

#set colour
color_scheme_set("mix-blue-red")

#trace plot
traceplot = mcmc_trace(mcmc, pars = c("R_0", "eta"),
           facet_args = list(nrow=2))

#density plot

densityplot = mcmc_dens_overlay(mcmc, pars = c("R_0", "eta"),
                  facet_args = list(nrow=2))

mcmc_acf(mcmc, pars = c("R_0", "gamma", "eta"))

ggarrange(traceplot,densityplot,
          nrow = 1, ncol = 2)


autocorr(mcmc)






polio_posterior_chain1 = data.frame(mcmc[[1]])

admission_df = data.frame(n=1:10000, admission = NA)
for(row in 1:10000){
  

  
   value= sum(polio_posterior_chain1[row, 1:365])
  

  admission_df[row,2] = value
  
 
  
}


ggplot(admission_df, aes(x=admission))+geom_density()

write.csv(polio_posterior_chain1, "polio_posterior_BSOL.csv")



quantile(polio_posterior_chain1$R_0, c(0.025,0.5,0.975))


quantile(polio_posterior_chain1$eta, c(0.025,0.5,0.975))


quantile(polio_posterior_chain1$gamma, c(0.025,0.5,0.975))














