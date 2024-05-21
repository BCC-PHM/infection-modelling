library(rjags)

# Prepare data
time=30             # Time points
N <- 1144900        #population size 
case = runif(time ,0,10)  # Observed cases data
 # Initial number of infected individuals

data <- list(
  time = time,
  N = N,
  case = round(runif(time ,0,1),0)
)



# Write model specification
model_code <- "
model {
  # Prior distributions for parameters
  sigma ~ dgamma(7.617034, 53)     # Rate of transition from exposed to infectious
  gamma ~ dgamma(14, 201)     # Rate of recovery/removal
  R_0 ~ dunif(1,10)
  
  #contact rate
   beta = R_0*gamma # Transmission rate
  
  #death from infection
  alpha = 0
  
  #fixed birth and death
  birth <- 0.000035
  death <- 0.000034
  
  # Immunity loss rate
  omega ~ dgamma(2,360)        # Rate of loss of immunity
  
  #fixed vacinnation rate
  vacc_rate <- 0

  # Initial conditions
  S[1] <- N-1
  E[1] <- 0
  I[1] <- 1
  R[1] <- 0

  # SEIR dynamics
  for (t in 2:time) {
    dS[t] <- birth * N - beta * S[t-1] * I[t-1] / N - death * S[t-1]+ omega * R[t-1] - vacc_rate * S[t-1]
    dE[t] <- beta * S[t-1] * I[t-1] / N - sigma * E[t-1] - death * E[t-1]
    dI[t] <- sigma * E[t-1] - gamma * I[t-1] - (death+alpha) * I[t-1]
    dR[t] <- gamma * I[t-1] - death * R[t-1]- omega * R[t-1] + vacc_rate * S[t-1]

    S[t] <- S[t-1] + dS[t]
    E[t] <- E[t-1] + dE[t]
    I[t] <- I[t-1] + dI[t]
    R[t] <- R[t-1] + dR[t]}

    for (t in 2:time) {
    case[t] ~ dpois(I[t])
  }
}"

# Compile model
model <- jags.model(textConnection(model_code), data = data,  n.chains = 3, n.adapt = 1000)


##sample from the posterior
mcmc = coda.samples(model = model, variable.names = c("R_0", "sigma", "gamma", "omega"), n.iter = 50000, thin = 5)

head(mcmc)

plot(mcmc, density = TRUE)


# par(mfrow=c(4,1),mar=c(4,4,1,1))
# plot(data.frame(mcmc[["beta"]]), ,type="l")
# plot(data.frame(mcmc[["sigma"]]), ,type="l")
# plot(data.frame(mcmc[["gamma"]]), ,type="l")
# plot(data.frame(mcmc[["case"]]), ,type="l")


x