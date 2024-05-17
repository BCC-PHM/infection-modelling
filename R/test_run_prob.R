source("infection-modelling2.R")
set.seed(1)

# simuation Time / Day
tsteps = 365

# Load in parameters
params_df <- readxl::read_excel(
  "../data/parameters.xlsx", 
  sheet = "test")

# Define initial conditions
inits <- list(
  # susceptible ratio
  s = 0.99,
  # exposed ratio
  e = 0.001,
  # infected ratio
  i = 0,
  # recovered ratio
  r = 0
)



# Create new parameter sample with 100 iterations
params= param_sample_prob(params_df, 
                          iter = 1000, 
                          vac = 0.5)  #set a vaccine coverage rate here 


# Run SEIR simulation
SEIR_df <- SEIR_model_prob(params, inits, tsteps) 


#plot the data for s for now
plot_SEIR_prob(SEIR_df)











total_infected <- c()
















for (vac_i in vac) {
  # Create new parameter sample with 100 iterations
  params= param_sample(params_df, iter = 100)
  
  # replace vaccination rate with current list value
  for(i in seq_along(params)) {
  params[[i]]$vac <- vac_i
  }
  
  # Run SEIR simulation
  SEIR_df <- SEIR_model(sampled_params, inits, tsteps) 
  
  # Extract last infection sum value
  total_infected <- c(total_infected, tail(SEIR_df$infected.sum, n=1))
}




params= param_sample(params_df, iter = 100)



SEIR_df <- SEIR_model(params =params, inits, tsteps)



