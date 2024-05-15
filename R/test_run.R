source("infection-modelling.R")
set.seed(1)

# simuation Time / Day
tsteps = 365
# Test random vaccination %'s from 0 to 1% (per day)
vac <- runif(200, 0, 1e-2)

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

total_infected <- c()
for (vac_i in vac) {
  # Create new parameter sample
  sampled_params <- param_sample(params_df)
  
  # replace vaccination rate with current list value
  sampled_params$vac <- vac_i
  
  # Run SEIR simulation
  SEIR_df <- SEIR_model(sampled_params, inits, tsteps) 
  
  # Extract last infection sum value
  total_infected <- c(total_infected, tail(SEIR_df$infected.sum, n=1))
}

# Calculate vaccination and infection percentages
test_df <- data.frame(
  "vaccination" = 100*vac,
  "infected" = 100*total_infected
)

# Plot infection % as function of vaccination %
plt <- ggplot(test_df, aes(x = vaccination, y = infected)) +
  geom_point(size = 0.5) +
  theme_bw() +
  # Add 95% confidence band
  geom_smooth(formula = 'y ~ x', method=lm, 
              level=0.95, color='#1f77b4', fill='lightblue') + 
  xlab("Daily Vaccination %") +
  ylab("Infection % (1 year)")


# Save plot
ggsave("../output/figures/SEIR-sampling-test.png", plt, 
       width = 5, height = 3)
