set.seed(2)
source("infection-modelling.R")

nsims = 100#0

# Values from polio_R_estimate.R
R_init_vals = sn::rsn(n=nsims, xi=0.8207131, omega=0.0386484, alpha=2.0402655)

# Load in parameters
params_df <- readxl::read_excel(
  "../data/parameters.xlsx", 
  sheet = "Polio")

# Checking impact of smaller R0
params_df[params_df$Parameter == "R0",]$Value = 1
params_df[params_df$Parameter == "R0",]$lowerCI = 0.8
params_df[params_df$Parameter == "R0",]$upperCI = 1.2

# Check if any R0 value violates the bounds
testit::assert("R0 violated bounds of (0 to 1)", {
  !(any(R_init_vals <= 0) | any(R_init_vals > 1))
})

R_init = R_init_vals[0]

# This may seem off since England is "worse" than Birmingham
# however, this is a reflection of the fact that Birmingham
# has a younger population and therefore a larger percentage
# of the *total* population are vaccinated each year.
vac_vals <- list(
  "No vac" = 0,
  #"National" = 2.722e-05, # 567057 / 57106400 / 365,
  "No change" = 3.186e-05,
  "100%" = 4.21e-5
)

output_df_list <- list()

# Loop over all vaccination values
for (j in 1:length(vac_vals)) {
  # create empty list for run dfs
  run_dfs <- list()
  # loop over all R0 values
  for (i in 1:nsims) {
    r_init = R_init_vals[[i]]
    # exposed ratio
    e_init = 3e-4
    # infected ratio
    i_init = 7e-4
    # susceptible ratio
    s_init = 1-r_init-e_init-i_init
    
    # Define initial conditions
    inits <- list(s = s_init, e = e_init, i = i_init, r = r_init)
    
    sampled_params <- param_sample(params_df)
    # Change vac with current value being tested
    sampled_params$vac = vac_vals[[j]]
    sampled_params$R0 = 1
    
    SEIR_df <- SEIR_model(sampled_params, inits, timesteps = 365*5, t0 = 1) 
    
    run_dfs[[i]] <- SEIR_df %>%
      # extract last data for last date
      filter(Time == max(Time))
  }
  
  output_df_list[[j]] = data.table::rbindlist(run_dfs)
  output_df_list[[j]]$vac_scenario = names(vac_vals)[j]
  
}

full_sim_df <- data.table::rbindlist(output_df_list) %>%
  select(c(vac_scenario, Infected_sum))

infected_sum_results <- full_sim_df %>%
  group_by(vac_scenario) %>%
  summarise(
    median = median(Infected_sum),
    CI_lower = quantile(Infected_sum, 0.025),
    CI_upper = quantile(Infected_sum, 0.975)
  )

polio_result<-ggplot(infected_sum_results, aes(x = vac_scenario, y = median, 
                                 fill = vac_scenario)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper,), width=.2,
                position=position_dodge(.9))+
  theme_bw() +  
  theme(
    legend.position = "None"
  ) +
  labs(
    y ='Median Proportion of Population Infected\nwith Polio (5 Years)', 
    x = 'Vaccination Scenario'
    ) +
  viridis::scale_fill_viridis(discrete = TRUE, option = "G", 
                              begin = 0.3, end = 0.8)
polio_result
# ggsave("../output/figures/polio_result.pdf", polio_result, 
#        width = 5, height = 4, dpi = 300)

polio_result_dist <- ggplot(full_sim_df, aes(x = Infected_sum, 
                        linetype = vac_scenario, color = vac_scenario)) +
  geom_freqpoly(bins = 50, linewidth = 1) +
  theme_bw() +
  #facet_wrap(~vac_scenario) +
  theme(strip.background = element_rect(fill="white")) +
  #scale_x_continuous(trans='log10') +
  #theme(legend.position="none") +
  xlab("Total Proportion of Population Infected with Polio (5 Years)") +
  viridis::scale_color_viridis(discrete = TRUE, option = "G", 
                               begin = 0.3, end = 0.8) +
  ylab("Number of Rimulation Runs") +
  labs(color='Vaccination Scenario', linetype = 'Vaccination Scenario') +
  theme(
    legend.position = "inside", 
    legend.position.inside =  c(.2, .83)
    )
polio_result_dist
# ggsave("../output/figures/polio_result_dist.pdf", polio_result_dist, 
#        width = 5, height = 4, dpi = 300)
  
