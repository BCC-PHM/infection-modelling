set.seed(2)
source("infection-modelling.R")


# Load in parameters
params_df <- readxl::read_excel(
  "../data/parameters.xlsx", 
  sheet = "Polio")

# convert parameter df to list
params <- setNames(as.list(params_df$Value),    
                   params_df$Parameter)

vac_vals <- list(
  "No new vaccinations" = 0,
  "No change" = 3.186e-05,
  "100%" = 4.21e-5
)

models <- data.frame(
  R_init = c(0.5, 0.75),
  R0 = c(2, 4),
  model = c("R(today)=0.5, R0 = 2", "R(today)=0.75, R0 = 4")
)

num_models <- nrow(models)


output_df_list <- list()

# Loop over all vaccination values
for (j in 1:length(vac_vals)) {
  # create empty list for run dfs
  run_dfs <- list()
  # loop over all R0 values
  for (i in 1:num_models) {
    
    r_init = models$R_init[[i]]
    # exposed ratio
    e_init = 1e-5
    # infected ratio
    i_init = 1e-5
    # susceptible ratio
    s_init = 1-r_init-e_init-i_init
    
    # Define initial conditions
    inits <- list(s = s_init, e = e_init, i = i_init, r = r_init)
    
    # convert parameter df to list
    params$R0 = models$R0[[i]]
  
    # Change vac with current value being tested
    params$vac = vac_vals[[j]]
    
    SEIR_df <- SEIR_model(params, inits, timesteps = 365*5, t0 = 1) 
    
    run_dfs[[i]] <- SEIR_df %>%
      # extract last data for last date
      filter(Time == max(Time)) %>%
      mutate(model = models$model[[i]])
  }
  
  output_df_list[[j]] = data.table::rbindlist(run_dfs)
  output_df_list[[j]]$vac_scenario = names(vac_vals)[j]
  
}

full_sim_df <- data.table::rbindlist(output_df_list) %>%
  select(c(vac_scenario, Infected_sum, model)) %>%
  mutate(
    Infection_difference = Infected_sum - mean(full_sim_df$Infected_sum)
  )


polio_result<-ggplot(full_sim_df, aes(x = model, y = Infected_sum, fill = vac_scenario)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  theme_bw() +
  theme(legend.position = "inside", legend.position.inside = c(0.2, 0.8)) +
  labs(
    y ='Proportion of Population\nInfected with Polio (5 years)', 
    x = 'Model',
    fill = "Vaccination Scenario"
    ) +
  viridis::scale_fill_viridis(discrete = TRUE, option = "G", 
                              begin = 0.3, end = 0.8) + 
  scale_y_continuous(labels = scales::percent) 

polio_result
ggsave("../output/figures/polio_fixed_param_result.png", polio_result, 
       width = 6, height = 4, dpi = 300)
