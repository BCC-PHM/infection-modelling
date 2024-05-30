set.seed(2)
source("infection-modelling.R")

num_param_samples <- 100

pop_size <- 1144900

vac_vals <- list(
  "No new vaccinations" = 0,
  "No change" = 3.186e-05,
  "100%" = 4.21e-5
)

set.seed(0)

post_chain <- read.csv("../data/Polio/polio_posterior_final 1.csv") %>%
  sample_n(num_param_samples)

params <- list(
  "1/sigma" = 4,
  "1/omega" = 1e99,
  "vac" = 3.374e-5,
  "alpha" = 0,
  "mu_b" = 3.464e-5,
  "mu_d" = 3.449e-5
)

output_df_list <- list()

# Loop over all vaccination values
for (j in 1:length(vac_vals)) {
  # create empty list for run dfs
  run_dfs <- list()
  # loop over all R0 values
  for (i in 1:num_param_samples) {
    
    inits <- list(
      r = 0.8,
      e = 8.0e-5,
      s = 1 - 0.8 - 8.0e-5 - 3.0e-4,
      i = 3.0e-5
    )
    
    # convert parameter df to list
    params$R0 = post_chain$R_0[[i]]
    params$hosp = post_chain$eta[[i]]
    params$`1/gamma` = 1/post_chain$gamma[[i]]
    # Change vac with current value being tested
    params$vac = vac_vals[[j]]
    
    SEIR_df <- SEIR_model(params, inits, timesteps = 365*5, t0 = 1) 
    
    run_dfs[[i]] <- SEIR_df %>%
      # extract last data for last date
      filter(Time == max(Time)) %>%
      mutate(run = i)
  }
  
  output_df_list[[j]] = data.table::rbindlist(run_dfs)
  output_df_list[[j]]$vac_scenario = names(vac_vals)[j]
  
}

full_sim_df <- data.table::rbindlist(output_df_list) %>%
  select(c(vac_scenario, Infected_sum, Admissions_sum, run)) %>%
  mutate(
    Infected_sum = Infected_sum * pop_size,
    Admissions_sum =Admissions_sum * pop_size
    ) %>%
  group_by(vac_scenario) %>%
  pivot_longer(
    cols = c("Infected_sum", "Admissions_sum"),
    names_to='Count_Type',
    values_to='Count') %>%
  group_by(vac_scenario, Count_Type) %>%
  summarise(
    median_count = quantile(Count, 0.5),
    lowerCI = quantile(Count, 0.25),
    upperCI = quantile(Count, 0.75)
  ) 


polio_result<-ggplot(full_sim_df, aes(x = vac_scenario, y = median_count, 
                                      fill = vac_scenario)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin = lowerCI, ymax = upperCI), position = "dodge", 
                width = 0.25) +
  theme_bw() +
  #theme(legend.position = "inside", legend.position.inside = c(0.1, 0.9)) +
  labs(
    y ='Estimated Number (5 years)', 
    x = 'Model',
    fill = "Vaccination Scenario"
    ) +
  viridis::scale_fill_viridis(discrete = TRUE, option = "G", 
                              begin = 0.3, end = 0.8) +
  facet_wrap(~Count_Type, scales ="free_y", ncol = 1)

polio_result
ggsave("../output/figures/polio_fixed_param_result.png", polio_result, 
       width = 6, height = 4, dpi = 300)
