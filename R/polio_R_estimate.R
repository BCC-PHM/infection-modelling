# Estimating Birmingham's polio resistance
library("latex2exp")
source("infection-modelling.R")
set.seed(1)

# simuation Time / Day
tsteps = 365
# Test random vaccination %'s from 0 to 1% (per day)
#vac <- runif(200, 0, 1e-2)

# Load in parameters
params_df <- readxl::read_excel(
  "../data/parameters.xlsx", 
  sheet = "Polio")

vac_data <- read.csv(
  "../data/Polio/IPV_age_5_booster.csv",
  check.names = FALSE
) %>%
  filter(
    AreaName == "Birmingham"
  ) %>%
  mutate(
    Year = `Time period`,
    vac = Count/`Total Population`/365,
    ) %>%
  select(
    c(Year, vac, Count, `Total Population`)
  )

output_df_list <- list()

for (run in 1:1000){
  r = 0.85 + runif(1, -0.15, 0.1)
  e = 10^runif(1, -4, -1)
  s = 1-r-e
  
  # Define initial conditions
  inits <- list(
    # susceptible ratio
    s = s,
    # exposed ratio
    e = e,
    # infected ratio
    i = 0,
    # recovered ratio
    r = r
  )
  
  sampled_params <- param_sample(params_df)
  
  yearly_dfs <- list()
  
  for (i in 1:nrow(vac_data)) {
    vac_i <- vac_data$vac[i]
    sampled_params$vac <- vac_i
    
    SEIR_df <- SEIR_model(sampled_params, inits, tsteps, t0 = 1+365*i) 
    
    SEIR_df$Year = vac_data$Year[i]
    
    yearly_dfs[[i]] <- SEIR_df
    
    # Get new inits
    inits <- tail(SEIR_df, n=1) %>%
      rename(s = susceptible, e = exposed, i = infectious, r = recovered) %>%
      select(c(s, e, i, r)) %>%
      as.list(SEIR_df)
    
  }
  
  run_df = data.table::rbindlist(yearly_dfs)
  run_df$run = run
  
  output_df_list[[run]] = run_df
}

full_sim_df = data.table::rbindlist(output_df_list) %>%
  mutate(
    date = as.Date("31/03/2015", format = "%d/%m/%Y") + time
  )



summary_stats_r <- full_sim_df %>%
  group_by(date) %>%
  summarise(
    median = median(recovered),
    CI_lower = quantile(recovered, 0.025),
    CI_upper = quantile(recovered, 0.975)
  )



# Plot infection % as function of vaccination %
plt <- ggplot() +
  geom_line(data = full_sim_df, 
            aes(x = date, y = recovered, group = as.factor(run)), 
            color = "grey", alpha = 0.3, linewidth = 0.3) +
  geom_line(data =summary_stats_r, aes(x = date, y = median), color = "purple", 
            size =0.9) +
  geom_line(data =summary_stats_r,aes(x = date, y = CI_lower), 
            linetype = "dotted", color = "purple", size = 0.9) +
  geom_line(data = summary_stats_r, aes(x = date, y = CI_upper), 
            linetype = "dotted", color = "purple", size = 0.9) +
  geom_ribbon(data =summary_stats_r, aes(x = date, ymax = CI_upper, ymin = CI_lower), 
              fill = "purple", alpha=0.05)+
  theme_bw() +
  ylab("Percentage of Birmingham with\nPolio Resistance") +
  xlab("") +   
  scale_y_continuous(
    # Make y-axis percentages
    labels = scales::percent, 
    # Force y-axis to start at zero
    expand = c(0, 0), limits = c(0.7, 1)
  ) +
  annotate("text", x = as.Date("01/12/2021", format = "%d/%m/%Y"), y=0.98, size = 3.5,
           label = TeX("$e_0\\in $ (0.0001, 0.1)"))+
  annotate("text", x = as.Date("01/12/2021", format = "%d/%m/%Y"), y=0.96, size = 3.5,
           label = TeX("$n = 1000"))+  
  annotate("text", x = as.Date("01/12/2021", format = "%d/%m/%Y"), y=0.94, size = 3.5,
           label = TeX("$r_t = $0.85 (0.80, 0.91)"))

plt

estimate_today <- tail(summary_stats_r, 1)
print(estimate_today)


# Save plot
ggsave("../output/figures/polio_estimate.png", plt, 
       width = 5, height = 3.5)