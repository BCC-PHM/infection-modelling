# Estimating Birmingham's polio resistance
library("latex2exp")
source("../infection-modelling.R")

# Set random seed
set.seed(1)

# Load in parameters
params_df <- readxl::read_excel(
  "../../data/parameters.xlsx", 
  sheet = "Polio")

vac_data <- read.csv(
  "../../data/Polio/IPV_age_5_booster.csv",
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
    
    SEIR_df <- SEIR_model(sampled_params, inits, timesteps = 365, t0 = 1+365*i) 
    
    SEIR_df$Year = vac_data$Year[i]
    
    yearly_dfs[[i]] <- SEIR_df
    
    # Get new inits
    inits <- tail(SEIR_df, n=1) %>%
      rename(s = Susceptible, e = Exposed, i = Infectious, r = Recovered) %>%
      select(c(s, e, i, r)) %>%
      as.list(SEIR_df)
    
  }
  
  run_df = data.table::rbindlist(yearly_dfs)
  run_df$run = run
  
  output_df_list[[run]] = run_df
}

full_sim_df = data.table::rbindlist(output_df_list) %>%
  mutate(
    date = as.Date("31/03/2015", format = "%d/%m/%Y") + Time
  ) 

full_sim_long <- full_sim_df  %>%
  select(-c("Died", "Infected_sum")) %>%
  reshape2::melt(id.vars = c("run", "Time", "Year", "date"), 
                 variable.name = 'Group',
                 value.name = 'Fraction')

summary_stats_r <- full_sim_long %>%
  group_by(date, Group) %>%
  summarise(
    median = median(Fraction),
    CI_lower = quantile(Fraction, 0.025),
    CI_upper = quantile(Fraction, 0.975)
  )



# Plot infection % as function of vaccination %
plt <- ggplot() +
  geom_line(data = full_sim_long, 
            aes(x = date, y = Fraction, group = as.factor(run)), 
            color = "grey", alpha = 0.3, linewidth = 0.3) +
  geom_line(data =summary_stats_r, aes(x = date, y = median, color = Group), 
            size =1) +
  geom_line(data =summary_stats_r,aes(x = date, y = CI_lower, color = Group), 
            linetype = "dashed", size = 0.9) +
  geom_line(data = summary_stats_r, aes(x = date, y = CI_upper, color = Group), 
            linetype = "dashed", size = 0.9) +
  theme_bw() +
  ylab("Simulated SEIR Population Percentages") +
  xlab("") +   
  scale_y_continuous(
    # Make y-axis percentages
    labels = scales::percent, 
  ) +
  facet_wrap(.~Group, scales = "free_y") +
  theme(strip.background = element_rect(fill="white")) +
  theme(legend.position="none")

plt

# Save plot
ggsave("../../output/figures/polio_estimate.png", plt, 
       width = 7, height = 5, dpi = 300)

# Plot distribution of values for today
all_vals_today <- full_sim_long %>% 
  filter(date == max(full_sim_long$date))
plt_hist <- ggplot(all_vals_today, aes(x=Fraction, fill = Group)) +
  geom_histogram(color = "black", linewidth = 0.2) +
  facet_wrap(.~Group, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white")) +
  theme(legend.position="none") +
  ylab("Number of Rimulation Runs")

plt_hist
ggsave("../../output/figures/polio_hist.png", plt_hist, 
       width = 7, height = 5, dpi = 300)

# Fit recovered distribution to skewed normal distribution
R_today <- all_vals_today %>%
  filter(Group == "Recovered") %>%
  select(Fraction)
mod <- sn::selm(Fraction ~ 1, data=R_today)
R_dist_params <- extractSECdistr(mod)
print(R_dist_params)

# Extract group fractions for last time step and save
estimate_today <- summary_stats_r[summary_stats_r$date == max(summary_stats_r$date),]
writexl::write_xlsx(estimate_today, "../../output/data/polio_inits.xslx")
print(estimate_today)