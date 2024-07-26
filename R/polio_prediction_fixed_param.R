library(dplyr)
library(ggplot2)
source("infection-modelling.R")

set.seed(0)

num_param_samples <- 1000

pop_size <- 1361159

vac_vals <- list(
  "50% vaccination rate" = 4.21e-5/2,
  "No change" = 3.186e-05,
  "100% vaccination rate" = 4.21e-5
)

post_chain <- read.csv("polio_posterior_BSOL.csv") %>%
  sample_n(num_param_samples)

params <- list(
  "1/sigma" = 4,
  "1/omega" = 1e99,
  "1/gamma" = 1/0.09,
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
    set.seed(0)
    inits <- list(
      r = 0.803,
      e = 0.0000837,
      s = 1 - 0.803 - 0.0000837 - 0.000232,
      i = 0.000232
    )
    
    # convert parameter df to list
    params$R0 = post_chain$R_0[[i]]
    params$hosp = post_chain$eta[[i]]
  
    # Change vac with current value being tested
    params$vac = vac_vals[[j]]
    
    SEIR_df <- SEIR_surge_model(params, inits, timesteps = 365*10, t0 = 1, t_surge = 365*5, surge_factor = 100) #change to surge model
    
    run_dfs[[i]] <- SEIR_df %>%
      mutate(
        run = i,
        R0 = post_chain$R_0[[i]],
        eta = post_chain$eta[[i]]
        )
  }
  
  output_df_list[[j]] = data.table::rbindlist(run_dfs)
  output_df_list[[j]]$vac_scenario = names(vac_vals)[j]
  
}

full_sim_df <- data.table::rbindlist(output_df_list) %>%
  select(c(vac_scenario, Infectious, Admissions, Infected_sum, Admissions_sum, Time, run, R0, eta)) %>%
  mutate(
    `Polio Infections` = Infected_sum * pop_size,
    `Acute Polio Admissions` = Admissions_sum * pop_size,
     Infectious =  Infectious* pop_size,
     Admissions = Admissions* pop_size
    ) %>%
  group_by(vac_scenario, Time) %>%
  tidyr::pivot_longer(
    cols = c("Polio Infections", "Acute Polio Admissions", "Infectious", "Admissions"),
    names_to='Count_Type',
    values_to='Count') %>%
  ungroup()

result_snapshots <- full_sim_df %>%
  filter(Time %in% c(365*5, 365*10)) %>%
  mutate(
    Year = case_when(
      Time == 365*5 ~ "5 Years",
      Time == 365*10 ~ "10 Years")
    ,
    Year = factor(Year, 
                  level = c("5 Years", "10 Years")),
    vac_scenario = factor(vac_scenario, 
                          level = c("50% vaccination rate", "No change","100% vaccination rate"))
  )

# extract last data for last date


sim_quantiles <- result_snapshots  %>%
  filter(Count_Type != "Infectious" & Count_Type != "Admissions") %>% 
  group_by(vac_scenario, Count_Type, Year) %>%
  summarise(
    median_count = quantile(Count, 0.5),
    lowerCI = quantile(Count, 0.25),
    upperCI = quantile(Count, 0.75)
  ) 


########################################
##########   Quantile plot   ###########
########################################


polio_quantiles <- ggplot(sim_quantiles, aes(x = Year, y = median_count, 
                                             fill = vac_scenario)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge(), width = 0.75) +
  geom_errorbar(aes(ymin = lowerCI, ymax = upperCI), 
                position = position_dodge(0.75), 
                width = 0.15) +
  theme_bw() +
  #theme(legend.position = "inside", legend.position.inside = c(0.15, 0.85)) +
  labs(
    y ='Estimated Number (Cummulative)', 
    x = '',
    fill = "Vaccination Scenario"
  ) +
  viridis::scale_fill_viridis(discrete = TRUE, option = "G", 
                              begin = 0.3, end = 0.8) +
  facet_wrap(~Count_Type, scales ="free_y", ncol = 1) +
  theme(strip.background = element_rect(fill="white"))

polio_quantiles
ggsave("../output/figures/polio_fixed_param_result.pdf", polio_quantiles, 
       width = 8, height = 5, dpi = 300)

########################################
########## Distribution plot ###########
########################################

polio_dists <- ggplot(result_snapshots , aes(x = Count, color = vac_scenario)) +
  geom_density(linewidth = 1.2) +
  theme_bw() +
  #theme(legend.position = "inside", legend.position.inside = c(0.1, 0.9)) +
  labs(
    x ='Estimated Number (Cummulative)', 
    color = "Vaccination Scenario"
  ) +
  viridis::scale_fill_viridis(discrete = TRUE, option = "G", 
                              begin = 0.3, end = 0.8) +
  facet_wrap(~Count_Type*Year, scales ="free", ncol = 2)+
  theme(strip.background = element_rect(fill="white"))

polio_dists

###############################################
##########   Count over time plot   ###########
###############################################

over_time <- full_sim_df %>%
  group_by(Time, vac_scenario, Count_Type) %>%
  summarise(
    median_count = quantile(Count, 0.5),
    lowerCI = quantile(Count, 0.25),
    upperCI = quantile(Count, 0.75)
  ) %>%
  mutate(
    Years = Time/365,
    vac_scenario = factor(vac_scenario, 
                          level = c("50% vaccination rate", "No change","100% vaccination rate"))
  )

time_plt <- ggplot(over_time, aes(x = Years, y = median_count, color = vac_scenario)) +
  geom_ribbon(aes(ymin = lowerCI,
                  ymax = upperCI,
                  fill = vac_scenario),
              alpha = 0.3,
              colour = NA) +
  geom_line(linewidth = 1.2) +
  theme_bw() +
  facet_wrap(~Count_Type, scale = "free_y", ncol = 1) +
  #theme(legend.position = "inside", legend.position.inside = c(0.1, 0.9)) +
  labs(
    x ='Years from today', 
    y = "Estimated Number (Cummulative)",
    color = "Vaccination Scenario"
  ) +
  viridis::scale_color_viridis(discrete = TRUE, option = "G", 
                              begin = 0.3, end = 0.8) +
  viridis::scale_fill_viridis(discrete = TRUE, option = "G", 
                               begin = 0.3, end = 0.8) +
  theme(strip.background = element_rect(fill="white")) +
  guides(fill = "none")

time_plt
ggsave("../output/figures/polio_time_plot.pdf", time_plt, 
       width = 7, height = 5, dpi = 300)

bifurcation <- ggplot(result_snapshots, aes(x = R0, y = Infected_sum * pop_size, color = vac_scenario)) +
  geom_point(size = 0.6) +
  theme_bw() +
  viridis::scale_color_viridis(discrete = TRUE, option = "G", 
                             begin = 0.3, end = 0.8) +
  labs(
    y = "Total Infection Estimate (10 Years)", 
    x = latex2exp::TeX("Basic Reproduction Number $R_0$"), 
    color = "Vaccination Scenario"
  ) +
  scale_y_log10() +
  theme(legend.position = "inside", legend.position.inside = c(0.2, 0.75)) 

ggsave("../output/figures/polio_bifurcation.pdf", bifurcation, 
       width = 6, height = 3, dpi = 300)



########################################################################
##########   cumulative and non-cumulative  over time plot   ###########
########################################################################

##by 1000

over_time_surge_1000 <- full_sim_df %>%
  group_by(Time, vac_scenario, Count_Type) %>%
  summarise(
    median_count = quantile(Count, 0.5),
    lowerCI = quantile(Count, 0.25),
    upperCI = quantile(Count, 0.75)
  ) %>%
  mutate(
    Years = Time/365,
    vac_scenario = factor(vac_scenario, 
                          level = c("50% vaccination rate", "No change","100% vaccination rate"))
  )

over_time_surge_1000$Count_Type = ifelse(over_time_surge_1000$Count_Type == "Admissions", "Acute Polio Admissions ",
                                    ifelse(over_time_surge_1000$Count_Type == "Infectious", "Polio Infections ", over_time_surge_1000$Count_Type))


time_surge_plt_1000 <- ggplot(over_time_surge_1000, aes(x = Years, y = median_count, color = vac_scenario)) +
  geom_ribbon(aes(ymin = lowerCI,
                  ymax = upperCI,
                  fill = vac_scenario),
              alpha = 0.3,
              colour = NA) +
  geom_line(linewidth = 1.2) +
  theme_bw() +
  facet_wrap(~factor(Count_Type, c("Acute Polio Admissions", "Polio Infections",  "Acute Polio Admissions ", "Polio Infections ")), 
             scale = "free_y", ncol = 2) +
  ylab("  Estimated Number (Non-cumulative)             Estimated Number (Cumulative)\n")+
  theme(axis.title.y = element_text(hjust = 0.25, size = 14),
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+
  #theme(legend.position = "inside", legend.position.inside = c(0.1, 0.9)) +
  labs(
    x ='Years from today', 
    color = "Vaccination Scenario",
    title = "Year 5 Surge with 1000 times Increased Infection Rate"
  ) +
  viridis::scale_color_viridis(discrete = TRUE, option = "G", 
                               begin = 0.3, end = 0.8) +
  viridis::scale_fill_viridis(discrete = TRUE, option = "G", 
                              begin = 0.3, end = 0.8) +
  theme(strip.background = element_rect(fill="white")) +
  guides(fill = "none")

time_surge_plt_1000


##by 100

over_time_surge_100 <- full_sim_df %>%
  group_by(Time, vac_scenario, Count_Type) %>%
  summarise(
    median_count = quantile(Count, 0.5),
    lowerCI = quantile(Count, 0.25),
    upperCI = quantile(Count, 0.75)
  ) %>%
  mutate(
    Years = Time/365,
    vac_scenario = factor(vac_scenario, 
                          level = c("50% vaccination rate", "No change","100% vaccination rate"))
  )

over_time_surge_100$Count_Type = ifelse(over_time_surge_100$Count_Type == "Admissions", "Acute Polio Admissions ",
                                    ifelse(over_time_surge_100$Count_Type == "Infectious", "Polio Infections ", over_time_surge_100$Count_Type))


time_surge_plt_100 <- ggplot(over_time_surge_100, aes(x = Years, y = median_count, color = vac_scenario)) +
  geom_ribbon(aes(ymin = lowerCI,
                  ymax = upperCI,
                  fill = vac_scenario),
              alpha = 0.3,
              colour = NA) +
  geom_line(linewidth = 1.2) +
  theme_bw() +
  facet_wrap(~factor(Count_Type, c("Acute Polio Admissions", "Polio Infections",  "Acute Polio Admissions ", "Polio Infections ")), 
             scale = "free_y", ncol = 2) +
  ylab("  Estimated Number (Non-cumulative)             Estimated Number (Cumulative)\n")+
  theme(axis.title.y = element_text(hjust = 0.25, size = 14),
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+
  #theme(legend.position = "inside", legend.position.inside = c(0.1, 0.9)) +
  labs(
    x ='Years from today', 
    color = "Vaccination Scenario",
    title = "Year 5 Surge with 100 times Increased Infection Rate"
  ) +
  viridis::scale_color_viridis(discrete = TRUE, option = "G", 
                               begin = 0.3, end = 0.8) +
  viridis::scale_fill_viridis(discrete = TRUE, option = "G", 
                              begin = 0.3, end = 0.8) +
  theme(strip.background = element_rect(fill="white")) +
  guides(fill = "none")

time_surge_plt_100





ggplot(over_time_surge_100, aes(x = Years, y = median_count, color = vac_scenario)) +
  geom_ribbon(aes(ymin = lowerCI,
                  ymax = upperCI,
                  fill = vac_scenario),
              alpha = 0.3,
              colour = NA) +
  geom_line(linewidth = 1.2) +
  theme_bw() +
  facet_wrap(~factor(Count_Type, c("Acute Polio Admissions ", "Polio Infections ")), 
             scale = "free_y", ncol = 2) +
  ylab("  Estimated Number (Non-cumulative)\n")+
  theme(axis.title.y = element_text(hjust = 0.25, size = 14),
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+
  #theme(legend.position = "inside", legend.position.inside = c(0.1, 0.9)) +
  labs(
    x ='Years from today', 
    color = "Vaccination Scenario",
    title = "Year 5 Surge with 100 times Increased Infection Rate"
  ) +
  viridis::scale_color_viridis(discrete = TRUE, option = "G", 
                               begin = 0.3, end = 0.8) +
  viridis::scale_fill_viridis(discrete = TRUE, option = "G", 
                              begin = 0.3, end = 0.8) +
  theme(strip.background = element_rect(fill="white")) +
  guides(fill = "none")





over_time_surge_100 %>% filter(Count_Type == "Polio Infections ") %>% 
  ggplot(aes(x = Years, y = median_count, color = vac_scenario)) +
  geom_ribbon(aes(ymin = lowerCI,
                  ymax = upperCI,
                  fill = vac_scenario),
              alpha = 0.3,
              colour = NA) +
  geom_line(linewidth = 1.2) +
  theme_bw()+
  facet_wrap(~factor(Count_Type, c("Polio Infections ")), 
              ncol = 1)+
  theme(axis.title.y = element_text(hjust = 0.25, size = 14),
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+
  labs(
    x ='Years from today', 
    color = "Vaccination Scenario",
    title = "Year 5 Surge with 100 times Increased Infection Rate"
  ) +
  viridis::scale_color_viridis(discrete = TRUE, option = "G", 
                               begin = 0.3, end = 0.8) +
  viridis::scale_fill_viridis(discrete = TRUE, option = "G", 
                              begin = 0.3, end = 0.8) +
  theme(strip.background = element_rect(fill="white")) +
  guides(fill = "none")+
  scale_y_continuous(limits = c(0,1000), breaks = seq(0,1000, by=250))


  
scales = list(
  scale_y_continuous(limits = c(0,250), breaks = seq(0,250,by=50)),
  scale_y_continuous(limits = c(0,100000), breaks = seq(0,100000,by=25000)),
  scale_y_continuous(limits = c(0.00,0.20), breaks = seq(0.00,0.20, by=0.05)),
  scale_y_continuous(limits = c(0,800), breaks = seq(0,800,by=100))
)




time_surge_plt_100+facetted_pos_scales(y=scales)


time_surge_plt_1000+facetted_pos_scales(y=scales)









ggsave("../output/figures/polio_time_plot.pdf", time_plt, 
       width = 7, height = 5, dpi = 300)

bifurcation <- ggplot(result_snapshots, aes(x = R0, y = Infected_sum * pop_size, color = vac_scenario)) +
  geom_point(size = 0.6) +
  theme_bw() +
  viridis::scale_color_viridis(discrete = TRUE, option = "G", 
                               begin = 0.3, end = 0.8) +
  labs(
    y = "Total Infection Estimate (10 Years)", 
    x = latex2exp::TeX("Basic Reproduction Number $R_0$"), 
    color = "Vaccination Scenario"
  ) +
  scale_y_log10() +
  theme(legend.position = "inside", legend.position.inside = c(0.2, 0.75)) 

ggsave("../output/figures/polio_bifurcation.pdf", bifurcation, 
       width = 6, height = 3, dpi = 300)














