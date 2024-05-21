# Assume fixed values for 1/gamma, 1/sigma. Hence calculate R0 
# based on hospitalisation data assuming total number of infections is less than
# 100 times the number of hospitalisations

library(dplyr)
library(ggplot2)
library(GGally)

source("infection-modelling.R")

set.seed(0)

nsims = 10000

# Load in parameters
params_df <- readxl::read_excel(
  "../data/parameters.xlsx", 
  sheet = "Polio") 
# convert parameter df to list
params <- setNames(as.list(params_df$Value),    
         params_df$Parameter)

R0 = runif(n=nsims, min = 0, max = 6)
R_inits = runif(n=nsims, min = 0, max = 1)

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

# Load hospitalisation data
path = "//svwvap1126.addm.ads.brm.pri/PHSensitive$/Intelligence/2. Requests/REQ3280 - HES Inpatients data on polio, flu and rubella for Justin/REQ3280_Inpatients_data_10years_on_polio_flu_and_rubella_BirminghamResidents.xlsx"
hosp_data <- readxl::read_excel(path)%>%
  filter(
    Condition == "Polio" &
    Year %in% vac_data$Year
  ) %>%
  group_by(Year) %>%
  summarise(
    `No of Admissions` = sum(`No of Admissions`)
  )

total_admissions <- sum(hosp_data$`No of Admissions`)

total_infections <- list()
R_today <- list()
last_year_infections <- c()

for (i in 1:nsims) {
  inits <- list(
    # susceptible ratio
    s = 1-R_inits[[i]]-2e-5,
    # exposed ratio
    e = 1e-5,
    # infected ratio
    i = 1e-5,
    # recovered ratio
    r = R_inits[[i]]
  )
  
  params$R0 = R0[[i]]
  
  infections_i = 0
  
  for (j in 1:nrow(vac_data)) {
    vac_i <- vac_data$vac[j]
    params$vac <- vac_i
    
    SEIR_df <- SEIR_model(params, inits, timesteps = 365, t0 = 1+365*i) 
    
    # Get new inits
    inits <- tail(SEIR_df, n=1) %>%
      rename(s = Susceptible, e = Exposed, i = Infectious, r = Recovered) %>%
      select(c(s, e, i, r)) %>%
      as.list(SEIR_df)
    
    yearly_infections <- vac_data$`Total Population`[j] * tail(SEIR_df$Infected_sum, 1)
    # print(yearly_infections)
    # infections_i <- infections_i + yearly_infections
    # incompatible_year_i <- incompatible_year_i + 
    #   # too few
    #   (yearly_infections < hosp_data$`No of Admissions`[j]) +
    #   # too many
    #   (yearly_infections > 100*hosp_data$`No of Admissions`[j]) 
  }
  
  total_infections[[i]] <- infections_i
  R_today[[i]] <- tail(SEIR_df$Recovered, 1)
  last_year_infections[[i]] <- yearly_infections

}

result <- data.frame(
  R0 = R0,
  R_init = R_inits,
  total_infections = unlist(total_infections),
  R_today = unlist(R_today),
  last_year_infections = unlist(last_year_infections)
) %>%
  mutate(
    infections = case_when(
      last_year_infections < tail(hosp_data$`No of Admissions`, 1) ~ "Too few",
      last_year_infections > 100*tail(hosp_data$`No of Admissions`, 1) ~ "Too many",
      TRUE ~ "Possible",
    )
  )

# ggplot(result, aes(x = R0, y = R_today, color = total_infections)) +
#   geom_point() +
#   theme_bw() +
#   scale_color_gradient(name = "Total Infectons", trans = "log",
#                        breaks = 10^(-2:5))
# 
# ggplot(result, aes(x = R_init, y = R_today, color = R0)) +
#   geom_point() +
#   theme_bw() +
#   scale_color_gradient(name = "R0")


gm <- ggpairs(result, columns = c("R0", "R_init", "R_today"), 
        mapping = aes(color = infections, alpha = 0.7)) + 
  theme_bw()
gm

plt <- result %>%
  filter(infections == "Possible") %>%
  ggplot(aes(x = R0, y = R_today)) +
  geom_point() +
  theme_bw()
plt