source("infection-modelling.R")

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

# simuation Time / Day
tsteps = 365*5

params_df <- readxl::read_excel("../data/parameters.xlsx", sheet = "test")

params_df$Distribution = "Fixed"
output <- param_sample(params_df)

df <- SEIR_model(output, inits, tsteps) %>%
  reshape2::melt(id.vars = "time", 
                 variable.name = 'Group',
                 value.name = 'Fraction') %>%
  filter(
    Group %in% c("susceptible","exposed","infectious","recovered")
  ) %>%
  mutate(
    Years = time/365
  )

plt <- ggplot(df, aes(x = Years, y = Fraction, col = Group)) +
  geom_line(linewidth = 1.1) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave("../output/figures/SEIR-test.png", plt, 
       width = 5, height = 3)