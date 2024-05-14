library(dplyr)
library(ggplot2)

#https://github.com/XuelongSun/Dynamic-Model-of-Infectious-Diseases/blob/master/SIR.ipynb

# population
N = 1e7 + 10 + 5
# simuation Time / Day
T = 170
# susceptiable ratio
s = rep(0, T)
# exposed ratio
e = rep(0, T)
# infective ratio
i = rep(0, T)
# remove ratio
r = rep(0, T)

# contact rate
lamda = 0.5
# recover rate
gamma = 0.0821
# exposed period
sigma = 1 / 4
# vaccination rate
vac = 0.006

# initial infective people
i[1] = 10.0 / N
s[1] = 1e7 / N
e[1] = 40.0 / N

for (t in 1:(T-1)){
  s[t + 1] = s[t] - lamda * s[t] * i[t] - vac * s[t]
  e[t + 1] = e[t] + lamda * s[t] * i[t] - sigma * e[t]
  i[t + 1] = i[t] + sigma * e[t] - gamma * i[t]
  r[t + 1] = r[t] + gamma * i[t] + vac * s[t]
}

df <- data.frame(
  "time" = 1:T,
  "susceptible" = s,
  "exposed" = e,
  "infectious" = i,
  "recovered" = r
) %>%
  reshape2::melt(id.vars = "time", 
       variable.name = 'Group',
       value.name = 'Number')

plt <- ggplot(df, aes(x = time, y = Number, col = Group)) +
  geom_line(linewidth = 1.2) +
  theme_bw() +
  ggtitle("Vaccination rate = 0.002") + theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave("../output/figures/SEIR-test-with-vac.png", plt, 
       width = 5, height = 3)