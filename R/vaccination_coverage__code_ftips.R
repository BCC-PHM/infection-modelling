library(ggplot2)
library(dplyr)

X1_year <- read.csv("../data/Polio/vaccines/1year.csv", check.names = F)
X2_year <- read.csv("../data/Polio/vaccines/2year.csv", check.names = F)
X5_year <- read.csv("../data/Polio/vaccines/5year.csv", check.names = F)

combined_year = rbind(X1_year, X2_year, X5_year)

combined_year$`Time period` <- gsub("^20", "", combined_year$`Time period`)

combined_year$`Time period` = factor(combined_year$`Time period`)

plt <- combined_year %>% ggplot(aes(x=`Time period`, y=Value, color = AreaName, group = AreaName))+
  geom_point()+
  geom_line()+
  geom_ribbon(aes(ymin=`Lower CI 95.0 limit`, ymax= `Upper CI 95.0 limit`, fill=AreaName),alpha=0.5,colour=NA)+
  geom_text(data=combined_year %>% filter(`Time period` == "2022/23") %>% select(`Time period`, `Indicator Name`,AreaName,Value), 
            aes(label=paste0(round(Value,1), "%")), vjust=-0.5, hjust=0.6,size=4,show.legend = FALSE)+
  labs(x="Year",
       y="Population IPV Coverage (%)")+
  scale_color_manual(values = c("#3488a6", "#404040")) +
  scale_fill_manual(values = c("#3488a6", "#404040")) +
  theme_bw() +
  theme(legend.title = element_blank(),
         strip.background = element_rect(fill="white")) +
  facet_wrap(~factor(`Indicator Name`, levels = c("Dtap IPV  Hib HepB (1 year old)",
                                                  "Dtap  IPV  Hib HepB (2 years old)",
                                                  "DTaP and IPV booster (5 years)" )), nrow = 3, ncol = 1)
  


combined_year$`Indicator Name`


ggsave("../output/IPV-vaccination.svg", plt, width = 6.5, height = 5)
ggsave("../output/IPV-vaccination.pdf", plt, width = 6.5, height = 5)