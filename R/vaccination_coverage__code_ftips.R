X1_year <- read_csv("~/R projects/infection-modelling/data/1 year.csv")
X2_year <- read_csv("~/R projects/infection-modelling/data/2year.csv")
X5_year <- read_csv("~/R projects/infection-modelling/data/5year.csv")

combined_year = rbind(X1_year, X2_year, X5_year)

combined_year$`Time period` = factor(combined_year$`Time period`)

combined_year %>% ggplot(aes(x=`Time period`, y=Value, color = AreaName, group = AreaName))+
                  geom_point()+
                  geom_line()+
                  geom_ribbon(aes(ymin=`Lower CI 95.0 limit`, ymax= `Upper CI 95.0 limit`, fill=AreaName),alpha=0.5,colour=NA)+
                  geom_text(data=combined_year %>% filter(`Time period` == "2022/23") %>% select(`Time period`, `Indicator Name`,AreaName,Value), 
                            aes(label=paste0(round(Value,2), "%")), vjust=-0.4, hjust=0.5,size=5,show.legend = FALSE)+
                  labs(x="Year",
                      y="Percentage(%)")+
                  theme_bw() +
                  theme(legend.title = element_blank(),
                        strip.background = element_rect(fill="white"),
                        axis.title.y = element_text(hjust = 0.5, size = 16),
                        plot.title = element_text(hjust = 0.5, size = 18),
                        axis.title.x = element_text(size = 16),
                        axis.text = element_text(size = 14),
                        strip.text = element_text(size = 14),
                        legend.text = element_text(size = 14))+
                  facet_wrap(~factor(`Indicator Name`, levels = c("Population vaccination coverage: Dtap IPV  Hib HepB (1 year old)",
                                                                  "Population vaccination coverage: Dtap  IPV  Hib HepB (2 years old)",
                                                                  "Population vaccination coverage: DTaP and IPV booster (5 years)" )), nrow = 3, ncol = 1, scales = "free_x")+
                  facetted_pos_scales(y=y_scales)
                  
                  
y_scales = list(
          scale_y_continuous(limits = c(60,100), breaks = seq(60,100, by=10)),
          scale_y_continuous(limits = c(60,100), breaks = seq(60,100, by=10)),
          scale_y_continuous(limits = c(60,100), breaks = seq(60,100, by=10))
)


combined_year$`Indicator Name`