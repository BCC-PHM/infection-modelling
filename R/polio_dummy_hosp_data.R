# Creating dummy H(t) data for polio
library("dplyr")

# hospital_data = readxl::read_excel(
#   paste(
#     "//svwvap1126.addm.ads.brm.pri/PHSensitive$/Intelligence/2. Requests/REQ3280 - HES Inpatients data on polio, flu and rubella for Justin",
#     "/REQ3280_Inpatients_data_10years_on_polio_flu_and_rubella_BirminghamResidents.xlsx",
#     sep = ""
#   )
# )

#filter Polio only and only data in 2022/23
Poli_data  = hospital_data %>% 
  filter(Condition == "Polio" & Year == "2022/2023") %>%
  mutate(
    Date =  as.Date(Date_of_Admission, format = "%d%m%Y"),
    Day = as.numeric(Date - as.Date("01/04/2022", format = "%d/%m/%Y"))
  ) %>%
  select(Day) %>%
  arrange(Day) %>%
  as.list()

H_data <- rep(0, 365)

for (item in Poli_data$Day) {
  # Random length between 7 and 14 days
  hospital_stay = sample(7:14, 1)
  
  # From day in data
  d1 = item
  
  # To day + stay length or 365 
  d2 = min(c(item + hospital_stay, 365))
  
  # Add one to these elements
  H_data[d1:d2] = H_data[d1:d2] + 1
}

# create new data frame
dummy_data <- data.frame(
  Date = as.Date("01/04/2022", format = "%d/%m/%Y") + 0:364,
  hospitalised = H_data
)

# Save data
writexl::write_xlsx(dummy_data, "../data/Polio/dummy_data.xlsx")
