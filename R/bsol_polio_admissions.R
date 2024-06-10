library(ggplot2)
library(dplyr)
library(lubridate)
library(tidyr)
library("scales")

#load data
PolioAdmissions <- readxl::read_excel("FILEPATH")

# view data types
str(PolioAdmissions)

#creation of year column
PolioAdmissions <- PolioAdmissions %>%
  mutate(AdmissionYear = year(AdmissionDate))

#Creation of Aggregated Diagnosis Description Column
PolioAdmissions$AggregatedDiagnosisDescription <- ifelse(PolioAdmissions$DiagnosisDescription %in% c("Acute nonparalytic poliomyelitis", "Acute paralytic poliomyelitis, other and unspecified", "Acute poliomyelitis, unspecified"), "Acute", PolioAdmissions$DiagnosisDescription)

#################### Total costs #################
Cost <- sum(PolioAdmissions$Cost)
print(Cost)

######################number of admissions #######################################
# whole dataset total number
WholeDataSetTotalAdmissions <- PolioAdmissions %>%
  summarise(
    Total = n()
  )

DiagnosisWholeDataSetAdmissions <-  PolioAdmissions %>%
  group_by(AggregatedDiagnosisDescription) %>%
  summarise(
    Count = n()
    )

# Summarize total admissions per year
TotalAdmissions <- PolioAdmissions %>%
  filter(AdmissionYear != 2024) %>% #do not include 2024 as not a complete year
  group_by(AdmissionYear) %>%
  summarise(
    Total = n()
  )
# Summarize admissions per year for each diagnosis
DiagnosisAdmissions <- PolioAdmissions %>%
  filter(AdmissionYear != 2024) %>%
  group_by(AdmissionYear, AggregatedDiagnosisDescription) %>%
  summarise(
    Count = n()
  ) %>%
  pivot_wider(names_from = AggregatedDiagnosisDescription, values_from = Count, values_fill = list(Count = 0))

# Combine the total admissions with the diagnosis-specific admissions
YearNumberAdmissions <- TotalAdmissions %>%
  left_join(DiagnosisAdmissions, by = "AdmissionYear")

# Convert the summary dataframe to long format for plotting
LongYearNumberAdmissions <- YearNumberAdmissions %>%
  pivot_longer(cols = -AdmissionYear, names_to = "Diagnosis", values_to = "Count") %>%
  mutate(
    plot_value =
      case_when(Count > 0 & Count < 5 ~ 2.5, TRUE ~ Count)
    )

# Add a column to mark suppressed points
LongYearNumberAdmissions <- LongYearNumberAdmissions %>%
  mutate(Suppressed = ifelse(plot_value == 2.5, TRUE, FALSE))

#Plot the line graph
AdmissionsNumberPlot <- ggplot(data=LongYearNumberAdmissions, aes(x=AdmissionYear, y=plot_value, color=Diagnosis)) +
  geom_line() +
  geom_point(aes(shape = Suppressed)) +
  theme_bw() +
  scale_y_continuous(
    # Force y-axis to start at zero
    expand = c(0, 0), limits = c(0, 70)) +
  labs(title = "Polio Hospital Admissions per Year",
       x = "Year",
       y = "Number of Admissions") +
  annotate("text", x = 2014, y = 60, size = 3,hjust = 0,
           label = "Note:Values relating to less than 5 people\nhave been supressed for data security.") +
  scale_shape_manual(values = c(16,17),  # Only keeping the triangle shape
                     labels = c("Not Suppressed", "Small number suppressed")) +
  labs(shape = NULL)

AdmissionsNumberPlot

ggsave("figures/AdmissionsNumberPlot.png",
       AdmissionsNumberPlot, width = 8, height = 4)

#mean number of  admissions
mean(YearNumberAdmissions$Acute)
mean(YearNumberAdmissions$`Postpolio syndrome`)
mean(YearNumberAdmissions$`Sequelae of poliomyelitis`)
####################### cost per year ##################

# Summarize total cost per year
TotalCost <- PolioAdmissions %>%
  filter(AdmissionYear != 2024) %>%
  group_by(AdmissionYear) %>%
  summarise(TotalCost =sum(Cost))

# Summarize cost per year for each diagnosis
DiagnosisCost <- PolioAdmissions %>%
  filter(AdmissionYear != 2024) %>% #do not include 2024 as not a complete year
  group_by(AdmissionYear, AggregatedDiagnosisDescription) %>%
  summarise(TotalCost =sum(Cost)) %>%
  pivot_wider(names_from = AggregatedDiagnosisDescription, values_from = TotalCost, values_fill = list(TotalCost = 0))

# Combine the total admissions with the diagnosis-specific admissions
YearCostAdmissions <- TotalCost %>%
  left_join(DiagnosisCost, by = "AdmissionYear")

# Convert the summary dataframe to long format for plotting
LongYearCostAdmissions <- YearCostAdmissions %>%
  pivot_longer(cols = -AdmissionYear, names_to = "Diagnosis", values_to = "TotalCost")

#Plot the line graph
AdmissionsCostPlot <- ggplot(data=LongYearCostAdmissions, aes(x=AdmissionYear, y=TotalCost, color=Diagnosis)) +
  geom_line() +
  geom_point(aes(shape = Diagnosis)) +
  theme_bw() +
  scale_y_continuous(labels = comma,
                     # Force y-axis to start at zero
                     expand = c(0, 0), limits = c(0, 300000)) +
  labs(title = "Cost of Polio Hospital Admissions per Year, Birmingham and Solihull",
       x = "Year",
       y = "Cost of Admissions (£)") +
  theme(plot.title = element_text(size=11), axis.title=element_text(size=9))

AdmissionsCostPlot

ggsave("figures/AdmissionsCostPlot.png",
       AdmissionsCostPlot, width = 8, height = 4)

########################lengths of admission ############################
# Calculate the number of days between the two dates
PolioAdmissions <- PolioAdmissions %>%
  mutate(days_diff = as.numeric(difftime(DischargeDate, AdmissionDate, units = "days")))

#Median admission length Acute, PPS and Sequale
DiagnosesAdmissionLengthMedian <- PolioAdmissions %>%
  group_by(AggregatedDiagnosisDescription) %>%
  summarise(MedianAdmissionLength =median(days_diff))

#histogram of admission length
AdmissionLengthHist <- ggplot(data= PolioAdmissions, aes(x = days_diff)) +
  geom_histogram(binwidth=1, colour="black", fill="#B8CCE4") +
  #add vertical line showing the median
  geom_vline(aes(xintercept=median(days_diff)),
             color="black", linetype="dashed", lwd=1) +
  annotate("text", x=2.1, y=165, label="Median", angle=0, size=2.5) +
  #label y axis
  labs( y = "Count", x = "Length of admissions") +
  theme_bw() +
  scale_y_continuous(
    # Force y-axis to start at zero
    expand = c(0, 0), limits = c(0, 200)) +
  labs(title = "")
#view histogram
print(AdmissionLengthHist)


#################################### age ##################
# all admissions (including individuals admitted multiple times )
AgeHist <- ggplot(data= PolioAdmissions, aes(x = AgeOnAdmission)) +
  geom_histogram(binwidth=1, colour="black", fill="#B8CCE4", linewidth=0.2) +
  #add vertical line showing the median
  geom_vline(aes(xintercept=median(AgeOnAdmission)),
             color="black", linetype="dashed", lwd=0.5) +
  annotate("text", x=79, y=15, label="Polio Median", angle=0, size=2.5) +
  #label y axis
  labs( y = "Count", x = "Age on Admission") +
  theme_bw() +
  scale_y_continuous(
    # Force y-axis to start at zero
    expand = c(0, 0), limits = c(0, 20)) +
  labs(title = "Histogram showing age on admission, Birmingham and Solihull 2014 to 2024") +
  facet_wrap(~AggregatedDiagnosisDescription) +
  theme(plot.title = element_text(size=10), axis.text=element_text(size=8), axis.title=element_text(size=8))
#view histogram
print(AgeHist)
#save histogram
ggsave("figures/AgeHist.png",
       AgeHist, width = 12, height = 6)

DiagnosesAgeMedian <- PolioAdmissions %>%
  group_by(AggregatedDiagnosisDescription) %>%
  summarise(MedianAge =median(AgeOnAdmission))

median(PolioAdmissions$AgeOnAdmission)

#IQR(PolioAdmissions$AgeOnAdmission)

#quartileAge <- quantile(PolioAdmissions$AgeOnAdmission, probs = c(0.25,0.5,0.75))
#### number of individuals ##############
n_distinct(PolioAdmissions$NHSNumber)

UniqueIndividuals <- PolioAdmissions %>%
  filter(NHSNumber != "NULL") %>%
  group_by(NHSNumber) %>%
  summarise(
    Count = n()
  )


