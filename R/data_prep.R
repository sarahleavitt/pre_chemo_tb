#Sarah V. Leavitt
#TB Duration of Infectiousness

## Data Prep ##

setwd("~/Boston University/Duration_of_Infectiousness/duration_code")
rm(list = ls())
options(scipen=999)
options(digits = 10)

source("duration_functions.R")
reload_source()

## Reading in all of the datasets
dataList <- read_excel_allsheets("../Duration data_20Oct2020.xlsx")

## Cleaning list of datasets to only include sheets I am interested in
dataList$Documentation <- NULL
dataList$`Data dictionary` <- NULL
dataList$`12` <- NULL
dataList$`45` <- NULL
dataList$`67` <- NULL


#### Mortality Data ####

## Converting the study data to individual mortality data
indAll <- map_dfr(dataList, studyToInd, outcome = "mortality")

#Removing people who are censored at the study start (no observations)
indAll <- indAll %>% filter(!(eventIC == 0 & interval_l == 0))

## Formatting individual dataset
metaData <- indAll %>%
  mutate(start_type = ifelse(start_type == "Unknown", "Entry", start_type),
         #Adding overall average length of stay to bring start_type = 'Exit' start_type = 'Entry'
         interval_l = ifelse(start_type == "Exit", interval_l + 165/365, interval_l),
         interval_r = ifelse(start_type == "Exit", interval_r + 165/365, interval_r),
         study_id = ifelse(grepl("1029", paper_id), "1029", paper_id),
         severity = ifelse(stratified_var %in% c("Incipient disease",
                                                 "Minimal disease",
                                                 "Stage I",
                                                 "TB plus I (early)"), "Minimal",
                           ifelse(stratified_var %in% c("Moderately advanced disease",
                                                        "Stage II",
                                                        "TB plus II (intermediate)"), "Moderate",
                                  ifelse(stratified_var %in% c("Far advanced disease",
                                                               "Stage III",
                                                               "TB plus III (advanced)"), "Advanced",
                                         "Unknown"))),
         severity = factor(severity, levels = c("Minimal", "Moderate", "Advanced", "Unknown"))) %>%
  unite(study_sev, study_id, severity, remove = FALSE) %>%
  #Removing the cohorts that have duplicated data
  filter(!cohort_id == "75_1") %>%
  arrange(study_sev)

write.csv(metaData, "../metaData.csv", row.names = FALSE)



#### Cure Data ####

cureList <- list(dataList$`1029_1055`, dataList$`1029_1056`, dataList$`48_1000_1029`,
                 dataList$`45_2`, dataList$`67_2`)

## Converting the study data to individual mortality data
indCure <- map_dfr(cureList, studyToInd, outcome = "cure")

write.csv(indCure, "../cureData.csv", row.names = FALSE)


