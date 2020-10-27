#Sarah V. Leavitt
#Boston University
#Pre-chemotherapy TB Analysis

##############################################################################
# This program transforms the life-table study data extracted from the 
# publications and transforms it into individual-level data
##############################################################################

options(scipen=999)
options(digits = 10)

rm(list = ls())
source("R/utils.R")
reload_source()

## Reading in all of the datasets
dataList <- read_excel_allsheets("data/pre_chemo_data.xlsx")
dataList$`Data dictionary` <- NULL


#### Mortality Data ####

#Converting the study data to individual mortality data
indAll <- map_dfr(dataList, studyToInd, outcome = "mortality")

#Removing people who are censored at the study start
indAll <- indAll %>% filter(!(death == 0 & interval_l == 0))

#Formatting individual dataset
mortalityData <- indAll %>%
  mutate(start_type = ifelse(start_type == "Unknown", "Entry", start_type),
         #Adding overall average length of stay to bring start_type = 'Exit' start_type = 'Entry'
         interval_l = ifelse(start_type == "Exit", interval_l + 165/365, interval_l),
         interval_r = ifelse(start_type == "Exit", interval_r + 165/365, interval_r),
         severity = ifelse(severity %in% c("Incipient disease",
                                                 "Minimal disease",
                                                 "Stage I",
                                                 "TB plus I (early)"), "Minimal",
                           ifelse(severity %in% c("Moderately advanced disease",
                                                        "Stage II",
                                                        "TB plus II (intermediate)"), "Moderate",
                                  ifelse(severity %in% c("Far advanced disease",
                                                               "Stage III",
                                                               "TB plus III (advanced)"), "Advanced",
                                         "Unknown"))),
         severity = factor(severity, levels = c("Minimal", "Moderate", "Advanced", "Unknown"))) %>%
  unite(study_sev, study_id, severity, remove = FALSE) %>%
  arrange(study_sev)

write.csv(mortalityData, "data/mortality_data.csv", row.names = FALSE)



#### Cure Data ####

cureList <- list(dataList$`1029_1055`, dataList$`1029_1056`, dataList$`48_1000_1029`,
                 dataList$`45`, dataList$`67`)

## Converting the study data to individual mortality data
cureData <- map_dfr(cureList, studyToInd, outcome = "cure", timepoints = 3)

write.csv(cureData, "data/cure_data.csv", row.names = FALSE)


