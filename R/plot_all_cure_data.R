#Laura F. White
#Boston University
#Pre-chemotherapy TB Analysis

##############################################################################
# This program generates the plots used to show all the cure data reported in 
# the supplemental materials
##############################################################################

options(scipen=999)
options(digits = 10)

rm(list = ls())
source("R/utils.R")
reload_source()

# Reading in the study_id correspondence table
studyid <- read.csv("data/study_id.csv")

# Reading in cure data
cureData <- read.csv("data/cure_data_all.csv")

# add in author and date of study
cureData2 <- inner_join(cureData,studyid,by="study_id")

# plot only those without severity
ggplot(cureData2[cureData2$severity=="None",],
       aes(x=interval_r,y=cureRate,group=cohort_id,color=first_author))+
  geom_point()+geom_line()+xlim(0,10)+ylim(0,1)+
  labs(title="Studies without severity reported",
       x="Year since diagnosis",y="Probability of Natural Recovery",color="Author")

# plot studies with severity reported
ggplot(cureData2[cureData2$severity!="None",],
       aes(x=interval_r,y=cureRate,group=cohort_id,color=severity))+
  geom_point()+geom_line()+xlim(0,10)+ylim(0,1)+
  labs(title="Studies with severity reported",
       x="Year since diagnosis",y="Probability of Natural Recovery",color="Severity")
