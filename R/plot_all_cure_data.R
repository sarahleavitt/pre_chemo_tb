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

### Create plots of the data by severity --------------------------------------------
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

ggplot(cureData2[cureData2$severity!="None",],
       aes(x=interval_r,y=cureRate,group=cohort_id,color=first_author))+
  geom_point(aes(shape=severity))+geom_line()+xlim(0,10)+ylim(0,1)+
  labs(title="Studies with severity reported",
       x="Year since diagnosis",y="Probability of Natural Recovery",color="Author",shape="Severity")

### Create a data frame with all the data -------------------------------------------
cureData3 <- cureData2[,c('study_id','severity','interval_r','n','cureRate','first_author','year')]

# round interval_r so we are looking at cureRate within a one year interval of the integer time
# also remove interval data that is greater than 10
cureData3$interval_r <- round(cureData3$interval_r)

cureData4 <- subset(cureData3,cureData3$interval_r<=10)
cureSummary <- spread(cureData4,key="interval_r",value="cureRate")

write.csv(cureSummary,file="data/cureDataSummary.csv")


