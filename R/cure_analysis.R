#Sarah V. Leavitt
#Boston University
#Pre-chemotherapy TB Analysis

##############################################################################
# This program performs the cure Bayesian meta-analysis and creates tables
# to summarize the results
##############################################################################

options(scipen=999)
options(digits = 10)
set.seed(150183)

rm(list = ls())
source("R/utils.R")
reload_source()

#Reading in the study_id correspondence table
studyid <- read.csv("data/study_id.csv")

#Reading in cure data
cureData <- read.csv("data/cure_data.csv")

cureData <- cureData %>%
  mutate(severity = factor(severity, levels = c("Far advanced disease", 
                                                "Moderately advanced disease",
                                                "Minimal disease"),
                           labels = c("Advanced", "Moderate", "Minimal")),
         severityAdv = ifelse(severity != "Advanced", "Min/Mod", severity))



#### Bayesian Logistic Model------------------------------------------------------------------------

#Model
m_cure <- function(){
  
  for (i in 1:nStudy){
    
    #Distribution of cure counts for each severity level
    cAdv[i] ~ dbin(p_0[i], nAdv[i])
    cMod[i] ~ dbin(p_1[i], nMod[i])
    cMin[i] ~ dbin(p_2[i], nMin[i])
    
    #Logit models
    logit(p_0[i]) <- u[i]
    logit(p_1[i]) <- u[i] + bmod
    logit(p_2[i]) <- u[i] + bmin
    
    #Study-level odds for each severity level
    oddsAdv[i] <- exp(u[i])
    oddsMod[i] <- exp(u[i] + bmod)
    oddsMin[i] <- exp(u[i] + bmin)
    
    #Random effects
    u[i] ~ dnorm(alpha, tau)
  }
  
  #Overall odds ratios
  ORmod <- exp(bmod)
  ORmin <- exp(bmin)
  
  #Priors
  bmod ~ dnorm(0, 0.001)
  bmin ~ dnorm(0, 0.001)
  alpha ~ dnorm(0,0.001)
  tau ~ dgamma(1,1)
  
  #Random effect variance
  theta = 1/tau
}

#Parameters to track
par_cure <- c("alpha", "bmod", "bmin", "theta",
              "ORmod", "ORmin", "oddsAdv", "oddsMod", "oddsMin")



#### Three year only analysis, all studies ---------------------------------------------------------

## Running the model
cureData_3 <- cureData %>% filter(study_id != "79_1023")

cureAggregate_3 <- cureData_3 %>%
  group_by(study_id) %>%
  summarize(nMin = sum(severity == "Minimal"),
            nMod = sum(severity == "Moderate"),
            nAdv = sum(severity == "Advanced"),
            cMin = sum(severity == "Minimal" & cure == 1),
            cMod = sum(severity == "Moderate" & cure == 1),
            cAdv = sum(severity == "Advanced" & cure == 1),
            .groups = "drop")

#Data
dt_3 <- list(nStudy = length(unique(cureData_3$study_id)),
             nMin = cureAggregate_3$nMin,
             nMod = cureAggregate_3$nMod,
             nAdv = cureAggregate_3$nAdv,
             cMin = cureAggregate_3$cMin,
             cMod = cureAggregate_3$cMod,
             cAdv = cureAggregate_3$cAdv
)

#Fitting the model
fit_3 <- jags(data = dt_3, model.file = m_cure,
                parameters.to.save = par_cure,
                n.iter = 11000, n.burnin = 1000,
                n.chains = 1, n.thin = 20)

#Extracting data
mcmc_3 <- as.mcmc(fit_3)
eval_3 <- mcmc_3
summary(eval_3)$quantiles

png("Figures/xyplot_cure_3year.png")
xyplot(eval_3[, c("alpha", "bmod", "bmin", "theta")])
dev.off()
png("Figures/autocorr_cure_3year.png")
autocorr.plot(eval_3[, c("alpha", "bmod", "bmin", "theta")])
dev.off()




#### Three and four year analysis, all studies -----------------------------------------------------

## Running the model

cureAggregate_4 <- cureData %>%
  group_by(study_id) %>%
  summarize(nMin = sum(severity == "Minimal"),
            nMod = sum(severity == "Moderate"),
            nAdv = sum(severity == "Advanced"),
            cMin = sum(severity == "Minimal" & cure == 1),
            cMod = sum(severity == "Moderate" & cure == 1),
            cAdv = sum(severity == "Advanced" & cure == 1),
            .groups = "drop")

#Data
dt_4 <- list(nStudy = length(unique(cureData$study_id)),
             nMin = cureAggregate_4$nMin,
             nMod = cureAggregate_4$nMod,
             nAdv = cureAggregate_4$nAdv,
             cMin = cureAggregate_4$cMin,
             cMod = cureAggregate_4$cMod,
             cAdv = cureAggregate_4$cAdv
)

#Fitting the model
fit_4 <- jags(data = dt_4, model.file = m_cure,
              parameters.to.save = par_cure,
              n.iter = 11000, n.burnin = 1000,
              n.chains = 1, n.thin = 20)

#Extracting data
mcmc_4 <- as.mcmc(fit_4)
eval_4 <- mcmc_4
summary(eval_4)$quantiles

png("Figures/xyplot_cure_4year.png")
xyplot(eval_4[, c("alpha", "bmod", "bmin", "theta")])
dev.off()
png("Figures/autocorr_cure_4year.png")
autocorr.plot(eval_4[, c("alpha", "bmod", "bmin", "theta")])
dev.off()



#### Three year only analysis, all studies ---------------------------------------------------------

## Running the model
cureData_3sub <- cureData %>% filter(study_id %in% c("1029", "45"))

cureAggregate_3sub <- cureData_3sub %>%
  group_by(study_id) %>%
  summarize(nMin = sum(severity == "Minimal"),
            nMod = sum(severity == "Moderate"),
            nAdv = sum(severity == "Advanced"),
            cMin = sum(severity == "Minimal" & cure == 1),
            cMod = sum(severity == "Moderate" & cure == 1),
            cAdv = sum(severity == "Advanced" & cure == 1),
            .groups = "drop")

#Data
dt_3sub <- list(nStudy = length(unique(cureData_3sub$study_id)),
             nMin = cureAggregate_3sub$nMin,
             nMod = cureAggregate_3sub$nMod,
             nAdv = cureAggregate_3sub$nAdv,
             cMin = cureAggregate_3sub$cMin,
             cMod = cureAggregate_3sub$cMod,
             cAdv = cureAggregate_3sub$cAdv
)

#Fitting the model
fit_3sub <- jags(data = dt_3sub, model.file = m_cure,
              parameters.to.save = par_cure,
              n.iter = 11000, n.burnin = 1000,
              n.chains = 1, n.thin = 20)

#Extracting data
mcmc_3sub <- as.mcmc(fit_3sub)
eval_3sub <- mcmc_3sub
summary(eval_3sub)$quantiles

png("Figures/xyplot_cure_3subyear.png")
xyplot(eval_3sub[, c("alpha", "bmod", "bmin", "theta")])
dev.off()
png("Figures/autocorr_cure_3subyear.png")
autocorr.plot(eval_3sub[, c("alpha", "bmod", "bmin", "theta")])
dev.off()





#### Tables of results -----------------------------------------------------------------------------

make_cure_tab <- function(eval_tab, aggregate_tab){
  
  #Odds ratios for minimal and moderate vs. advanced
  or <- as.data.frame(summary(eval_tab[, c("ORmin", "ORmod")])$quantiles) %>%
    mutate(rownames = row.names(.),
           severity = ifelse(grepl("min", rownames), "Minimal", "Moderate"),
           OR_CI = paste0(round(`50%`, 2), " (", round(`2.5%`, 2), ", ", round(`97.5%`, 2), ")")) %>%
    select(severity, OR_CI)
  
  #Table of counts per study
  cureTab <- aggregate_tab %>%
    mutate(study_id = as.character(study_id)) %>%
    left_join(studyid, by = "study_id") %>%
    mutate(pMin = 100 * round(cMin / nMin, 2),
           pMod = 100 * round(cMod / nMod, 2),
           pAdv = 100 * round(cAdv / nAdv, 2),
           Min_cure = ifelse(nMin == 0, "-", paste0(cMin, " (", pMin, "%)")),
           Mod_cure = paste0(cMod, " (", pMod, "%)"),
           Adv_cure = paste0(cAdv, " (", pAdv, "%)")) %>%
    select(first_author, Min_total = nMin, Min_cure, Mod_total = nMod, Mod_cure,
           Adv_total = nAdv, Adv_cure) %>%
    arrange(first_author)
  
  return(list(or, cureTab))
}


## Three year, all studies
tabs_3 <- make_cure_tab(eval_3, cureAggregate_3)
or_3 <- tabs_3[[1]]
cureTab_3 <- tabs_3[[2]]

## Three and four year, all studies
tabs_4 <- make_cure_tab(eval_4, cureAggregate_4)
or_4 <- tabs_4[[1]]
cureTab_4 <- tabs_4[[2]]

## Three year only, US post-1930s
tabs_3sub <- make_cure_tab(eval_3sub, cureAggregate_3sub)
or_3sub <- tabs_3sub[[1]]
cureTab_3sub <- tabs_3sub[[2]]


