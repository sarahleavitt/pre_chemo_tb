#Sarah V. Leavitt
#TB Duration of Infectiousness

## Mortality Analysis ##

setwd("~/Boston University/Duration_of_Infectiousness/duration_code")
rm(list = ls())
options(scipen=999)
options(digits = 10)

source("duration_functions.R")
reload_source()

library(R2jags)
library(lattice)
library(mcmcplots)

metaData <- read.csv("../metaData.csv")



########################## All-cause Mortality Analysis #############################

#interval = 1 implies x1 < t <= x2
#For some reason dinterval() isn't working as it should so make all obs have interval = 1
#and set x2 to be 10000 (close enough to infinity) for right censored
metaData <- metaData %>% 
  mutate(study_sev_num = as.numeric(factor(study_sev)),
         study_id_num = as.numeric(factor(study_id)),
         time = ifelse(event == 0, NA, interval_l),
         interval = 1,
         x1 = ifelse(eventIC == 2, 0, interval_l),
         x2 = ifelse(eventIC == 3, interval_r,
                     ifelse(eventIC == 0, 10000, interval_l)),
         sev_mod = as.numeric(severity == "Moderate"),
         sev_adv = as.numeric(severity == "Advanced"),
         sev_unk = as.numeric(severity == "Unknown"))



######## Full model ########

dt_all <- list(N = nrow(metaData),
                interval = metaData$interval,
                lim = cbind(metaData$x1, metaData$x2),
                time = rep(NA, nrow(metaData)),
                n_frail = length(unique(metaData$study_sev_num)),
                frail = metaData$study_sev_num
                
)

par_all <- c("mu", "theta", "sdlog", "med_all", "meanlog", "med_ind",
             "pred_all", "pred1", "pred5", "pred10")

m_all <- function(){
  
  for(j in 1:n_frail){
    
    #Frailty for each cohort
    ui[j] ~ dnorm(mu, tau)
    #Meanlog of lognormal distribution for each cohort
    meanlog[j] <- ui[j]
    #Median of lognormal distribution for each cohort
    med_ind[j] <- exp(meanlog[j])
    
    #Prediction of 1/5/10 year survival
    pred1[j] <- 1 - plnorm(1, meanlog[j], taulog)
    pred5[j] <- 1 - plnorm(5, meanlog[j], taulog)
    pred10[j] <- 1 - plnorm(10, meanlog[j], taulog)
  }
  
  for(i in 1:N){
    #Setting up the survival model
    interval[i] ~ dinterval(time[i], lim[i, ])
    time[i] ~ dlnorm(meanlog[frail[i]], taulog)
  }
  
  #Priors
  mu ~ dnorm(0, 0.0001)
  taulog ~ dgamma(1, 1)
  tau ~ dgamma(1, 1)
  
  #Prediction for confidence bands
  for(t in 1:30){
    pred_all[t] <- 1 - plnorm(t, mu, taulog)
  }
  
  ## Derived parameters ##
  
  #Variance of the frailty distribution
  theta <- 1/tau
  #sdlog of all survival densities (cohort-specific and overall)
  sdlog <- sqrt(1/taulog)
  #Overall median
  med_all <- exp(mu)
}


fit_all <- jags(data = dt_all, model.file = m_all,
                 parameters.to.save = par_all,
                n.iter = 11000, n.burnin = 1000,
                 n.chains = 1, n.thin = 20)

mcmc_all <- as.mcmc(fit_all)
eval_all <- mcmc_all[, c("mu", "theta", "sdlog")]
#summary(eval_all)
#xyplot(eval_all)
#autocorr.plot(eval_all)
#densityplot(eval_all)

png("../Figures/xyplot_all.png")
xyplot(eval_all)
dev.off()
png("../Figures/autocorr_all.png")
autocorr.plot(eval_all)
dev.off()

res_all <- as.data.frame(summary(mcmc_all)$quantiles)


###### Fixed Effect for Severity #######

metaData_sev <- metaData %>%
  filter(severity != "Unknown") %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)),
         study_id_num = as.numeric(factor(study_id)))

cohort_data <- metaData_sev %>%
  group_by(study_sev_num) %>%
  summarize(study_id_num = first(study_id_num),
            sev_mod = first(sev_mod),
            sev_adv = first(sev_adv),
            sev_unk = first(sev_unk),
            .groups = "drop")

dt_sev <- list(N = nrow(metaData_sev),
               interval = metaData_sev$interval,
               lim = cbind(metaData_sev$x1, metaData_sev$x2),
               time = rep(NA, nrow(metaData_sev)),
               n_frail = length(unique(metaData_sev$study_id_num)),
               frail = metaData_sev$study_id_num,
               sev_mod = metaData_sev$sev_mod,
               sev_adv = metaData_sev$sev_adv,
               n_study_sev = nrow(cohort_data),
               frail2 = cohort_data$study_id_num,
               study_sev_mod = cohort_data$sev_mod,
               study_sev_adv = cohort_data$sev_adv
)

par_sev <- c("theta", "sdlog", "alpha", "bmod", "badv",
             "meanlog_ind", "meanlog_min", "meanlog_mod", "meanlog_adv",
             "med_ind", "med_min", "med_mod", "med_adv",
             "pred1", "pred5", "pred10", "pred_min", "pred_mod", "pred_adv")

m_sev <- function(){
  
  #Frailty for each cohort
  for(j in 1:n_frail){
    ui[j] ~ dnorm(alpha, tau)
  }
  
  for(i in 1:N){
    #Meanlog of lognormal distribution for each individual
    meanlog[i] <- ui[frail[i]] + bmod*sev_mod[i] + badv*sev_adv[i]
    #Setting up the survival model
    interval[i] ~ dinterval(time[i], lim[i, ])
    time[i] ~ dlnorm(meanlog[i], taulog)
  }
  
  #Priors
  alpha ~ dnorm(0, 0.0001)
  bmod ~ dnorm(0, 0.0001)
  badv ~ dnorm(0, 0.0001)
  taulog ~ dgamma(1, 1)
  tau ~ dgamma(1, 1)
  
  #Variance of the frailty distribution
  theta <- 1/tau
  #sdlog of all survival densities (cohort-specific and overall)
  sdlog <- sqrt(1/taulog)
  #Overall median survival for each severity
  meanlog_min <- alpha
  meanlog_mod <- alpha + bmod
  meanlog_adv <- alpha + badv
  med_min <- exp(alpha)
  med_mod <- exp(alpha + bmod)
  med_adv <- exp(alpha + badv)
  
  #Prediction for confidence bands
  for(t in 1:30){
    pred_min[t] <- 1 - plnorm(t, meanlog_min, taulog)
    pred_mod[t] <- 1 - plnorm(t, meanlog_mod, taulog)
    pred_adv[t] <- 1 - plnorm(t, meanlog_adv, taulog)
  }
  
  #Cohort-specific meanlog of the lognormal density and median survival
  for(k in 1:n_study_sev){
    meanlog_ind[k] <- ui[frail2[k]] + bmod*study_sev_mod[k] + badv*study_sev_adv[k]
    med_ind[k] <- exp(meanlog_ind[k])
    
    #Prediction of 1/5/10 year survival
    pred1[k] <- 1 - plnorm(1, meanlog_ind[k], taulog)
    pred5[k] <- 1 - plnorm(5, meanlog_ind[k], taulog)
    pred10[k] <- 1 - plnorm(10, meanlog_ind[k], taulog)
  }
}

fit_sev <- jags(data = dt_sev, model.file = m_sev,
                 parameters.to.save = par_sev,
                n.iter = 11000, n.burnin = 1000,
                 n.chains = 1, n.thin = 20)

mcmc_sev <- as.mcmc(fit_sev)
eval_sev <- mcmc_sev[, c("alpha", "bmod", "badv", "theta", "sdlog")]
#summary(eval_sev)
#xyplot(eval_sev)
#autocorr.plot(eval_sev)

png("../Figures/xyplot_sev.png")
xyplot(eval_sev)
dev.off()
png("../Figures/autocorr_sev.png")
autocorr.plot(eval_sev)
dev.off()

res_sev <- as.data.frame(summary(mcmc_sev)$quantiles)




######## Sanatorium only model ########

san <- metaData %>%
  filter(sanatorium == 1) %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)))

dt_san <- list(N = nrow(san),
               interval = san$interval,
               lim = cbind(san$x1, san$x2),
               time = rep(NA, nrow(san)),
               n_frail = length(unique(san$study_sev_num)),
               frail = san$study_sev_num
               
)

fit_san <- jags(data = dt_san, model.file = m_all,
                parameters.to.save = par_all,
                n.iter = 11000, n.burnin = 100,
                n.chains = 1, n.thin = 20)

mcmc_san <- as.mcmc(fit_san)
eval_san <- mcmc_san[, c("mu", "theta", "sdlog")]
#summary(eval_san)
#xyplot(eval_san)
#autocorr.plot(eval_san)

png("../Figures/xyplot_san.png")
xyplot(eval_san)
dev.off()
png("../Figures/autocorr_san.png")
autocorr.plot(eval_san)
dev.off()

res_san <- as.data.frame(summary(mcmc_san)$quantiles)



######## Non-sanatorium only model ########

nosan <- metaData %>%
  filter(sanatorium == 0) %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)))

dt_nosan <- list(N = nrow(nosan),
               interval = nosan$interval,
               lim = cbind(nosan$x1, nosan$x2),
               time = rep(NA, nrow(nosan)),
               n_frail = length(unique(nosan$study_sev_num)),
               frail = nosan$study_sev_num
               
)

fit_nosan <- jags(data = dt_nosan, model.file = m_all,
                parameters.to.save = par_all,
                n.iter = 11000, n.burnin = 1000,
                n.chains = 1, n.thin = 20)

mcmc_nosan <- as.mcmc(fit_nosan)
eval_nosan <- mcmc_nosan[, c("mu", "theta", "sdlog")]
#summary(eval_nosan)
#xyplot(eval_nosan)
#autocorr.plot(eval_nosan)

png("../Figures/xyplot_nosan.png")
xyplot(eval_nosan)
dev.off()
png("../Figures/autocorr_nosan.png")
autocorr.plot(eval_nosan)
dev.off()

res_nosan <- as.data.frame(summary(mcmc_nosan)$quantiles)




############################### TB Mortailty Analysis ###############################

#interval = 1 implies x1 < t <= x2
#For some reason dinterval() isn't working as it should so make all obs have interval = 1
#and set x2 to be 10000 (close enough to infinity) for right censored
metaData_tb <- metaData %>% 
  filter(!paper_id %in% c("12", "65", "91")) %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)),
         study_id_num = as.numeric(factor(study_id)),
         time = ifelse(death_tb == 0, NA, interval_l),
         interval = 1,
         x1 = ifelse(death_tbIC == 2, 0, interval_l),
         x2 = ifelse(death_tbIC == 3, interval_r,
                     ifelse(death_tbIC == 0, 10000, interval_l)),
         sev_mod = as.numeric(severity == "Moderate"),
         sev_adv = as.numeric(severity == "Advanced"),
         sev_unk = as.numeric(severity == "Unknown"))



######## Full model ########

dt_all_tb <- list(N = nrow(metaData_tb),
               interval = metaData_tb$interval,
               lim = cbind(metaData_tb$x1, metaData_tb$x2),
               time = rep(NA, nrow(metaData_tb)),
               n_frail = length(unique(metaData_tb$study_sev_num)),
               frail = metaData_tb$study_sev_num
)

fit_all_tb <- jags(data = dt_all_tb, model.file = m_all,
                parameters.to.save = par_all,
                n.iter = 11000, n.burnin = 1000,
                n.chains = 1, n.thin = 20)

mcmc_all_tb <- as.mcmc(fit_all_tb)
eval_all_tb <- mcmc_all_tb[, c("mu", "theta", "sdlog")]
#summary(eval_all_tb)
#xyplot(eval_all_tb)
#autocorr.plot(eval_all_tb)

png("../Figures/xyplot_all_tb.png")
xyplot(eval_all_tb)
dev.off()
png("../Figures/autocorr_all_tb.png")
autocorr.plot(eval_all_tb)
dev.off()

res_all_tb <- as.data.frame(summary(mcmc_all_tb)$quantiles)


###### Fixed Effect for Severity #######

metaData_sev_tb <- metaData_tb %>%
  filter(severity != "Unknown") %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)),
         study_id_num = as.numeric(factor(study_id)))

cohort_data_tb <- metaData_sev_tb %>%
  group_by(study_sev_num) %>%
  summarize(study_id_num = first(study_id_num),
            sev_mod = first(sev_mod),
            sev_adv = first(sev_adv),
            sev_unk = first(sev_unk))


dt_sev_tb <- list(N = nrow(metaData_sev_tb),
               interval = metaData_sev_tb$interval,
               lim = cbind(metaData_sev_tb$x1, metaData_sev_tb$x2),
               time = rep(NA, nrow(metaData_sev_tb)),
               n_frail = length(unique(metaData_sev_tb$study_id_num)),
               frail = metaData_sev_tb$study_id_num,
               sev_mod = metaData_sev_tb$sev_mod,
               sev_adv = metaData_sev_tb$sev_adv,
               n_study_sev = nrow(cohort_data_tb),
               frail2 = cohort_data_tb$study_id_num,
               study_sev_mod = cohort_data_tb$sev_mod,
               study_sev_adv = cohort_data_tb$sev_adv
)

fit_sev_tb <- jags(data = dt_sev_tb, model.file = m_sev,
                parameters.to.save = par_sev,
                n.iter = 11000, n.burnin = 1000,
                n.chains = 1, n.thin = 20)

mcmc_sev_tb <- as.mcmc(fit_sev_tb)
eval_sev_tb <- mcmc_sev_tb[, c("alpha", "bmod", "badv", "theta", "sdlog")]
#summary(eval_sev_tb)
#xyplot(eval_sev_tb)
#autocorr.plot(eval_sev_tb)
#densityplot(eval_sev_tb)

png("../Figures/xyplot_sev_tb.png")
xyplot(eval_sev_tb)
dev.off()
png("../Figures/autocorr_sev_tb.png")
autocorr.plot(eval_sev_tb)
dev.off()

res_sev_tb <- as.data.frame(summary(mcmc_sev_tb)$quantiles)


######## Sanatorium only model ########

san_tb <- metaData_tb %>%
  filter(sanatorium == 1) %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)))

dt_san_tb <- list(N = nrow(san_tb),
               interval = san_tb$interval,
               lim = cbind(san_tb$x1, san_tb$x2),
               time = rep(NA, nrow(san_tb)),
               n_frail = length(unique(san_tb$study_sev_num)),
               frail = san_tb$study_sev_num
               
)

fit_san_tb <- jags(data = dt_san_tb, model.file = m_all,
                parameters.to.save = par_all,
                n.iter = 11000, n.burnin = 1000,
                n.chains = 1, n.thin = 20)

mcmc_san_tb <- as.mcmc(fit_san_tb)
eval_san_tb <- mcmc_san_tb[, c("mu", "theta", "sdlog")]
#summary(eval_san_tb)
#xyplot(eval_san_tb)
#autocorr.plot(eval_san_tb)

png("../Figures/xyplot_san_tb.png")
xyplot(eval_san_tb)
dev.off()
png("../Figures/autocorr_san_tb.png")
autocorr.plot(eval_san_tb)
dev.off()

res_san_tb <- as.data.frame(summary(mcmc_san_tb)$quantiles)



######## Non-sanatorium only model ########

nosan_tb <- metaData_tb %>%
  filter(sanatorium == 0) %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)))

dt_nosan <- list(N = nrow(nosan_tb),
                 interval = nosan_tb$interval,
                 lim = cbind(nosan_tb$x1, nosan_tb$x2),
                 time = rep(NA, nrow(nosan_tb)),
                 n_frail = length(unique(nosan_tb$study_sev_num)),
                 frail = nosan_tb$study_sev_num
                 
)

fit_nosan_tb <- jags(data = dt_nosan, model.file = m_all,
                  parameters.to.save = par_all,
                  n.iter = 5000, n.burnin = 500,
                  n.chains = 1, n.thin = 10)

mcmc_nosan_tb <- as.mcmc(fit_nosan_tb)
eval_nosan_tb <- mcmc_nosan_tb[, c("mu", "theta", "sdlog")]
#summary(eval_nosan_tb)
#xyplot(eval_nosan_tb)
#autocorr.plot(eval_nosan_tb)

png("../Figures/xyplot_nosan_tb.png")
xyplot(eval_nosan_tb)
dev.off()
png("../Figures/autocorr_nosan_tb.png")
autocorr.plot(eval_nosan_tb)
dev.off()

res_nosan_tb <- as.data.frame(summary(mcmc_nosan_tb)$quantiles)


############################### Compiling Results ###############################


#### Function to get information about the dataset for each run ####

getData <- function(data){

  #Finding concordance to frailty IDs
  tab <- data %>%
    select(study_id, study_id_num, study_sev, study_sev_num) %>%
    filter(!duplicated(study_sev))
  
  #Finding n studies, severity groups, cohorts, individuals
  counts <- c("nStudies" = length(unique(data$study_id_num)),
                 "nSeverity" = length(unique(data$study_sev)),
                 "nCohorts" = length(unique(data$cohort_id)),
                 "nIndividuals" = nrow(data))
  
  return(list(tab, counts))
}

data_all <- getData(metaData)
data_sev <- getData(metaData_sev)
data_san <- getData(san)
data_nosan <- getData(nosan)

data_all_tb <- getData(metaData_tb)
data_sev_tb <- getData(metaData_sev_tb)
data_san_tb <- getData(san_tb)
data_nosan_tb <- getData(nosan_tb)


save(res_all, res_sev, res_san, res_nosan, data_all, data_sev, data_san, data_nosan,
     res_all_tb, res_sev_tb, res_san_tb, res_nosan_tb, data_all_tb, data_sev_tb, data_san_tb, data_nosan_tb,
     file = "bayesian_mortality.RData")

