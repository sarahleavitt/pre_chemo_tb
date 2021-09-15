#Sarah V. Leavitt
#Boston University
#Pre-chemotherapy TB Analysis

##############################################################################
# This program performs the mortality Bayesian meta-analysis and saves the
# results as an R object. This takes around 30 minutes to run.
##############################################################################

options(scipen=999)
options(digits = 10)
set.seed(150183)

rm(list = ls())
source("R/utils.R")
reload_source()

#Reading in individual mortality data
mortalityData <- read.csv("data/mortality_data.csv")


#### Models ----------------------------------------------------------------------------------------

#### Complete model (no fixed effect for severity) ####

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

#Parameters to track
par_all <- c("mu", "theta", "sdlog", "med_all", "meanlog", "med_ind",
             "pred_all", "pred1", "pred5", "pred10")



#### Stratified model (fixed effect for severity) ####

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

#Parameters to track
par_sev <- c("theta", "sdlog", "alpha", "bmod", "badv",
             "meanlog_ind", "meanlog_min", "meanlog_mod", "meanlog_adv",
             "med_ind", "med_min", "med_mod", "med_adv",
             "pred1", "pred5", "pred10", "pred_min", "pred_mod", "pred_adv")





#### TB Mortality Analysis--------------------------------------------------------------------------

#interval = 1 implies x1 < t <= x2
#dinterval() isn't working as it should so make all obs have interval = 1
#and set x2 to be 10000 (close enough to infinity) for right censored
mortalityData_tb <- mortalityData %>% 
  mutate(study_sev_num = as.numeric(factor(study_sev)),
         study_id_num = as.numeric(factor(study_id)),
         time = ifelse(death_tb == 0, NA, interval_l),
         interval = 1,
         x1 = interval_l,
         x2 = ifelse(death_tb == 0, 10000, interval_r),
         sev_mod = as.numeric(severity == "Moderate"),
         sev_adv = as.numeric(severity == "Advanced"),
         sev_unk = as.numeric(severity == "Unknown"))



#### TB mortality: complete model ####

#Data
dt_all_tb <- list(N = nrow(mortalityData_tb),
               interval = mortalityData_tb$interval,
               lim = cbind(mortalityData_tb$x1, mortalityData_tb$x2),
               time = rep(NA, nrow(mortalityData_tb)),
               n_frail = length(unique(mortalityData_tb$study_sev_num)),
               frail = mortalityData_tb$study_sev_num
)

#Fitting the model
fit_all_tb <- jags(data = dt_all_tb, model.file = m_all,
                parameters.to.save = par_all,
                n.iter = 11000, n.burnin = 1000,
                n.chains = 1, n.thin = 20)

#Extracting results
mcmc_all_tb <- as.mcmc(fit_all_tb)
eval_all_tb <- mcmc_all_tb[, c("mu", "theta", "sdlog")]
res_all_tb <- as.data.frame(summary(mcmc_all_tb)$quantiles)

png("Figures/xyplot_all_tb.png")
xyplot(eval_all_tb)
dev.off()
png("Figures/autocorr_all_tb.png")
autocorr.plot(eval_all_tb)
dev.off()



#### TB mortality: stratified model ####

#Subsetting and formatting data
mortalityData_sev_tb <- mortalityData_tb %>%
  filter(severity != "Unknown") %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)),
         study_id_num = as.numeric(factor(study_id)))

cohort_data_tb <- mortalityData_sev_tb %>%
  group_by(study_sev_num) %>%
  summarize(study_id_num = first(study_id_num),
            sev_mod = first(sev_mod),
            sev_adv = first(sev_adv),
            sev_unk = first(sev_unk))

#Data
dt_sev_tb <- list(N = nrow(mortalityData_sev_tb),
               interval = mortalityData_sev_tb$interval,
               lim = cbind(mortalityData_sev_tb$x1, mortalityData_sev_tb$x2),
               time = rep(NA, nrow(mortalityData_sev_tb)),
               n_frail = length(unique(mortalityData_sev_tb$study_id_num)),
               frail = mortalityData_sev_tb$study_id_num,
               sev_mod = mortalityData_sev_tb$sev_mod,
               sev_adv = mortalityData_sev_tb$sev_adv,
               n_study_sev = nrow(cohort_data_tb),
               frail2 = cohort_data_tb$study_id_num,
               study_sev_mod = cohort_data_tb$sev_mod,
               study_sev_adv = cohort_data_tb$sev_adv
)

#Fitting model
fit_sev_tb <- jags(data = dt_sev_tb, model.file = m_sev,
                parameters.to.save = par_sev,
                n.iter = 11000, n.burnin = 1000,
                n.chains = 1, n.thin = 20)

#Extracting results
mcmc_sev_tb <- as.mcmc(fit_sev_tb)
eval_sev_tb <- mcmc_sev_tb[, c("alpha", "bmod", "badv", "theta", "sdlog")]
res_sev_tb <- as.data.frame(summary(mcmc_sev_tb)$quantiles)

png("Figures/xyplot_sev_tb.png")
xyplot(eval_sev_tb)
dev.off()
png("Figures/autocorr_sev_tb.png")
autocorr.plot(eval_sev_tb)
dev.off()





#### Sanatorium Sensitivity Analysis----------------------------------------------------------------


#### TB mortality: sanatorium/hospital ####

#Subsetting data
san_tb <- mortalityData_tb %>%
  filter(sanatorium == "Yes") %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)))

#Data
dt_san_tb <- list(N = nrow(san_tb),
                  interval = san_tb$interval,
                  lim = cbind(san_tb$x1, san_tb$x2),
                  time = rep(NA, nrow(san_tb)),
                  n_frail = length(unique(san_tb$study_sev_num)),
                  frail = san_tb$study_sev_num
                  
)

#Fitting the model
fit_san_tb <- jags(data = dt_san_tb, model.file = m_all,
                   parameters.to.save = par_all,
                   n.iter = 11000, n.burnin = 1000,
                   n.chains = 1, n.thin = 20)

#Extracting results
mcmc_san_tb <- as.mcmc(fit_san_tb)
eval_san_tb <- mcmc_san_tb[, c("mu", "theta", "sdlog")]
res_san_tb <- as.data.frame(summary(mcmc_san_tb)$quantiles)

png("Figures/xyplot_san_tb.png")
xyplot(eval_san_tb)
dev.off()
png("Figures/autocorr_san_tb.png")
autocorr.plot(eval_san_tb)
dev.off()



#### TB mortality: non-sanatorium ####

#Subsetting data
nosan_tb <- mortalityData_tb %>%
  filter(sanatorium == "No") %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)))

#Data
dt_nosan <- list(N = nrow(nosan_tb),
                 interval = nosan_tb$interval,
                 lim = cbind(nosan_tb$x1, nosan_tb$x2),
                 time = rep(NA, nrow(nosan_tb)),
                 n_frail = length(unique(nosan_tb$study_sev_num)),
                 frail = nosan_tb$study_sev_num
                 
)

#Fitting model
fit_nosan_tb <- jags(data = dt_nosan, model.file = m_all,
                     parameters.to.save = par_all,
                     n.iter = 5000, n.burnin = 500,
                     n.chains = 1, n.thin = 10)

#Extracting results
mcmc_nosan_tb <- as.mcmc(fit_nosan_tb)
eval_nosan_tb <- mcmc_nosan_tb[, c("mu", "theta", "sdlog")]
res_nosan_tb <- as.data.frame(summary(mcmc_nosan_tb)$quantiles)

png("Figures/xyplot_nosan_tb.png")
xyplot(eval_nosan_tb)
dev.off()
png("Figures/autocorr_nosan_tb.png")
autocorr.plot(eval_nosan_tb)
dev.off()




#### Compiling Results------------------------------------------------------------------------------


## Function to get information about the dataset for each run
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

data_all_tb <- getData(mortalityData_tb)
data_sev_tb <- getData(mortalityData_sev_tb)
data_san_tb <- getData(san_tb)
data_nosan_tb <- getData(nosan_tb)


save(res_all_tb, res_sev_tb, res_san_tb, res_nosan_tb,
     data_all_tb, data_sev_tb, data_san_tb, data_nosan_tb,
     file = "R/bayesian_mortality.RData")

