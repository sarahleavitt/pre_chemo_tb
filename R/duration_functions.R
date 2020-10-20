
#Loading packages
reload_source <- function(){
  if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
  if (!require('tidyr')) install.packages('tidyr'); library('tidyr')
  if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
  if (!require('purrr')) install.packages('purrr'); library('purrr')
  if (!require('naniar')) install.packages('naniar'); library('naniar')
  if (!require('markdown')) install.packages('markdown'); library('markdown')
  if (!require('readxl')) install.packages('readxl'); library('readxl')
  if (!require('survminer')) install.packages('survminer'); library('survminer')
  if (!require('survival')) install.packages('survival'); library('survival')
  if (!require('flexsurv')) install.packages('flexsurv'); library('flexsurv')
  if (!require('parfm')) install.packages('parfm'); library('parfm')
  if (!require('lme4')) install.packages('lme4'); library('lme4')
  if (!require('gridExtra')) install.packages('gridExtra'); library('gridExtra')
  if (!require('knitr')) install.packages('knitr'); library('knitr')
  if (!require('stringr')) install.packages('stringr'); library('stringr')
}


## Function to read all sheets from an excel file
read_excel_allsheets <- function(filename, tibble = TRUE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}


## Function to turn one sheet of study data into individual-level data
studyToInd <- function(DF, outcome = c("mortality", "cure")){
  
  DF <- DF %>% filter(!is.na(cohort_id))
  DF <- DF[, !grepl("\\.\\.", names(DF))]
  
  DF <- DF %>% mutate(n = round(n, 0),
                      c1a = round(c1a, 0),
                      c1b = round(c1b, 0),
                      c2 = round(c2, 0),
                      l = round(l, 0))
  
  #Finding all-cause mortality
  if(is.na(DF$c1a_plus_b)[1]){
    DF <- DF %>% mutate(c1a_plus_b = c1a + c1b) 
  }
  
  if(outcome == "mortality"){
    cohorts <- DF %>%
      group_by(cohort_id) %>%
      group_modify(~ cohortToInd(.x)) %>%
      ungroup()
  }

  if(outcome == "cure"){
    cohorts <- DF %>%
      group_by(cohort_id) %>%
      group_modify(~ cohortToIndCure(.x)) %>%
      ungroup()
  }
  
  cohorts <- as.data.frame(cohorts)
}


## Function to turn one cohort into individual cure data (run by studyToInd)
cohortToIndCure <- function(DFc, timepoints = c(3, 5, 10)){
  
  nTotal <- as.numeric(DFc[1, "n"])
  
  ind <- NULL
  for(t in timepoints){
    nCure <- as.numeric(DFc[DFc$interval_r == t, "c2"])
    nLost <- sum(DFc[DFc$interval_r <= t, "l"])
    n <- nTotal - nLost
    
    cure <- cbind.data.frame("time" = t, "cure" = rep(1, nCure))
    censor <- cbind.data.frame("time" = t, "cure" = rep(0, n - nCure))
    ind <- bind_rows(ind, cure, censor)
  }
  
  ind <- ind %>% mutate(start_type = DFc$start_type[1],
                        sanatorium = DFc$sanatorium[1],
                        severity = DFc$stratified_var[1],
                        paper_id = as.character(DFc$paper_id[1]),
                        study_id = ifelse(grepl("1029", paper_id), "1029", paper_id))
}



## Function to turn one cohort into individual mortality data (run by studyToInd)
cohortToInd <- function(DFc){
  
  #Creating individual data where event = 1 means that the subject died and event = 0 means they were censored
  ind <- NULL
  for(i in 1:nrow(DFc)){
    row <- DFc[i, ]
    interval_l <- row$interval_l
    interval_r <- row$interval_r
    
    #If there is TB verses all-cause mortality, split the two, if not combine
    if(!is.na(DFc$c1a[1])){
      
      death_tb <- cbind.data.frame(rep(interval_l, row$c1a),
                                   rep(interval_r, row$c1a),
                                   rep(1, row$c1a))
      names(death_tb) <- c("interval_l", "interval_r", "death_tb")
      death_all <- cbind.data.frame(rep(interval_l, row$c1b),
                                    rep(interval_r, row$c1b),
                                    rep(1, row$c1b))
      names(death_all) <- c("interval_l", "interval_r", "death_all")
      
      event <- death_tb %>%
        bind_rows(death_all) %>%
        replace_na(list(death_tb = 0, death_all = 0)) %>%
        mutate(event = death_all + death_tb)
    }else{
      event <- cbind.data.frame(rep(interval_l, row$c1a_plus_b),
                                rep(interval_r, row$c1a_plus_b),
                                rep(1, row$c1a_plus_b))
      names(event) <- c("interval_l", "interval_r", "event")
    }
    
    #Changing interval_l values of 0 to be left censored at interval_r and setting interval_r to missing
    event <- event %>%
      mutate(interval_l_new = ifelse(interval_l == 0, interval_r, interval_l),
             interval_r = ifelse(interval_l == 0, NA, interval_r)) %>%
      select(-interval_l) %>%
      rename(interval_l = interval_l_new)
    
    #Finding those who were LTFU
    censor <- cbind.data.frame(rep(interval_l, row$l),
                               rep(interval_r, row$l),
                               rep(0, row$l),
                               rep(0, row$l),
                               rep(0, row$l))
    names(censor) <- c("interval_l", "interval_r", "death_tb", "death_all", "event")
    
    indTemp <- bind_rows(event, censor)
    ind <- bind_rows(ind, indTemp)
  }
  
  #Adding individuals who are censored at the end of the follow-up
  total <- DFc$n[1]
  haveOutcome <- nrow(ind)
  endTime <- DFc$interval_r[nrow(DFc)]
  censor <- cbind.data.frame(rep(endTime, total - haveOutcome), 
                             rep(0, total - haveOutcome),
                             rep(0, total - haveOutcome),
                             rep(0, total - haveOutcome))
  names(censor) <- c("interval_l", "death_tb", "death_all", "event")
  
  ind <- bind_rows(ind, censor)
  
  #Adding interval censor info
  #For interval censoring 0=right censored, 1=event at time, 2=left censored, 3=interval censored
  ind2 <- ind %>%
    mutate(eventIC = ifelse(event == 0, 0,
                            ifelse(is.na(interval_r), 2, 3)),
           death_tbIC = ifelse(death_tb == 0, 0,
                               ifelse(is.na(interval_r), 2, 3)),
           start_type = DFc$start_type[1],
           sanatorium = DFc$sanatorium[1],
           stratified_var = DFc$stratified_var[1],
           paper_id = as.character(DFc$paper_id[1]))
  
  return(ind2)
}




## Function to format results of Bayesian analysis
formatBayesian <- function(res, data, label, fixed = FALSE){
  
  #Parameter values
  if(fixed == FALSE){
    
    param <- res[row.names(res) %in% c("mu", "sdlog", "theta", "med_all",
                                       "pred_all[1]", "pred_all[5]", "pred_all[10]"), ]
    names(param) <- c("cilb", "lowerquant", "est", "upperquant", "ciub")
    param <- param %>% mutate(value = row.names(param),
                              value = ifelse(value == "mu", "meanlog",
                                      ifelse(value == "med_all", "median", gsub("_all\\[|\\]", "", value))),
                              label = label)
    
    #Overall survival and density curves
    sdlog <- param %>% filter(value == "sdlog") %>% pull(est)
    meanlog <- param %>% filter(value == "meanlog") %>% pull(est)
    x <- seq(0, 30, 0.01)
    dens <- dlnorm(x, meanlog, sdlog)
    surv <- plnorm(x, meanlog, sdlog, lower.tail = FALSE)
    
    #Credible bounds for survival curves
    credint <- res[grepl("pred_all", row.names(res)), ]
    credint <- credint %>%
      mutate(x = as.numeric(str_extract(row.names(.), "[0-9]+"))) %>%
      select(x, surv_est = `50%`, cilb = `2.5%`, ciub = `97.5%`) %>%
      bind_rows(c(x = 0, surv_est = 1, cilb = 1, ciub = 1))
    
    surv_dens <- cbind.data.frame(x, dens, surv, "label" = label) %>%
      full_join(credint, by = "x")
    
  } else{
    
    param <- res[row.names(res) %in% c("meanlog_min", "meanlog_mod", "meanlog_adv",
                                       "sdlog", "theta", "med_min", "med_mod", "med_adv",
                                       "pred_min[1]", "pred_min[5]", "pred_min[10]",
                                       "pred_mod[1]", "pred_mod[5]", "pred_mod[10]",
                                       "pred_adv[1]", "pred_adv[5]", "pred_adv[10]"), ]
    names(param) <- c("cilb", "lowerquant", "est", "upperquant", "ciub")
    param <- param %>% mutate(label = label,
                              severity = ifelse(grepl("min", row.names(.)), "Minimal",
                                         ifelse(grepl("mod", row.names(.)), "Moderate",
                                         ifelse(grepl("adv", row.names(.)), "Advanced", NA))),
                              value = gsub("_[a-z]*$|_[a-z]*\\[|\\]", "", row.names(.)),
                              value = ifelse(value == "med", "median", value))
    
    #Overall survival and density curves
    sdlog <- param %>% filter(value == "sdlog") %>% pull(est)
    surv_dens <- NULL
    x <- seq(0, 30, 0.1)
    for(sev in c("Minimal", "Moderate", "Advanced")){
      meanlog <- param %>% filter(severity == sev, value == "meanlog") %>% pull(est)
      
      dens <- dlnorm(x, meanlog, sdlog)
      surv <- plnorm(x, meanlog, sdlog, lower.tail = FALSE)
      densTemp <- cbind.data.frame(x, "severity" = sev, meanlog, sdlog, dens, surv)
      
      surv_dens <- bind_rows(surv_dens, densTemp)
    }
    
    #Credible bounds for survival curves
    credint <- res[grepl("pred_[a-z]*", row.names(res)), ]
    credint <- credint %>%
      mutate(x = as.numeric(str_extract(row.names(.), "[0-9]+")),
             severity = ifelse(grepl("min", row.names(.)), "Minimal",
                               ifelse(grepl("mod", row.names(.)), "Moderate",
                                      ifelse(grepl("adv", row.names(.)), "Advanced", NA)))) %>%
      select(x, severity, surv_est = `50%`, cilb = `2.5%`, ciub = `97.5%`) %>%
      bind_rows(cbind.data.frame(x = 0, severity = "Minimal", surv_est = 1, cilb = 1, ciub = 1),
                cbind.data.frame(x = 0, severity = "Moderate", surv_est = 1, cilb = 1, ciub = 1),
                cbind.data.frame(x = 0, severity = "Advanced", surv_est = 1, cilb = 1, ciub = 1))
    
    surv_dens <- surv_dens %>% 
      mutate(label = label) %>%
      full_join(credint, by = c("x", "severity")) %>%
      mutate(severity = factor(severity, levels = c("Minimal", "Moderate", "Advanced", "Unknown")))
  }
  
  
  
  #Individual study median
  med_ind <- as.data.frame(res[grepl("med_ind", row.names(res)),])
  med_ind <- med_ind %>%
    mutate(study_sev_num = as.numeric(gsub("med_ind\\[|\\]", "", row.names(med_ind))),
           value = "median")
  
  #Individual study predictions
  pred_ind <- as.data.frame(res[grepl("pred[0-9]*\\[", row.names(res)),])
  pred_ind <- pred_ind %>%
    mutate(rown = row.names(.),
           study_sev_num = as.numeric(gsub("pred[0-9]*\\[|\\]", "", rown)),
           value = str_extract(rown, "pred[0-9]*")) %>%
    select(-rown)
  
  #Individual study meanlog
  if(fixed == FALSE){
    mean_ind <- as.data.frame(res[grepl("meanlog", row.names(res)),])
    mean_ind <- mean_ind %>%
      mutate(study_sev_num = as.numeric(gsub("meanlog\\[|\\]", "", row.names(mean_ind))),
             value = "meanlog")
  }else{
    mean_ind <- as.data.frame(res[grepl("meanlog_ind", row.names(res)),])
    mean_ind <- mean_ind %>%
      mutate(study_sev_num = as.numeric(gsub("meanlog_ind\\[|\\]", "", row.names(mean_ind))),
             value = "meanlog")
  }
  
  #Combining the above
  ind_est <- bind_rows(med_ind, mean_ind, pred_ind)
  names(ind_est) <- c("cilb", "lowerquant", "est", "upperquant", "ciub", "study_sev_num", "value")
  ind_est <- ind_est %>%
    full_join(data[[1]], by = "study_sev_num") %>%
    left_join(covar, by = "study_sev") %>%
    mutate(label = label) %>%
    select(-study_sev_num, -study_id_num)
  
  #Finding density and survival estimates for each study
  sdlog <- param %>% filter(value == "sdlog") %>% pull(est)
  par_data <- ind_est %>% filter(value == "meanlog")
  ind_surv <- NULL
  x <- seq(0, 30, 0.1)
  for(i in 1:nrow(par_data)){
    row <- par_data[i,]
    
    dens <- dlnorm(x, row$est, sdlog)
    surv <- plnorm(x, row$est, sdlog, lower.tail = FALSE)
    densTemp <- cbind.data.frame(x, "study_sev" = row$study_sev, "meanlog" = row$est, sdlog, dens, surv)
    
    ind_surv <- bind_rows(ind_surv, densTemp)
  }
  ind_surv <- ind_surv %>%
    left_join(covar, by = "study_sev") %>%
    mutate(label = label,
           severity = factor(severity, levels = c("Minimal", "Moderate", "Advanced", "Unknown")))
  
  #Combining study-specific and overall median and predictions
  pred_comb <- ind_est %>% 
    filter(value == "median" | grepl("pred", value)) %>%
    bind_rows(param %>% filter(value == "med" | grepl("pred", value))) %>%
    mutate(shape = ifelse(is.na(study_sev), "Overall", "Individual"))

  
  if(fixed == FALSE){
    pred_comb <- pred_comb %>%
      replace_na(list(severity = "",
                      study_sev = "Overall",
                      first_author = "Overall")) %>%
      mutate(severity = factor(severity, levels = c("Minimal", "Moderate", "Advanced",
                                                    "Unknown", "")))
  }else{
    pred_comb <- pred_comb %>%
      mutate(study_sev = ifelse(is.na(study_sev), paste("Overall", severity, sep = "_"),
                                study_sev),
             first_author = ifelse(is.na(first_author), "Overall", first_author),
             severity = ifelse(grepl("Overall", study_sev) & severity == "Minimal", "",
                               ifelse(grepl("Overall", study_sev) & severity == "Moderate", " ",
                                      ifelse(grepl("Overall", study_sev) & severity == "Advanced", "  ",
                                             severity))),
             severity = factor(severity, levels = c("Minimal", "", "Moderate", " ", "Advanced", "  ")))
    
  }
  
  pred_comb <- pred_comb %>% 
    mutate(pred_label = factor(value, levels = c("pred1", "pred5", "pred10", "median"),
                               labels = c("1-Year", "5-Year", "10-Year", "Median")))
  
  return(list("surv_dens" = surv_dens, "param" = param,
              "pred_comb" = pred_comb, "ind_surv" = ind_surv))
}



## Log likelihood for interval censoring
logl <- function(params, DF, dist, interval_l = "interval_l", interval_r = "interval_r",
                 event = "eventIC"){
  
  #0=right censored, 1=event at interval_r, 2=left censored, 3=interval censored
  if(length(params) == 1){
    DF$f <- ifelse(DF[, event] == 1,
                   distL[[dist]][[1]](DF[, interval_r], params),
                   ifelse(DF[, event] == 2,
                          distL[[dist]][[2]](DF[, interval_r], params),
                          ifelse(DF[, event] == 3,
                                 distL[[dist]][[2]](DF[, interval_r], params) -  
                                   distL[[dist]][[2]](DF[, interval_l], params),
                                 distL[[dist]][[2]](DF[, interval_l], params, lower.tail = FALSE))))
  }else if(length(params == 2)){
    DF$f <- ifelse(DF[, event] == 1,
                   distL[[dist]][[1]](DF[, interval_r], params[1], params[2]),
                   ifelse(DF[, event] == 2,
                          distL[[dist]][[2]](DF[, interval_l], params[1], params[2]),
                          ifelse(DF[, event] == 3,
                                 distL[[dist]][[2]](DF[, interval_r], params[1], params[2]) -
                                   distL[[dist]][[2]](DF[, interval_l], params[1], params[2]),
                                 distL[[dist]][[2]](DF[, interval_l], params[1], params[2], lower.tail = FALSE))))
  }
  
  
  DF$logf <- ifelse(DF$f == 0, 0, log(DF$f))
  -sum(DF$logf, na.rm = TRUE)
}


## Function to perform interval censoring mortality analysis
survIC <- function(DF, dist, ipars, cl = 0.95){
  
  # ## Manual approach ##
  # if(dist == "weibull"){
  #   opt2 <- optim(par = ipars, logl, DF = DF, dist = dist, method = "L-BFGS-B",
  #                lower = c(0.005, 0.005), hessian = TRUE)
  # }else{
  #   opt <- optim(par = ipars, logl, DF = DF, dist = dist, method = "L-BFGS-B",
  #                lower = c(0.00001, 0.00001), hessian = TRUE)
  # }
  
  # est <- opt$par
  # cov <- solve(opt$hessian)
  # se <- sqrt(diag(cov))
  # lcl <- est - qnorm(1 - (1 - cl)/2) * se
  # ucl <- est + qnorm(1 - (1 - cl)/2) * se
  
  
  # ## survival Package ##
  opt <- survreg(Surv(event = eventIC, time = interval_l, time2 = interval_r,
                          type = "interval") ~ 1, data = DF, dist = dist)
  
  # Converting results
  if(dist == "exp"){
    est <- c("rate" = as.numeric(exp(-opt$coefficients)))
  }else if(dist == "weibull"){
    est <- c("shape" = as.numeric(1/opt$scale), scale = as.numeric(exp(opt$coefficients)))
  }else if(dist == "lognormal"){
    est <- c("meanlog" = as.numeric(opt$coefficients), "sdlog" = as.numeric(opt$scale))
  }
  
  x <- seq(0, 40, 0.1)
  if(dist == "exp"){
    median <- distL[[dist]][[3]](0.5, est)
    dens <- distL[[dist]][[1]](x, est)
    surv <- distL[[dist]][[2]](x, est, lower.tail = FALSE)
    DF <- 1
  }else{
    median <- distL[[dist]][[3]](0.5, est[1], est[2])
    dens <- distL[[dist]][[1]](x, est[1], est[2])
    surv <- distL[[dist]][[2]](x, est[1], est[2], lower.tail = FALSE)
    DF <- 2
  }
  
  #Returning list of results
  res <- list(pars = est,
              # cilb = lcl,
              # ciub = ucl,
              median = median,
              loglik = opt$loglik[1],
              DF = DF,
              aic = -2*opt$loglik[1] + 2*DF,
              x = x,
              dens = dens,
              surv = surv)
  
  return(res) 
}




## Function to return plots and results from mortality survival analysis
survAnalysis <- function(DF, dists = c("exp", "weibull", "lognormal")){
  
  parDF <- NULL
  resDF <- NULL
  aicDF <- NULL
  medDF <- NULL
  
  for(d in dists){
    
    # if(d == "exp"){
    #   m <- survIC(DF, d, 1)
    # }else{
    #   m <- survIC(DF, d, c(1,1))
    # }
    
    m <- survIC(DF, d)
    
    DF1 <- cbind.data.frame("x" = m$x, "surv" = m$surv, "dens" = m$dens,
                            "dist" = d)
    DF2 <- cbind.data.frame("dist" = d, "aic" = m$aic)
    DF3 <- cbind.data.frame("dist" = d, "est" = m$par)
                            #"cilb" = m$cilb, "ciub" = m$ciub)
    DF4 <- cbind.data.frame("dist" = d, "median" = m$median)
    
    resDF <- bind_rows(resDF, DF1)
    aicDF <- bind_rows(aicDF, DF2)
    parDF <- bind_rows(parDF, DF3)
    medDF <- bind_rows(medDF, DF4)
  }
  
  return(list(data = DF, aic = aicDF, par = parDF, median = medDF, density = resDF))
}

plotSurv <- function(survRes){

  n_dists <- length(unique(survRes$density$dist)) 
  
  p1 <- ggsurvplot(
    fit = survfit(Surv(interval_l, event) ~ 1, data = survRes$data),
    data = DF,
    xlab = "Years",
    ylab = "Overall survival probability",
    legend = "none")
  
  p2 <- ggplot(survRes$density, aes(x = x, y = surv, color = dist, linetype = dist)) +
    geom_line(size = 1) +
    #geom_smooth(aes(ymin = survl, ymax = survu, fill = dist, color = dist),
    #            stat = "identity") +
    xlab("Years") + ylab("Survival, 1 - F(t)") +
    scale_y_continuous(limits = c(0, 1)) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title = element_blank()) +
    guides(color = guide_legend(nrow=2, byrow=TRUE),
           linetype = guide_legend(nrow=2, byrow=TRUE))
  
  p3 <- ggplot(survRes$density, aes(x = x, y = dens, color = dist, linetype = dist)) +
    geom_line(size = 1) +
    #geom_smooth(aes(ymin = densu, ymax = densl, fill = dist, color = dist),
    #            stat = "identity") +
    xlab("Years") + ylab("Density, f(t)") +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title = element_blank()) +
    guides(color = guide_legend(nrow=2, byrow=TRUE),
           linetype = guide_legend(nrow=2, byrow=TRUE))
  

  grid.arrange(p1$plot, p2, p3, nrow = 1)
}


metaIndividual <- function(data, meta, fixedEffect = FALSE){
  
  data2 <- data[!duplicated(data$study_sev), c("study_id", "study_sev", "severity", "sanatorium")]
  
  if(fixedEffect == FALSE){
    
    coeff <- meta$coefficients[1]
    
    if(nrow(data2) <= 5){
      frail <- as.numeric(meta$coefficients[-1])
    }else{
      frail <- meta$frail
    }
    
    frailEst <- cbind.data.frame("study_sev" = data2$study_sev,
                                 "frail" = frail,
                                 "fixed" = as.numeric(coeff),
                                 "est" = coeff + frail)
  }
  
  if(fixedEffect == TRUE){
    
    coeff <- cbind.data.frame("severity" = c("Minimal", "Moderate", "Advanced", "Unknown"),
                              "fixed" = meta$coefficients)
    coeff <- coeff %>% mutate(fixed = ifelse(severity != "Minimal", fixed + meta$coefficients[1], fixed))
    
    frailEst <- cbind.data.frame("study_id" = unique(data2$study_id),
                                 "frail" = meta$frail)
    
    frailEst <- frailEst %>%
      full_join(data2, by = "study_id") %>%
      full_join(coeff, by = "severity") %>%
      mutate(est = frail + fixed)
  }
  
  dist <- ifelse(meta$dist == "exponential", "exp", meta$dist)
  
  # Finding individual parameters
  if(dist == "exp"){
    frailEst$par1 <- exp(-frailEst$est) #rate
  }else if(dist == "weibull"){
    frailEst$par1 <- 1 / meta$scale #shape
    frailEst$par2 <- exp(frailEst$est)
  }else if(dist == "lognormal"){
    frailEst$par1 <- frailEst$est #meanlog
    frailEst$par2 <- meta$scale #sdlog
  }
  
  #Finding density and survival estimates
  densData <- NULL
  x <- seq(0, 40, 0.1)
  for(i in 1:nrow(frailEst)){
    row <- frailEst[i,]
    
    if(dist == "exp"){
      dens <- distL[[dist]][[1]](x, row$par1)
      surv <- distL[[dist]][[2]](x, row$par1, lower.tail = FALSE)
      densTemp <- cbind.data.frame(dist, x, row$study_sev, row$par1, NA, dens, surv)
    }else{
      dens <- distL[[dist]][[1]](x, row$par1, row$par2)
      surv <- distL[[dist]][[2]](x, row$par1, row$par2, lower.tail = FALSE)
      densTemp <- cbind.data.frame(dist, x, row$study_sev, row$par1, row$par2, dens, surv)
    }
    densData <- bind_rows(densData, densTemp)
  }
  
  names(densData) <- c("dist", "x", "study_sev", "par1", "par2", "dens", "surv")
  densData2 <- densData %>% full_join(data2, by = "study_sev")
  
  return(densData2)
}


performMeta <- function(data, model, fixedEffect = FALSE,
                        dists = c("exp", "weibull", "lognormal")){

  form <- as.formula(paste0("Surv(event = eventIC, time = interval_l, time2 = interval_r, type = 'interval') ~ ",
         model))
  
  #form <- as.formula(paste0("Surv(event = event, time = interval_l) ~ ", model))
  
  indRes <- NULL
  aic <- NULL
  for(d in dists){
    meta <- survreg(form, data = data, dist = d)
    
    indRes_temp <- metaIndividual(data, meta, fixedEffect)
    aic_temp <- bind_cols("dist" = d, "aic" = -2*meta$loglik[2] + 2*meta$df[2])
    
    indRes <- bind_rows(indRes, indRes_temp)
    aic <- bind_rows(aic, aic_temp)
  }
  
  return(list("indRes" = indRes, "aic" = aic))
}


