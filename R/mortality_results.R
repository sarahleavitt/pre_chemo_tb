#Sarah V. Leavitt
#Boston University
#Pre-chemotherapy TB Analysis

##############################################################################
# This program creates figures and tables of the mortality analysis results
##############################################################################

options(scipen=999)
options(digits = 10)

rm(list = ls())
source("R/utils.R")
reload_source()

#Reading in the study_id correspondence table
studyid <- read.csv("data/study_id.csv")

#Reading in individual mortality data and analysis results
mortality <- read.csv("data/mortality_data.csv")
load('R/bayesian_mortality.RData')

#Formatting results
form_comp_all <- formatBayesian(mortality, res_comp_all, data_comp_all, "Combined_all")
form_sev_all <- formatBayesian(mortality, res_sev_all, data_sev_all, "Severity_all", fixed = TRUE)
form_comp_us <- formatBayesian(mortality, res_comp_us, data_comp_us, "Combined_us")
form_sev_us <- formatBayesian(mortality, res_sev_us, data_sev_us, "Severity_us", fixed = TRUE)
form_comp_post <- formatBayesian(mortality, res_comp_post, data_comp_post, "Combined_post")
form_sev_post <- formatBayesian(mortality, res_sev_post, data_sev_post, "Severity_post", fixed = TRUE)
form_san <- formatBayesian(mortality, res_san, data_san, "Sanatorium/hospital")
form_nosan <- formatBayesian(mortality, res_nosan, data_nosan, "Non-Sanatorium")



#### Table of Main Results -------------------------------------------------------------------------

#Combining raw tables from results lists
raw_tab <- bind_rows(form_comp_all$param, form_sev_all$param,
                     form_comp_us$param, form_sev_us$param,
                     form_comp_post$param, form_sev_post$param,
                     form_san$param, form_nosan$param)

#Adding formatted, ordered labels
raw_tab <- raw_tab %>%
  mutate(Severity = ifelse(is.na(severity), label,
                           ifelse(label == "Severity_all", paste0(severity, "_all"),
                                  ifelse(label == "Severity_us", paste0(severity, "_us"),
                                         ifelse(label == "Severity_post", paste0(severity, "_post"),
                                                label)))),
         Severity = factor(Severity, level = c("Minimal_all", "Moderate_all", "Advanced_all", "Combined_all",
                                               "Minimal_us", "Moderate_us", "Advanced_us", "Combined_us",
                                               "Minimal_post", "Moderate_post", "Advanced_post", "Combined_post",
                                               "Sanatorium/hospital", "Non-Sanatorium"))) %>%
  arrange(Severity)

#Extracting the survival probabilities
pred1_tab <- raw_tab %>%
  filter(value == "pred1") %>%
  mutate(`1-Year Survival (95% CI)` = paste0(round(est, 2), " (",
                                             round(cilb, 2), ", ",
                                             round(ciub, 2), ")")) %>%
  select(Severity, `1-Year Survival (95% CI)`)

pred5_tab <- raw_tab %>%
  filter(value == "pred5") %>%
  mutate(`5-Year Survival (95% CI)` = paste0(round(est, 2), " (",
                                             round(cilb, 2), ", ",
                                             round(ciub, 2), ")")) %>%
  select(Severity, `5-Year Survival (95% CI)`)

pred10_tab <- raw_tab %>%
  filter(value == "pred10") %>%
  mutate(`10-Year Survival (95% CI)` = paste0(round(est, 2), " (",
                                              round(cilb, 2), ", ",
                                              round(ciub, 2), ")")) %>%
  select(Severity, `10-Year Survival (95% CI)`)


#Extracting the distribution parameters
sdlog <- raw_tab %>%
  filter(value == "sdlog") %>%
  select(label, sdlog = est)

dist_tab <- raw_tab %>%
  filter(value == "meanlog") %>%
  full_join(sdlog, by = "label") %>%
  mutate(`Survival Distribution` = paste0("lognormal(", round(est, 2), ", ",
                                          round(sdlog, 2), ")")) %>%
  select(Severity, `Survival Distribution`)


#Finding number of papers, cohorts, individual for each analysis
data_comp_all2 <- as.data.frame(t(data_comp_all[[2]]))
data_comp_all2$Severity <- "Combined_all"
data_comp_us2 <- as.data.frame(t(data_comp_us[[2]]))
data_comp_us2$Severity <- "Combined_us"
data_comp_post2 <- as.data.frame(t(data_comp_post[[2]]))
data_comp_post2$Severity <- "Combined_post"
data_san2 <- as.data.frame(t(data_san[[2]]))
data_san2$Severity <- "Sanatorium/hospital"
data_nosan2 <- as.data.frame(t(data_nosan[[2]]))
data_nosan2$Severity <- "Non-Sanatorium"

counts_comp <- bind_rows(data_comp_all2, data_comp_us2, data_comp_post2, data_san2, data_nosan2)

#Counts stratified by severity for all studies
mortality_sev <- mortality %>% filter(severity != "Unknown")

counts_all <- mortality_sev %>%
  group_by(severity) %>%
  summarize(nStudies = length(unique(study_id)),
            nCohorts = length(unique(cohort_id)),
            nIndividuals = n(),
            .groups = "drop") %>%
  mutate(Severity = paste0(severity, "_all"))


#Counts stratified by severity for US studies
mortality_us <- mortality %>% filter(severity != "Unknown",
                                     study_id %in% c("1029", "93", "45", "63", "67", "90_1016"))

counts_us <- mortality_us %>%
  group_by(severity) %>%
  summarize(nStudies = length(unique(study_id)),
            nCohorts = length(unique(cohort_id)),
            nIndividuals = n(),
            .groups = "drop") %>%
  mutate(Severity = paste0(severity, "_us"))


#Counts stratified by severity for US post-1930s subset
mortality_post <- mortality %>% filter(severity != "Unknown", study_id %in% c("1029", "93", "45"))

counts_post <- mortality_post %>%
  group_by(severity) %>%
  summarize(nStudies = length(unique(study_id)),
            nCohorts = length(unique(cohort_id)),
            nIndividuals = n(),
            .groups = "drop") %>%
  mutate(Severity = paste0(severity, "_post"))


counts_initial <- bind_rows(counts_comp, counts_all, counts_us, counts_post) %>%
  select(Severity, `Number of Studies` = nStudies,
         `Number of Cohorts` = nCohorts,
         `Number of Individuals` = nIndividuals)



#Combining the tables
final_tab <- dist_tab %>%
  full_join(pred1_tab, by = c("Severity")) %>%
  full_join(pred5_tab, by = c("Severity")) %>%
  full_join(pred10_tab, by = c("Severity")) %>%
  full_join(counts_initial, by = c("Severity"))

#Separating three tables for paper
all_tab <- final_tab %>% filter(grepl("all", Severity))
us_tab <- final_tab %>% filter(grepl("us", Severity))
post_tab <- final_tab %>% filter(grepl("post", Severity))
san_tab <- final_tab %>% filter(grepl("San", Severity))

#Variance of frailty terms
theta <- raw_tab %>%
  filter(value == "theta") %>%
  mutate(theta = round(est, 2)) %>%
  select(label, theta)



#### Summary Survival Curves -----------------------------------------------------------------------

plot_sum_curves <- function(complete_model, stratified_model){
  
  #Survival curve for complete model
  p1 <- ggplot(complete_model$surv_dens) +
    geom_line(aes(x = x, y = surv),
              color = "black", size = 1, linetype = "solid") +
    geom_smooth(aes(x = x, y = surv_est, ymin = cilb, ymax = ciub),
                stat = "identity", linetype = 0, alpha = 0.25, na.rm = TRUE) +
    scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
    scale_x_continuous(name = "Years", limits = c(0, 30)) +
    theme_bw()
  
  
  #TB survival for fixed effect model
  p2 <- ggplot(stratified_model$surv_dens) +
    geom_line(aes(x = x, y = surv, color = severity),
              size = 1, linetype = "solid") +
    geom_smooth(aes(x = x, y = surv_est, ymin = cilb, ymax = ciub, fill = severity),
                stat = "identity", linetype = 0, alpha = 0.15, na.rm = TRUE) +
    scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
    scale_x_continuous(name = "Years", limits = c(0, 30)) +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_color_manual("", values = c("Minimal" = "seagreen", "Moderately advanced" = "goldenrod1",
                                      "Far advanced" = "firebrick2", "Unknown" = "grey50")) +
    scale_fill_manual("",
                      values = c("Minimal" = "seagreen", "Moderately advanced" = "goldenrod1",
                                 "Far advanced" = "firebrick2", "Unknown" = "grey50"))
  
  p_comb <- arrangeGrob(p1, p2, nrow = 2)
  return(p_comb)
}

#All studies
p_sum_all <- plot_sum_curves(form_comp_all, form_sev_all)
ggsave("Figures/summary_curves_all.png", p_sum_all, width = 5.5, height = 8)

#US studies
p_sum_us <- plot_sum_curves(form_comp_us, form_sev_us)
ggsave("Figures/summary_curves_us.png", p_sum_us, width = 5.5, height = 8)

#US post-1930s studies
p_sum_post <- plot_sum_curves(form_comp_post, form_sev_post)
ggsave("Figures/summary_curves_post.png", p_sum_post, width = 5.5, height = 8)



#### Individual Survival Curves --------------------------------------------------------------------

plot_ind_curves <- function(complete_model, stratified_model){
  
  #Survival curves for complete model
  p1 <- ggplot(complete_model$ind_surv) +
    geom_line(aes(x = x, y = surv, group = study_sev, color = severity),
              size = 0.7, alpha = 0.3) +
    geom_line(data = complete_model$surv_dens, aes(x = x, y = surv),
              color = "black", size = 1, linetype = "longdash") +
    scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
    scale_x_continuous(name = "Years", limits = c(0, 30)) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_color_manual("Disease Severity",
                       values = c("Minimal" = "seagreen", "Moderately advanced" = "goldenrod1",
                                  "Far advanced" = "firebrick2", "Unknown" = "grey50"))
  
  #Survival curves for stratified model
  p2 <- ggplot(stratified_model$ind_surv) +
    geom_line(aes(x = x, y = surv, group = study_sev, color = severity),
              size = 0.7, alpha = 0.3) +
    geom_line(data = stratified_model$surv_dens,
              aes(x = x, y = surv, color = severity),
              linetype = "longdash", size = 1) +
    scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
    scale_x_continuous(name = "Years", limits = c(0, 30)) +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_color_manual("", drop = FALSE,
                       values = c("Minimal" = "seagreen", "Moderately advanced" = "goldenrod1",
                                  "Far advanced" = "firebrick2", "Unknown" = "grey50"))
  
  p_comb <- arrangeGrob(p1, p2, nrow = 2)
  return(p_comb)
}

#All studies
p_ind_all <- plot_ind_curves(form_comp_all, form_sev_all)
ggsave("Figures/individual_curves_all.png", p_ind_all, width = 5.5, height = 8)

#US studies
p_ind_us <- plot_ind_curves(form_comp_us, form_sev_us)
ggsave("Figures/individual_curves_us.png", p_ind_us, width = 5.5, height = 8)

#US post-1930s studies
p_ind_post <- plot_ind_curves(form_comp_post, form_sev_post)
ggsave("Figures/individual_curves_post.png", p_ind_post, width = 5.5, height = 8)




#### All studies survival curves by category -------------------------------------------------------

ggplot(form_sev_all$ind_surv) +
  geom_line(aes(x = x, y = surv, group = study_sev, color = category),
            size = 0.7, alpha = 0.3) +
  geom_line(data = form_sev_all$surv_dens,
            aes(x = x, y = surv, group = severity),
            linetype = "longdash", size = 1, color = "grey50") +
  facet_wrap(~severity) +
  scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
  scale_x_continuous(name = "Years", limits = c(0, 30)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual("Type of Study",
                     values = c("Non-US" = "mediumturquoise", "US post-1930" = "mediumvioletred",
                                "US pre-1930" = "royalblue3")) +
  ggsave("Figures/survival_curves_category.png", width = 8, height = 4.5)



##### Forest plots --------------------------------------------------------------------------------

pred_plot_all <- form_comp_all$pred_comb %>%
  mutate(category = ifelse(is.na(category), "Overall", category),
         category = factor(category, levels = c("US post-1930", "US pre-1930",
                                                "Non-US", "Overall"))) %>%
  arrange(desc(category), desc(first_author)) %>%
  mutate(first_author = factor(first_author, levels=unique(first_author)))

pred_plot_sev <- form_sev_all$pred_comb %>%
  mutate(category = ifelse(is.na(category), "Overall", category),
         category = factor(category, levels = c("US post-1930", "US pre-1930", 
                                                "Non-US", "Overall"))) %>%
  arrange(desc(category), desc(first_author)) %>%
  mutate(first_author = factor(first_author, levels=unique(first_author)))


#TB survival for full model
ggplot(pred_plot_all %>% filter(value != "median"),
             aes(x = est, y = first_author, xmin = cilb, xmax = ciub,
                 shape = shape, color = category)) +
  facet_grid(severity ~ pred_label, scales = "free", space = "free") +
  geom_point(aes(color = category)) +
  geom_point(data = pred_plot_all %>% filter(shape == "Overall" & value != "median"),
             color = 'black', shape = 18, size = 3) +
  geom_errorbar(width = 0.5) +
  scale_x_continuous(name = "1-Year Survival Probability", limits = c(0, 1),
                     breaks = c(0, 0.5, 1)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom") +
  scale_color_manual("Type of Study",
                     values = c("Non-US" = "mediumturquoise", "US post-1930" = "mediumvioletred",
                                "US pre-1930" = "royalblue3", "Overall" = "black")) +
  scale_shape_discrete(guide = "none") +
  ggsave("Figures/forest_full.png", width = 7, height = 6)

#TB survival for stratified model
ggplot(pred_plot_sev %>% filter(value != "median"),
              aes(x = est, y = first_author, xmin = cilb, xmax = ciub, shape = shape,
                  color = category)) +
  geom_point(aes(color = category)) +
  geom_point(data = pred_plot_sev %>% filter(shape == "Overall" & value != "median"),
             color = 'black', shape = 18, size = 3) +
  geom_errorbar(width = 0.5) +
  scale_x_continuous(name = "1-Year Survival Probability", limits = c(0, 1),
                     breaks = c(0, 0.5, 1)) +
  facet_grid(severity ~ pred_label, scales = "free", space = "free") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom") +
  scale_color_manual("Type of Study",
                     values = c("Non-US" = "mediumturquoise", "US post-1930" = "mediumvioletred",
                                "US pre-1930" = "royalblue3", "Overall" = "black")) +
  scale_shape_discrete(guide = "none") +
  ggsave("Figures/forest_stratified.png", width = 7, height = 6)



#### Sanatorium Sensitivity Analysis----------------------------------------------------------------

ggplot(bind_rows(form_san$ind_surv, form_nosan$ind_surv)) +
  facet_wrap(~label, nrow = 2) +
  geom_line(aes(x = x, y = surv, group = study_sev, color = severity),
            size = 1, alpha = 0.3) +
  geom_line(data = bind_rows(form_san$surv_dens, form_nosan$surv_dens),
            aes(x = x, y = surv),
            color = "black", size = 1, linetype = "longdash") +
  geom_smooth(data = bind_rows(form_san$surv_dens, form_nosan$surv_dens),
              aes(x = x, y = surv_est, ymin = cilb, ymax = ciub),
              stat = "identity", linetype = 0, alpha = 0.4) +
  scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
  scale_x_continuous(name = "Years", limits = c(0, 30)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual("", values = c("Minimal" = "seagreen", "Moderately advanced" = "goldenrod1",
                                "Far advanced" = "firebrick2", "Unknown" = "grey50")) +
  ggsave("Figures/sanatorium_curves.png", width = 5.5, height = 8)




