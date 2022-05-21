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
mortalityData <- read.csv("data/mortality_data.csv")
load('R/bayesian_mortality.RData')

#Formatting results
form_all_tb <- formatBayesian(mortalityData, res_all_tb, data_all_tb, "Combined")
form_sev_tb <- formatBayesian(mortalityData, res_sev_tb, data_sev_tb, "Severity", fixed = TRUE)
form_san_tb <- formatBayesian(mortalityData, res_san_tb, data_san_tb, "Sanatorium/hospital")
form_nosan_tb <- formatBayesian(mortalityData, res_nosan_tb, data_nosan_tb, "Non-Sanatorium")
form_all_sub <- formatBayesian(mortalityData, res_all_sub, data_all_sub, "Combined_Sub")
form_sev_sub <- formatBayesian(mortalityData, res_sev_sub, data_sev_sub, "Severity_Sub", fixed = TRUE)


#### All studies summary survival curves -----------------------------------------------------------

#TB survival for full model
p1 <- ggplot(form_all_tb$surv_dens) +
  geom_line(aes(x = x, y = surv),
            color = "black", size = 1, linetype = "solid") +
  geom_smooth(aes(x = x, y = surv_est, ymin = cilb, ymax = ciub),
              stat = "identity", linetype = 0, alpha = 0.25, na.rm = TRUE) +
  scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
  scale_x_continuous(name = "Years", limits = c(0, 30)) +
  theme_bw()


#TB survival for fixed effect model
p2 <- ggplot(form_sev_tb$surv_dens) +
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

grid.arrange(p1, p2, nrow = 2)
ggsave("Figures/survival_curves.png",
       arrangeGrob(p1, p2, nrow = 2),
       width = 5.5, height = 8)



#### All studies individual survival curves --------------------------------------------------------

#TB survival for full model
s1 <- ggplot(form_all_tb$ind_surv) +
  geom_line(aes(x = x, y = surv, group = study_sev, color = severity),
            size = 0.7, alpha = 0.3) +
  geom_line(data = form_all_tb$surv_dens, aes(x = x, y = surv),
            color = "black", size = 1, linetype = "longdash") +
  scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
  scale_x_continuous(name = "Years", limits = c(0, 30)) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual("Disease Severity",
                     values = c("Minimal" = "seagreen", "Moderately advanced" = "goldenrod1",
                                "Far advanced" = "firebrick2", "Unknown" = "grey50"))

#TB survival for stratified model
s2 <- ggplot(form_sev_tb$ind_surv) +
  geom_line(aes(x = x, y = surv, group = study_sev, color = severity),
            size = 0.7, alpha = 0.3) +
  geom_line(data = form_sev_tb$surv_dens,
            aes(x = x, y = surv, color = severity),
            linetype = "longdash", size = 1) +
  scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
  scale_x_continuous(name = "Years", limits = c(0, 30)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual("", drop = FALSE,
                     values = c("Minimal" = "seagreen", "Moderately advanced" = "goldenrod1",
                                "Far advanced" = "firebrick2", "Unknown" = "grey50"))

grid.arrange(s1, s2, nrow = 2)
ggsave("Figures/survival_curves_supp.png",
       arrangeGrob(s1, s2, nrow = 2),
       width = 5.5, height = 8)



#### US post-1930s summary survival curves ---------------------------------------------------------

#TB survival for full model
p1_sub <- ggplot(form_all_sub$surv_dens) +
  geom_line(aes(x = x, y = surv),
            color = "black", size = 1, linetype = "solid") +
  geom_smooth(aes(x = x, y = surv_est, ymin = cilb, ymax = ciub),
              stat = "identity", linetype = 0, alpha = 0.25, na.rm = TRUE) +
  scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
  scale_x_continuous(name = "Years", limits = c(0, 30)) +
  theme_bw()


#TB survival for fixed effect model
p2_sub <- ggplot(form_sev_sub$surv_dens) +
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

grid.arrange(p1_sub, p2_sub, nrow = 2)
ggsave("Figures/survival_curves_subset.png",
       arrangeGrob(p1_sub, p2_sub, nrow = 2),
       width = 5.5, height = 8)



#### US post-1930s individual survival curves ------------------------------------------------------

#TB survival for full model
s1_sub <- ggplot(form_all_sub$ind_surv) +
  geom_line(aes(x = x, y = surv, group = study_sev, color = severity),
            size = 0.7, alpha = 0.3) +
  geom_line(data = form_all_tb$surv_dens, aes(x = x, y = surv),
            color = "black", size = 1, linetype = "longdash") +
  scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
  scale_x_continuous(name = "Years", limits = c(0, 30)) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual("Disease Severity",
                     values = c("Minimal" = "seagreen", "Moderately advanced" = "goldenrod1",
                                "Far advanced" = "firebrick2", "Unknown" = "grey50"))

#TB survival for stratified model
s2_sub <- ggplot(form_sev_sub$ind_surv) +
  geom_line(aes(x = x, y = surv, group = study_sev, color = severity),
            size = 0.7, alpha = 0.3) +
  geom_line(data = form_sev_sub$surv_dens,
            aes(x = x, y = surv, color = severity),
            linetype = "longdash", size = 1) +
  scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
  scale_x_continuous(name = "Years", limits = c(0, 30)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual("", drop = TRUE,
                     values = c("Minimal" = "seagreen", "Moderately advanced" = "goldenrod1",
                                "Far advanced" = "firebrick2", "Unknown" = "grey50"))

grid.arrange(s1_sub, s2_sub, nrow = 2)
ggsave("Figures/survival_curves_supp_subset.png",
       arrangeGrob(s1_sub, s2_sub, nrow = 2),
       width = 5.5, height = 8)


#### All studies survival curves by category -------------------------------------------------------

ggplot(form_sev_tb$ind_surv) +
  geom_line(aes(x = x, y = surv, group = study_sev, color = category),
            size = 0.7, alpha = 0.3) +
  geom_line(data = form_sev_tb$surv_dens,
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

pred_plot_all <- form_all_tb$pred_comb %>%
  mutate(category = ifelse(is.na(category), "Overall", category),
         category = factor(category, levels = c("US post-1930", "US pre-1930",
                                                "Non-US", "Overall"))) %>%
  arrange(desc(category), desc(first_author)) %>%
  mutate(first_author = factor(first_author, levels=unique(first_author)))

pred_plot_sev <- form_sev_tb$pred_comb %>%
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

#TB survival for full model
s1s <- ggplot(bind_rows(form_san_tb$ind_surv, form_nosan_tb$ind_surv)) +
  facet_wrap(~label, nrow = 2) +
  geom_line(aes(x = x, y = surv, group = study_sev, color = severity),
            size = 1, alpha = 0.3) +
  geom_line(data = bind_rows(form_san_tb$surv_dens, form_nosan_tb$surv_dens),
            aes(x = x, y = surv),
            color = "black", size = 1, linetype = "longdash") +
  geom_smooth(data = bind_rows(form_san_tb$surv_dens, form_nosan_tb$surv_dens),
              aes(x = x, y = surv_est, ymin = cilb, ymax = ciub),
              stat = "identity", linetype = 0, alpha = 0.4) +
  scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
  scale_x_continuous(name = "Years", limits = c(0, 30)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual("", values = c("Minimal" = "seagreen", "Moderately advanced" = "goldenrod1",
                                "Far advanced" = "firebrick2", "Unknown" = "grey50")) +
  ggsave("Figures/sanatorium_curves.png", width = 5.5, height = 8)



#### Table of Main Results -------------------------------------------------------------------------

#Combining raw tables from results lists
raw_tab <- bind_rows(form_all_tb$param, form_sev_tb$param,
                     form_san_tb$param, form_nosan_tb$param,
                     form_all_sub$param, form_sev_sub$param)

#Adding formatted, ordered labels
raw_tab <- raw_tab %>%
  mutate(Severity = ifelse(label == "Severity", severity,
                           ifelse(label == "Severity_Sub", paste0(severity, "_Sub"),
                                  label)),
         Severity = factor(Severity, level = c("Minimal", "Moderate", "Advanced", "Combined",
                                               "Minimal_Sub", "Moderate_Sub", "Advanced_Sub",
                                               "Combined_Sub", "Sanatorium/hospital",
                                               "Non-Sanatorium"))) %>%
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


#Finding number of papers, cohorts, individuals for each analysis
data_all_tb2 <- as.data.frame(t(data_all_tb[[2]]))
data_san_tb2 <- as.data.frame(t(data_san_tb[[2]]))
data_nosan_tb2 <- as.data.frame(t(data_nosan_tb[[2]]))
data_all_sub2 <- as.data.frame(t(data_all_sub[[2]]))
data_all_tb2$Severity <- "Combined"
data_san_tb2$Severity <- "Sanatorium/hospital"
data_nosan_tb2$Severity <- "Non-Sanatorium"
data_all_sub2$Severity <- "Combined_Sub"
counts_tb <- bind_rows(data_all_tb2, data_san_tb2, data_nosan_tb2, data_all_sub2)

#Counts stratified by severity for all studies
mortalityData_sev <- mortalityData %>% filter(severity != "Unknown")
counts_sev <- mortalityData_sev %>%
  group_by(severity) %>%
  summarize(nStudies = length(unique(study_id)),
            nCohorts = length(unique(cohort_id)),
            nIndividuals = n(),
            .groups = "drop") %>%
  rename(Severity = severity)

#Counts stratified by severity for US post-1930s subset
mortalityData_sub <- mortalityData %>% filter(severity != "Unknown",
                                              study_id %in% c("1029", "93", "45"))
counts_sub <- mortalityData_sub %>%
  group_by(severity) %>%
  summarize(nStudies = length(unique(study_id)),
            nCohorts = length(unique(cohort_id)),
            nIndividuals = n(),
            .groups = "drop") %>%
  mutate(Severity = paste0(severity, "_Sub"))


counts_initial <- bind_rows(counts_tb, counts_sev, counts_sub) %>%
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
main_tab <- final_tab %>% filter(!grepl("San|Sub", Severity))
sub_tab <- final_tab %>% filter(grepl("Sub", Severity))
san_tab <- final_tab %>% filter(grepl("San", Severity))

#Variance of frailty terms
theta <- raw_tab %>%
  filter(value == "theta") %>%
  mutate(theta = round(est, 2)) %>%
  select(label, theta)
