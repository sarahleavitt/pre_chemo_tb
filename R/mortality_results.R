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

# Formatting results
form_all <- formatBayesian(mortalityData, res_all, data_all, "All-cause mortality: Combined")
form_sev <- formatBayesian(mortalityData, res_sev, data_sev, "All-cause mortality: Severity", fixed = TRUE)
form_all_tb <- formatBayesian(mortalityData, res_all_tb, data_all_tb, "TB-specific mortality: Combined")
form_sev_tb <- formatBayesian(mortalityData, res_sev_tb, data_sev_tb, "TB-specific mortality: Severity", fixed = TRUE)

form_san <- formatBayesian(mortalityData, res_san, data_san, "All-cause mortality: Sanatorium")
form_nosan <- formatBayesian(mortalityData, res_nosan, data_nosan, "All-cause mortality: Non-Sanatorium")
form_san_tb <- formatBayesian(mortalityData, res_san_tb, data_san_tb, "TB-specific mortality: Sanatorium")
form_nosan_tb <- formatBayesian(mortalityData, res_nosan_tb, data_nosan_tb, "TB-specific mortality: Non-Sanatorium")


#### Main text survival curves --------------------------------------------------------------------

#All-cause survival for full model
s1 <- ggplot(data = form_all$surv_dens) +
  geom_line(aes(x = x, y = surv),
            color = "black", size = 1, linetype = "solid") +
  geom_smooth(aes(x = x, y = surv_est, ymin = cilb, ymax = ciub),
              stat = "identity", linetype = 0, alpha = 0.25, na.rm = TRUE) +
  scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
  scale_x_continuous(name = "Years", limits = c(0, 30)) +
  theme_bw() +
  ggtitle("All-cause mortality")

#TB survival for full model
s2 <- ggplot(form_all_tb$surv_dens) +
  geom_line(aes(x = x, y = surv),
            color = "black", size = 1, linetype = "solid") +
  geom_smooth(aes(x = x, y = surv_est, ymin = cilb, ymax = ciub),
              stat = "identity", linetype = 0, alpha = 0.25, na.rm = TRUE) +
  scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
  scale_x_continuous(name = "Years", limits = c(0, 30)) +
  theme_bw() +
  ggtitle("TB-specific mortality")


#Overall survival for fixed effects model
s1f <- ggplot(form_sev$surv_dens) +
  geom_line(aes(x = x, y = surv, color = severity),
            size = 1, linetype = "solid") +
  geom_smooth(aes(x = x, y = surv_est, ymin = cilb, ymax = ciub, fill = severity),
              stat = "identity", linetype = 0, alpha = 0.15, na.rm = TRUE) +
  scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
  scale_x_continuous(name = "Years", limits = c(0, 30)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual("", drop = TRUE,
                     values = c("Minimal" = "seagreen", "Moderately advanced" = "goldenrod1",
                                "Far advanced" = "firebrick2", "Unknown" = "grey50")) +
  scale_fill_manual("", guide = FALSE,
                    values = c("Minimal" = "seagreen", "Moderately advanced" = "goldenrod1",
                               "Far advanced" = "firebrick2", "Unknown" = "grey50"))

#TB survival for fixed effect model
s2f <- ggplot(form_sev_tb$surv_dens) +
  geom_line(aes(x = x, y = surv, color = severity),
            size = 1, linetype = "solid") +
  geom_smooth(aes(x = x, y = surv_est, ymin = cilb, ymax = ciub, fill = severity),
              stat = "identity", linetype = 0, alpha = 0.15, na.rm = TRUE) +
  scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
  scale_x_continuous(name = "Years", limits = c(0, 30)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual("", guide = FALSE,
                     values = c("Minimal" = "seagreen", "Moderately advanced" = "goldenrod1",
                                "Far advanced" = "firebrick2", "Unknown" = "grey50")) +
  scale_fill_manual("",
                    values = c("Minimal" = "seagreen", "Moderately advanced" = "goldenrod1",
                               "Far advanced" = "firebrick2", "Unknown" = "grey50"))

grid.arrange(s1, s2, s1f, s2f, nrow = 2)
ggsave("Figures/survival_curves.png",
       arrangeGrob(s1, s2, s1f, s2f, nrow = 2),
       width = 9, height = 8.5)



#### Supplement survival curves -------------------------------------------------------------------

#All-cause survival for full model
s1 <- ggplot(form_all$ind_surv) +
  geom_line(aes(x = x, y = surv, group = study_sev, color = severity),
            size = 0.7, alpha = 0.3) +
  geom_line(data = form_all$surv_dens, aes(x = x, y = surv),
            color = "black", size = 1, linetype = "longdash") +
  scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
  scale_x_continuous(name = "Years", limits = c(0, 30)) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual("Disease Severity",
                     values = c("Minimal" = "seagreen", "Moderately advanced" = "goldenrod1",
                                "Far advanced" = "firebrick2", "Unknown" = "grey50")) +
  ggtitle("All-cause mortality")

#TB survival for full model
s2 <- ggplot(form_all_tb$ind_surv) +
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
                                "Far advanced" = "firebrick2", "Unknown" = "grey50")) +
  ggtitle("TB-specific mortality")


#Overall survival for fixed effects model
s1f <- ggplot(form_sev$ind_surv) +
  geom_line(aes(x = x, y = surv, group = study_sev, color = severity),
            size = 0.7, alpha = 0.3) +
  geom_line(data = form_sev$surv_dens, aes(x = x, y = surv, color = severity),
            size = 1, linetype = "longdash") +
  scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
  scale_x_continuous(name = "Years", limits = c(0, 30)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual("", drop = FALSE,
                     values = c("Minimal" = "seagreen", "Moderately advanced" = "goldenrod1",
                                "Far advanced" = "firebrick2", "Unknown" = "grey50"),
                     labels = c("Minimal", "Moderately\nadvanced", "Far\nadvanced", "Unknown"))

#TB survival for fixed effect model
s2f <- ggplot(form_sev_tb$ind_surv) +
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
                                "Far advanced" = "firebrick2", "Unknown" = "grey50"),
                     labels = c("Minimal", "Moderately\nadvanced", "Far\nadvanced", "Unknown"))

grid.arrange(s1, s2, s1f, s2f, nrow = 2)
ggsave("Figures/survival_curves_supp.png",
       arrangeGrob(s1, s2, s1f, s2f, nrow = 2),
       width = 9, height = 8.5)



##### Forest plots --------------------------------------------------------------------------------

#All-cause survival for full model
f1 <- ggplot(form_all$pred_comb %>% filter(value != "median"),
       aes(x = est, y = first_author, xmin = cilb, xmax = ciub, shape = shape)) +
  geom_point(color = "black") +
  geom_point(data = form_all$pred_comb %>% filter(shape == "Overall" & value != "median"),
             color = 'black', shape = 18, size = 3) +
  geom_errorbar(width = 0.5) +
  scale_x_continuous(name = "Survival Probability", limits = c(0, 1),
                     breaks = c(0, 0.5, 1)) +
  facet_grid(severity ~ pred_label, scales = "free", space = "free") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_blank())

#TB survival for full model
f2 <- ggplot(form_all_tb$pred_comb %>% filter(value != "median"),
             aes(x = est, y = first_author, xmin = cilb, xmax = ciub, shape = shape)) +
  geom_point(color = "black") +
  geom_point(data = form_all_tb$pred_comb %>% filter(shape == "Overall" & value != "median"),
             color = 'black', shape = 18, size = 3) +
  geom_errorbar(width = 0.5) +
  scale_x_continuous(name = "1-Year Survival Probability", limits = c(0, 1),
                     breaks = c(0, 0.5, 1)) +
  facet_grid(severity ~ pred_label, scales = "free", space = "free") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_blank())

#All-cause survival for fixed effect model
f1f <- ggplot(form_sev$pred_comb %>% filter(value != "median"),
       aes(x = est, y = first_author, xmin = cilb, xmax = ciub, shape = shape)) +
  geom_point(color = "black") +
  geom_point(data = form_sev$pred_comb %>% filter(shape == "Overall" & value != "median"),
             color = 'black', shape = 18, size = 3) +
  geom_errorbar(width = 0.5) +
  scale_x_continuous(name = "Survival Probability", limits = c(0, 1),
                     breaks = c(0, 0.5, 1)) +
  facet_grid(severity ~ pred_label, scales = "free", space = "free") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_blank())

#TB survival for fixed effect model
f2f <- ggplot(form_sev_tb$pred_comb %>% filter(value != "median"),
              aes(x = est, y = first_author, xmin = cilb, xmax = ciub, shape = shape)) +
  geom_point(color = "black") +
  geom_point(data = form_sev_tb$pred_comb %>% filter(shape == "Overall" & value != "median"),
             color = 'black', shape = 18, size = 3) +
  geom_errorbar(width = 0.5) +
  scale_x_continuous(name = "1-Year Survival Probability", limits = c(0, 1),
                     breaks = c(0, 0.5, 1)) +
  facet_grid(severity ~ pred_label, scales = "free", space = "free") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_blank())

grid.arrange(f1, f2, f1f, f2f, nrow = 2)
ggsave("Figures/forest_plots.png",
       arrangeGrob(f1, f2, f1f, f2f, nrow = 2),
       height = 9, width = 11)



#### Sanatorium Sensitivity Analysis----------------------------------------------------------------

#Overall survival for full model
s1s <- ggplot(bind_rows(form_san$ind_surv, form_nosan$ind_surv)) +
  facet_wrap(~label) +
  geom_line(aes(x = x, y = surv, group = study_sev, color = severity),
            size = 1, alpha = 0.3) +
  geom_line(data = bind_rows(form_san$surv_dens, form_nosan$surv_dens), aes(x = x, y = surv),
            color = "black", size = 1, linetype = "longdash") +
  geom_smooth(data = bind_rows(form_san$surv_dens, form_nosan$surv_dens),
              aes(x = x, y = surv_est, ymin = cilb, ymax = ciub),
              stat = "identity", linetype = 0, alpha = 0.4) +
  scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
  scale_x_continuous(name = "Years", limits = c(0, 30)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual("Disease Severity",
                     values = c("Minimal" = "seagreen", "Moderately advanced" = "goldenrod1",
                                "Far advanced" = "firebrick2", "Unknown" = "grey50")) 

#TB survival for full model
s2s <- ggplot(bind_rows(form_san_tb$ind_surv, form_nosan_tb$ind_surv)) +
  facet_wrap(~label) +
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
  scale_color_manual("Disease Severity",
                     values = c("Minimal" = "seagreen", "Moderately advanced" = "goldenrod1",
                                "Far advanced" = "firebrick2", "Unknown" = "grey50"))

grid.arrange(s1s, s2s, nrow = 2)
ggsave("Figures/sanatorium_curves.png", arrangeGrob(s1s, s2s, nrow = 2),
       height = 7, width = 7)



#### Table of Main Results -------------------------------------------------------------------------

#Combining raw tables from results lists
raw_tab <- bind_rows(form_all$param, form_sev$param, form_san$param, form_nosan$param,
                     form_all_tb$param, form_sev_tb$param, form_san_tb$param, form_nosan_tb$param)

#Adding formatted, ordered labels
raw_tab <- raw_tab %>%
  mutate(Severity = ifelse(is.na(severity) & grepl("Combined", label), "Combined",
                           ifelse(grepl("Non-", label), "Non-Sanatorium",
                                  ifelse(grepl("San", label), "Sanatorium", severity))),
         Severity = factor(Severity, level = c("Minimal", "Moderate", "Advanced", "Combined",
                                               "Sanatorium", "Non-Sanatorium")),
         Outcome = gsub(": [A-z]*-*[A-z]*", "", label)) %>%
  arrange(Outcome, Severity)

#Extracting the survival probabilities
pred1_tab <- raw_tab %>%
  filter(value == "pred1") %>%
  mutate(`1-Year Survival (95% CI)` = paste0(round(est, 2), " (",
                            round(cilb, 2), ", ",
                            round(ciub, 2), ")")) %>%
    select(Outcome, Severity, `1-Year Survival (95% CI)`)

pred5_tab <- raw_tab %>%
  filter(value == "pred5") %>%
  mutate(`5-Year Survival (95% CI)` = paste0(round(est, 2), " (",
                                            round(cilb, 2), ", ",
                                            round(ciub, 2), ")")) %>%
  select(Outcome, Severity, `5-Year Survival (95% CI)`)

pred10_tab <- raw_tab %>%
  filter(value == "pred10") %>%
  mutate(`10-Year Survival (95% CI)` = paste0(round(est, 2), " (",
                                            round(cilb, 2), ", ",
                                            round(ciub, 2), ")")) %>%
  select(Outcome, Severity, `10-Year Survival (95% CI)`)


#Extracting the distribution parameters
sdlog <- raw_tab %>%
  filter(value == "sdlog") %>%
  select(label, sdlog = est)

dist_tab <- raw_tab %>%
  filter(value == "meanlog") %>%
  full_join(sdlog, by = "label") %>%
  mutate(`Survival Distribution` = paste0("lognormal(", round(est, 2), ", ",
                                          round(sdlog, 2), ")")) %>%
  select(Outcome, Severity, `Survival Distribution`)


#Finding number of papers, cohorts, individuals for each analysis
data_all <- as.data.frame(t(data_all[[2]]))
data_san <- as.data.frame(t(data_san[[2]]))
data_nosan <- as.data.frame(t(data_nosan[[2]]))
data_all_tb <- as.data.frame(t(data_all_tb[[2]]))
data_san_tb <- as.data.frame(t(data_san_tb[[2]]))
data_nosan_tb <- as.data.frame(t(data_nosan_tb[[2]]))
data_all$Severity <- "Combined"
data_san$Severity <- "Sanatorium"
data_nosan$Severity <- "Non-Sanatorium"
data_all_tb$Severity <- "Combined"
data_san_tb$Severity <- "Sanatorium"
data_nosan_tb$Severity <- "Non-Sanatorium"

counts_all <- bind_rows(data_all, data_san, data_nosan)
counts_all$Outcome <- "All-cause mortality"
counts_tb <- bind_rows(data_all_tb, data_san_tb, data_nosan_tb)
counts_tb$Outcome <- "TB-specific mortality"

mortalityData_sev <- mortalityData %>% filter(severity != "Unknown") #same for both all-cause and tb-specific
counts_sev <- mortalityData_sev %>%
  group_by(severity) %>%
  summarize(nStudies = length(unique(study_id)),
            nCohorts = length(unique(cohort_id)),
            nIndividuals = n(),
            .groups = "drop") %>%
  rename(Severity = severity)

counts_initial <- bind_rows(counts_all, counts_tb, counts_sev) %>%
  replace_na(list(Outcome = "All-cause mortality")) %>%
  bind_rows(counts_sev) %>%
  replace_na(list(Outcome = "TB-specific mortality")) %>%
  select(Outcome, Severity, `Number of Studies` = nStudies, `Number of Cohorts` = nCohorts,
         `Number of Individuals` = nIndividuals)
  

#Combining the tables
final_tab <- dist_tab %>%
  full_join(pred1_tab, by = c("Outcome", "Severity")) %>%
  full_join(pred5_tab, by = c("Outcome", "Severity")) %>%
  full_join(pred10_tab, by = c("Outcome", "Severity")) %>%
  full_join(counts_initial, by = c("Outcome", "Severity"))

#Separating sanatorium
main_tab <- final_tab %>% filter(!grepl("San", Severity))
san_tab <- final_tab %>% filter(grepl("San", Severity))

#Variance of frailty terms
theta <- raw_tab %>%
  filter(value == "theta") %>%
  mutate(theta = round(est, 2)) %>%
  select(label, theta)
