# Data-centric SDM Step 6) Analyze and visualize differences between models
## Sara Hansen
## Modified July 25, 2024

library(conflicted)
library(tidyverse)
library(tools)
library(janitor)
library(rstatix)
library(plotrix)
library(ggtext)
library(ggpubr)
library(ggrepel)
library(ggpattern)
library(ggh4x)
library(mgcv)
library(xgboost)
library(FedData)
library(caret) 
library(epm)
library(USA.state.boundaries)
library(terra)
library(sf)
library(tidyterra)
library(FedData)

options(max.print = 1e9)
options(dplyr.print_max = 1e9)

############################ READ IN MODEL METRICS #############################
# Read in metrics calculated in Step 5
metricsList <- list.files(path = "output", pattern = "metrics.txt", full.names = TRUE)
names(metricsList) <- tools::file_path_sans_ext(basename(metricsList))
# Including .txt makes sure this excludes any of the "_alt" ones from Step 9

metrics <- lapply(metricsList, read.delim, header = TRUE, sep = "\t",
                  fileEncoding = "UTF-8", row.names = NULL) %>%
  bind_rows(.id = "id") %>%
  rename(evalThreshold = threshold) %>%
  mutate(varRemovalThreshold = case_when(str_detect(id, "_B") ~ 0.5,
                                         TRUE ~ 0.7))

rm(metricsList)
################################################################################

# The 5 questions are: How are model results and performance affected by 
# 1) the method of removing missing values in explanatory data prior to modeling (by row or by column)?
# 2) the threshold for excluding correlated explanatory variables prior to modeling (0.7 or 0.5)?
# 3) the scale of response data (Michigan or Saginaw Bay)?
# 4) the source of response data (specimens, observations, or both)? 
# 5) the type of response data (presence-absence or presence-background)?

# While evaluating each question, we will control for the others

########################## DATA PROCESSING QUESTIONS ###########################

# Question 1: Method of removing missing data
# Variable of interest is removalType

# First exclude any groups where a comparison can't be made
# So any case where we can't compare colrm and rowrm groups
metrics %>%
  count(scale, modelType, dataSource, dataType, varRemovalThreshold, removalType)
# In Question 1 below, we'll exclude:
# mich, ANN, "all" data at 0.7 threshold (the col remove group did not build)
# mich, ANN, "observation" at 0.5 threshold
# mich, BRT, "specimen" at 0.5 threshold
# sag, BRT, presence-absence
# sag, ANN, presence-background
# sag, GAM, "presence-absence"

# Note that Saginaw Bay had only one varRemovalThreshold, labeled as 0.7

question1 <- metrics %>%
  mutate(remove = case_when(scale == "mich" & modelType == "ANN" & dataSource == "all" &
                              varRemovalThreshold == 0.7 ~ "remove",
                            scale == "mich" & modelType == "ANN" & dataSource == "observation" &
                              varRemovalThreshold == 0.5 ~ "remove",
                            scale == "mich" & modelType == "BRT" & dataSource == "specimen" &
                              varRemovalThreshold == 0.5 ~ "remove",
                            scale == "sag" & modelType == "BRT" & dataType == "presence-absence" ~ "remove",
                            scale == "sag" & modelType == "ANN" & dataType == "presence-background" ~ "remove",
                            scale == "sag" & modelType == "GAM" & dataType == "presence-absence" ~ "remove",
                            TRUE ~ "keep")) %>% dplyr::filter(remove != "remove") %>%
  group_by(modelClass, modelType, scale, dataSource, dataType, varRemovalThreshold, removalType) %>%
  mutate(count = n()) %>% dplyr::filter(count >=3) %>% dplyr::select(-c(remove, count)) %>%
  ungroup() %>%
  group_by(modelClass, modelType, scale, dataSource, dataType, varRemovalThreshold) %>%
  rstatix::wilcox_test(mae ~ removalType, paired = FALSE) %>%
  mutate(sig = case_when(p > 0.05 ~ "ns", p <= 0.05 ~ "s")) %>%
  left_join(metrics %>%
              group_by(modelClass, modelType, scale, dataSource, dataType, 
                       removalType, varRemovalThreshold) %>%
              summarize(mean = mean(mae)) %>%
              pivot_wider(names_from = removalType,
                          values_from = mean) %>%
              dplyr::rename(colrmMean = 'col remove',
                            rowrmMean = 'row remove'),
            by = c("modelClass", "modelType", "scale", "dataSource", 
                   "dataType", "varRemovalThreshold"))

# Optional AUC ROC evaluation
#question1auc <- metrics %>%
  #mutate(remove = case_when(scale == "mich" & modelType == "ANN" & dataSource == "all" &
                              #varRemovalThreshold == 0.7 ~ "remove",
                            #scale == "mich" & modelType == "ANN" & dataSource == "observation" &
                              #varRemovalThreshold == 0.5 ~ "remove",
                            #scale == "mich" & modelType == "BRT" & dataSource == "specimen" &
                              #varRemovalThreshold == 0.5 ~ "remove",
                            #scale == "sag" & modelType == "BRT" & dataType == "presence-absence" ~ "remove",
                            #scale == "sag" & modelType == "ANN" & dataType == "presence-background" ~ "remove",
                            #scale == "sag" & modelType == "GAM" & dataType == "presence-absence" ~ "remove",
                            #TRUE ~ "keep")) %>% dplyr::filter(remove != "remove" & modelClass != "Envelope") %>%
  #group_by(modelClass, modelType, scale, dataSource, dataType, varRemovalThreshold, removalType) %>%
  #mutate(count = n()) %>% dplyr::filter(count >=3) %>% dplyr::select(-c(remove, count)) %>%
  #ungroup() %>%
  #group_by(modelClass, modelType, scale, dataSource, dataType, varRemovalThreshold) %>%
  #rstatix::wilcox_test(aucROC ~ removalType, paired = FALSE) %>%
  #mutate(sig = case_when(p > 0.05 ~ "ns", p <= 0.05 ~ "s")) %>%
  #left_join(metrics %>%
              #group_by(modelClass, modelType, scale, dataSource, dataType, 
                       #removalType, varRemovalThreshold) %>%
              #summarize(mean = mean(aucROC)) %>%
              #pivot_wider(names_from = removalType,
                          #values_from = mean) %>%
              #dplyr::rename(colrmMean = 'col remove',
                            #rowrmMean = 'row remove'),
            #by = c("modelClass", "modelType", "scale", "dataSource", 
                   #"dataType", "varRemovalThreshold"))

#question1maeauc <- question1 %>%
  #dplyr::select(modelClass, modelType, scale, dataSource, dataType,
                #varRemovalThreshold, mae_P = p) %>%
  #left_join(question1auc %>%
              #dplyr::select(modelClass, modelType, scale, dataSource, dataType,
                            #varRemovalThreshold, auc_P = p),
            #by = c("modelClass", "modelType", "scale", "dataSource", "dataType",
                   #"varRemovalThreshold"))

#question1maeauc %>%
  #dplyr::filter(mae_P > 0.05 & auc_P <= 0.05)
# AUC was significant in 5 cases when MAE wasn't

#question1maeauc %>%
  #dplyr::filter(mae_P <= 0.05 & auc_P > 0.05)
# MAE was significant 10 times when AUC wasn't significant
# not including the envelope models where AUC couldn't be used at all

# Supplementary table for Question 1 (Appendix S7 in manuscript)
suppTable1 <- question1 %>%
  mutate(across(.cols = c(colrmMean, rowrmMean),
                ~ round(., 4))) %>%
  mutate(modelClass = replace(modelClass, modelClass == "Machine learning", "ML")) %>%
  mutate(across(.cols = c(scale, dataSource, dataType, sig),
                ~ str_to_sentence(.))) %>%
  dplyr::select("Model class" = modelClass, "Model type" = modelType,
                "Scale" = scale, "Data source" = dataSource, "Data type" = dataType,
                "Correlation threshold for removing explanatory variables" = varRemovalThreshold,
                "Group 1 n" = n1, "Group 2 n" = n2,
                "Wilcox statistic" = statistic, "P" = p, "Significance" = sig,
                "Mean MAE for variable-removed models" = colrmMean,
                "Mean MAE for observation-removed models" = rowrmMean)

#write_delim(suppTable1, "tables/suppTable1.txt", delim = "\t", quote = "none")

# Results table for Question 1 (Table 3 in manuscript)
# "Var" refers to variable, meaning the models with missing values removed by variable (column)
# "Obs" refer to observation, meaning the models with missing values removed by observation (row)
table1top <- question1 %>%
  dplyr::filter(dataType != "presence-only") %>%
  mutate(sigValue = case_when(colrmMean < rowrmMean & sig == "s" ~ "Var *",
                              colrmMean < rowrmMean & sig == "ns" ~ "Var",
                              rowrmMean < colrmMean & sig == "s" ~ "Obs *",
                              rowrmMean < colrmMean & sig == "ns" ~ "Obs")) %>%  
  dplyr::select(modelClass, modelType, scale, dataSource, dataType, varRemovalThreshold, sigValue) %>%
  pivot_wider(id_cols = c(modelClass, modelType),
              names_from = c(varRemovalThreshold, scale, dataType, dataSource), values_from = sigValue) %>%
  arrange(factor(modelClass, levels = c("Machine learning", "Regression", "Envelope"))) %>%
  dplyr::select("Model class" = modelClass, "Model type" = modelType,
                "0.7_mich_presence-background_all",
                "0.7_mich_presence-background_specimen",
                "0.7_mich_presence-background_observation",
                "0.7_sag_presence-background_all",
                "0.7_sag_presence-absence_all",
                "0.5_mich_presence-background_all",
                "0.5_mich_presence-background_specimen",
                "0.5_mich_presence-background_observation")

table1bottom <- question1 %>%
  dplyr::filter(dataType == "presence-only") %>%
  mutate(sigValue = case_when(colrmMean < rowrmMean & sig == "s" ~ "Var *",
                              colrmMean < rowrmMean & sig == "ns" ~ "Var",
                              rowrmMean < colrmMean & sig == "s" ~ "Obs *",
                              rowrmMean < colrmMean & sig == "ns" ~ "Obs")) %>%  
  dplyr::select(modelClass, modelType, scale, dataSource, varRemovalThreshold, sigValue) %>%
  pivot_wider(id_cols = c(modelClass, modelType),
              names_from = c(varRemovalThreshold, scale, dataSource), values_from = sigValue) %>%
  arrange(factor(modelClass, levels = c("Machine learning", "Regression", "Envelope"))) %>%
  dplyr::select("Model class" = modelClass, "Model type" = modelType,
                "0.7_mich_all",
                "0.7_mich_specimen",
                "0.7_mich_observation",
                "0.7_sag_all",
                "0.5_mich_all",
                "0.5_mich_specimen",
                "0.5_mich_observation")

# table1top and table1bottom will be combined for the paper
#write_delim(table1top, "tables/table1top.txt", delim = "\t", quote = "none")
#write_delim(table1bottom, "tables/table1bottom.txt", delim = "\t", quote = "none")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# To control for the method of removing missing values for future questions,
# retain only the better performing group of 4 in each case

betterMAEgroup <- metrics %>%
  group_by(modelClass, modelType, scale, dataType, dataSource, removalType, varRemovalThreshold) %>%
  mutate(count = n()) %>%
  #exclude any where one of the groups has fewer than 3
  dplyr::filter(count >= 3) %>%
  dplyr::select(-count) %>%
  mutate(removalType = case_when(removalType == "row remove" ~ "rowrm",
                                 removalType == "col remove" ~ "colrm")) %>%
  summarize(mean = mean(mae)) %>%
  pivot_wider(names_from = removalType, values_from = mean) %>%
  mutate(betterMAE = case_when(rowrm < colrm | (is.na(colrm) & !is.na(rowrm)) ~ "row remove",
                               rowrm > colrm | (is.na(rowrm) & !is.na(colrm)) ~ "col remove")) %>%
  dplyr::select(-c(rowrm, colrm))

metrics2 <- metrics %>%
  right_join(betterMAEgroup, by = c("modelClass", "modelType", "scale", "dataType", "dataSource",
                                    "varRemovalThreshold", "removalType" = "betterMAE"))

rm(betterMAEgroup)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Question 2: Threshold for excluding correlated explanatory variables prior to modeling
# Variable of interest is varRemovalThreshold

# Note that Saginaw Bay had the same explanatory variable set at both thresholds
# so no comparisons can be made at the Saginaw Bay scale

# Identify any groups that cannot be evaluated
metrics2 %>%
  dplyr::filter(scale != "sag") %>%
  count(scale, modelType, dataSource, varRemovalThreshold)
# mich, GAM, "specimen" only built at 0.5 threshold, so exclude

question2 <- metrics2 %>%
  mutate(remove = case_when(scale == "sag" ~ "remove",
                            scale == "mich" & modelType == "GAM" & dataSource == "specimen" ~ "remove",
                            TRUE ~ "keep")) %>%
  dplyr::filter(remove != "remove") %>%
  dplyr::select(-remove) %>%
  arrange(modelClass, modelType) %>%
  group_by(modelClass, scale, modelType, dataSource, dataType) %>%
  rstatix::wilcox_test(mae ~ varRemovalThreshold, paired = FALSE, p.adjust.method = "bonferroni") %>%
  mutate(sig = case_when(p >= 0.05 ~ "ns",
                         p < 0.05 ~ "s")) %>%
  left_join(metrics2 %>%
              group_by(modelClass, scale, modelType, dataSource, dataType, varRemovalThreshold) %>%
              summarize(mean = mean(mae)) %>%
              pivot_wider(names_from = varRemovalThreshold,
                          values_from = mean),
            by = c("modelClass", "modelType", "scale", "dataType", "dataSource"))

# Optional AUC ROC evaluation
#question2auc <- metrics2 %>%
  #mutate(remove = case_when(scale == "sag" ~ "remove",
                            #scale == "mich" & modelType == "GAM" & dataSource == "specimen" ~ "remove",
                            #TRUE ~ "keep")) %>%
  #dplyr::filter(remove != "remove" & modelClass != "Envelope") %>%
  #dplyr::select(-remove) %>%
  #arrange(modelClass, modelType) %>%
  #group_by(modelClass, scale, modelType, dataSource, dataType) %>%
  #rstatix::wilcox_test(aucROC ~ varRemovalThreshold, paired = FALSE) %>%
  #mutate(sig = case_when(p >= 0.05 ~ "ns",
                         #p < 0.05 ~ "s")) %>%
  #left_join(metrics2 %>%
              #group_by(modelClass, scale, modelType, dataSource, dataType, varRemovalThreshold) %>%
              #summarize(mean = mean(aucROC)) %>%
              #pivot_wider(names_from = varRemovalThreshold,
                          #values_from = mean),
            #by = c("modelClass", "modelType", "scale", "dataType", "dataSource"))

#question2maeauc <- question2 %>%
  #dplyr::select(modelClass, modelType, scale, dataSource, dataType, mae_P = p) %>%
  #left_join(question2auc %>%
              #dplyr::select(modelClass, modelType, scale, dataSource, dataType, auc_P = p),
            #by = c("modelClass", "modelType", "scale", "dataSource", "dataType"))

#question2maeauc %>%
  #dplyr::filter(mae_P <= 0.05 |auc_P <= 0.05)
# MAE was only significant for Envelope models, so AUC was never significant

# Supplementary table for Question 2 (Appendix S8 in manuscript)
suppTable2 <- question2 %>%
  mutate(across(.cols = c("0.5", "0.7"),
                ~ round(., 4))) %>%
  mutate(modelClass = replace(modelClass, modelClass == "Machine learning", "ML")) %>%
  mutate(across(.cols = c(scale, dataSource, dataType, sig),
                ~ str_to_sentence(.))) %>%
  dplyr::select("Model class" = modelClass, "Model type" = modelType, "Scale" = scale, "Data source" = dataSource, 
                "Data type" = dataType, "Group 1 n" = n1, "Group 2 n" = n2,
                "Wilcox statistic" = statistic, "P" = p, "Significance" = sig,
                "Mean MAE for 0.5 variable removal threshold" = "0.5", 
                "Mean MAE for 0.7 variable removal threshold" = "0.7")

#write_delim(suppTable2, "tables/suppTable2.txt", delim = "\t", quote = "none")

# Results table for Question 2 (Table 4 in manuscript)
table2top <- question2 %>%
  dplyr::filter(modelClass != "Envelope") %>%
  mutate(sigValue = case_when(`0.5` < `0.7` & sig == "s" ~ "0.5 *",
                              `0.5` < `0.7` & sig == "ns" ~ "0.5",
                              `0.7` < `0.5` & sig == "s" ~ "0.7 *",
                              `0.7` < `0.5` & sig == "ns" ~ "0.7")) %>%
  dplyr::select(modelClass, modelType, scale, dataSource, dataType, sigValue) %>%
  pivot_wider(id_cols = c(modelClass, modelType),
              names_from = c(scale, dataSource, dataType), values_from = sigValue) %>%
  dplyr::select("Model class" = modelClass, "Model type" = modelType,
                "mich_all_presence-background", "mich_specimen_presence-background",
                "mich_observation_presence-background")

table2bottom <- question2 %>%
  dplyr::filter(modelClass == "Envelope") %>%
  mutate(sigValue = case_when(`0.5` < `0.7` & sig == "s" ~ "0.5 *",
                              `0.5` < `0.7` & sig == "ns" ~ "0.5",
                              `0.7` < `0.5` & sig == "s" ~ "0.7 *",
                              `0.7` < `0.5` & sig == "ns" ~ "0.7")) %>%
  dplyr::select(modelClass, modelType, scale, dataSource, sigValue) %>%
  pivot_wider(id_cols = c(modelClass, modelType),
              names_from = c(scale, dataSource), values_from = sigValue) %>%
  dplyr::select("Model class" = modelClass, "Model type" = modelType,
                "mich_all", "mich_specimen", "mich_observation")

#write_delim(table2top, "tables/table2top.txt", delim = "\t", quote = "none")
#write_delim(table2bottom, "tables/table2bottom.txt", delim = "\t", quote = "none")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# To control for the threshold for excluding correlated variables,
# retain only the better performing group of 4 in each case (even though most aren't significant)
betterMAEgroup2 <-  metrics2 %>%
  group_by(modelClass, modelType, scale, dataType, dataSource, varRemovalThreshold) %>%
  mutate(count = n()) %>%
  #exclude any where one of the groups has fewer than 3
  dplyr::filter(count >= 3) %>%
  dplyr::select(-count) %>%
  summarize(mean = mean(mae)) %>%
  pivot_wider(names_from = varRemovalThreshold, values_from = mean) %>%
  mutate(betterMAE = case_when(`0.5` < `0.7` | (is.na(`0.7`) & !is.na(`0.5`)) ~ 0.5,
                               `0.5` > `0.7` | (is.na(`0.5`) & !is.na(`0.7`)) ~ 0.7)) %>%
  dplyr::select(-c(`0.5`, `0.7`))

metrics3 <- metrics2 %>%
  right_join(betterMAEgroup2, by = c("modelClass", "modelType", "scale", "dataType", "dataSource",
                                    "varRemovalThreshold" = "betterMAE"))

rm(betterMAEgroup2)

# Optional test to see sensitivity and specificity distribution
#summary(metrics3$sensitivity) # median = 0.874 (wihout Envelope models it is 0.889)
#summary(metrics3$specificity) # median = 0.8375

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

################################################################################

################# RELATIVE MODEL PERFORMANCE BY TYPE AND CLASS #################

# These are not part of the 5 questions
# How did each model perform overall?

# Order by Mean Absolute Error
modelTypeOrder <- metrics3 %>%
  group_by(modelType) %>%
  summarise(meanMAE = mean(mae)) %>%
  arrange(meanMAE)

# Question 0a: Mean Absolute Error based on model type and class
question0a <- metrics3 %>%
  left_join(modelTypeOrder,
            by = "modelType") %>%
  mutate(modelType = fct_reorder(modelType, meanMAE)) %>%
  rstatix::wilcox_test(mae ~ modelType, paired = FALSE, p.adjust.method = "bonferroni") %>%
  bind_rows(data.frame(group1 = modelTypeOrder$modelType,
                       group2 = modelTypeOrder$modelType,
                       p = NA)) %>%
  left_join(modelTypeOrder %>%
              rename(group1MeanMAE = meanMAE),
            by = c("group1" = "modelType")) %>%
  left_join(modelTypeOrder %>%
              rename(group2MeanMAE = meanMAE),
            by = c("group2" = "modelType")) %>%
  mutate(group1ModelClass = 
           case_when(str_detect(group1, "BIOCLIM|DOMAIN") ~ "Envelope",
                     str_detect(group1, "GLM|GAM|MARS|MaxEnt|MaxNet") ~ "Regression",
                     str_detect(group1, "BRT|RF|Cforest|ANN|XGBoost") ~ "Machine learning"),
         group2ModelClass =
           case_when(str_detect(group2, "BIOCLIM|DOMAIN") ~ "Envelope",
                     str_detect(group2, "GLM|GAM|MARS|MaxEnt|MaxNet") ~ "Regression",
                     str_detect(group2, "BRT|RF|Cforest|ANN|XGBoost") ~ "Machine learning"))

# Visualize summary statistics for each model type
plot0a1 <- metrics3 %>%
  ggplot(aes(x = mae, y = reorder(modelType, -mae), color = modelClass,
             fill = modelClass)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 4) +
  geom_vline(xintercept = mean(metrics3$mae),
             linetype = 2, color = "#96939b") +
  theme_classic() +
  scale_color_manual(values = c("#F9BFC0", "#e3a800", "#3A7D1C")) +
  scale_fill_manual(values = alpha(c("#F9BFC0", "#e3a800", "#3A7D1C"), 0.25)) +
  labs(x = "Mean absolute error (MAE)", y = "Model type",
       fill = "Model class", color = "Model class")+
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(margin = ggplot2::margin(t = 0.25, r = 0, b = 0, l = 0, unit = "cm")),
        axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 0.25, b = 0, l = 0, unit = "cm")),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10),
        legend.position = c(0.8, 0.7),
        legend.key.height = unit(0.6, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.spacing.y = unit(0.15, "cm"),
        plot.margin = ggplot2::margin(t = 1, r = 0.5, b = 1, l = 0.5, unit = "cm")) +
  guides(fill = guide_legend(byrow = TRUE),
         color = guide_legend(byrow = TRUE))

plot0a1

#ggsave("figures/plot0a1.pdf", width = 7.25, height = 5, units = "in")
#ggsave("figures/plot0a1.png", width = 7.25, height = 5, units = "in", dpi = 600)

# Visualize statistically significant pairwise differences between model types
question0aForPlot <- question0a %>%
  mutate(pfill = case_when(is.na(p) ~ "Not applicable",
                           p <=0.05 ~ "Significant (P < 0.05)",
                           p > 0.05 ~ "Not significant"),
         p = round(p, 3),
         borderColor = case_when(group1ModelClass == "Envelope" &
                                   group2ModelClass == "Envelope" &
                                   group1 != group2 ~ "Envelope",
                                 group1ModelClass == "Regression" &
                                   group2ModelClass == "Regression" &
                                   group1 != group2 ~ "Regression",
                                 group1ModelClass == "Machine learning" &
                                   group2ModelClass == "Machine learning" &
                                   group1 != group2 ~ "Machine learning",
                                 TRUE ~ "Across classes"),
         borderWidth = case_when(borderColor == "Across classes" ~ 0.5,
                                 TRUE ~ 1))

plot0a2 <- ggplot(question0aForPlot) +
  geom_tile(aes(x = fct_reorder(group1, group1MeanMAE), 
                y = fct_reorder(group2, -group2MeanMAE),
                fill = pfill),
            height = 0.82, width = 0.82) +
  geom_tile(aes(x = fct_reorder(group1, group1MeanMAE), 
                y = fct_reorder(group2, -group2MeanMAE),
                color = borderColor), alpha = 0,
            linewidth = question0aForPlot$borderWidth,
            height = 0.82, width = 0.82, lineend = "square") +
  coord_equal() +
  theme_classic() +
  scale_color_manual(values = c("#000000", "#F9BFC0", "#e3a800", "#3A7D1C"),
                     name = "Group comparison") +
  scale_fill_manual(values = c("#000000", "#F2F1F3", "#7F7B85"),
                    name = "Significance",
                    labels = c("*Not applicable*", "Not significant",
                               "Significant (*P* < 0.05)")) +  
  theme(legend.position = c(0.95, 0.7),
        legend.spacing.y = unit(0.15, "cm"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1),
        legend.text = element_markdown(size = 9),
        legend.key.height = unit(0.6, "cm"),
        legend.key.width = unit(0.6, "cm"),
        plot.margin = margin(t = 0.0761035, r = 1.415, b = 0.0761035, l = 1.415, unit = "in")) +
  guides(color = guide_legend(override.aes = list(linewidth = c(0.5, 1, 1, 1)),
                              byrow = TRUE),
         fill = guide_legend(byrow = TRUE))
plot0a2

# Note on margin: Default is 5.5 points, and there are 72.27 points per inch,
# so the default in inches is 0.0761035 inches, which we'll keep for the top and bottom.
# By setting the left and right margins larger,
# we fill in the empty space that would come from combining with plot0a1.

#ggsave("figures/plot0a2.pdf", width = 7.25, height = 5, units = "in")
#ggsave("figures/plot0a2.png", width = 7.25, height = 5, units = "in", dpi = 600)

# Comparisons between machine learning models and GAM were usually not significant
# but comparisons between GAM and regression models, and sometimes other regression models were
# Machine learning models are not all alike, but they tend to perform similarly

# Combined plot0a1 and plot0a2 into one figure, plot0a (Figure 1 in manuscript)
plot0a <- ggarrange(plot0a1, plot0a2, nrow = 2, labels = c("A", "B"))
plot0a
#ggsave("figures/plot0a.pdf", width = 7.25, height = 9, units = "in")
#ggsave("figures/plot0a.png", width = 7.25, height = 9, units = "in", dpi = 600)

# Are differences in MAE bewteen the three model classes overall significant?
metrics3 %>%
  rstatix::wilcox_test(mae ~ modelClass, paired = FALSE, p.adjust.method = "bonferroni")
# Machine learning and regression are not significantly differently, likely because of GAM

metrics3 %>%
  group_by(modelClass) %>%
  summarize(mean(mae))
# Machine learning is best (0.246), then Regression (0.281), then Envelope (0.557; not all significant)

# Question 0b: Sensitivity and Specificity based on model type and class
question0b <- metrics3 %>%
  dplyr::filter(modelClass != "Envelope") %>%
  group_by(modelType) %>%
  summarize(meanSens = mean(sensitivity),
            meanSpec = mean(specificity),
            seSens = std.error(sensitivity),
            seSpec = std.error(specificity)) %>%
  dplyr::mutate(upperSens = meanSens + seSens,
                lowerSens = meanSens - seSens,
                upperSpec = meanSpec + seSpec,
                lowerSpec = meanSpec - seSpec,
                modelClass = case_when(
                  str_detect(modelType, "GLM|GAM|MARS|MaxEnt|MaxNet") ~ "Regression",
                  str_detect(modelType, "BRT|RF|Cforest|ANN|XGBoost") ~ "Machine learning"))

# Build convex hull around error bars to visualize ranges in each model calss
outerPoints <- bind_rows(question0b %>% dplyr::select(modelClass, x = lowerSens, y = meanSpec),
                         question0b %>% dplyr::select(modelClass, x = upperSens, y = meanSpec),
                         question0b %>% dplyr::select(modelClass, x = meanSens, y = lowerSpec),
                         question0b %>% dplyr::select(modelClass, x = meanSens, y = upperSpec))

find_hull <- function(outerPoints) outerPoints[chull(outerPoints$x, outerPoints$y), ]
hulls <- plyr::ddply(outerPoints, "modelClass", find_hull)

# Table to show means and significance for each metric in the plot
table0b <- metrics3 %>%
  dplyr::filter(modelClass != "Envelope") %>%
  group_by(modelClass) %>%
  summarize(Sensitivity = round(mean(sensitivity), 3),
            Specificity = round(mean(specificity), 3)) %>%
  t() %>%
  janitor::row_to_names(1) %>%
  as.data.frame() %>%
  mutate("P-value" = c(
    metrics3 %>% dplyr::filter(modelClass != "Envelope") %>%
      rstatix::wilcox_test(sensitivity ~ modelClass, paired = FALSE) %>% 
      pull(p) %>% round(3),
    metrics3 %>% dplyr::filter(modelClass != "Envelope") %>%
      rstatix::wilcox_test(specificity ~ modelClass, paired = FALSE) %>% 
      pull(p) %>% round(3))) %>%
  ggpubr::ggtexttable(theme = ttheme("classic", base_size = 9,
                                     colnames.style = colnames_style(linewidth = 1,
                                                                     linecolor = "#000000",
                                                                     fill = "#F2F1F3",
                                                                     color = "#000000",
                                                                     size = 8))) %>%
  ggpubr::table_cell_bg(row = 2:3, column = 2, fill = "#e3a800", 
                        color = "#000000", alpha = 0.2) %>%
  ggpubr::table_cell_bg(row = 2:3, column = 3, fill = "#28560b", 
                        color = "#000000", alpha = 0.2) %>%
  tab_add_border(from.col = 2, linewidth = 1.5) %>%
  tab_add_hline(from.col = 2, linewidth = 1) %>%
  tab_add_vline(from.row = 2, linewidth = 1) %>%
  thead_add_border(from.col = 2, linewidth = 1.5) %>%
  tab_add_title(text = "    Mean values by model class", 
                face = "italic",
                size = 9, color = "#96939b",
                padding = unit(0.2, "in"))
table0b

# Visualize metrics by model type and class (Figure 2 in manuscript)
plot0b <- ggplot(data = question0b,
                 aes(x = meanSens, y = meanSpec, 
                     color = modelClass, label = modelType)) +
  geom_pointrange(aes(ymin = lowerSpec, ymax = upperSpec),
                  linewidth = 1/.pt, linetype = 2) +
  geom_pointrange(aes(xmin = lowerSens, xmax = upperSens),
                  linewidth = 1/.pt, linetype = 2) +
  ggrepel::geom_label_repel(size = 8/.pt,
                           box.padding = 0.4) +
  theme_classic() +
  scale_color_manual(values = c("#e3a800", "#28560b")) +
  geom_polygon(data = hulls, aes(x = x, y = y, fill = modelClass,
                                 label = NULL),
               alpha = 0.1, linetype = 0) +
  scale_fill_manual(values = c("#e3a800", "#28560b")) +
  labs(x = "Sensitivity", y = "Specificity",
       color = "Model class", fill = "Model class") +
  annotation_custom(ggplotGrob(table0b), # adding the table created before
                    xmin = 0.74, xmax = 0.82,
                    ymin = 0.57, ymax = 0.67) + 
  theme(legend.position = c(0.89, 0.85),
        legend.background = element_blank(),
        legend.spacing.y = unit(0.15, "cm"),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(margin = ggplot2::margin(t = 0.25, r = 0, b = 0, l = 0, unit = "cm")),
        axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 0.25, b = 0, l = 0, unit = "cm")),
        axis.text = element_text(size = 9),
        legend.title = element_text(size = 10),
        legend.text = element_markdown(size = 9),
        legend.key.height = unit(0.25, "in"),
        legend.key.width = unit(0.25, "in")) +
  guides(fill = guide_legend(byrow = TRUE),
         color = guide_legend(byrow = TRUE))

plot0b
#ggsave("figures/plot0b.pdf", width = 7.25, height = 6, units = "in")
#ggsave("figures/plot0b.png", width = 7.25, height = 6, units = "in", dpi = 600)

################################################################################

################### DATA INPUT QUESTIONS: MODEL PERFORMANCE ####################

# Question 3: Scale of response data
# Variable of interest is scale
# Use only presence-background or presence-only data, because Michigan has no absences

question3 <- bind_rows(
  metrics3 %>%
    dplyr::filter(dataType != "presence-absence" &
                    dataSource == "all") %>%
    group_by(modelType) %>%
    rstatix::wilcox_test(mae ~ scale, paired = FALSE, p.adjust.method = "bonferroni") %>%
    left_join(metrics3 %>%
                group_by(modelType, scale) %>%
                summarize(meanMAE = mean(mae)) %>%
                pivot_wider(names_from = "scale",
                            values_from = "meanMAE"),
              by = "modelType"),
  metrics3 %>%
    dplyr::filter(dataType != "presence-absence" &
                    dataSource == "all") %>%
    group_by(modelType) %>%
    rstatix::wilcox_test(sensitivity ~ scale, paired = FALSE, p.adjust.method = "bonferroni") %>%
    left_join(metrics3 %>%
                group_by(modelType, scale) %>%
                summarize(meanSens = mean(sensitivity)) %>%
                pivot_wider(names_from = "scale",
                            values_from = "meanSens"),
              by = "modelType"),
  metrics3 %>%
    dplyr::filter(dataType == "presence-background" &
                    dataSource == "all") %>%
    group_by(modelType) %>%
    rstatix::wilcox_test(specificity ~ scale, paired = FALSE) %>%
    left_join(metrics3 %>%
                group_by(modelType, scale) %>%
                summarize(meanSpec = mean(specificity)) %>%
                pivot_wider(names_from = "scale",
                            values_from = "meanSpec"),
              by = "modelType")) %>%
  mutate(sig = case_when(p <= 0.05 ~ "s",
                         p > 0.05 ~ "ns")) %>%
  dplyr::rename(metric = .y.) %>%
  mutate(metric = case_when(metric == "mae" ~ "Mean absolute error (MAE)",
                            metric == "sensitivity" ~ "Sensitivity (true positive rate)",
                            metric == "specificity" ~ "Specificity (true negative rate)"))

# Optional AUC evaluation
#question3auc <- metrics3 %>%
  #dplyr::filter(dataType != "presence-absence" & dataSource == "all" &
                  #modelClass != "Envelope") %>%
  #group_by(modelType) %>%
  #rstatix::wilcox_test(aucROC ~ scale, paired = FALSE) %>%
  #left_join(metrics3 %>%
              #group_by(modelType, scale) %>%
              #summarize(meanAUC = mean(aucROC)) %>%
              #pivot_wider(names_from = "scale",
                         #values_from = "meanAUC"),
            #by = "modelType")

#question3maeauc <- question3auc %>%
  #dplyr::select(modelType, auc_P = p) %>%
  #left_join(question3 %>% dplyr::filter(metric == "Mean absolute error (MAE)") %>%
              #dplyr::select(modelType, mae_P = p),
            #by = "modelType")

#question3maeauc %>%
  #dplyr::filter(mae_P <= 0.05 | auc_P <= 0.05)
# AUC doesn't provide new information and misses significance for ANN

# Supplementary table for Question 3 (Appendix S9 in manuscript)
suppTable3 <- question3 %>%
  mutate(across(.cols = c(mich, sag),
                ~ round(., 4))) %>%
  mutate(sig = str_to_sentence(sig)) %>%
  dplyr::select("Model type" = modelType, "Metric" = metric, 
                "Group 1 n" = n1, "Group 2 n" = n2,
                "Wilcox statistic" = statistic, "P" = p, "Significance" = sig,
                "Mean value for Michigan" = mich, 
                "Mean value for Saginaw Bay" = sag)

#write_delim(suppTable3, "figures/suppTable3.txt", delim = "\t", quote = "none")

# Visualize differences in each metric between scales
question3Long <- question3 %>%
  pivot_longer(cols = c(mich, sag),
               names_to = "scale",
               values_to = "meanValue") %>%
  mutate(meanValue = case_when(metric == "Mean absolute error (MAE)" ~ 1 - meanValue,
                               TRUE ~ meanValue),
         metric = case_when(metric == "Mean absolute error (MAE)" ~ "1 - mean absolute error (MAE)",
                            TRUE ~ metric))

# Plot for Question 3 (Figure 3 in manuscript)
plot3 <- ggplot(question3Long, aes(x = meanValue, 
                                   y = reorder(modelType, desc(modelType)))) +
  geom_line() +
  geom_point(aes(color = scale, group = metric),
             size = 4) +
  geom_text(data = question3Long %>%
              mutate(sig = case_when(sig == "s" ~ "*",
                                     sig == "ns" ~ "")) %>%
              arrange(metric, modelType, meanValue) %>%
              distinct(modelType, metric, sig, .keep_all = TRUE),
            aes(label = sig), position = position_nudge(x = -0.085, y = -0.1),
            size = 6) +
  facet_wrap(vars(metric)) +
  theme_classic() +
  scale_x_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00),
                     expand = c(0.04, 0.04)) +
  labs(x = "Mean value", y = "Model type",
       color = "Scale") +
  scale_color_manual(values = c("#37406B", "#96BBDA"),
                     labels = c("Michigan", "Saginaw Bay")) +
  theme(axis.line = element_blank(),
        panel.border = element_rect(color = "#000000", fill = NA,
                                    linewidth = 1),
        panel.grid.major.y = element_line(color = "#CDCBCF"),
        panel.spacing = unit(0.49, "cm"),
        strip.background = element_rect(color = "#000000", fill = NA,
                                        linewidth = 1),
        legend.position = "bottom",
        axis.title = element_text(size = 10),
        axis.title.x = element_text(margin = ggplot2::margin(t = 0.25, r = 0, b = 0, l = 0, unit = "cm")),
        axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 0.25, b = 0, l = 0, unit = "cm")),
        axis.text = element_text(size = 9),
        legend.title = element_text(size = 10),
        legend.key.height = unit(0.25, "in"),
        legend.key.width = unit(0.25, "in"),
        plot.margin = ggplot2::margin(t = 9, r = 9, b = 9, l = 9, unit = "pt"))

plot3

#ggsave("figures/plot3.pdf", width = 7.25, height = 5, units = "in")
#ggsave("figures/plot3.png", width = 7.25, height = 5, units = "in", dpi = 600)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Question 4: Source of response data
# Variable of interest is dataSource
# Michigan scale only (Saginaw Bay is one year of research data, no sources)
# Note that the observation and specimen datasets are independent from one another,
# but not from the "all" dataset, so we will not be comparing with all data

question4 <- bind_rows(
  metrics3 %>%
    dplyr::filter(dataType != "presence-absence" & scale == "mich" &
                    dataSource != "all") %>%
    group_by(modelType) %>%
    rstatix::wilcox_test(mae ~ dataSource, paired = FALSE) %>%
    left_join(metrics3 %>%
                group_by(modelType, dataSource) %>%
                summarize(meanMAE = mean(mae)) %>%
                pivot_wider(names_from = "dataSource",
                            values_from = "meanMAE"),
              by = "modelType"),
  metrics3 %>%
    dplyr::filter(dataType != "presence-absence" & scale == "mich" &
                    dataSource != "all") %>%
    group_by(modelType) %>%
    rstatix::wilcox_test(sensitivity ~ dataSource, paired = FALSE) %>%
    left_join(metrics3 %>%
                group_by(modelType, dataSource) %>%
                summarize(meanSens = mean(sensitivity)) %>%
                pivot_wider(names_from = "dataSource",
                            values_from = "meanSens"),
              by = "modelType"),
  metrics3 %>%
    dplyr::filter(dataType == "presence-background" & scale == "mich" &
                    dataSource != "all") %>%
    group_by(modelType) %>%
    rstatix::wilcox_test(specificity ~ dataSource, paired = FALSE) %>%
    left_join(metrics3 %>%
                group_by(modelType, dataSource) %>%
                summarize(meanSpec = mean(specificity)) %>%
                pivot_wider(names_from = "dataSource",
                            values_from = "meanSpec"),
              by = "modelType")) %>%
  mutate(sig = case_when(p <= 0.05 ~ "s",
                         p > 0.05 ~ "ns")) %>%
  dplyr::rename(metric = .y.) %>%
  mutate(metric = case_when(metric == "mae" ~ "Mean absolute error (MAE)",
                            metric == "sensitivity" ~ "Sensitivity (true positive rate)",
                            metric == "specificity" ~ "Specificity (true negative rate)"))

# Optional AUC evaluation:
#question4auc <- metrics3 %>%
  #dplyr::filter(dataType == "presence-background" & scale == "mich" & 
                  #dataSource != "all" & modelClass != "Envelope") %>%
  #group_by(modelType) %>%
  #rstatix::wilcox_test(aucROC ~ dataSource, paired = FALSE) %>%
  #left_join(metrics3 %>%
              #dplyr::filter(dataSource != "all") %>%
              #group_by(modelType, dataSource) %>%
              #summarize(meanAUC = mean(aucROC)) %>%
              #pivot_wider(names_from = "dataSource",
                          #values_from = "meanAUC"),
            #by = "modelType")

#question4maeauc <- question4auc %>%
  #dplyr::select(modelType, auc_P = p) %>%
  #left_join(question4 %>% dplyr::filter(metric == "Mean absolute error (MAE)") %>%
              #dplyr::select(modelType, mae_P = p),
            #by = "modelType")

#question4maeauc %>%
  #dplyr::filter(mae_P <= 0.05 | auc_P <= 0.05)
# AUC doesn't provide new information and misses several comparisons' significance

# Supplementary table for Question 4 (Appendix S10 in manuscript)
suppTable4 <- question4 %>%
  mutate(across(.cols = c(specimen, observation),
                ~ round(., 4))) %>%
  mutate(across(.cols = c(group1, group2, sig),
                ~ str_to_sentence(.))) %>%
  dplyr::select("Model type" = modelType, "Metric" = metric, 
                "Group 1 n" = n1, "Group 2 n" = n2,
                "Wilcox statistic" = statistic, "P" = p, 
                "Significance" = sig,
                "Mean value for specimens only" = specimen,
                "Mean value for observations only" = observation)

#write_delim(suppTable4, "tables/suppTable4.txt", delim = "\t", quote = "none")

# Visualize differences in each metric between data sources
question4Long <- question4 %>%
  pivot_longer(cols = c(observation, specimen),
               names_to = "dataSource",
               values_to = "meanValue") %>%
  mutate(meanValue = case_when(metric == "Mean absolute error (MAE)" ~ 1 - meanValue,
                               TRUE ~ meanValue),
         metric = case_when(metric == "Mean absolute error (MAE)" ~ "1 - mean absolute error (MAE)",
                            TRUE ~ metric))

# Plot for Question 4 (Figure 4 in manuscript)
plot4 <- ggplot(question4Long, aes(x = meanValue, 
                                   y = reorder(modelType, desc(modelType)))) +
  geom_line() +
  geom_point(aes(color = dataSource, group = metric),
             size = 4) +
  geom_text(data = question4Long %>%
              mutate(sig = case_when(sig == "s" ~ "*",
                                     sig == "ns" ~ "")) %>%
              arrange(metric, modelType, meanValue) %>%
              distinct(modelType, metric, sig, .keep_all = TRUE),
            aes(label = sig), position = position_nudge(x = -0.085, y = -0.1),
            size = 6) +
  facet_wrap(vars(metric)) +
  theme_classic() +
  scale_x_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00),
                     expand = c(0.04, 0.04)) +
  labs(x = "Mean value", y = "Model type",
       color = "Data source") +
  scale_color_manual(values = c("#37406B", "#96BBDA"),
                     labels = c("Observations only", "Specimens only")) +
  theme(axis.line = element_blank(),
        panel.border = element_rect(color = "#000000", fill = NA,
                                    linewidth = 1),
        panel.grid.major.y = element_line(color = "#CDCBCF"),
        panel.spacing = unit(0.49, "cm"),
        strip.background = element_rect(color = "#000000", fill = NA,
                                        linewidth = 1),
        legend.position = "bottom",
        axis.title = element_text(size = 10),
        axis.title.x = element_text(margin = ggplot2::margin(t = 0.25, r = 0, b = 0, l = 0, unit = "cm")),
        axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 0.25, b = 0, l = 0, unit = "cm")),
        legend.title = element_text(size = 10),
        legend.key.height = unit(0.25, "in"),
        legend.key.width = unit(0.25, "in"),
        plot.margin = ggplot2::margin(t = 9, r = 9, b = 9, l = 9, unit = "pt"))

plot4

#ggsave("figures/plot4.pdf", width = 7.25, height = 5, units = "in")
#ggsave("figures/plot4.png", width = 7.25, height = 5, units = "in", dpi = 600)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Question 5: Type of response data
# Variable of interest is dataType
# Saginaw Bay only, presence-background vs. presence-absence

question5 <- bind_rows(
  metrics3 %>%
    dplyr::filter(scale == "sag" & 
                    dataType != "presence-only") %>%
    group_by(modelType) %>%
    rstatix::wilcox_test(mae ~ dataType, paired = FALSE) %>%
    left_join(metrics3 %>%
                group_by(modelType, dataType) %>%
                summarize(meanMAE = mean(mae)) %>%
                pivot_wider(names_from = "dataType",
                            values_from = "meanMAE"),
              by = "modelType"),
  metrics3 %>%
    dplyr::filter(scale == "sag" & 
                    dataType != "presence-only") %>%
    group_by(modelType) %>%
    rstatix::wilcox_test(sensitivity ~ dataType, paired = FALSE) %>%
    left_join(metrics3 %>%
                group_by(modelType, dataType) %>%
                summarize(meanSens = mean(sensitivity)) %>%
                pivot_wider(names_from = "dataType",
                            values_from = "meanSens"),
              by = "modelType"),
  metrics3 %>%
    dplyr::filter(scale == "sag" & 
                    dataType != "presence-only") %>%
    group_by(modelType) %>%
    rstatix::wilcox_test(specificity ~ dataType, paired = FALSE) %>%
    left_join(metrics3 %>%
                group_by(modelType, dataType) %>%
                summarize(meanSpec = mean(specificity)) %>%
                pivot_wider(names_from = "dataType",
                            values_from = "meanSpec"),
              by = "modelType")) %>%
  mutate(sig = case_when(p <= 0.05 ~ "s",
                         p > 0.05 ~ "ns")) %>%
  dplyr::rename(metric = .y.) %>%
  mutate(metric = case_when(metric == "mae" ~ "Mean absolute error (MAE)",
                            metric == "sensitivity" ~ "Sensitivity (true positive rate)",
                            metric == "specificity" ~ "Specificity (true negative rate)"))

# Optional AUC evaluation:
#question5auc <- metrics3 %>%
  #dplyr::filter(scale == "sag" & dataType != "presence-only") %>%
  #group_by(modelType) %>%
  #rstatix::wilcox_test(aucROC ~ dataType, paired = FALSE, p.adjust.method = "bonferroni") %>%
  #left_join(metrics3 %>%
              #group_by(modelType, dataType) %>%
              #summarize(meanAUC = mean(aucROC)) %>%
              #pivot_wider(names_from = "dataType",
                          #values_from = "meanAUC"),
            #by = "modelType")

#question5maeauc <- question5auc %>%
  #dplyr::select(modelType, auc_P = p) %>%
  #left_join(question5 %>% dplyr::filter(metric == "Mean absolute error (MAE)") %>%
              #dplyr::select(modelType, mae_P = p),
            #by = "modelType")

#question5maeauc %>%
  #dplyr::filter(mae_P <= 0.05 | auc_P <= 0.05)
# AUC and MAE are each significant once when the other is not

# Supplementary table for Question 5 (Appendix S11 in manuscript)
suppTable5 <- question5 %>%
  mutate(across(.cols = c(`presence-absence`, `presence-background`),
                ~ round(., 4))) %>%
  mutate(sig = str_to_sentence(sig)) %>%
  dplyr::select("Model type" = modelType, "Metric" = metric, 
                "Group 1 n" = n1, "Group 2 n" = n2,
                "Wilcox statistic" = statistic, "P" = p, "Significance" = sig,
                "Mean value for presence-absence data" = `presence-absence`, 
                "Mean value for presence-background data" = `presence-background`)

#write_delim(suppTable5, "tables/suppTable5.txt", delim = "\t", quote = "none")

# Visualize differences in each metric between data types
question5Long <- question5 %>%
  pivot_longer(cols = c(`presence-absence`, `presence-background`),
               names_to = "dataType",
               values_to = "meanValue") %>%
  mutate(meanValue = case_when(metric == "Mean absolute error (MAE)" ~ 1 - meanValue,
                               TRUE ~ meanValue),
         metric = case_when(metric == "Mean absolute error (MAE)" ~ "1 - mean absolute error (MAE)",
                            TRUE ~ metric))

# Plot for Question 5 (Figure 5 in manuscript)
plot5 <- ggplot(question5Long, aes(x = meanValue, 
                                   y = reorder(modelType, desc(modelType)))) +
  geom_line() +
  geom_point(aes(color = dataType, group = metric),
             size = 4) +
  geom_text(data = question5Long %>%
              mutate(sig = case_when(sig == "s" ~ "*",
                                     sig == "ns" ~ "")) %>%
              arrange(metric, modelType, meanValue) %>%
              distinct(modelType, metric, sig, .keep_all = TRUE),
            aes(label = sig), position = position_nudge(x = -0.085, y = -0.1),
            size = 6) +
  facet_wrap(vars(metric)) +
  theme_classic() +
  scale_x_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00),
                     expand = c(0.04, 0.04)) +
  labs(x = "Mean value", y = "Model type",
       color = "Data type") +
  scale_color_manual(values = c("#37406B", "#96BBDA"),
                     labels = c("Presence-absence", "Presence-background")) +
  theme(axis.line = element_blank(),
        panel.border = element_rect(color = "#000000", fill = NA,
                                    linewidth = 1),
        panel.grid.major.y = element_line(color = "#CDCBCF"),
        panel.spacing = unit(0.49, "cm"),
        strip.background = element_rect(color = "#000000", fill = NA,
                                        linewidth = 1),
        legend.position = "bottom",
        axis.title = element_text(size = 10),
        axis.title.x = element_text(margin = ggplot2::margin(t = 0.25, r = 0, b = 0, l = 0, unit = "cm")),
        axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 0.25, b = 0, l = 0, unit = "cm")),
        axis.text = element_text(size = 9),
        legend.title = element_text(size = 10),
        legend.key.height = unit(0.25, "in"),
        legend.key.width = unit(0.25, "in"),
        plot.margin = ggplot2::margin(t = 9, r = 9, b = 9, l = 9, unit = "pt"))

plot5

#ggsave("figures/plot5.pdf", width = 7.25, height = 4.5, units = "in")
#ggsave("figures/plot5.png", width = 7.25, height = 4.5, units = "in", dpi = 600)

################################################################################

##################### DATA INPUT QUESTIONS: MODEL OUTCOMES #####################
# To assess model outcomes based on our three data input variables, we will:
## 1) Find the model type that performs the best most often
## 2) Rebuild that type of model, using each entire dataset
## 3) Assess the outcomes of those models (variable contribution and predictions)

topModelMetrics <- metrics3 %>%
  dplyr::filter(dataType != "presence-only") %>%
  group_by(scale, dataSource, dataType, modelType) %>%
  summarize(mean = mean(mae)) %>%
  arrange(scale, dataSource, dataType, mean) 

topModelMetrics %>%
  group_by(scale, dataSource, dataType) %>%
  dplyr::slice(1:3)

metrics %>%
  group_by(scale, dataSource, dataType, modelType) %>%
  count()

# XGBoost and ANN were each the top model twice
# but XGBoost was always in the top three, unlike ANN
# (GAM was also always in top three, but only the top for Michigan specimens)
# XGBoost was also the one one of those three that never failed in the entire model set

# So XGBoost is a good choice for comparing across scales and data sources/types
# We'll produce a final "best" model for each entire dataset

# We want to remove missing values in the best way for each dataset
# and use the best variable removal threshold for the XGBoost models
metrics3 %>%
  dplyr::filter(modelType == "XGBoost") %>%
  count(scale, dataSource, dataType, modelType, removalType, varRemovalThreshold)
# 0.7 is always the better varRemovalThreshold, use the individual removalTypes indicated

# Ready to build models
# Read in the five datasets -> all the 0.7 variable removal threshold
michPB <- read.delim("output/michPB.txt", header = TRUE,
                     sep = "\t", quote = "", fileEncoding = "UTF-8-BOM")
michSpecPB <- read.delim("output/michSpecPB.txt", header = TRUE,
                         sep = "\t", quote = "", fileEncoding = "UTF-8-BOM")
michObsPB <- read.delim("output/michObsPB.txt", header = TRUE,
                        sep = "\t", quote = "", fileEncoding = "UTF-8-BOM")
sagPB <- read.delim("output/sagPB.txt", header = TRUE,
                    sep = "\t", quote = "", fileEncoding = "UTF-8-BOM")
sagPA <- read.delim("output/sagPA.txt", header = TRUE,
                    sep = "\t", quote = "", fileEncoding = "UTF-8-BOM")

# This function will build models, show importance of each variable, 
# predict on test data, and predict across study area
buildTopModels <- function() {
  
  #rowrm
  michPB_train0 <- michPB %>%
    dplyr::mutate(land = as.factor(land)) %>%
    drop_na() %>%
    dplyr::select(-c(x, y, id))
  michPB_train <- model.matrix(~ ., michPB_train0) %>% as.data.frame() %>% 
    dplyr::select(-1) %>% as.matrix()
  
  #colrm
  michSpecPB_train0 <- michSpecPB %>%
    dplyr::mutate(land = as.factor(land)) %>%
    dplyr::select(where(~ !any(is.na(.)))) %>% 
    dplyr::select(-c(x, y, id))
  michSpecPB_train <- model.matrix(~ ., michSpecPB_train0) %>% as.data.frame() %>% 
    dplyr::select(-1) %>% as.matrix()
  
  #rowrm
  michObsPB_train0 <- michObsPB %>%
    dplyr::mutate(land = as.factor(land)) %>%
    drop_na() %>%
    dplyr::select(-c(x, y, id))
  michObsPB_train <- model.matrix(~ ., michObsPB_train0) %>% as.data.frame() %>% 
    dplyr::select(-1) %>% as.matrix()
  
  #rowrm
  set.seed(789)
  sagPB_train0 <- sagPB %>%
    dplyr::mutate(land = as.factor(land)) %>%
    drop_na() %>%
    dplyr::select(-c(x, y, id))
  sagPB_train <- model.matrix(~ ., sagPB_train0) %>% as.data.frame() %>% 
    dplyr::select(-1) %>% as.matrix()
  
  #colrm
  sagPA_train0 <- sagPA %>%
    dplyr::mutate(land = as.factor(land)) %>%
    dplyr::select(where(~ !any(is.na(.)))) %>% 
    dplyr::select(-c(x, y, id))
  sagPA_train <- model.matrix(~ ., sagPA_train0) %>% as.data.frame() %>% 
    dplyr::select(-1) %>% as.matrix()
  
  q1 <- ncol(michPB_train)
  
  q2 <- ncol(michSpecPB_train)
  
  q3 <- ncol(michObsPB_train)
  
  q4 <- ncol(sagPB_train)
  
  q5 <- ncol(sagPA_train)
  
  # Build models (keeping parameters the same as before)
  
  set.seed(512)
  michPB_XGBoost = xgboost::xgboost(data = michPB_train[, 2:q1],
                                    label = michPB_train[, 1],
                                    nrounds = 5000,
                                    max_depth = 10,
                                    eta = 0.8,
                                    objective = "binary:logistic",
                                    verbose = 0)
  
  michPB_varImportance = xgboost::xgb.importance(colnames(michPB_XGBoost),
                                                 model = michPB_XGBoost)
    
  michPB_XGBoost_n = nrow(michPB_train)
  
  set.seed(512)
  michSpecPB_XGBoost = xgboost::xgboost(data = michSpecPB_train[, 2:q2],
                                        label = michSpecPB_train[, 1],
                                        nrounds = 5000,
                                        max_depth = 10,
                                        eta = 0.8,
                                        objective = "binary:logistic",
                                        verbose = 0)
  
  michSpecPB_varImportance = xgboost::xgb.importance(colnames(michSpecPB_XGBoost),
                                                     model = michSpecPB_XGBoost)
    
  michSpecPB_XGBoost_n = nrow(michSpecPB_train)
    
  set.seed(512)
  michObsPB_XGBoost = xgboost::xgboost(data = michObsPB_train[, 2:q3],
                                       label = michObsPB_train[, 1],
                                       nrounds = 5000,
                                       max_depth = 10,
                                       eta = 0.8,
                                       objective = "binary:logistic",
                                       verbose = 0)
    
  michObsPB_varImportance = xgboost::xgb.importance(colnames(michObsPB_XGBoost),
                                                    model = michObsPB_XGBoost)
    
  michObsPB_XGBoost_n = nrow(michObsPB_train)
    
  set.seed(512)
  sagPB_XGBoost = xgboost::xgboost(data = sagPB_train[, 2:q4],
                                   label = sagPB_train[, 1],
                                   nrounds = 5000,
                                   max_depth = 10,
                                   eta = 0.8,
                                   objective = "binary:logistic",
                                   verbose = 0)
    
  sagPB_varImportance = xgboost::xgb.importance(colnames(sagPB_XGBoost),
                                                model = sagPB_XGBoost)
    
  sagPB_XGBoost_n = nrow(sagPB_train)
    
    
  set.seed(512)
  sagPA_XGBoost = xgboost::xgboost(data = sagPA_train[, 2:q5],
                                   label = sagPA_train[, 1],
                                   nrounds = 5000,
                                   max_depth = 10,
                                   eta = 0.8,
                                   objective = "binary:logistic",
                                   verbose = 0)
    
  sagPA_varImportance = xgboost::xgb.importance(colnames(sagPA_XGBoost),
                                                model = sagPA_XGBoost)
    
  sagPA_XGBoost_n = nrow(sagPA_train)
  
  list(michPB_XGBoost = michPB_XGBoost, michPB_varImportance = michPB_varImportance,
       michPB_XGBoost_n = michPB_XGBoost_n, michSpecPB_XGBoost = michSpecPB_XGBoost,
       michSpecPB_varImportance = michSpecPB_varImportance,
       michSpecPB_XGBoost_n = michSpecPB_XGBoost_n, michObsPB_XGBoost = michObsPB_XGBoost,
       michObsPB_varImportance = michObsPB_varImportance,
       michObsPB_XGBoost_n = michObsPB_XGBoost_n, sagPB_XGBoost = sagPB_XGBoost,
       sagPB_varImportance = sagPB_varImportance, sagPB_XGBoost_n = sagPB_XGBoost_n,
       sagPA_XGBoost = sagPA_XGBoost, sagPA_varImportance = sagPA_varImportance,
       sagPA_XGBoost_n = sagPA_XGBoost_n)
}

topModels <- buildTopModels()

# See importance of each variable and number of observations (as one data frame)
# with all land cover classes merged into one
topModelsVarImportance <- topModels[c(2, 5, 8, 11, 14)] %>%
  bind_rows(.id = "id") %>%
  mutate(Feature = case_when(Feature == "boat" ~ "Distance from nearest boat launch",
                             Feature == "elev" ~ "Elevation",
                             str_detect(Feature, "land") ~ "Land cover (all classes)",
                             Feature == "nit" ~ "Nitrogen",
                             Feature == "phos" ~ "Phosphorus application",
                             Feature == "pop" ~ "Population density",
                             Feature == "solar" ~ "Solar radiation",
                             Feature == "wc10" ~ "WorldClim BIO10",
                             Feature == "wc18" ~ "WorldClim BIO18",
                             Feature == "wind" ~ "Wind speed"),
         Feature = str_wrap(Feature, width = 20)) %>%
  group_by(id, Feature) %>%
  mutate(Gain2 = sum(Gain)) %>%
  ungroup() %>%
  distinct(id, Feature, .keep_all = TRUE)

# Check to see only Land cover was affected by summarizing:
#topModelsVarImportance %>% dplyr::filter(Gain2 != Gain) %>% count(id, Feature)

numberObservations <- topModels[c(3, 6, 9, 12, 15)] %>%
  bind_rows() %>%
  pivot_longer(cols = everything(),
               names_to = "id", values_to = "n") %>%
  mutate(id = str_remove(id, "_n"),
         scale = case_when(str_detect(id, "sag") ~ "Saginaw Bay",
                           str_detect(id, "mich") ~ "Michigan"),
         data = case_when(str_detect(id, "michPB") ~ "All data",
                          str_detect(id, "Spec") ~ "Specimens",
                          str_detect(id, "Obs") ~ "Observations",
                          str_detect(id, "sagPB") ~ "Presence-background",
                          str_detect(id, "sagPA") ~ "Presence-absence"),
         n = str_pad(n, width = 4, side = "right", pad = " ")) %>%
  dplyr::select(-id)

topModelsForPlot <- topModelsVarImportance %>%
  dplyr::filter(Gain2 >= 0.05) %>%
  group_by(id) %>% 
  mutate(Feature = fct_reorder(Feature, Gain2)) %>%
  bind_rows(topModelsVarImportance %>%
              mutate(extra = case_when(Gain2 < 0.05 ~ Gain2,
                                       Gain2 >= 0.05 ~ 0)) %>%
              group_by(id) %>%
              summarize(Gain2 = sum(extra)) %>%
              mutate(Feature = "All others")) %>%
  mutate(scale = case_when(str_detect(id, "sag") ~ "Saginaw Bay",
                           str_detect(id, "mich") ~ "Michigan"),
         data = case_when(str_detect(id, "michPB") ~ "All data",
                          str_detect(id, "Spec") ~ "Specimens",
                          str_detect(id, "Obs") ~ "Observations",
                          str_detect(id, "sagPB") ~ "Presence-background",
                          str_detect(id, "sagPA") ~ "Presence-absence"),
         Gain2 = Gain2*100,
         order = case_when(Feature == "All others" ~ 2,
                           TRUE ~ 1)) %>%
  left_join(numberObservations, by = c("scale", "data")) %>%
  ungroup()

# Annotation layer for plot
annot <- data.frame(scale = c("Michigan", "Saginaw Bay"),
                    x = c(2, 1.5),
                    y = c(-0.055, -0.055),
                    label = c("Michigan data source", "Saginaw Bay data type"))

# Visualize relative importance of each variable in the final models
# (any variables with <5% contribution is collapsed into one group)
# (Figure 6 in manuscript)
plot6 <- ggplot(topModelsForPlot, aes(x = data, y = Gain2)) +
  ggpattern::geom_bar_pattern(aes(fill = reorder(Feature, order), group = id,
                                  pattern = reorder(Feature, order)),
                              position = "fill", stat = "identity",
                              pattern_size = 0.4, pattern_color = "#FFFFFF",
                              pattern_density = 0.05, pattern_fill = "#FFFFFF",
                              pattern_spacing = 0.02) +
  geom_text(aes(y = 0.015, label = str_c("n = ", as.character(n))),
             size = 8/.pt, color = "#000000", hjust = "left") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  labs(fill = "Explanatory variable", pattern = "Explanatory variable",
       y = "Contribution to final model", x = NULL) +
  coord_flip(ylim = c(0, 1), clip = "off") +
  scale_fill_manual(values = c("#3A7D1C", "#E3A800", "#3A7D1C", "#E3A800",
                               "#37406B", "#96BBDA", "#37406B", "#96BBDA",
                               "#FAC7C8", "#FAC7C8",
                               "#000000")) +
  scale_pattern_manual(values = c("none", "none", "stripe", "stripe", 
                                  "none", "none", "stripe", "stripe",
                                  "none", "stripe",
                                  "none")) +
  facet_wrap(~ scale, nrow = 2, scales = "free") +
  force_panelsizes(rows = c(1, 2/3)) +
  geom_text(data = annot, aes(x = x, y = y, label = label),
            angle = 90, size = 10/.pt) +
  theme_classic() +
  theme(panel.border = element_rect(color = "#000000", fill = NA,
                                    linewidth = 1),
        panel.grid.major.y = element_line(color = "#CDCBCF"),
        panel.spacing = unit(0.49, "cm"),
        strip.background = element_rect(color = "#000000", fill = NA,
                                        linewidth = 1),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(margin = ggplot2::margin(t = 0.25, r = 0, b = 0, l = 0, unit = "cm")),
        axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 0.25, b = 0, l = 0, unit = "cm")),
        axis.text = element_text(size = 8),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.key.height = unit(0.25, "in"),
        legend.key.width = unit(0.25, "in"),
        legend.position = "bottom",
        legend.justification = 0.4,
        legend.spacing.y = unit(0.15, "cm"),
        plot.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.8, unit = "cm")) +
  guides(fill = guide_legend(byrow = TRUE, nrow = 3),
         pattern = guide_legend(byrow = TRUE))

plot6

# If the plot viewer is too small error will be:
# Error in seq.default(from, to, by) : invalid '(to - from)/by

#ggsave("figures/plot6.pdf", width = 7.25, height = 9, units = "in")
#ggsave("figures/plot6.png", width = 7.25, height = 9, units = "in", dpi = 600)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Generate predictive maps from each model

# Read in necessary data (created in Step 4)
# with land cover as characters, not integers

landcoverClasses <- FedData::pal_nlcd() %>% dplyr::select(ID, Class)

# Michigan
fullMichScaled <- read.delim("output/fullMichScaled.txt", header = TRUE,
                              sep = "\t", quote = "", fileEncoding = "UTF-8") %>%
    left_join(landcoverClasses, by = c("land" = "ID")) %>%
    dplyr::select(-land) %>%
    rename(land = Class)

data("state_boundaries_wgs84")
baseMich <- subset(state_boundaries_wgs84, NAME == "Michigan" & TYPE == "Land",
                   select = Shape) %>% 
  terra::vect()
crs(baseMich) <- "+proj=longlat"
  
michBoundingBox <- terra::ext(c(-85, -82.2, 41.5, 46.5)) %>% terra::vect()
crs(michBoundingBox) <- crs(baseMich)

# Saginaw Bay
fullSagScaled <- read.delim("output/fullSagScaled.txt", header = TRUE,
                             sep = "\t", quote = "", fileEncoding = "UTF-8") %>%
  left_join(landcoverClasses, by = c("land" = "ID")) %>%
  dplyr::select(-land) %>%
  rename(land = Class)

sagBoundingBox <- terra::ext(c(-84, -83.35, 43.55, 44.05)) %>% terra::vect()
crs(sagBoundingBox) <- crs(baseMich)

# Mapping functin generates predicted probability of suitable habitat
predMap <- function(i, model, baseMap, boundingBox, title) {
  
  if("land" %in% colnames(i)) {
  
  i2 <- model.matrix(~., i %>% dplyr::select(-c(x, y, id))) %>%
    as.data.frame() %>%
    dplyr::select(all_of(model$feature_names)) %>%
    drop_na() %>% as.matrix() } else {
      i2 <- i %>% dplyr::select(-c(x, y, id)) %>%
        drop_na() %>% as.matrix()
    }
  
  predXY <- stats::predict(object = model,
                           newdata = i2,
                           type = "response") %>%
    bind_cols(i %>% drop_na() %>% dplyr::select(c(x, y))) %>%
    relocate(1, .after = y) %>%
    rename("pred" = 3) %>%
    terra::rast(type = "xyz", crs = "+proj=longlat +ellps=WGS84 +datum=WGS84")
  
  if(str_detect(deparse(substitute(boundingBox)), "mich")) {
    a = -85.0
    b = -82.5
    d = 42.0
    e = 46.0 }
  
  else if(str_detect(deparse(substitute(boundingBox)), "sag")) {
    a = -83.9
    b = -83.4
    d = 43.6
    e = 44.0 }
  
  ggplot() +
    geom_sf(data = baseMap %>% terra::crop(boundingBox), fill = "#FFFFFF") +
    geom_sf(data = boundingBox, fill = NA, linewidth = 2/.pt) +
    geom_spatraster(data = predXY) +
    scale_x_continuous(breaks = c(a, b), expand = c(0, 0)) +
    scale_y_continuous(breaks = c(d, e), expand = c(0, 0)) +
    scale_fill_gradient(low = "#DAE7F2", high = "#293051", na.value = NA) +
    labs(fill = "Probability of suitable habitat", 
         title = title) +
    theme_classic() +
    theme(axis.text = element_text(size = 8),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 8.5),
          legend.position = "bottom",
          legend.title = element_text(size = 8.5),
          legend.text = element_text(size = 8))
}

# Scale = Michigan, source = specimens + observations, type = presence-background
michPBmap <- predMap(fullMichScaled, topModels[["michPB_XGBoost"]], 
                     baseMich, michBoundingBox,
                     "All data")

# Scale = Michigan, source = specimens, type = presence-background
michSpecPBmap <- predMap(fullMichScaled, topModels[["michSpecPB_XGBoost"]],
                         baseMich, michBoundingBox,
                         "Specimens only")

# Scale = Michigan, source = observations, type = presence-background
michObsPBmap <- predMap(fullMichScaled, topModels[["michObsPB_XGBoost"]],
                        baseMich, michBoundingBox,
                        "Observations only")

# Combine Michigan predictive maps (Appendix S12 in manuscript)
michPredMaps <- ggpubr::ggarrange(michPBmap, michSpecPBmap, michObsPBmap,
                                  nrow = 1, ncol = 3, common.legend = TRUE,
                                  legend = "bottom")
michPredMaps
#ggsave("figures/michPredMaps.pdf", width = 7.25, height = 6, units = "in")
#ggsave("figures/michPredMaps.png", width = 7.25, height = 6, units = "in", dpi = 600)


# Scale = Saginaw Bay, type = presence-background
sagPBmap <- predMap(fullSagScaled, topModels[["sagPB_XGBoost"]],
                    baseMich, sagBoundingBox,
                    "Presence-background")

# Scale = Saginaw Bay, type = presence-absence
sagPAmap <- predMap(fullSagScaled, topModels[["sagPA_XGBoost"]],
                    baseMich, sagBoundingBox,
                    "Presence-absence")

# Combine Saginaw Bay predictive maps (Appendix S14 in manuscript)
sagPredMaps <- ggpubr::ggarrange(sagPBmap, sagPAmap, nrow = 1, ncol = 2, 
                                 common.legend = TRUE,  legend = "bottom")

sagPredMaps
#ggsave("figures/sagPredMaps.pdf", width = 7.25, height = 4.5, units = "in")
#ggsave("figures/sagPredMaps.png", width = 7.25, height = 4.5, units = "in", dpi = 600)

# Compare the predicted outcomes within each scale, across data source or type

# Function generates values only
predValues <- function(i, model) {
  
  if("land" %in% colnames(i)) {
    
    i2 <- model.matrix(~., i %>% dplyr::select(-c(x, y, id))) %>%
      as.data.frame() %>%
      dplyr::select(all_of(model$feature_names)) %>%
      drop_na() %>% as.matrix() } else {
        i2 <- i %>% dplyr::select(-c(x, y, id)) %>%
          drop_na() %>% as.matrix()
      }
  
  stats::predict(object = model,
                 newdata = i2,
                 type = "response") %>%
    bind_cols(i %>% drop_na() %>% dplyr::select(c(x, y))) %>%
    relocate(1, .after = y) %>%
    rename("pred" = 3)
}

# Predicted values
michPB_predValues <- predValues(fullMichScaled, topModels[["michPB_XGBoost"]])
michSpecPB_predValues <- predValues(fullMichScaled, topModels[["michSpecPB_XGBoost"]])
michObsPB_predValues <- predValues(fullMichScaled, topModels[["michObsPB_XGBoost"]])
sagPB_predValues <- predValues(fullSagScaled, topModels[["sagPB_XGBoost"]])
sagPA_predValues <- predValues(fullSagScaled, topModels[["sagPA_XGBoost"]])

# Function maps differences in prediction values
comparisonMap <- function(pred1, pred2, i, baseMap, boundingBox, title) {
  diffXY <- data.frame(diff = pred1$pred - pred2$pred) %>%
    bind_cols(i %>% drop_na() %>% dplyr::select(c(x, y))) %>%
    relocate(1, .after = y) %>%
    rename("pred" = 3) %>%
    terra::rast(type = "xyz", crs = "+proj=longlat +ellps=WGS84 +datum=WGS84")
  
  if(str_detect(deparse(substitute(boundingBox)), "mich")) {
    a = -85
    b = -82.5
    d = 42
    e = 46 }
  
  else if(str_detect(deparse(substitute(boundingBox)), "sag")) {
    a = -83.9
    b = -83.4
    d = 43.6
    e = 44.0 }
  
  ggplot() +
    geom_sf(data = baseMap %>% terra::crop(boundingBox), fill = "#FFFFFF") +
    geom_sf(data = boundingBox, fill = NA, linewidth = 2/.pt) +
    geom_spatraster(data = diffXY) +
    scale_x_continuous(breaks = c(a, b), expand = c(0, 0)) +
    scale_y_continuous(breaks = c(d, e), expand = c(0, 0)) +
    scale_fill_gradient2(low = "#f6999c", mid = "#FFFFFF",
                         high = "#28560b", na.value = NA) +
    labs(fill = "Difference in predicted probability of suitable habitat", 
         title = title) +
    theme_classic() +
    theme(axis.text = element_text(size = 8),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 8.5),
          legend.position = "bottom",
          legend.title = element_text(size = 8.5),
          legend.text = element_text(size = 8))
}

# Michigan comparisons (based on data source)
# Mean difference all - observation
mean(michPB_predValues$pred - michObsPB_predValues$pred)
# -0.0237 -> on average, the addition of specimens lowers probabilities by 2%

# Mean difference all - specimen
mean(michPB_predValues$pred - michSpecPB_predValues$pred)
# -0.096 -> on average, the addition of observations lowers probabilities by 10%

# Mean difference observation - specimen
mean(michObsPB_predValues$pred - michSpecPB_predValues$pred)
# -0.120 -> on average, observations predict 12% lower probabilities than specimens

# Specimens tend to predict higher probability, which we can see in the means and maps

michObstoSpec <- comparisonMap(michObsPB_predValues, michSpecPB_predValues, 
                               fullMichScaled, baseMich, michBoundingBox, 
                               "Observations only - Specimens only")

michAlltoObs <- comparisonMap(michPB_predValues, michObsPB_predValues,
                               fullMichScaled, baseMich, michBoundingBox,
                               "All data - Observations only")

michAlltoSpec <- comparisonMap(michPB_predValues, michSpecPB_predValues,
                                fullMichScaled, baseMich, michBoundingBox,
                                "All data - Specimens only")

# Combine Michigan comparison maps (Appendix S13 in manuscript)
michCompMaps <- ggpubr::ggarrange(michObstoSpec, michAlltoObs, michAlltoSpec,
                                  nrow = 1, ncol = 3, common.legend = TRUE,
                                  legend = "bottom")
michCompMaps
#ggsave("figures/michCompMaps.pdf", width = 7.25, height = 6, units = "in")
#ggsave("figures/michCompMaps.png", width = 7.25, height = 6, units = "in", dpi = 600)

# Saginaw Bay comparisons (based on data type)
# Mean difference presence-background - presence-absence
mean(sagPB_predValues$pred - sagPA_predValues$pred)
# -0.252 -> on average, presence-background predicts 25% lower probability than presence-absence

# Saginaw Bay comparison map (Appendix S15 in manuscript)
sagCompMap <- comparisonMap(sagPB_predValues, sagPA_predValues, 
                           fullSagScaled, baseMich, sagBoundingBox, 
                           "Presence-background - presence-absence")
sagCompMap
#ggsave("figures/sagCompMap.pdf", width = 7.25, height = 6, units = "in")
#ggsave("figures/sagCompMap.png", width = 7.25, height = 6, units = "in", dpi = 600)

################################################################################

rm(list = ls())
