####################################
# 1_data_process
####################################

# Packages
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, forcats)

# Load data
data  <-  read.csv("./data/Learning.csv")

# Remove individuals who did not participate (more than 15 NAs), remove trials [36-40] (Only in associative) and split treatment into Temp and Cort
data_asso <- data %>%
  group_by(lizard_id) %>%
    filter(sum(is.na(FC_associative)) <= 15) %>%
    filter(trial_associative <= 35) %>%
  ungroup()  %>%
    mutate(group = factor(group,
     levels = c("R_B", "B_R"),
     labels=c("R_B"="Red", "B_R"="Blue"))) %>%
    mutate(temp = gsub("[AB]_", "", trt),
          cort = gsub("_[2][38]", "", trt))  %>%
    mutate(temp = factor(temp,
      levels = c("23", "28"),
      labels = c("23" = "Cold", "28" = "Hot"))) %>%
    mutate(cort = factor(cort,
      levels = c("B", "A"),
      labels = c("B" = "CORT", "A" = "Control"))) %>%
    mutate( # Standarize data by trial (i.e. make the first trial where each individual participated their trial 1)
      first_non_na = min(which(!is.na(FC_associative))),
      Associative_Trial = ifelse(!is.na(first_non_na),trial_associative - first_non_na + 1, trial_associative))%>% 
      filter(Associative_Trial >= 1) %>%
  ungroup() %>% 
data.frame()


data_rev <- data %>%
  group_by(lizard_id) %>%
    filter(sum(is.na(FC_associative)) <= 15) %>%
    filter(sum(is.na(FC_reversal)) <= 15) %>%
   ungroup()  %>%
    mutate(group = factor(group,
     levels = c("R_B", "B_R"),
     labels=c("R_B"="Red", "B_R"="Blue"))) %>%
    mutate(temp = gsub("[AB]_", "", trt),
          cort = gsub("_[2][38]", "", trt))  %>%
    mutate(temp = factor(temp,
      levels = c("23", "28"),
      labels = c("23" = "Cold", "28" = "Hot"))) %>%
    mutate(cort = factor(cort,
      levels = c("B", "A"),
      labels = c("B" = "CORT", "A" = "Control"))) %>%
    mutate(trial_reversal=as.numeric(trial_reversal)) %>%
  data.frame()


