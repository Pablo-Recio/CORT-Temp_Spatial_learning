####################################
# 1_data_process
####################################
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes, forcats)
# A) We extract the lizards we are using from the main data_base
meta_data <- read.csv("./data/lizard_database.csv") %>%
  filter(lizard_id %in% c("LD3158_23", "LD3177_23", "LD3104_23", "LD3100_23", "LD3027_23", "LD3196_23", "LD3025_23",
"LD2930_23", "LD2953_23", "LD2904_23", "LD2981_23", "LD3078_23", "LD3102_23", "LD2997_23", "LD3090_23", "LD3098_23",
"LD3041_23", "LD2875_23", "LD2950_23", "LD3023_23", "LD2993_23", "LD2958_23", "LD3031_23", "LD3065_23", "LD3038_23",
"LD3198_23", "LD2851_23", "LD3193_23", "LD2952_23", "LD2995_23", "LD3033_23", "LD2928_23", "LD3081_23", "LD2982_23",
"LD3103_23", "LD2924_23", "LD3195_23", "LD3028_23", "LD3176_23", "LD3099_23", "LD3079_23", "LD3159_23", "LD3082_23",
"LD3093_23", "LD2983_23", "LD3039_23", "LD2960_23", "LD3101_23", "LD2926_23", "LD2959_23", "LD3042_23", "LD2908_23",
"LD3086_23", "LD2925_23", "LD3014_23", "LD3051_23", "LD2984_23", "LD2968_23", "LD2967_23", "LD3009_23", "LD2922_23",
"LD3034_23", "LD3095_23", "LD2971_23", "LD3191_23", "LD2873_23", "LD3109_23", "LD3178_23", "LD3188_23", "LD3187_23",
"LD3154_23", "LD3157_23", "LD3179_23", "LD3084_23", "LD2911_23", "LD3201_23")) %>%
data.frame()
write.csv(meta_data, "./data/spaL.csv") 

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


