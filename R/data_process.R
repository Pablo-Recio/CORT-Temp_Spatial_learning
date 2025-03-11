####################################
# 1_data_process
####################################
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes, forcats)
source(here("R", "func.R"))
#
# A) PROCESSING THE LEARNING DF 
#### A.1) Modify the learning df
learning_df <- read.csv(here("./data/Spatial_learn.csv")) %>%
  group_by(lizard_id) %>%
    filter(sum(is.na(choice)) <= 15) %>%
  ungroup()  %>%
    mutate(temp = gsub("[AB]_", "", trt),
          cort = gsub("_[2][38]", "", trt))  %>%
    mutate(temp = factor(temp,
      levels = c("23", "28"),
      labels = c("23" = "Cold", "28" = "Hot"))) %>%
    mutate(cort = factor(cort,
      levels = c("A", "B"),
      labels = c("A" = "Control", "B" = "CORT"))) %>%
  ungroup() %>%
data.frame()
write.csv(learning_df, here("./output/databases_clean/learning_df.csv"))
#
#### A.2) Extract the slope of each individual
#
## Fit the model only if it has not been fit yet (if refit = TRUE)
if(refit){
  model_errors <- brm(errors ~ day + (1 + day|lizard_id),
                data = learning_df,
                negbinomial(link = "log"),
                chains = 4, cores = 4, iter = 8000, warmup = 2000,
                control = list(adapt_delta = 0.99, max_treedepth = 12))
  # Write the model to rds files
  saveRDS(model_errors, file = paste0(here("output/models/model_errors.rds")))
  } else {
    # Read the model from a file
    model_errors <- readRDS(file = paste0(here("output/models/model_errors.rds")))
} 
#
#### A.3) Modify posteriors df to get the mean slope and intercept of each individual (the latter would be useful for the figures)
posteriors_errors <- as_draws_df(model_errors) %>%
  dplyr::select(starts_with("r_")) %>%
  rename_with(~ str_replace(.x, "r_lizard_id\\[(\\d+),\\s*(Intercept|day)\\]", "\\1_\\2")) %>%
  summarise(across(everything(), list(
    mean = ~ mean(.x, na.rm = TRUE),
    sd = ~ sd(.x, na.rm = TRUE)
    ))) %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
  mutate(
    lizard_id = gsub("_.*", "", variable),
    type = gsub("^\\d+_", "", variable)
  ) %>%
  distinct(lizard_id, type, value) %>%
  pivot_wider(
    names_from = type,  # Create columns for each type
    values_from = value  # Fill columns with the values
  )
#
#### A.5) Merging initial df with posteriors to get the a second learning df with individual learning slopes and treatments
mod_learning <- learning_df %>%
  group_by(lizard_id) %>%
    filter(day == 1) %>%
  dplyr::select(lizard_id, clutch, age, temp, cort) %>% # We modify here the spal df to merge it with the estiamates per invidual
data.frame()
learning_slopes <- merge(mod_learning, posteriors_errors, by = "lizard_id")
#
#
# B) PROCESSING THE PHYSIOLOGY DFs
#### B.1) Modify the original dfs
mit_df <- read.csv(here("./data/mit_data.csv"))
euth <- read.csv(here("./data/Euthanasia.csv"))
mit_damage <- read.csv(here("./data/mit_damage.csv")) %>%
  filter(startsWith(X, "Specimen_")) %>%
  mutate(tube.dam = sapply(strsplit(X, "_"), function(x) x[4])) %>%
data.frame()
#
#### B.2) Merge mito_df and mit_damage and modify the resultant df to make it clearer
physio_df <- merge(mit_df, euth, by = c("plate", "tube.no")) %>%
  filter(type != "control", ch == "pe" | ch == "percp.cy5.5" | ch == "apc") %>%
  mutate(ch = recode(ch, "pe" = "mit_potential", "percp.cy5.5" = "ROS", "apc" = "mit_density")) %>%
  dplyr::select(lizard_id, sex, ch, geo.mean) %>%
  group_by(lizard_id, ch, sex) %>%
  summarize(geo.mean = mean(geo.mean, na.rm = TRUE),
            .groups = 'drop') %>%
  pivot_wider(
    names_from = ch,
    values_from = geo.mean
  ) %>%
data.frame()
#
#### B.3) Modify slightly Euthanasia.csv and merge with the damage df
euth_mod <- euth %>%
  filter(type != "control", plate.dam == 1)
damage_df <- merge(mit_damage, euth_mod, by = "tube.dam") %>%
  dplyr::select(lizard_id, mean_DNAdamage, mean_peroxidation) %>%
  group_by(lizard_id) %>%
  summarize(DNAdamage = mean(mean_DNAdamage, na.rm = TRUE),
            peroxidation = mean(mean_peroxidation, na.rm = TRUE),
            .groups = 'drop') %>%
data.frame()
#
#### B.4) Get the mito_df
final_mito_df <- merge(physio_df, damage_df, by = "lizard_id",
                  all.x = TRUE, all.y = FALSE)
write.csv(final_mito_df, "./output/Checking/final_mit.csv")
#
# C) MERGING THE DF TO GET THE FINAL ONE
clean_df <- merge(learning_slopes, final_mito_df, by = "lizard_id") %>%
  rename("int_mean" = "Intercept_mean", "int_sd" = "Intercept_sd",
        "slope_mean" = "day_mean", "slope_sd" = "day_sd") %>%
  mutate(DNAdage = log(DNAdamage),
        peroxidation = log(peroxidation)) %>%
  dplyr::select(lizard_id, clutch, age, sex, temp, cort,
                int_mean, int_sd, slope_mean, slope_sd,
                mit_density, mit_potential, ROS, DNAdamage, peroxidation) %>%
  mutate(across(where(~ is.double(.) & !is.integer(.)) & !any_of("slope_sd"),
                ~ (. - mean(., na.rm = TRUE)) / (2 * sd(., na.rm = TRUE)))) %>%
data.frame()
write.csv(clean_df, file = "./output/databases_clean/clean_df.csv")
