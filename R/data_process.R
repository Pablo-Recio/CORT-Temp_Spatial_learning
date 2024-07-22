####################################
# 1_data_process
####################################
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes, forcats)
#
# A) PROCESSING THE LEARNING DF
data <- read.csv(here("./data/Spatial_learn.csv"))
#### A.1) Remove individuals who did not participate (more than 15 NAs), remove trials [36-40] (Only in associative) and split treatment into Temp and Cort
data_spal <- data %>%
  group_by(lizard_id) %>%
    filter(sum(is.na(choice)) <= 15) %>%
  ungroup()  %>%
    mutate(temp = gsub("[AB]_", "", trt),
          cort = gsub("_[2][38]", "", trt))  %>%
    mutate(temp = factor(temp,
      levels = c("23", "28"),
      labels = c("23" = "Cold", "28" = "Hot"))) %>%
    mutate(cort = factor(cort,
      levels = c("B", "A"),
      labels = c("B" = "CORT", "A" = "Control"))) %>%
  ungroup() %>% 
data.frame()
#
#
#### A.2) Extract the slope of each individual
#
## Fit the models only if they have not been fit yet (if refit = TRUE)
if(refit==TRUE){
  # Choice
  model_choice <- brm(choice ~ 1 + (1 + day|lizard_id) + (1|clutch),
                data = data_spal,
                family = bernoulli(link = "logit"),
                chains = 4, cores = 4, iter = 8000, warmup = 2000, 
                control = list(adapt_delta = 0.99, max_treedepth = 12),
                prior = custom_priors)
  # Errors
  model_errors <- brm(errors ~ 1 + (1 + day|lizard_id) + (1|clutch),
                data = data_spal,
                negbinomial(link = "log"),
                chains = 4, cores = 4, iter = 8000, warmup = 2000,
                control = list(adapt_delta = 0.99, max_treedepth = 12))
  # Write the models to two files
  saveRDS(model_choice, file = paste0(here("output/models/model_choice.rds")))
  saveRDS(model_errors, file = paste0(here("output/models/model_errors.rds")))
  } else {
    # Read the model from a file
    model_choice <- readRDS(file = paste0(here("output/models/model_choice.rds")))
    model_errors <- readRDS(file = paste0(here("output/models/model_errors.rds")))
  } 
#
## Extract posteriors
posteriors_choice <- as_draws_df(model_choice)
posteriors_errors <- as_draws_df(model_errors)
#
#
#### A.3) Modify posteriors df to get the mean slope and intercept of each individual (the latter would be useful for the figures)
#
# Getting slopes
source(here("R", "func.R"))
## Choice
post_choice_slope <- tidy_post_df("choice", "slope")
## Errors
post_errors_slope <- tidy_post_df("errors", "slope")
#
#### A.4) Merge posteriors of the slopes for the analyses with other variables
learning_posteriors <- merge(post_choice_slope, post_errors_slope, by = "lizard_id", all = TRUE)
#
#
#### A.5) Merging initial df with posteriors to get the final learning df with individual learning slopes and treatments
  mod_spal <- data_spal %>%
    group_by(lizard_id) %>%
      filter(day == 1) %>%
    select(lizard_id, clutch, age, temp, cort) %>% # We modify here the spal df to merge it with the estiamates per invidual
  data.frame()
learning_df <- merge(mod_spal, learning_posteriors, by = "lizard_id")
#
#
#
# B) 