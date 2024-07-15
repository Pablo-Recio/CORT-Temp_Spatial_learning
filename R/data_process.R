####################################
# 1_data_process
####################################
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes, forcats)
# A) Processing the learning df
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
    mutate( # Standarize data by trial (i.e. make the first trial where each individual participated their trial 1)
      first_non_na = min(which(!is.na(choice))),
      day = ifelse(!is.na(first_non_na),choice - first_non_na + 1, day))%>% 
      filter(day >= 1) %>%
  ungroup() %>% 
data.frame()
#### A.2) Extract the slope of each individual
  #Fit the models only if they have not been fit yet (if refit=TRUE)
  if(refit){
    # Formula models
    formula <- choice ~ 1 + (1 + day|lizard_id) + (1|clutch)
    # Fit the models
    # Choice
    model_choice <- brm(formula,
                data = data_spal,
                family = bernoulli(link = "logit"),
                chains = 4, cores = 4, iter = 3000, warmup = 1000, control = list(adapt_delta = 0.99))
    # Errors
    model_errors <- brm(formula,
                data = data_spal,
                negbinomial(link = "log"),
                chains = 4, cores = 4, iter = 3000, warmup = 1000, control = list(adapt_delta = 0.99))
    # Write the models to two files
    saveRDS(model_choice, file = paste0(here("output/models/model_choice.rds")))
    saveRDS(model_errors, file = paste0(here("output/models/model_errors.rds")))

  } else {
      # Read the model from a file
      model_choice <- readRDS(file = paste0(here("output/models/model_choice.rds")))
      model_errors <- readRDS(file = paste0(here("output/models/model_errors.rds")))
  } 
  # Extract posteriors
  posteriors_choice <- as_draws_df(model_choice)
  posteriors_errors <- as_draws_df(model_errors)

