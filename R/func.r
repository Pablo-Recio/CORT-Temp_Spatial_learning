#################### 
####################
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes, ggh4x)
#
####################
####################
# Function to make the first cleaning for the initial posteriors to estimate individual learning slopes or intercepts
#' @title tidy_post_df
#' @param df to select the posteriors to use. Options: choice or errors
#' @param effect to select whether we want the intercept or the slopes. Options: slope or intercept
tidy_post_df <- function(df, effect) {
  if(df == "choice") {
    data <- posteriors_choice
    label_df <- "choice_"
  } else if(df == "errors") {
    data <- posteriors_errors
    label_df <- "errors_"
  } else {
    stop("df not valid")
  }
  if(effect == "slope") {
    label_effect <- "day"
    label_df <- paste0(label_df, "slope_")
    pattern <- "^r_lizard_id\\[\\d+,day\\]$"
  } else if(effect == "intercept") {
    label_effect <- "Intercept"
    label_df <- paste0(label_df, "intercept_")
    pattern <- "^r_lizard_id\\.\\d+\\.Intercept\\."
    data <- data %>%
      select(b_Intercept, matches("^r_lizard.*Intercept]$")) %>%
      mutate(across(-b_Intercept, ~ . + b_Intercept)) %>%
    data.frame()
  } else {
    stop("effect not valid")
  }
  tidy_post_1 <- data %>%
    select(matches(pattern)) %>%
    summarise(across(everything(), list(
      mean = ~ mean(.x, na.rm = TRUE),
      sd = ~ sd(.x, na.rm = TRUE),
      se = ~ sd(.x, na.rm = TRUE)/sqrt(length(.x))
    ))) %>%
  data.frame()
  tidy_post_2 <- tidy_post_1 %>%
    pivot_longer(cols = everything(),
              names_to = "lizard_id", 
              values_to = "value") %>% # Extract the relevant columns and reshape them
    mutate(lizard_id = gsub("r_lizard_id\\.", "", lizard_id),  # Remove prefix
           lizard_id = gsub(paste0(label_effect, "\\."), "", lizard_id)) %>%      # Replace all remaining dots with underscores
    separate(lizard_id, into = c("lizard_id", "stat"),
            sep = "\\.", extra = "merge", fill = "right") %>% # Separate the parameter names to get lizard_id, statistic type, and effect type
    pivot_wider(names_from = stat, values_from = value, names_prefix = label_df) %>% # Split the colum stats into mean, sd, and se
    data.frame()
  return(tidy_post_2)
}
####################
####################
# Extract sample size per group
#' @title sample
#' @description Estract sample size per group
#' @param df To select the df
#' @param sp To select the species of interest ("deli"/"guich")
#' @param corti To select cort treatment ("CORT"/"Control")
#' @param therm To select temp treatment ("Cold"/"Hot")
sample <- function(df, sp, corti, therm){
  sample_size <- df %>%
                filter(species == sp, cort == corti, temp == therm, trial_reversal == 1) %>%
                group_by(lizard_id) %>%
                summarise(n = n()) %>%
                summarise(total_count = sum(n)) %>%
                pull(total_count)
  return(sample_size)
}
####################
####################
# Fit models and extract posteriors both per each treatment:
#' @title fit_m Function
#' @description Fit brm models for different the associative task
#' @param df To select the df
#' @param sp To select the species of interest ("deli"/"guich")
#' @param com To select whether the analyses are done for the full dataset ("complete") or for only those lizards 
# that met the learning criterion ("suppl")
#' @param refit To choose whether to refit the models (TRUE, default) or use the ones already made (FALSE)
#' @return Raw posteriors of fitted brm model for each treatment, species, and group (df)
fit_m <- function(df, sp, com, refit = TRUE) {
  formula <- (FC_reversal ~ trial_reversal*cort*temp + (1 + trial_reversal|lizard_id) + (1|clutch)) 
  #Specify species
    if (sp == "deli"){
      sp_data <- df %>%
            group_by(lizard_id) %>%
            filter(species == "delicata") %>%
            ungroup() %>%
      data.frame() 
    } else {
      if(sp == "guich"){
        sp_data <- df %>%
              group_by(lizard_id) %>%
              filter(species == "guichenoti") %>%
              ungroup() %>%
        data.frame()
      } else {
        stop("Species not valid")
      }
    }
  #Fit the model only if it has not been fit yet (if refit=TRUE)
  if(refit){
    # Fit the model
    model <- brm(formula,
                data = sp_data,
                family = bernoulli(link = "logit"),
                chains = 4, cores = 4, iter = 3000, warmup = 1000, control = list(adapt_delta = 0.99))
    # Write the model to a file
    saveRDS(model, file = paste0(here("output/models/"), sp, com, ".rds"))
  } else {
      # Read the model from a file
      model <- readRDS(file = paste0(here("output/models/"), sp, com, ".rds"))
  } 
  # Extract posteriors
  posteriors <- as_draws_df(model)
  return(posteriors)
}
###################
###################
# Tidy estimates from fit_m
#' @title tidy_post
#' @description The df obtained by fit_m contains the estimates of the Intercept and fixed effects of "Choice ~ trial (1 + lizard)"
#' for each treatment, species, and group,but it extracts the estimates of the four chains, each chain in a column, here we want to
#' tidy up to two columns (Intercept and trial), plus the ones indicating species, group, tand treatment
#' @param df Dataframe used
#' @return Same df but tidy
tidy_post <- function(df) {
  # Select data
  data <- df
  # Split original df into four by chain, and give the columns a common name
  res1 <- data%>%select(matches("^X1.b_"),treatment_level, species_level,bias_level)%>%rename_all(~sub("^X1.b_", "", .))
  res2 <- data%>%select(matches("^X2.b_"),treatment_level, species_level,bias_level)%>%rename_all(~sub("^X2.b_", "", .))
  res3 <- data%>%select(matches("^X3.b_"),treatment_level, species_level,bias_level)%>%rename_all(~sub("^X3.b_", "", .))
  res4 <- data%>%select(matches("^X4.b_"),treatment_level, species_level,bias_level)%>%rename_all(~sub("^X4.b_", "", .))
  # Bind data again
  new_df <- bind_rows(res1, res2,res3,res4)
return(new_df) 
}
####################
####################
# Estimate pmcm
#' @title pMCMC Function
#' @param x The vector for the posterior distribution. Note that this will test the null hypothesis that the parameter of interest is significantly different from 0. 
#' @param null A numeric value decsribing what the null hypothesis should be
#' @param twotail Whether to conduct a one-tailed hypothesis or a two-tailed hypotheses. Default = true indicating a two-tailed test will be done.
#' @param dir The direction of the one-tail test (<) or (>)
pmcmc <- function(x, null = 0, twotail = TRUE, dir){
  if(twotail){
    2*(1 - max(table(x<=null) / length(x)))
  } else{
    if(dir == "<"){
    (max(sum(x>=null) / length(x)))
    } else{
      if(dir == ">"){
        (max(sum(x<=null) / length(x)))
      } else{
        stop("dir not valid")
      }
    }
  }
}
####################
####################
# Function to format numbers with n decimal places
#' @title format_dec
#' @param x The object
#' @param n The number of decimals
format_dec <- function(x, n) {
  z <- sprintf(paste0("%.",n,"f"), x)
  return(z)
}
####################
####################
# Function to format p_values with n decimal places
#' @title format_p
#' @param x The object
#' @param n The number of decimals
format_p <- function(x, n) {
  z <- sprintf(paste0("%.",n,"f"), x)
  tmp <- ifelse(as.numeric(z) <= 0.001, "< 0.001",
         ifelse(as.numeric(z) <= 0.05 & as.numeric(z) > 0.001, "< 0.05",
                paste0("= ", as.character(z))))
  return(tmp)
}
####################
####################
# Function to create the plot for the slopes
#' @title plot_slopes
#' @param type to select the df for the violin plot
plot_slopes <- function(type){
  if(type == "choice") {
    data <- posteriors_choice
    label_df <- "choice_"
  } else if(type == "errors") {
    data <- posteriors_errors
    label_df <- "errors_"
  } else {
    stop("df not valid")
  }
  label_effect <- "day"
  pattern <- "^r_lizard_id\\[\\d+,day\\]$"
  tidy_post <- data %>%
    select(matches(pattern)) %>%
    pivot_longer(cols = everything(),
              names_to = "lizard_id", 
              values_to = "value") %>% # Extract the relevant columns and reshape them
    mutate(lizard_id = gsub("r_lizard_id\\[|,day\\]", "", lizard_id)) %>%  # Remove prefix
  data.frame()
  mod_spal <- data_spal %>%
    group_by(lizard_id) %>%
      filter(day == 1) %>%
    select(lizard_id, clutch, age, temp, cort) %>% # We modify here the spal df to merge it with the estiamates per invidual
    mutate(lizard_id = as.character(lizard_id)) %>%
  data.frame()
  tidy_post_fig <- merge(mod_spal, tidy_post, by = "lizard_id")
###
  data_fig <- tidy_post_fig %>%
    mutate(Treatment = paste(cort, temp, sep = "-")) %>%
    mutate(Treatment = factor(Treatment,
      levels = c("CORT-Cold", "Control-Cold", "CORT-Hot", "Control-Hot"),
      labels = c("CORT-Cold" = "CORT-Cold (n=20)",
                "Control-Cold" = "Control-Cold (n=20)",
                "CORT-Hot" = "CORT-Hot (n=19)",
                "Control-Hot" = "Control-Hot (n=20)")
    )) %>%
  data.frame()
  data_plot <- data_fig %>%
    group_by(Treatment) %>%
    summarize(
      Mean = mean(value),
      SD = sd(value),
      SE = sd(value)/sqrt(length(value))
    ) %>%
    ungroup() %>%
  data.frame()
write.csv(data_plot, here("output/Checking/data_plot.csv"))
write.csv(data_fig, here("output/Checking/data_fig.csv"))
#
# Make the plot
  plot2 <- ggplot(data_fig, aes(x = Treatment, y = value, fill = Treatment)) +
  geom_flat_violin(alpha = 0.5) +
  scale_fill_manual(values = c("CORT-Cold (n=20)"="#00008B", "Control-Cold (n=20)"="#68bde1", 
                    "CORT-Hot (n=19)"="#b50101", "Control-Hot (n=20)"="#fa927d")) +
  geom_point(data = data_plot, aes(y = Mean, x = Treatment), position = position_dodge(width = 0.75), color = "black", fill = "black", size = 3) +
  geom_segment(data = data_plot, aes(y = Mean - SD, yend = Mean + SD, x = Treatment, xend = Treatment), size = 1.5, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  ylim(min(data_fig$value), max(data_fig$value)) +
  coord_flip() +
  theme_classic() +
  labs(y = "Slope estimates", x = "Treatments") +
  theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5, "mm"),
    axis.text.y = element_blank(),   # Remove y-axis labels
    axis.title = element_text(size = 12, family = "Times"),
    axis.text = element_text(size = 10, family = "Times"),
    legend.position = "right",
    legend.title = element_text(size = 12, family = "Times"),
    legend.text = element_text(size = 11, family = "Times")
    )
  return(plot2)
}
####################
####################
# Function to create the plot per species (A-B or C-D) for fig-results
#' @title plotting
#' @param sp to select the species for the labels
#' @param df_prob to select the df for probability 
#' @param df_violin to select the df for the violin plot
#' @param df_points to select the df for the points and geom_bars
plotting <- function(sp, df_prob, df_violin, df_points){
  # Specify labels depending on species and relevel the factor treatment for the legend
    df_violin$Treatment <- factor(df_violin$Treatment,
      levels = c("Control-Hot (n =  20 )", "CORT-Hot (n =  20 )", "Control-Cold (n = 19 )", "CORT-Cold (n = 20 )"))
    img <- readPNG("./Others/Deli.png")
    note <- paste("L. delicata")
  # First part of the plot, the probabilities of choosing right over trial
  plot1 <- ggplot(df_prob, aes(x = Trial, y = Mean_Predicted_prob, color = Treatment)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("CORT-Cold"="#00008B", "Control-Cold"="#68bde1", "CORT-Hot"="#b50101", "Control-Hot"="#fa927d")) +
  geom_ribbon(aes(ymin = Mean_Predicted_prob - SE_Predicted_prob, ymax = Mean_Predicted_prob + SE_Predicted_prob, fill = Treatment), color = NA, alpha = 0.075) + 
  scale_fill_manual(values = c("CORT-Cold"="darkblue", "Control-Cold"="#68bde1", "CORT-Hot"="#b50101", "Control-Hot"="#fa927d")) +
  theme_classic() +
  labs(y = "Predicted probability of correct choice", x = "Trial") +
  theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5, "mm")) +
  theme(
    axis.title = element_text(size = 12, family = "Times"),
    axis.text = element_text(size = 10, family = "Times"),
    legend.position = "none"
  )
  #
  # Second part of the plot: the estimates of each treatment in violin plot
  plot2 <- ggplot(df_violin, aes(x = Treatment, y = Value, fill = Treatment)) +
  geom_flat_violin(alpha = 0.5) +
  scale_fill_manual(values = c("CORT-Cold"="#00008B", "Control-Cold"="#68bde1", "CORT-Hot"="#b50101", "Control-Hot"="#fa927d")) +
  geom_point(data = df_points, aes(y = Mean, x = Treatment), position = position_dodge(width = 0.75), color = "black", fill = "black", size = 3) +
  geom_segment(data = df_points, aes(y = Mean - SD, yend = Mean + SD, x = Treatment, xend = Treatment), size = 1.5, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  ylim(-0.015, max(df_violin$Value)) +
  coord_flip() +
  theme_classic() +
  labs(y = "Slope estimates", x = "Treatments") +
  theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5, "mm"),
    axis.text.y = element_blank(),   # Remove y-axis labels
    axis.title = element_text(size = 12, family = "Times"),
    axis.text = element_text(size = 10, family = "Times"),
    legend.position = "right",
    legend.title = element_text(size = 12, family = "Times"),
    legend.text = element_text(size = 11, family = "Times")
    )
  #
  # Combine them
  plot <- plot_grid(plot1, plot2, labels = lab, nrow = 1, rel_widths = c(0.45, 0.45)) +
    annotation_custom(rasterGrob(img), xmin = 0.73, xmax = 0.98, ymin = 0.73, ymax = 0.98) +
    annotate("text", x = 0.5, y = 0.05, label = note, size = 5, color = "black", fontface = "italic")
  return(plot)
}