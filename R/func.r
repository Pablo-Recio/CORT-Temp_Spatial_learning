#################### 
####################
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes, ggh4x)
#
####################
####################
# Extract sample size per group
#' @title sample
#' @description Estract sample size per group
#' @param df To select the df
#' @param corti To select cort treatment ("CORT"/"Control")
#' @param therm To select temp treatment ("Cold"/"Hot")
sample <- function(df, corti, therm){
  sample_size <- df %>%
                filter(cort == corti, temp == therm) %>%
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
#' @param cat To select whether the analyses are preliminary ("prel") or finel ("fin")
#' @param refit To choose whether to refit the models (TRUE, default) or use the ones already made (FALSE)
#' @param fam To specify the family of the model
#' @param var To specify the response variable
#' @param for To specify the formula
#' @return Raw posteriors of fitted brm model for each treatment, species, and group (df)
fit_m <- function(df, cat, fam, var, formula, refit = TRUE) {
  #Fit the model only if it has not been fit yet (if refit=TRUE)
  if(refit){
    # Fit the model
    model <- brm(formula,
                data = df,
                family = fam,
                chains = 4, cores = 4, iter = 8000, warmup = 2000, control = list(adapt_delta = 0.99))
    # Write the model to a file
    saveRDS(model, file = paste0(here("output/models/"), var, "_", cat, ".rds"))
  } else {
      # Read the model from a file
      model <- readRDS(file = paste0(here("output/models/"), var, "_", cat, ".rds"))
  } 
  # Extract posteriors
  posteriors <- as_draws_df(model)
  return(posteriors)
}
####################
####################
# Function to generate four Q-Q plots for one variable
#' @title qq_plots_single 
#' @description Generate four Q-Q plots for one variable
#' @param df Database
#' @param vle Variable chosen
#' @param label Label assigned to each variable
#' @return fig_qq
qq_plots_single <- function(df, vle, label) {
  # Remove NAs
  plot_df <- data.frame(value = na.omit(df[[vle]])) 

    # Normal Q-Q Plot
    p1 <- ggplot(plot_df, aes(sample = value)) +
      stat_qq() + 
      stat_qq_line(color = "red") +
      ggtitle(paste(label, "- Normal Q-Q Plot")) +
      theme_classic()

    # Log-Normal Q-Q Plot
    p2 <- ggplot(plot_df, aes(sample = log(value +1))) +
      stat_qq() + 
      stat_qq_line(color = "red") +
      ggtitle(paste(label, "- Log-Normal Q-Q Plot")) +
      theme_classic()

    # Exponential Q-Q Plot
    p3 <- ggplot(plot_df, aes(sample = value)) +
      stat_qq(distribution = qexp, dparams = list(rate = 1/mean(plot_df$value))) +
      stat_qq_line(distribution = qexp, dparams = list(rate = 1/mean(plot_df$value)), color = "red") +
      ggtitle(paste(label, "- Exponential Q-Q Plot")) +
      theme_classic()

    # Arrange plots in a grid
    fig_qq <- plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(0.33, 0.34, 0.33))
    
  return(fig_qq)
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
#' @param equal To write whether we want the equal or not
format_p <- function(x, n, equal) {
  if (equal == TRUE){
    label <- "="
  } else {
    label <- ""
  }
  z <- sprintf(paste0("%.",n,"f"), x)
  tmp <- ifelse(as.numeric(z) <= 0.001, "< 0.001",
         ifelse(as.numeric(z) <= 0.05 & as.numeric(z) > 0.001, "< 0.05",
                paste0(label, as.character(z))))
  return(tmp)
}
###################
###################
# Function to make the first cleaning for the posteriors and create a big df that includes
# all variables and their predictors
#' @title tidy_post
#' @param df to select the posteriors to use. Options: choice or errors
tidy_post <- function(df) {
  df <- as.data.frame(df)
  summary_df <- df %>%
    dplyr::select(starts_with("b_")) %>%  # Select only parameter estimates
    pivot_longer(everything(), names_to = "Predictor", values_to = "Estimate") %>%
    group_by(Predictor) %>%
    summarise(
      Mean = format_dec(mean(Estimate), 3),
      CI_lower = format_dec(quantile(Estimate, 0.05, na.rm = TRUE), 3),
      CI_upper = format_dec(quantile(Estimate, 0.95, na.rm = TRUE), 3),
      PMCMC = format_dec(pmcmc(Estimate, null = 0, twotail = TRUE), 3),
      .groups = "drop"
    )
  
  return(summary_df)
}
####################
####################
# Function to rename and refine the posterior databases for preparing the tables.
#' @title refine_post
#' @param df to select the databse
refine_post <- function(df) {
  variable_order <- c("Mit density",
                     "Mit potential",
                     "ROS",
                     "DNA damage",
                     "Lipid peroxidation")
  predictor_order <- c("b_Intercept",
                       "b_cortCORT",
                       "b_tempHot",
                       "b_cortCORT:tempHot",
                       "b_age")
  df <- df %>%
    mutate(
      Model = case_when(
        Model == "m_def_mit_density" ~ "Mit density",
        Model == "m_def_mit_potential" ~ "Mit potential",
        Model == "m_def_ROS" ~ "ROS",
        Model == "m_def_DNAdamage" ~ "DNA damage",
        Model == "m_def_peroxidation" ~ "Lipid peroxidation",
        TRUE ~ Model
      )
    ) %>%
    mutate(CI_95 = paste0("[", CI_lower, " , ", CI_upper, "]")) %>%
    rename(
      Variable = Model,
      Predictors = Predictor,
      'Estimate Mean' = Mean, 
      '95% CI' = CI_95, 
      PMCMC = PMCMC
    ) %>%
    dplyr::select(
      Variable,
      Predictors,
      'Estimate Mean',
      '95% CI',
      PMCMC) %>%
    mutate(Variable = factor(Variable, levels = variable_order)) %>%
    mutate(Predictors = factor(Predictors, levels = predictor_order))
  
  return(df)
}
####################
####################
# Function to get estimates to make the contrasts between treatments
#' @title post_values
#' @param df To select the df
#' @param fac To indicate if there are other parameters to take into account
post_values <- function(df, fac){
  df <- as.data.frame(df)
  #Getting the estimated values per each treatment
  if(fac == "none"){
    Control_Cold <- df$b_Intercept
    Control_Hot <- df$b_Intercept + df$b_tempHot
    CORT_Cold <- df$b_Intercept + df$b_cortCORT
    CORT_Hot <- df$b_Intercept + df$b_cortCORT + df$b_tempHot + df$`b_cortCORT:tempHot`
  } else if (fac == "sex"){
    # Males
    Control_Cold_males <- base::sample(df$b_Intercept, size = 12000, replace = FALSE)
    Control_Hot_males <- base::sample(df$b_Intercept + df$b_tempHot, size = 12000, replace = FALSE)
    CORT_Cold_males <- base::sample(df$b_Intercept + df$b_cortCORT, size = 12000, replace = FALSE)
    CORT_Hot_males <- base::sample(df$b_Intercept + df$b_cortCORT + df$b_tempHot +
                                      df$`b_cortCORT:tempHot`, size = 12000, replace = FALSE)
    # Females
    Control_Cold_fem <- base::sample(df$b_Intercept + df$b_sexFemale, size = 12000, replace = FALSE)
    Control_Hot_fem <- base::sample(df$b_Intercept + df$b_tempHot + df$b_sexFemale, size = 12000, replace = FALSE)
    CORT_Cold_fem <- base::sample(df$b_Intercept + df$b_cortCORT + df$b_sexFemale, size = 12000,
                                    replace = FALSE)
    CORT_Hot_fem <- base::sample(df$b_Intercept + df$b_cortCORT + df$b_tempHot +
                                    df$`b_cortCORT:tempHot` + df$b_sexFemale, size = 12000, replace = FALSE)
    # Pool
    CORT_Cold <- c(CORT_Cold_males, CORT_Cold_fem)
    CORT_Hot <- c(CORT_Hot_males, CORT_Hot_fem)
    Control_Cold <- c(Control_Cold_males, Control_Cold_fem)
    Control_Hot <- c(Control_Hot_males, Control_Hot_fem)
  }
  data_values <- data.frame(
    CORT_Cold = CORT_Cold,
    CORT_Hot = CORT_Hot,
    Control_Cold = Control_Cold,
    Control_Hot = Control_Hot
  )
  return(data_values)
}
####################
####################
# Function to create the plots for mit variables
#' @title plotting
#' @param df to select the df (only posteriors)
#' @param lab to select the appropriate label for axis Y
plotting <- function(df, lab){
  if(lab %in% c("DNA damage", "Lipid peroxidation")){
    n_CORTCold <- paste0("CORT-Cold (n=", n_damage$CORT_Cold, ")")
    n_CORTHot <- paste0("CORT-Hot (n=", n_damage$CORT_Hot, ")")
    n_ControlCold <- paste0("Control-Cold (n=", n_damage$Control_Cold, ")")
    n_ControlHot <- paste0("Control-Hot (n=", n_damage$Control_Hot, ")")
  } else {
    n_CORTCold <- paste0("CORT-Cold (n=", n_list$CORT_Cold, ")")
    n_CORTHot <- paste0("CORT-Hot (n=", n_list$CORT_Hot, ")")
    n_ControlCold <- paste0("Control-Cold (n=", n_list$Control_Cold, ")")
    n_ControlHot <- paste0("Control-Hot (n=", n_list$Control_Hot, ")")
  }
  # Build db for violin plots
  db_violin <- data.frame(
    Treatment = rep(c(n_CORTCold, n_CORTHot, n_ControlCold, n_ControlHot),
                      each = nrow(df)),
    Estimate = c(df$CORT_Cold, df$CORT_Hot, df$Control_Cold, df$Control_Hot)
  ) %>%
  mutate(Treatment = factor(Treatment,
                            levels = c(n_ControlHot, n_CORTHot, n_ControlCold, n_CORTCold)))
  # Build db for geom_bars
  db_bars <- db_violin %>%
    group_by(Treatment) %>%
    summarize(
      Mean = mean(Estimate),
      SD = sd(Estimate),
      SE = sd(Estimate)/sqrt(length(Estimate))
    )
  # Making the plot
  plot <- ggplot(db_violin, aes(x = Treatment, y = Estimate, fill = Treatment)) +
    geom_flat_violin(alpha = 0.5) +
    scale_fill_manual(values = set_names(
      c("#00008B", "#b50101", "#68bde1", "#fa927d"),
      unique(db_violin$Treatment))) +
    geom_point(data = db_bars, aes(y = Mean, x = Treatment), position = position_dodge(width = 0.75), color = "black", fill = "black", size = 3) +
    geom_segment(data = db_bars, aes(y = Mean - SD, yend = Mean + SD, x = Treatment, xend = Treatment), size = 1.5, color = "black") +
    ylim(min(db_violin$Estimate), max(db_violin$Estimate)) +
    coord_flip() +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme_classic() +
    labs(y = lab, x = "Treatments") +
    theme(plot.margin = margin(3, 3, 3, 3, "mm"),
      axis.text.y = element_blank(),   # Remove y-axis labels
      axis.title = element_text(size = 12, family = "Times"),
      axis.text = element_text(size = 10, family = "Times"),
     legend.position = "none",
      legend.title = element_text(size = 12, family = "Times"),
      legend.text = element_text(size = 11, family = "Times")
  )
  return(plot)
}
####################
####################
# Function to create the plot for the slopes
#' @title plot_slopes
#' @param df to select the df for the violin plot
plot_slopes <- function(df){
  data_fig <- df %>%
    mutate(treatment = factor(treatment,
      levels = c("CORT-Cold", "Control-Hot", "CORT-Hot", "Control-Cold"),
      labels = c("Control-Hot" = "Control-Hot (n=20)",
                "CORT-Hot" = "CORT-Hot (n=19)",
                "Control-Cold" = "Control-Cold (n=20)",
                "CORT-Cold" = "CORT-Cold (n=20)")
    )) %>%
  data.frame()
  data_plot <- data_fig %>%
    group_by(treatment) %>%
    summarize(
      Mean = mean(slopes),
      SD = sd(slopes),
      SE = sd(slopes)/sqrt(length(slopes))
    ) %>%
    ungroup() %>%
  data.frame()
#
# Make the plot
  plot2 <- ggplot(data_fig, aes(x = treatment, y = slopes, fill = treatment)) +
  geom_flat_violin(alpha = 0.5) +
  scale_fill_manual(values = c("CORT-Cold (n=20)"="#00008B",
                            "Control-Cold (n=20)"="#68bde1",
                            "CORT-Hot (n=19)"="#b50101",
                            "Control-Hot (n=20)"="#fa927d")) +
  geom_point(data = data_plot, aes(y = Mean, x = treatment), position = position_dodge(width = 0.75), color = "black", fill = "black", size = 3) +
  geom_segment(data = data_plot, aes(y = Mean - SD, yend = Mean + SD, x = treatment, xend = treatment), size = 1.5, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  ylim(min(data_fig$slopes), max(data_fig$slopes)) +
  coord_flip() +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_classic() +
  labs(y = "Slope estimates", x = "Treatments") +
  theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5, "mm"),
    axis.text.y = element_blank(),   # Remove y-axis labels
    axis.title = element_text(size = 12, family = "Times"),
    axis.text = element_text(size = 10, family = "Times"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 11, family = "Times")
    )
  return(plot2)
}
####################
####################
# Function to create the plot A for fig-results_learning
#' @title plot_errorsday
#' @param df to select the df
plot_errorsday <- function(df){
  df_prob <- df %>%
    mutate(treatment = factor(treatment,
          levels = c("Control-Hot", "CORT-Hot", "Control-Cold", "CORT-Cold"),
          labels = c("Control-Hot" = "Control-Hot (n=20)",
                    "CORT-Hot" = "CORT-Hot (n=19)",
                    "Control-Cold" = "Control-Cold (n=20)",
                    "CORT-Cold" = "CORT-Cold (n=20)")
          )) %>%
  data.frame()
  plot <- ggplot(df_prob, aes(x = trial, y = errors, color = treatment)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c("Control-Hot (n=20)"="#fa927d",
                                "CORT-Hot (n=19)"="#b50101",
                                "Control-Cold (n=20)"="#68bde1",
                                "CORT-Cold (n=20)"="#00008B")) +
    geom_ribbon(aes(ymin = errors - sd, ymax = errors + sd, fill = treatment), color = NA, alpha = 0.08) + 
    scale_fill_manual(values = c("Control-Hot (n=20)"="#fa927d",
                                "CORT-Hot (n=19)"="#b50101",
                                "Control-Cold (n=20)"="#68bde1",
                                "CORT-Cold (n=20)"="#00008B")) +
    theme_classic() +
    labs(y = "Predicted errors", x = "Trial") +
    theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5, "mm")) +
    theme(
      axis.title = element_text(size = 12, family = "Times"),
      axis.text = element_text(size = 10, family = "Times"),
      legend.position = "none"
    )
  return(plot)
}
