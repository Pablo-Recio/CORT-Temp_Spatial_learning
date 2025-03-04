####################################
# Exploratory analyses
####################################
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes,
ggh4x, cowplot, fitdistrplus, MASS, goftest)
#
source(here("R", "func.R"))
#
# 
data_expl <- clean_df %>%
  mutate(treatment = factor(str_c(cort, temp, sep = "-"),
                          levels = c("Control-Cold",
                                    "CORT-Cold",
                                    "Control-Hot",
                                    "CORT-Hot")))
#
# A) Make histograms for the main variables and figures factorizing main predictors
var <- c("slope_mean", "mit_density", "mit_potential", "ROS", "DNAdamage", "peroxidation")

hist_list <- list()
for (i in var) {
  bin_width <- (max(data_expl[[i]], na.rm = TRUE) - min(data_expl[[i]], na.rm = TRUE)) / sqrt(length(data_expl[[i]]))
  # Labels
  if (i == "slope_mean"){
     lab <- "Learning"
  } else if (i == "mit_density") {
     lab <- "Mitochondrial density"
  } else if (i == "mit_potential") {
     lab <- "Mitochondrial potential"
  } else if (i == "ROS") {
     lab <- "ROS production"
  } else if (i == "DNAdamage") {
     lab <- "DNA damage"
  } else if (i == "peroxidation") {
     lab <- "Lipid peroxidation"
  }
  # Histograms
  hist_plot <- ggplot(data_expl, aes(x = .data[[i]])) +
    geom_histogram(binwidth = bin_width, fill = "#062d00", color = "black", alpha = 0.5) +
    theme_classic() +
    labs(x = lab, y = "Counts") +
    theme(axis.title = element_text(size = 12, family = "Times"),
          axis.text = element_text(size = 10, family = "Times"))
  hist_list[[i]] <- hist_plot
}
fig_hist <- plot_grid(plotlist = hist_list, ncol = 3)
#
for (i in var) {
  if (i == "slope_mean"){
     lab <- "Learning"
  } else if (i == "mit_density") {
     lab <- "Mitochondrial density"
  } else if (i == "mit_potential") {
     lab <- "Mitochondrial potential"
  } else if (i == "ROS") {
     lab <- "ROS production"
  } else if (i == "DNAdamage") {
     lab <- "DNA damage"
  } else if (i == "peroxidation") {
     lab <- "Lipid peroxidation"
  }
  plot_age <- ggplot(data_expl, aes(x = age, y = .data[[i]], color = treatment)) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = c("CORT-Cold"="#00008B", "Control-Cold"="#68bde1", 
                                 "CORT-Hot"="#b50101", "Control-Hot"="#fa927d")) +
    theme_classic() +
    labs(y = lab, x = "Age (days)") + 
    theme(axis.title = element_text(size = 12, family = "Times"),
              axis.text = element_text(size = 10, family = "Times"))
    
  plot_sex <- ggplot(data_expl, aes(x = sex, y = .data[[i]], fill = treatment)) +
    geom_violin(alpha = 0.5, color = "black") +
    scale_fill_manual(values = c("CORT-Cold"="#00008B", "Control-Cold"="#68bde1", 
                                "CORT-Hot"="#b50101", "Control-Hot"="#fa927d")) +
    theme_classic() +
    labs(y = lab, x = "") + 
    theme(axis.title = element_text(size = 12, family = "Times"),
          axis.text = element_text(size = 10, family = "Times"))
    
  assign(paste0("fig_expl_", i), plot_grid(plot_age, plot_sex, nrow = 1, rel_widths = c(0.5, 0.5)))
  }
#
ggsave(here("./output/figures/exploratory/fig_hist.png"), fig_hist, width = 20, height = 20, dpi = 300)
ggsave(here("./output/figures/exploratory/fig_expl_learn.png"), fig_expl_slope_mean, width = 20, height = 10, dpi = 300)
ggsave(here("./output/figures/exploratory/fig_expl_mean_mitodensity.png"), fig_expl_mit_density, width = 20, height = 10, dpi = 300)
ggsave(here("./output/figures/exploratory/fig_expl_mean_potential.png"), fig_expl_mit_potential, width = 20, height = 10, dpi = 300)
ggsave(here("./output/figures/exploratory/fig_expl_mean_ros.png"), fig_expl_ROS, width = 20, height = 10, dpi = 300)
ggsave(here("./output/figures/exploratory/fig_expl_mean_dnadamage.png"), fig_expl_DNAdamage, width = 20, height = 10, dpi = 300)
ggsave(here("./output/figures/exploratory/fig_expl_mean_peroxidation.png"), fig_expl_peroxidation, width = 20, height = 10, dpi = 300)
#
# B) Get the distribution of the main independent variables
source(here("R", "func.R"))
#
qq_plots_list <- lapply(var, function(v) {
  label <- if (v == "slope_mean"){
    "Learning"
  } else if (v == "mit_density") {
    "Mitochondrial density"
  } else if (v == "mit_potential") {
    "Mitochondrial potential"
  } else if (v == "ROS") {
    "ROS production"
  } else if (v == "DNAdamage") {
    "DNA damage"
  } else if (v == "peroxidation") {
    "Lipid peroxidation"
  } else {
    v
  }
  qq_plots_single(data_expl, v, label)
})
#
final_qq_plot <- plot_grid(plotlist = qq_plots_list, ncol = 1)
ggsave(here("./output/figures/exploratory/fig_qqplots.png"), final_qq_plot, width = 20, height = 20, dpi = 300)
#
# Log-Normal: DNA Damage, Lipid Peroxidation
# Normal: Mitochondrial Density, Mitochondrial Potential, ROS Production, Learning slopes 