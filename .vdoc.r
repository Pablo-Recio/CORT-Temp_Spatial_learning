#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: setup
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes, ggh4x, cowplot, fitdistrplus, MASS, goftest, forcats, nortest, fitdistrplus, ggh4x, PupillometryR, png, grid, remotes, ggthemes, bayestestR, HDInterval, DiagrammeR, magick)
#
#
#
#| label: data_processing
# The result will be the final df with the data for the analysis. To do so, we first estimate the learning slope (choice and errors) for each individual and then merge the data with the mitochondrial data, extracted using the script in extract.R and extraction_finc.R (see R folder). The final df will be saved in (here("output/databases_clean/data_complete.csv") 
refit <- FALSE
source(here("R", "data_process.R"))
#
#
#
#| label: data_metrics
# Using the raw learning database, I calculate here some basic metrics mentioned later in the Methdos (see below).
data <- read.csv(here("./data/Spatial_learn.csv"))
#
# Total number of clutches
total_clutches <- data %>% 
  distinct(clutch) %>% 
  count() %>% 
  pull(n)
# Calculate total number of individuals
total_individuals <- data %>% 
  distinct(lizard_id) %>% 
  count() %>% 
  pull(n)
# Calculate number of individuals per clutch
individuals_per_clutch <- data %>%
  distinct(lizard_id, clutch) %>%
  group_by(clutch) %>%
  summarize(individual_count = n())
#
#
#
#| label: sampleSize
# List with the sample sizes from the main database.
source(here("R", "func.R"))
#
hormone <- c("CORT", "Control")
temperature <- c("Cold", "Hot")
#
n_list <- list()
#
for(k in 1:length(hormone)){
  for(l in 1:length(temperature)){
    list_name <- paste0(hormone[k], "_", temperature[l])
    n_list[[list_name]] <- sample(df = clean_df, corti = hormone[k], therm = temperature[l])
  }
}
#
df_damage <- clean_df %>%
  filter(DNAdamage != "NA")
n_damage <- list()
for(k in 1:length(hormone)){
  for(l in 1:length(temperature)){
    list_name <- paste0(hormone[k], "_", temperature[l])
    n_damage[[list_name]] <- sample(df = df_damage, corti = hormone[k], therm = temperature[l])
  }
}
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: models_mitochondrial
# Fitting the model and extraction of posteriors for both types of task and species using fit_m function (see func.r in R folder). The result everytime the function is used is a df with the posteriors of the model. The functions saves the model automatically in output/models; and when the parameter refit = FALSE then the posteriors are extracted from the model previously written instead of fitting the model again each time.
source(here("R", "func.R"))
# 
#
# Run models mitochondrial physiology 
var_m <- c("mit_density", "mit_potential", "ROS", "DNAdamage", "peroxidation")
for (p in var_m){ 
  if (p %in% c("mit_density", "mit_potential", "ROS", "DNAdamage")){
    formula <- paste0(p, "~ cort*temp + (1|clutch)")
  } else if (p == "peroxidation"){
    formula <- paste0(p, "~ cort*temp + age + (1|clutch)")
  } 
  pmodel_name <- paste0("m_def_", p)
  assign(pmodel_name, fit_m(df = clean_df,
                          cat = "def",
                          var = p,
                          formula = formula,
                          fam = gaussian(),
                          refit = FALSE),
        envir = .GlobalEnv)  # Assign to the global environment
} 
#
#
#
#| label: models_learning
#
formula_learn <- errors ~ day*cort*temp + (1 + day|lizard_id) + (1|clutch)
m_def_learn <- fit_m(df = learning_df,
                      cat = "def",
                      var = "learning",
                      formula = formula_learn,
                      fam = negbinomial(link = "log"),
                      refit = FALSE)
#
#
#
#| label: organise_posteriors_learning
# Rename some of the posteriors and make new estimates for the learning rate 
#
### Slopes
slope_ControlCold <- m_def_learn$b_day
slope_CORTCold <- m_def_learn$b_day + m_def_learn$`b_day:cortCORT`
slope_ControlHot <- m_def_learn$b_day + m_def_learn$`b_day:tempHot`
slope_CORTHot <- m_def_learn$b_day + m_def_learn$`b_day:cortCORT:tempHot` + m_def_learn$`b_day:cortCORT` + m_def_learn$`b_day:tempHot`
### Intercepts (for figure)
int_ControlCold <- m_def_learn$b_Intercept
int_CORTCold <- m_def_learn$b_Intercept + m_def_learn$`b_cortCORT`
int_ControlHot <- m_def_learn$b_Intercept + m_def_learn$`b_tempHot`
int_CORTHot <- m_def_learn$b_Intercept + m_def_learn$`b_cortCORT:tempHot` + m_def_learn$`b_cortCORT` + m_def_learn$`b_tempHot`
#
#
#
#| label: fig-learning
#| fig.cap: "Results for learning analyses. (a) the predicted number of errors over trials. The lines represent the mean predicted number of errors for each trial, and the shaded areas indicate the standard error of the mean; both were obtained using the slope and intercept estimates from the posterior distributions. (b) shows the distribution of the estimates of slopes per each treatment. The x-axis represents the slope estimate, and in the y-axis are the density of the estimates. Points and bars represent the mean and standard error of the estimated slopes, respectively. Dashed lines indicate value 0. The different colors in both panels indicate the different treatments. Asterisks in (b) indicate significant differences from 0."
source(here("R", "func.R"))
#
#
####### A) Create df
# Slopes df (employed in fig-learning B)
slope_list <- list(`Control-Cold` = slope_ControlCold,
                 `CORT-Cold` = slope_CORTCold,
                 `Control-Hot` = slope_ControlHot,
                 `CORT-Hot` = slope_CORTHot)
data_fig_learning_slopes <- do.call(rbind, lapply(names(slope_list), function(x) {
  data.frame(treatment = x, slopes = slope_list[[x]])
}))
# Intercepst df (only for fig-learning A)
intercept_list <- list(`Control-Cold` = int_ControlCold,
                 `CORT-Cold` = int_CORTCold,
                 `Control-Hot` = int_ControlHot,
                 `CORT-Hot` = int_CORTHot)
data_fig_learning_intercepts <- do.call(rbind, lapply(names(intercept_list), function(x) {
  data.frame(treatment = x, intercepts = intercept_list[[x]])
})) %>% 
  group_by(treatment) %>%
  summarize(int_sd = sd(intercepts),
            intercepts = mean(intercepts))
# df fig-learning A
data_fig_learning_slopesA <- data_fig_learning_slopes %>%
  group_by(treatment) %>%
  summarize(slope_sd = sd(slopes),
            slopes = mean(slopes),)
data_fig_learningA <- merge(data_fig_learning_slopesA, data_fig_learning_intercepts, by = "treatment")
treatment <- unique(data_fig_learningA$treatment)
fig_A_df <- data.frame()
for(t in treatment){ # Loop per treatment
  df <- data_fig_learningA %>%
    filter(treatment == t) %>%
    data.frame()
  # Variables selected
  m <- df$slopes
  u <- df$intercepts
  # Loop per treatment
  num_individuals <- length(u)
  for(x in 0:40){
    if (x == 0){
      sd <- df$int_sd
    } else {
      q25 <- df$slope_sd
    }
    value <- exp(u + m * x)
    temp_df <- data.frame(trial = rep(x, length(value)),
                          treatment = rep(t, length(value)),
                          errors = value,
                          sd = sd * value)
    fig_A_df <- rbind(fig_A_df, temp_df)
  }
}
# 
#### C) Plot the fig-learning A and fig-learning B separately 
fig_learningA <- plot_errorsday(fig_A_df)
fig_learningB <- plot_slopes(data_fig_learning_slopes)
#
#### D) Combine plots A and B to have fig-learning
fig_results_learning <- plot_grid(fig_learningA, fig_learningB, nrow = 1) +
  # Insert title for each plot
  annotate("text", x = 0.045, y = 0.935, label = "A", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("text", x = 0.545, y = 0.935, label = "B", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("text", x = 0.705, y = 0.89, label = "*", hjust = 1, vjust = 1, size = 7) +
  annotate("text", x = 0.68, y = 0.72, label = "*", hjust = 1, vjust = 1, size = 7) +
  annotate("text", x = 0.7, y = 0.55, label = "*", hjust = 1, vjust = 1, size = 7) +
  annotate("text", x = 0.69, y = 0.38, label = "*", hjust = 1, vjust = 1, size = 7)
#
ggsave("./output/figures/text/fig_results_errors.png", plot=fig_results_learning , width = 21, height = 7, units = "cm", dpi = 600)
knitr::include_graphics("./output/figures/text/fig_results_errors.png")
#
#
#
#
#
#
#
#| label: organise_posteriors_mit
#
# Organising the posteriors of the previous models for tables and figures
#
source(here("R", "func.R"))
#
#
post_mit <- data.frame()
# Names posteriors:
names <- c("m_def_mit_density", "m_def_mit_potential", "m_def_ROS", "m_def_DNAdamage", "m_def_peroxidation")
#
# Organising the results
for (pos in names) {
  model <- get(pos)      # Get the model from the global environment
  post_result <- tidy_post(model)        # Apply tidy_post to each model
  # Add a new column to the rmodel
  post_result$Model <- pos
  post_mit <- bind_rows(post_mit, post_result)
}
#
#
#
# Extracting the posteriors for the models with mitochondrial variables and the values of interest. Here, I am creating dfs for each of the variables with the values for all the prenatal conditions to make contrasts easier to write.
#
MD <- post_values(m_def_mit_density, "none")
MP <- post_values(m_def_mit_potential, "none")
ROS <- post_values(m_def_ROS, "none")
DNA <- post_values(m_def_DNAdamage, "none")
LP <- post_values(m_def_peroxidation, "none")
#
#
#
#| label: fig-results_mit
#| fig-cap: "Estimates of mitochondrial density (a), metabolic capacity (b), ROS (c), DNA damage (d), and lipid peroxidation (e) in the medial cortex as a function of the different prenatal conditions. Note that for (d) and (e) these analyses do not account for missing data resulting from a flow cytometer malfunction that impacted one plate, and so, sample sizes are lower for univariate analyses. The x-axis represents the estimated values and in the y-axis is the density of the estimates. Points and bars represent the mean and standard error of the estimated values, respectively. The different colors in both panels indicate the different treatments."
#| fig-name: "fig-results_oxidative"
#
source(here("R", "func.R"))
#
# A) Plotting the results for all variables
plot_mit_density <- plotting(MD, "Mit density")
plot_mit_potential <- plotting(MP, "Metabolic capacity")
plot_ros <- plotting(ROS, "ROS")
#
plot_legend_top <- plotting(ROS, "ROS") + theme(legend.position = "bottom", legend.title = element_blank())
gtable <- ggplot_gtable(ggplot_build(plot_legend_top))
legend_mit_top <- gtable$grobs[[which(sapply(gtable$grobs, function(x) x$name) == "guide-box")]]
#
plot_dnadamage <- plotting(DNA, "DNA damage")
plot_peroxidation <- plotting(LP, "Lipid peroxidation")
#
plot_legend_bottom <- plotting(DNA, "DNA damage") + theme(legend.position = "bottom", legend.title = element_blank())
gtable <- ggplot_gtable(ggplot_build(plot_legend_bottom))
legend_mit_bottom <- gtable$grobs[[which(sapply(gtable$grobs, function(x) x$name) == "guide-box")]]
#
# B) Organising the plots
fig_mit_top <- plot_grid(plot_mit_density, plot_mit_potential,
                        nrow = 1, rel_widths = c(1, 1))
fig_mit_mid <- plot_grid(plot_ros, NULL,
                        nrow = 1, rel_widths = c(1, 1))
fig_mit_bottom <- plot_grid(plot_dnadamage, plot_peroxidation,
                        nrow = 1, rel_widths = c(1, 1))
#
# C) Merging everything in final figure
# Create the figure grid with extra space for legend
fig_mit <- plot_grid(
  fig_mit_top, fig_mit_mid, NULL, fig_mit_bottom, NULL,
  nrow = 5, rel_heights = c(0.3, 0.3, 0.05, 0.3, 0.05))
# Final composition: Merge everything, adding images and legend
final_plot_mit <- ggdraw(fig_mit) +
  # Insert legend in the middle for ROS and energy variables
  draw_grob(legend_mit_top, x = 0.49, y = 0.39, width = 0.01, height = 0.01) +
  # Insert legend in the bottom for DNA damage and lipid peroxidation
  draw_grob(legend_mit_bottom, x = 0.49, y = 0.03, width = 0.01, height = 0.01) +
  # Insert title for each plot
  annotate("text", x = 0.0325, y = 0.985, label = "A", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("text", x = 0.531, y = 0.985, label = "B", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("text", x = 0.0325, y = 0.682, label = "C", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("text", x = 0.0325, y = 0.335, label = "D", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("text", x = 0.53, y = 0.335, label = "E", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  # Insert figure brain
  draw_image(here("./Others/brain_fig.png"), x = 0.53, y = 0.34, width = 0.45, height = 0.45)
#
# Print final plot
ggsave(here("./output/figures/text/results_mit.png"), plot = final_plot_mit, width = 21, height = 21, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/text/results_mit.png")
#
#
#
#| label: model_sem
# Making the SEM model by using a multivariate brms. The aim is to test the relationships between mitochondrial physiology and detection latency.
# Learning here was evaluated obtaining the learning slopes for each individual (see data_process.R).
# All continue variables were standardized (var/2SD) before running the models (see data_process.R). 
#
source(here("R", "func.R"))
#
refit <- FALSE
#
data_SEM <- clean_df %>%
  mutate(obs = 1:nrow(.),
        vec = rep(1, length(obs)))
if(refit){
  m_SEM <- brm(
    bf(slope_mean | se(slope_sd, sigma = FALSE) ~ mit_density + mit_potential + mi(DNAdamage) + mi(peroxidation) + (1|clutch) + (1|q|obs)) +
    bf(DNAdamage | mi() + se(vec, sigma = FALSE) ~ ROS + (1|clutch)+ (1|p|obs)) +
    bf(peroxidation | mi() + se(vec, sigma = FALSE) ~ ROS + (1|clutch) + (1|p|obs)) +
    bf(ROS | se(vec, sigma = FALSE) ~ mit_density + mit_potential + (1|clutch) + (1|t|obs)) +
  set_rescor(FALSE),
  family = gaussian(),
  data = data_SEM, 
  chains = 4, cores = 4, iter = 8000, warmup = 2000,
  control = list(adapt_delta = 0.99, max_treedepth = 15))
  # Save the model
  saveRDS(m_SEM, file = here("output/models/m_SEM.rds"))
} else {
  m_SEM <- readRDS(here("output/models/m_SEM.rds"))
}
#
#
#
#| label: sem_tidy
source(here("R", "func.R"))
#
# I am extracting here all the values for getting the total effects of each of the variables in the model. I am using the posterior values for each of the variables to get the total effects assuming that:
## total effect = direct effect + indirect effect + residual correlation
# In other words:
## total effect = 
# Extract the posteriors for the SEM model
post_sem <- as_draws_df(m_SEM) 
#
#### A) Get the direct paths / regression coefficients per each variable
# 
# A.1) Learning
coeff_mitodensity_learn <- post_sem$b_slopemean_mit_density
coeff_potential_learn <- post_sem$b_slopemean_mit_potential
coeff_dna_learn<- post_sem$bsp_slopemean_miDNAdamage
coeff_perox_learn <- post_sem$bsp_slopemean_miperoxidation
# 
# A.2) DNA damage
coeff_ros_dna <- post_sem$b_DNAdamage_ROS
#
# A.3) Lipid peroxidation
coeff_ros_perox <- post_sem$b_peroxidation_ROS
#
# A.4) ROS
coeff_mitodensity_ros <- post_sem$b_ROS_mit_density
coeff_potential_ros <- post_sem$b_ROS_mit_potential
#
#
#### B) Get the indirect paths for each variable
#
# B.1) Detection
undir_ros_learn <- coeff_ros_dna * coeff_dna_learn + coeff_ros_perox * coeff_perox_learn
undir_mitodensity_learn <- coeff_mitodensity_ros * coeff_ros_dna * coeff_dna_learn + coeff_mitodensity_ros * coeff_ros_perox * coeff_perox_learn
undir_potential_learn <- coeff_potential_ros * coeff_ros_dna * coeff_dna_learn + coeff_potential_ros * coeff_ros_perox * coeff_perox_learn
#
# B.2) DNA damage
undir_mitodensity_dna <- coeff_mitodensity_ros * coeff_ros_dna
undir_potential_dna <- coeff_potential_ros * coeff_ros_dna
#
# B.3) Lipid peroxidation
undir_mitodensity_perox <- coeff_mitodensity_ros * coeff_ros_perox
undir_potential_perox <- coeff_potential_ros * coeff_ros_perox
#
#
#### C) Get the total effects for each variable
#
# C.1) Detection
total_mitodensity_learn <- coeff_mitodensity_learn + undir_mitodensity_learn
total_potential_learn <- coeff_potential_learn + undir_potential_learn
total_ros_learn <- undir_ros_learn
total_dna_learn <- coeff_dna_learn
total_perox_learn <- coeff_perox_learn
#
# C.2) DNA damage
total_ros_dna <- coeff_ros_dna
total_mitodensity_dna <- undir_mitodensity_dna
total_potential_dna <- undir_potential_dna
#
# C.3) Lipid peroxidation
total_ros_perox <- coeff_ros_perox
total_mitodensity_perox <- undir_mitodensity_perox
total_potential_perox <- undir_potential_perox
#
# C.4) ROS
total_mitodensity_ros <- coeff_mitodensity_ros
total_potential_ros <- coeff_potential_ros
#
#
# D) Create a df with the values for each variable
#
# D.1) Detection (for example)
learn_sem <- data.frame(
  reference_variable = rep("learn", 5),
  predictor_modulator = c("mit_density",
                          "mit_potential",
                          "ROS",
                          "DNAdamage",
                          "peroxidation"),
  direct_effects = I(list(coeff_mitodensity_learn,
                          coeff_potential_learn,
                          NA,
                          coeff_dna_learn,
                          coeff_perox_learn)),
  indirect_effects = I(list(undir_mitodensity_learn,
                          undir_potential_learn,
                          undir_ros_learn,
                          NA,
                          NA)),
  total_effects = I(list(total_mitodensity_learn,
                        total_potential_learn,
                        total_ros_learn,
                        total_dna_learn,
                        total_perox_learn))
  ) %>% 
  group_by(predictor_modulator) %>%
  summarize(
    # Summarizing each list column (direct_effects, indirect_effects, etc.)
    mean_direct_effects = format_dec(mean(unlist(direct_effects), na.rm = TRUE), 3),
    mean_indirect_effects = format_dec(mean(unlist(indirect_effects), na.rm = TRUE), 3),
    mean_total_effects = format_dec(mean(unlist(total_effects), na.rm = TRUE), 3),
    q5_direct = format_dec(quantile(unlist(direct_effects), 0.05, na.rm = TRUE), 3),
    q95_direct = format_dec(quantile(unlist(direct_effects), 0.95, na.rm = TRUE), 3),
    q5_indirect = format_dec(quantile(unlist(indirect_effects), 0.05, na.rm = TRUE), 3),
    q95_indirect = format_dec(quantile(unlist(indirect_effects), 0.95, na.rm = TRUE), 3),
    q5_total = format_dec(quantile(unlist(total_effects), 0.05, na.rm = TRUE), 3),
    q95_total = format_dec(quantile(unlist(total_effects), 0.95, na.rm = TRUE), 3)) %>%
  mutate(source = "learn_sem")
#
# D.2) DNA damage
DNA_sem <- data.frame(
  reference_variable = rep("DNA damage", 3),
  predictor_modulator = c("mit_density",
                          "mit_potential",
                          "ROS"),
  direct_effects = I(list(NA,
                          NA,
                          coeff_ros_dna)),
  indirect_effects = I(list(undir_mitodensity_dna,
                          undir_potential_dna,
                          NA)),
  total_effects = I(list(total_mitodensity_dna,
                        total_potential_dna,
                        total_ros_dna))) %>% 
  group_by(predictor_modulator) %>%
  summarize(
    # Summarizing each list column (direct_effects, indirect_effects, etc.)
    mean_direct_effects = format_dec(mean(unlist(direct_effects), na.rm = TRUE), 3),
    mean_indirect_effects = format_dec(mean(unlist(indirect_effects), na.rm = TRUE), 3),
    mean_total_effects = format_dec(mean(unlist(total_effects), na.rm = TRUE), 3),
    q5_direct = format_dec(quantile(unlist(direct_effects), 0.05, na.rm = TRUE), 3),
    q95_direct = format_dec(quantile(unlist(direct_effects), 0.95, na.rm = TRUE), 3),
    q5_indirect = format_dec(quantile(unlist(indirect_effects), 0.05, na.rm = TRUE), 3),
    q95_indirect = format_dec(quantile(unlist(indirect_effects), 0.95, na.rm = TRUE), 3),
    q5_total = format_dec(quantile(unlist(total_effects), 0.05, na.rm = TRUE), 3),
    q95_total = format_dec(quantile(unlist(total_effects), 0.95, na.rm = TRUE), 3)) %>%
  mutate(source = "DNA_sem")
#
# D.3) Lipid peroxidation
perox_sem <- data.frame(
  reference_variable = rep("lipid peroxidation", 3),
  predictor_modulator = c("mit_density",
                          "mit_potential",
                          "ROS"),
  direct_effects = I(list(NA,
                          NA,
                          coeff_ros_perox)),
  indirect_effects = I(list(undir_mitodensity_perox,
                          undir_potential_perox,
                          NA)),
  total_effects = I(list(total_mitodensity_perox,
                        total_potential_perox,
                        total_ros_perox))) %>% 
  group_by(predictor_modulator) %>%
  summarize(
    # Summarizing each list column (direct_effects, indirect_effects, etc.)
    mean_direct_effects = format_dec(mean(unlist(direct_effects), na.rm = TRUE), 3),
    mean_indirect_effects = format_dec(mean(unlist(indirect_effects), na.rm = TRUE), 3),
    mean_total_effects = format_dec(mean(unlist(total_effects), na.rm = TRUE), 3),
    q5_direct = format_dec(quantile(unlist(direct_effects), 0.05, na.rm = TRUE), 3),
    q95_direct = format_dec(quantile(unlist(direct_effects), 0.95, na.rm = TRUE), 3),
    q5_indirect = format_dec(quantile(unlist(indirect_effects), 0.05, na.rm = TRUE), 3),
    q95_indirect = format_dec(quantile(unlist(indirect_effects), 0.95, na.rm = TRUE), 3),
    q5_total = format_dec(quantile(unlist(total_effects), 0.05, na.rm = TRUE), 3),
    q95_total = format_dec(quantile(unlist(total_effects), 0.95, na.rm = TRUE), 3)) %>%
  mutate(source = "perox_sem")
#
# D.4) ROS
ROS_sem <- data.frame(
  reference_variable = rep("ROS", 2),
  predictor_modulator = c("mit_density",
                          "mit_potential"),
  direct_effects = I(list(coeff_mitodensity_ros,
                          coeff_potential_ros)),
  indirect_effects = I(list(NA,
                          NA)),
  total_effects = I(list(total_mitodensity_ros,
                        total_potential_ros))) %>% 
  group_by(predictor_modulator) %>%
  summarize(
    # Summarizing each list column (direct_effects, indirect_effects, etc.)
    mean_direct_effects = format_dec(mean(unlist(direct_effects), na.rm = TRUE), 3),
    mean_indirect_effects = format_dec(mean(unlist(indirect_effects), na.rm = TRUE), 3),
    mean_total_effects = format_dec(mean(unlist(total_effects), na.rm = TRUE), 3),
    q5_direct = format_dec(quantile(unlist(direct_effects), 0.05, na.rm = TRUE), 3),
    q95_direct = format_dec(quantile(unlist(direct_effects), 0.95, na.rm = TRUE), 3),
    q5_indirect = format_dec(quantile(unlist(indirect_effects), 0.05, na.rm = TRUE), 3),
    q95_indirect = format_dec(quantile(unlist(indirect_effects), 0.95, na.rm = TRUE), 3),
    q5_total = format_dec(quantile(unlist(total_effects), 0.05, na.rm = TRUE), 3),
    q95_total = format_dec(quantile(unlist(total_effects), 0.95, na.rm = TRUE), 3)) %>%
  mutate(source = "ROS_sem")
#
#
# E) Merge everything into a single df
#
sem_results <- bind_rows(learn_sem, DNA_sem, perox_sem, ROS_sem)
#
#
#
#
#
#| label: fig-sem
#| fig-cap: "Structural equation model testing hypothesized direct and indirect effects of physiology on learning. Arrows indicate the directionality of the estimate, and the values show the mean predicted coefficient and its 95% confidence interval."
#| fig-name: "fig-sem"
# 
# A) Getting all the direct coefficients
density_learn <- paste0(format_dec(mean(coeff_mitodensity_learn), 3),
                  " [", format_dec(quantile(coeff_mitodensity_learn, 0.05), 3),
                  ", ", format_dec(quantile(coeff_mitodensity_learn, 0.95), 3), "]")
potential_learn <- paste0(format_dec(mean(coeff_potential_learn), 3),
                  " [", format_dec(quantile(coeff_potential_learn, 0.05), 3),
                  ", ", format_dec(quantile(coeff_potential_learn, 0.95), 3), "]")
dna_learn <- paste0(format_dec(mean(coeff_dna_learn), 3),
                  " [", format_dec(quantile(coeff_dna_learn, 0.05), 3),
                  ", ", format_dec(quantile(coeff_dna_learn, 0.95), 3), "]")
perox_learn <- paste0(format_dec(mean(coeff_perox_learn), 3),
                  " [", format_dec(quantile(coeff_perox_learn, 0.05), 3),
                  ", ", format_dec(quantile(coeff_perox_learn, 0.95), 3), "]")
#
ros_dna <- paste0(format_dec(mean(coeff_ros_dna), 3),
                  " [", format_dec(quantile(coeff_ros_dna, 0.05), 3),
                  ", ", format_dec(quantile(coeff_ros_dna, 0.95), 3), "]")
#
ros_perox <- paste0(format_dec(mean(coeff_ros_perox), 3),
                  " [", format_dec(quantile(coeff_ros_perox, 0.05), 3),
                  ", ", format_dec(quantile(coeff_ros_perox, 0.95), 3), "]")
#
density_ros <- paste0(format_dec(mean(coeff_mitodensity_ros), 3),
                  " [", format_dec(quantile(coeff_mitodensity_ros, 0.05), 3),
                  ", ", format_dec(quantile(coeff_mitodensity_ros, 0.95), 3), "]")
potential_ros <- paste0(format_dec(mean(coeff_potential_ros), 3),
                  " [", format_dec(quantile(coeff_potential_ros, 0.05), 3),
                  ", ", format_dec(quantile(coeff_potential_ros, 0.95), 3), "]")
#
imgSEM <- readPNG(here("Others", "SEM.png"))
plot_SEM <- rasterGrob(imgSEM, interpolate = TRUE)
#
fig_SEM <- ggdraw(plot_SEM) +
  annotate("text", x = 0.433, y = 0.19, label = density_learn, hjust = 1, vjust = 1, size = 3.5, family = "Times") +
  annotate("text", x = 0.433, y = 0.85, label = potential_learn, hjust = 1, vjust = 1, size = 3.5, family = "Times") +
  annotate("text", x = 0.877, y = 0.675, label = dna_learn, hjust = 1, vjust = 1, size = 3.5, family = "Times") +
  annotate("text", x = 0.877, y = 0.387, label = perox_learn, hjust = 1, vjust = 1, size = 3.5, family = "Times") +
  annotate("text", x = 0.547, y = 0.745, label = ros_dna, hjust = 1, vjust = 1, size = 3.5, family = "Times") +
  annotate("text", x = 0.547, y = 0.295, label = ros_perox, hjust = 1, vjust = 1, size = 3.5, family = "Times") +
  annotate("text", x = 0.325, y = 0.284, label = density_ros, hjust = 1, vjust = 1, size = 3.5, family = "Times") +
  annotate("text", x = 0.324, y = 0.757, label = potential_ros, hjust = 1, vjust = 1, size = 3.5, family = "Times")
ggsave(here("./output/figures/text/SEM.png"), plot = fig_SEM, width = 21, height = 10, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/text/SEM.png")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#
#
#| label: fig_learning_raw
#
source(here("R", "func.R"))
# Modify the df to plot the raw data
learning_df_plot <- learning_df %>%
  mutate(treatment = factor(trt, levels = c("B_23", "A_23", "B_28", "A_28"),
                      labels = c("B_23" = "CORT-Cold (n=20)",
                                 "A_23" = "Control-Cold (n=20)",
                                 "B_28" = "CORT-Hot (n=19)",
                                 "A_28" = "Control-Hot (n=20)"))) %>%
  group_by(treatment, day) %>%
  mutate(mean_raw = mean(na.omit(errors)),
         sd_raw = sd(na.omit(errors))) %>%
data.frame()
# Get the plot and combine it with the modeled data
fig_raw_learning <- plot_errorsday(fig_A_df) +
  geom_point(data = learning_df_plot, aes(x = day, y = mean_raw, color = treatment), size = 1.5) +
  scale_color_manual(values = c("CORT-Cold (n=20)" = "#00008B",
                                "Control-Cold (n=20)" = "#68bde1",
                                "CORT-Hot (n=19)" = "#b50101",
                                "Control-Hot (n=20)" = "#fa927d")) +
  facet_grid(~treatment) +
  ylim(c(min(na.omit(learning_df_plot$mean_raw)), max(na.omit(learning_df_plot$mean_raw)))) +  # Adjust y-axis limits
  theme(
    axis.title = element_text(size = 12, family = "Times"),
    axis.text = element_text(size = 10, family = "Times"),
    legend.position = "bottom",
    legend.title = element_blank()
  )
ggsave("./output/figures/suppl/Figure_S1_learning_raw.png", plot=fig_raw_learning , width = 21, height = 15, units = "cm", dpi = 600)
knitr::include_graphics("./output/figures/suppl/Figure_S1_learning_raw.png")
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#| label: tbl-learning
#| tbl-cap: "Learning tasks slopes per treatment"
source(here("R", "func.R"))
#
table_learn_df <- data_fig_learning_slopes %>%
  group_by(treatment) %>%
  summarize(mean_slope = format_dec(mean(slopes), 3),
            q025_slope = format_dec(quantile(slopes, 0.025), 3),
            q975_slope = format_dec(quantile(slopes, 0.975), 3),
            pMCMC = format_p(pmcmc(slopes), 3, equal = FALSE)) %>%
  mutate(`95 CI` = paste0("[", q025_slope, ", ", q975_slope, "]")) %>%
  dplyr::select(treatment, mean_slope, `95 CI`, pMCMC)
#
# Make Table
set_flextable_defaults(
 font.family = "Times New Roman",
 font.size = 10)
#
table_slopes <- flextable(table_learn_df) %>%
  align(align = "center", j = c(3,4), part = "body") %>%
  align(align = "center", j = c(1:4), part = "header") %>%
  flextable::compose(j = 1, value = as_paragraph("Treatment"), part = "header") %>%
  flextable::compose(j = 2, value = as_paragraph("Estimated slope"), part = "header") %>%
  bold(~`pMCMC` == "< 0.05", j = c("pMCMC", "mean_slope", "95 CI", "treatment"),
       bold = TRUE) %>%  # Bold when PMCMC is "<0.05"
  bold(~`pMCMC` == "< 0.001", j = c("pMCMC", "mean_slope", "95 CI", "treatment")) %>%  # Bold when PMCMC is "<0.001"
  autofit()
#
table_slopes
#
#
#
cat("\\newpage")
#
#
#
#
#| label: tbl-contrasts
#| tbl-cap: "Contrasts between prenatal conditions for mitochondrial physiology and learning."
#
source(here("R", "func.R"))
#
# A) Modify posteriors df for learning where we only get/compare the slopes
LEARN <- data.frame(Control_Cold = slope_ControlCold, 
                  CORT_Cold = slope_CORTCold,
                  Control_Hot = slope_ControlHot,
                  CORT_Hot = slope_CORTHot)

# B) Organise df for mit physiology & learning
var <- c("MD", "MP", "ROS", "DNA", "LP", "LEARN")
data_table <- data.frame()
for(x in var){
  df <- get(x)
  Temperature <- format_dec(mean(c(df$CORT_Hot, df$Control_Hot)) - mean(c(df$CORT_Cold, df$Control_Cold)), 3)
  pMCMC_temp <- format_p(pmcmc(c(df$CORT_Hot, df$Control_Hot) - c(df$CORT_Cold, df$Control_Cold)), 3, equal = FALSE)
  CORT <- format_dec(mean(c(df$Control_Hot, df$Control_Cold)) - mean(c(df$CORT_Hot, df$CORT_Cold)), 3)
  pMCMC_cort <- format_p(pmcmc(c(df$Control_Hot, df$Control_Cold) - c(df$CORT_Hot, df$CORT_Cold)), 3, equal = FALSE)
  Interaction <- format_dec((mean(df$Control_Hot) - mean(df$CORT_Hot)) - (mean(df$Control_Cold) - mean(df$CORT_Cold)), 3)
  pMCMC_int <- format_p(pmcmc((df$Control_Hot - df$CORT_Hot) - (df$Control_Cold - df$CORT_Cold)), 3, equal = FALSE)
  data_temp <- data.frame(Variable = x,
                          Temperature = as.numeric(Temperature),
                          pMCMC_temp = as.numeric(pMCMC_temp),
                          CORT = as.numeric(CORT),
                          pMCMC_cort = as.numeric(pMCMC_cort),
                          Interaction = as.numeric(Interaction),
                          pMCMC_int = as.numeric(pMCMC_int))
  data_table <- dplyr::bind_rows(data_table, data_temp)
}
# Modify the df
data_table_final <- data_table %>%
  pivot_longer(cols = c(Temperature, CORT, Interaction), 
               names_to = "Predictor", 
               values_to = "Contrast") %>%
  mutate(
    # Extract the pMCMC values from the corresponding columns
    `pMCMC contrast` = case_when(
      Predictor == "Temperature" ~ pMCMC_temp,
      Predictor == "CORT" ~ pMCMC_cort,
      Predictor == "Interaction" ~ pMCMC_int
    )
  ) %>%
  mutate(
    Variable = case_when(
      Variable == "MD" ~ "Mit density",
      Variable == "MP" ~ "Metabolic capacity",
      Variable == "ROS" ~ "ROS",
      Variable == "DNA" ~ "DNA damage",
      Variable == "LP" ~ "Lipid peroxidation",
      Variable == "LEARN" ~ "Learning slopes",
      TRUE ~ Variable
      )
    ) %>%
  dplyr::select(Variable, Predictor, Contrast, `pMCMC contrast`)
#
# C) Make the contrasts table:
#
set_flextable_defaults(
 font.family = "Times New Roman",
 font.size = 10)
#
contrast_table <- flextable(data_table_final) %>%
  align(align = "center", j = c(3,4), part = "body") %>%
  align(align = "center", j = c(1:4), part = "header") %>%
  bold(~`pMCMC contrast` < 0.05, j = c("pMCMC contrast", "Contrast", "Predictor"),
       bold = TRUE) %>%  # Bold when PMCMC is "<0.05"
  bold(~`pMCMC contrast` <0.001, j = c("pMCMC contrast", "Contrast", "Predictor")) %>%  # Bold when PMCMC is "<0.001"
  flextable::compose(i = c(2,3,5,6,8,9,11,12,14,15,17,18), j = 1, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the first column
  autofit()
#
contrast_table
#
#
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#
#
#
#| label: table_sem_results
#| tbl-cap: "Structural Equation Models"
#
# Table created from the df sem_results_OB (see above)
source(here("R", "func.R"))
#
# Modify database
table_sem_df <- sem_results %>%
  mutate(`Direct effects` = paste0(mean_direct_effects, " [", q5_direct, ", ", q95_direct, "]"),
         `Indirect effects` = paste0(mean_indirect_effects, " [", q5_indirect, ", ", q95_indirect, "]"),
         `Total effects` = paste0(mean_total_effects, " [", q5_total, ", ", q95_total, "]")) %>%
  dplyr::select(source, predictor_modulator, `Direct effects`, `Indirect effects`, `Total effects`) %>%
  mutate(
      predictor_modulator = case_when(
        predictor_modulator == "mit_density" ~ "Mitochondrial density",
        predictor_modulator == "mit_potential" ~ "Metabolic capacity",
        predictor_modulator == "DNAdamage" ~ "DNA damage",
        predictor_modulator == "peroxidation" ~ "Lipid peroxidation",
        predictor_modulator == "ROS" ~ "ROS",
        TRUE ~ predictor_modulator
      ),
      source = case_when(
        source == "learn_sem" ~ "Learning",
        source == "DNA_sem" ~ "DNA damage",
        source == "perox_sem" ~ "Lipid peroxidation",
        source == "ROS_sem" ~ "ROS",
        TRUE ~ source
      ),
      `Direct effects` = case_when(`Direct effects` == "NaN [NA, NA]" ~ " - ", TRUE ~ `Direct effects`),
      `Indirect effects` = case_when(`Indirect effects` == "NaN [NA, NA]" ~ " - ", TRUE ~ `Indirect effects`),
      `Total effects` = case_when(`Total effects` == "NaN [NA, NA]" ~ " - ", TRUE ~ `Total effects`)) %>%
  rename(Predictor = predictor_modulator, Response = source) %>%
  mutate(Response = factor(Response, levels = c("Learning", "DNA damage", "Lipid peroxidation", "ROS")),
         Predictor = factor(Predictor, levels = c("Mitochondrial density",
                                                "Metabolic capacity",
                                                "ROS",
                                                "DNA damage",
                                                "Lipid peroxidation"))) %>%
  arrange(Response, Predictor)
#
# Table format
set_flextable_defaults(
 font.family = "Times New Roman",
 font.size = 10)
#
table_sem <- flextable(table_sem_df) %>%
  align(align = "center", part = "header") %>%
  flextable::compose(i = c(2:5,7,8,10,11,13), j = 1, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the first column
  flextable::hline(i = c(5,8,11), part = "body") %>% 
  autofit()
#
table_sem
#
#
#
cat("\\newpage")
#
#
#
#
#
#
#
#| label: plotmod_mitdensity
#| caption: "Posterior predictive checks for the model of Mitochondrial Density in Olfactory Bulbs."
#
mod <- readRDS(here("output/models/mit_density_def.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_S2.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_S2.png")
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#| label: plotmod_potential
#| caption: "Posterior predictive checks for the model of metabolic capacity."
#
mod <- readRDS(here("output/models/mit_potential_def.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_S3.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_S3.png")
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#| label: plotmod_ros
#| caption: "Posterior predictive checks for the model of ROS."
#
mod <- readRDS(here("output/models/ROS_def.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_S4.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_S4.png")
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#| label: plotmod_dnadamage
#| caption: "Posterior predictive checks for the model of DNA damage."
#
mod <- readRDS(here("output/models/DNAdamage_def.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_S5.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_S5.png")
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#| label: plotmod_peroxidation
#| caption: "Posterior predictive checks for the model of lipid peroxidation."
#
mod <- readRDS(here("output/models/peroxidation_def.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_S6.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_S6.png")
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#| label: plotmod_learn
#| caption: "Posterior predictive checks for the model of learning."
#
mod <- readRDS(here("output/models/learning_def.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
plot_mod_3 <- plot_mod[[3]]
#
ggsave("./output/figures/suppl/Figure_S7A.png", plot = plot_mod_1,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_S7A.png")
#
ggsave("./output/figures/suppl/Figure_S7B.png", plot = plot_mod_2,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_S7B.png")
#
ggsave("./output/figures/suppl/Figure_S7C.png", plot = plot_mod_3,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_S7C.png")
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#
#
#| label: models_preliminary
# Fitting preliminary models to see if sex and age are relevant for our models
source(here("R", "func.R"))
#
#
var_m <- c("mit_density", "mit_potential", "ROS", "DNAdamage", "peroxidation")
formula_list_ <- list()
for (p in var_m){
  formula_list_[[p]] <- paste0(p, "~ cort*temp + age + sex + (1|clutch)")
  pmodel_name <- paste0("m_prel_", p)
  assign(pmodel_name, fit_m(df = clean_df,
                             cat = "prel",
                             var = p,
                             formula = formula_list_[[p]],
                             fam = gaussian(),
                             refit = FALSE),
          envir = .GlobalEnv)  # Assign to the global environment
}
#
#
# Run model learning including the effect of sex and interaction
formula_learn <- errors ~ day*cort*temp + sex + age + (1 + day|lizard_id) + (1|clutch)
m_prel_learn <- fit_m(df = learning_df,
                      cat = "prel",
                      var = "learning",
                      formula = formula_learn,
                      fam = negbinomial(link = "log"),
                      refit = FALSE)
#
#
#
#
#
#| tbl-cap: "Preliminary results of the models testing for Mitochondrial Density"
#| label: results_preliminary_mitdensity
#
sum_mitdensity_prel <- m_prel_mit_density %>%
  dplyr::select(starts_with("b_")) %>%  # Removes random effect terms
  summarise_draws()  %>%
  mutate(across(where(is.numeric), ~ as.numeric(format_dec(.x, 3))))
#
flextable(sum_mitdensity_prel)
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#| tbl-cap: "Preliminary results of the models testing for Metabolic capacity"
#| label: results_preliminary_potential
#
sum_potential_prel <- m_prel_mit_potential %>%
  dplyr::select(starts_with("b_")) %>%  # Removes random effect terms
  summarise_draws() %>%
  mutate(across(where(is.numeric), ~ as.numeric(format_dec(.x, 3))))
#
flextable(sum_potential_prel)
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#| tbl-cap: "Preliminary results of the models testing for ROS"
#| label: results_preliminary_ros
#
sum_m_ros_prel <- m_prel_ROS %>%
  dplyr::select(starts_with("b_")) %>%  # Removes random effect terms
  summarise_draws() %>%
  mutate(across(where(is.numeric), ~ as.numeric(format_dec(.x, 3))))
#
flextable(sum_m_ros_prel)
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#| tbl-cap: "Preliminary results of the models testing for DNA damage"
#| label: results_preliminary_dnadamage
#
sum_m_dnadamage_prel <- m_prel_DNAdamage %>%
  dplyr::select(starts_with("b_")) %>%  # Removes random effect terms
  summarise_draws() %>%
  mutate(across(where(is.numeric), ~ as.numeric(format_dec(.x, 3))))
#
flextable(sum_m_dnadamage_prel)
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#| tbl-cap: "Preliminary results of the models testing for lipid peroxidation"
#| label: results_preliminary_peroxidation
#
sum_m_peroxidation_prel <- m_prel_peroxidation %>%
  dplyr::select(starts_with("b_")) %>%  # Removes random effect terms
  summarise_draws() %>%
  mutate(across(where(is.numeric), ~ as.numeric(format_dec(.x, 3))))
#
flextable(sum_m_peroxidation_prel)
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#| tbl-cap: "Preliminary results of the models testing for lipid peroxidation"
#| label: results_preliminary_learning
#
sum_m_learn_prel <- m_prel_learn %>%
  dplyr::select(starts_with("b_")) %>%  # Removes random effect terms
  summarise_draws() %>%
  mutate(across(where(is.numeric), ~ as.numeric(format_dec(.x, 3))))
#
flextable(sum_m_learn_prel)
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#
#
