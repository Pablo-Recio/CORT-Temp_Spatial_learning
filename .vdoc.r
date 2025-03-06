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
#| label: fig-Methods
#| fig.cap: "Panel **A** shows the experimental design of the study. On top, the treatments applied to the eggs. On the bottom the learning task and the brain region extracted together with the physiological analyses performed. In panel **B** the details and measurements of the spatial learning maze."

knitr::include_graphics("./Others/SPAL_METH.svg")

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
#| label: learning_results_names
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
#| fig.cap: "Learning results"
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
  summarize(se_int = sd(intercepts)/sqrt(length(intercepts)),
            intercepts = mean(intercepts))
# df fig-learning A
data_fig_learning_slopesA <- data_fig_learning_slopes %>%
  group_by(treatment) %>%
  summarize(se_slope = sd(slopes)/sqrt(length(slopes)),
            slopes = mean(slopes),)
data_fig_learningB <- merge(data_fig_learning_slopesB, data_fig_learning_intercepts, by = "treatment")
treatment <- unique(data_fig_learningB$treatment)
fig_A_df <- data.frame(trial = integer(),
                      treatment = character(),
                      errors = numeric(),
                      se = numeric())
for(t in treatment){ # Loop per treatment
  df <- data_fig_learningB %>%
    filter(treatment == t) %>%
    data.frame()
  # Variables selected
  m <- df$slopes
  u <- df$intercepts
  # Loop per treatment
  num_individuals <- length(u)
  for(x in 0:40){
    if (x == 0){
      standard <- df$se_int
    } else {
      standard <- df$se_slope
    }
    value <- exp(u + m * x)
    temp_df <- data.frame(trial = rep(x, length(value)), treatment = rep(t, length(value)), errors = value, se = standard)
  }
  fig_A_df <- rbind(fig_A_df, temp_df)
}
str(fig_A_df)
# Transform the database for better use
fig_choice_df <- fig_df %>%
  mutate(Trial = gsub("X", "", Trial)) %>%
  mutate(Trial = as.numeric(Trial)) %>%
  group_by(Trial, Treatment) %>%
  summarize(
    Mean_Predicted_prob = mean(Value),
    SE_Predicted_prob = sd(Value)/sqrt(length(Value))
    ) %>%
  ungroup() %>%
  mutate(
    Treatment = factor(Treatment,
                        levels = c("CORT-Cold", "Control-Cold", "CORT-Hot", "Control-Hot")),
    ) %>%
data.frame()
#### C) Plot the violing plot 
violin_plot_learning <- plot_slopes(data_fig_learning_slopes)

#
## A.2) Make the plot for the probability of choosing the correct feeder first over trials
fig_prob_choice <- ggplot(fig_choice_df, aes(x = Trial, y = Mean_Predicted_prob, color = Treatment)) +
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
## A.3) Make the the plot for the mean slope per treatment using the probability of choosing the correct feeder
fig_slope_choice <- plot_slopes("choice")
#
#
## A.4) Combine plots of probability and slopes to have the "choice" figure
fig_results_choice <- plot_grid(fig_prob_choice, fig_slope_choice, nrow = 1) 
ggsave("./output/figures/fig_results_choice.png", plot=fig_results_choice, width = 25, height = 15, units = "cm", dpi = 600)
#
#
#
####### B) Errors figure
## B.1) Plot the errors using data_spal
errors_df <- data_spal %>%
  mutate(Treatment = paste(cort, temp, sep = "-")) %>%
    mutate(Treatment = factor(Treatment,
      levels = c("CORT-Cold", "Control-Cold", "CORT-Hot", "Control-Hot")
    )) %>%
  group_by(day, Treatment) %>%
  summarise(
    mean_errors = mean(errors, na.rm = TRUE),
    sd_errors = sd(errors, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  data.frame()
fig_errors <- ggplot(learn_df, aes(x = Trial, y = Value, color = Treatment, fill = Treatment)) +
  geom_smooth(se = TRUE, linewidth = 1, alpha = 0.075) +
  scale_color_manual(values = c("CORT-Cold" = "#00008B", "Control-Cold" = "#68bde1", "CORT-Hot" = "#b50101", "Control-Hot" = "#fa927d")) +
  scale_fill_manual(values = c("CORT-Cold" = "darkblue", "Control-Cold" = "#68bde1", "CORT-Hot" = "#b50101", "Control-Hot" = "#fa927d")) +
  theme_classic() +
  labs(y = "Mean amount of errors made", x = "Trial") +
  theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5, "mm")) +
  theme(
    axis.title = element_text(size = 12, family = "Times"),
    axis.text = element_text(size = 10, family = "Times"),
    legend.position = "none"
  )
## B.2) Make the the plot for the mean slope per treatment using the probability of choosing the correct feeder
fig_slope_errors <- plot_slopes("errors")
#
#
## B.3) Combine plots of probability and slopes to have the "choice" figure
fig_results_errors <- plot_grid(fig_errors, fig_slope_errors, nrow = 1) 
ggsave("./output/figures/fig_results_errors.png", plot=fig_results_errors, width = 25, height = 15, units = "cm", dpi = 600)
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
#| label: table_bayesR2
#| tbl-cap: "BayesR2 values of the final models"
#
data_bayes <- data.frame(
  Model = character(0),
  Mean = numeric(0),
  Error = numeric(0),
  Q2_5 = numeric(0),
  Q97_5 = numeric(0)
)
models <- c("mean_mitodensity_def_OB",
            "mean_potential_def_OB",
            "mean_ros_def_OB",
            "mean_dnadamage_def_OB",
            "mean_peroxidation_def_OB",
            "t_D_def_Chemical",
            "mean_mitodensity_def_OT",
            "mean_potential_def_OT",
            "mean_ros_def_OT",
            "mean_dnadamage_def_OT",
            "mean_peroxidation_def_OT",
            "t_D_def_Visual")
#
for (m in models){
  mod <- readRDS(here("output/models/", paste0(m, ".rds")))
  bayes <- bayes_R2(mod)
  data_bayes <- rbind(data_bayes, data.frame(
    Model = m,
    Mean = format_dec(bayes[1], 3),
    Error = format_dec(bayes[2], 3),
    Q2_5 = format_dec(bayes[3], 3),
    Q97_5 = format_dec(bayes[4], 3)
  ))
}
#
bayes_table_df <- data_bayes %>%
  mutate(Region = gsub(".*_", "", Model)    # Extract everything after the last "_"
  ) %>%
  mutate(Model = factor(Model,
                        levels = c("mean_mitodensity_def_OB",
                                  "mean_potential_def_OB",
                                  "mean_ros_def_OB",
                                  "mean_dnadamage_def_OB",
                                  "mean_peroxidation_def_OB",
                                  "t_D_def_Chemical",
                                  "mean_mitodensity_def_OT",
                                  "mean_potential_def_OT",
                                  "mean_ros_def_OT",
                                  "mean_dnadamage_def_OT",
                                  "mean_peroxidation_def_OT",
                                  "t_D_def_Visual"),
                        labels = c("m_def_mean_mitodensity_OB" = "Mit density",
                                  "m_def_mean_potential_OB" = "Mit potential",
                                  "m_def_mean_ros_OB" = "ROS",
                                  "m_def_mean_dnadamage_OB" = "DNA damage",
                                  "m_def_mean_peroxidation_OB" = "Peroxidation",
                                  "m_def_t_D_Chemical" = "Detection lat",
                                  "m_def_mean_mitodensity_OT" = "Mit density",
                                  "m_def_mean_potential_OT" = "Mit potential",
                                  "m_def_mean_ros_OT" = "ROS",
                                  "m_def_mean_dnadamage_OT" = "DNA damage",
                                  "m_def_mean_peroxidation_OT" = "Peroxidation",
                                  "m_def_t_D_Visual" = "Detection lat")),
      Region = factor(Region, levels = c("OB", "Chemical", "OT", "Visual"),
                      labels = c("OB" = "Olfactory bulbs",
                                "Chemical" = "Chemical",
                                "OT" = "Optic tecta",
                                "Visual" = "Visual"))) %>%
  dplyr::select(Region, Model, Mean, Error, Q2_5, Q97_5) %>%
  arrange(Region, Model)
#
#
# Create the table
#
## Table format
set_flextable_defaults(
 font.family = "Times New Roman",
 font.size = 10)
#
bayes_table <- flextable(bayes_table_df) %>%
  set_header_labels(
    Region = "Region/Stimulus",
    Mean = "Mean",
    Error = "Error",
    Q2_5 = "2.5%",
    Q97_5 = "97.5%") %>%
  align(align = "center", j = c(3:5), part = "body") %>%
  align(align = "center", j = c(1:5), part = "header") %>%
  flextable::compose(i = c(2:5,8:11), j = 1, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the first column
  flextable::hline(i = 6, part = "body") %>% 
  autofit()
#
bayes_table
#
#
#
cat("\\newpage")
#
#
#
#
#| label: results_OB_table
#| tbl-cap: "Results of the models testing for Olfactory Bulbs."
#| tbl-name: "results_OB"
#| tbl-label: "results_OB"
source(here("R", "func.R"))
# 
# A) Refining the df summarizing the posteriors for OB/Chemical stimulus (post_OB)
post_OB_refined <- refine_post(post_OB) %>%
  arrange(Variable, Predictors)
#
# B) Create table
## Table format
set_flextable_defaults(
 font.family = "Times New Roman",
 font.size = 10)
#
OB_table <- flextable(post_OB_refined) %>%
  align(align = "center", j = c(3:5), part = "body") %>%
  align(align = "center", j = c(1:5), part = "header") %>%
  bold(~`PMCMC` < 0.05, j = c("PMCMC", "Estimate Mean", "95% CI", "Predictors"),
       bold = TRUE) %>%  # Bold when PMCMC is "<0.05"
  bold(~`PMCMC` <0.001, j = c("PMCMC", "Estimate Mean", "95% CI", "Predictors")) %>%  # Bold when PMCMC is "<0.001"
  flextable::compose(i = c(2:4,6:8,10:12,14:18,20:23,25:27), j = 1, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the first column
  autofit()
#
OB_table
#
#
#
cat("\\newpage")
#
#
#
#| label: tbl-contrasts
#| tbl-cap: "Estimates of Associative learning slope for all the different treatments per each task, species and group. Mean shows the aritmetic mean of the estimates obtained from the posteriors of the model, and 95% CI indicates the 95% confidence interval of the mean. All p-values were obtained using pmcmc and test the hypothesis that the mean is equal to zero. In bold, those values that are significant (p-value <0.05)"
source(here("R", "func.R"))
#
############################## CREATING BIG DF FOR TABLE ##############################
# Building the vectors for titles of rows and columns
specie <- c("L. delicata", "L. guichenoti")
groups <- c("Red", "Blue")
test <- c("Associative", "Reversal")
treatments <- c("CORT-Cold", "Control-Cold", "CORT-Hot", "Control-Hot")
values <- c("Mean", "95% CI", "p-value")
# Building the vectors for estimated means, co.intervals(95%), and p-values for the slopes obtained from posteriors. p-values are obtained using pmcmc function (see func.R), assuming a two-tailed test that testes the hypothesis that the value (slopes in this case) is 0.
#First get estimates for both tasks
estimates_asso <- list(
  dar_CORTCold, dar_ControlCold, dar_CORTHot, dar_ControlHot, 
  dab_CORTCold, dab_ControlCold, dab_CORTHot, dab_ControlHot, 
  gar_CORTCold, gar_ControlCold, gar_CORTHot, gar_ControlHot, 
  gab_CORTCold, gab_ControlCold, gab_CORTHot, gab_ControlHot
)
#
estimates_rev <- list(
  drr_CORTCold, drr_ControlCold, drr_CORTHot, drr_ControlHot, 
  drb_CORTCold, drb_ControlCold, drb_CORTHot, drb_ControlHot, 
  grr_CORTCold, grr_ControlCold, grr_CORTHot, grr_ControlHot, 
  grb_CORTCold, grb_ControlCold, grb_CORTHot, grb_ControlHot
)
# Then get the mean, co.intervals(95%), and p-values
asso_mean <- format_dec(sapply(estimates_asso, mean), 3)
asso_interval_025 <- format_dec(sapply(estimates_asso, function(x) quantile(x,0.025)), 3)
asso_interval_975 <- format_dec(sapply(estimates_asso, function(x) quantile(x,0.975)), 3)
asso_intervals <- paste(asso_interval_025, asso_interval_975, sep = " , ")
asso_pvalue <- format_dec(sapply(estimates_asso, pmcmc), 3)
#
rev_mean <- format_dec(sapply(estimates_rev, mean), 3)
rev_interval_025 <- format_dec(sapply(estimates_rev, function(x) quantile(x,0.025)), 3)
rev_interval_975 <- format_dec(sapply(estimates_rev, function(x) quantile(x,0.975)), 3)
rev_intervals <- paste(rev_interval_025, rev_interval_975, sep = " , ")
rev_pvalue <- format_dec(sapply(estimates_rev, pmcmc), 3)
#
# Building the df for the associative
asso_df <- data.frame(
  Specie = rep(specie, each = length(groups) * length(treatments)),
  Group = rep(rep(groups, each = length(treatments)), times = length(specie)),
  Treatment = rep(rep(treatments, each = 1), times = length(groups) * length(specie)),
  Mean = rep(asso_mean, each = 1),
  CI = rep(asso_intervals, each = 1),
  PValue = rep(asso_pvalue, each = 1),
  Task = rep("Associative", length(asso_mean))
)
# Building the df for the reversal
rev_df <- data.frame(
  Specie = rep(specie, each = length(groups) * length(treatments)),
  Group = rep(rep(groups, each = length(treatments)), times = length(specie)),
  Treatment = rep(rep(treatments, each = 1), times = length(groups) * length(specie)),
  Mean = rep(rev_mean, each = 1),
  CI = rep(rev_intervals, each = 1),
  PValue = rep(rev_pvalue, each = 1),
  Task = rep("Reversal", length(rev_mean))
)
# Joining both dfs
table_data <- rbind(asso_df, rev_df)
table_data[, sapply(table_data, is.numeric)] <- lapply(table_data[, sapply(table_data, is.numeric)], function(x) format(x, scientific = FALSE))
#
############################## ADDING SAMPLE SIZE TO DF FOR TABLE ##############################
# Make n_list into a df
n_df <- as.data.frame(do.call(rbind, n_list)) %>%
  rename("n" = V1) %>%
  rownames_to_column("model") %>%
  separate(model, into = c("Specie", "Group", "cort", "temp"), sep = "_") %>%
  unite("Treatment", c("cort", "temp"), sep = "-") %>%
  mutate(Specie = factor(Specie,
                  labels = c(delicata = "L. delicata", guichenoti = "L. guichenoti")),
        Treatment = factor(Treatment,
                   levels = c("CORT-Cold", "Control-Cold", "CORT-Hot","Control-Hot")))
# Merge both dfs, put sample size together with the treatment, and organize the new df to make it look like the table
new_table_data <- merge(table_data, n_df) %>%
  rename('p-value' = 'PValue', '95% CI' = 'CI') %>% #Change the names of the columns for the table
  pivot_wider(names_from = Task, values_from = c(Mean, `95% CI`, `p-value`)) %>% # to split between Asociative and Reversal
  select(Specie, Group, Treatment, Mean_Associative, `95% CI_Associative`, `p-value_Associative`, Mean_Reversal, `95% CI_Reversal`, `p-value_Reversal`, n) %>% #To order the columns in the way I want for the table
  mutate(Specie = factor(Specie,
                  levels = c("L. delicata", "L. guichenoti")),
        Group = factor(Group,
                  levels = c("Red", "Blue")),
        Treatment = factor(Treatment, 
                  levels = c("CORT-Cold", "Control-Cold", "CORT-Hot", "Control-Hot")))%>%
  arrange(Specie, Group, Treatment) %>% # To arrange the rows the way I want
  unite("Treatment", c("Treatment", "n"), sep = " (n = ") %>%
  mutate(Treatment = paste0(Treatment, ")"))
write.csv(new_table_data, file= "./output/Checking/new_table_data.csv")
#
############################## MAKING THE TABLE ##############################
## Table format
set_flextable_defaults(
 font.family = "Times New Roman",
 fint.size = 10)
# Split the table_data df by task
real_table <- flextable(new_table_data) %>%
    bold(~ `p-value_Associative` < 0.05, ~ `p-value_Associative` + Mean_Associative + `95% CI_Associative`) %>%
    bold(~ `p-value_Reversal` < 0.05, ~ `p-value_Reversal` + Mean_Reversal + `95% CI_Reversal`) %>%
    set_table_properties(width = 1) %>%
    align(align="center", part="all") %>% 
    add_header_row(values = c("", "Associative task", "Reversal task"), colwidths = c(3, 3, 3)) %>%
    set_header_labels(Mean_Associative = "Mean",
                      `95% CI_Associative` = "95% CI",
                      `p-value_Associative` = "p-value",
                      Mean_Reversal = "Mean",
                      `95% CI_Reversal` = "95% CI",
                      `p-value_Reversal` = "p-value") %>%
    italic(j = 1, italic = TRUE, part = "body") %>% # To have names od species in italics
    flextable::compose(i = c(2:8,10:16), j = 1, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the first column
    flextable::compose(i = c(2:4,6:8,10:12,14:16), j = 2, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the second column
    hline(i = c(4,12), j = c(2:9), part = "body") %>% # To make some horizontal lines
    hline(i = c(8), j = c(1:9), part = "body") %>% # To make some horizontal lines
    vline(i = (1:16), j = c(3,6), part = "body") %>% # To make some vertical lines on body
    vline(j=c(3,6), part = "header") %>% # To make some vertical lines on header
    autofit() 
real_table
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
#| tbl-cap: "Preliminary results of the models testing for Mitochondrial Potential"
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
#| label: results_preliminary_peroxidation
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
ggsave("./output/figures/suppl/Figure_SX3.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX3.png")
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
#| label: plotmod_mitdensity_OT
#| caption: "Posterior predictive checks for the model of Mitochondrial Density in Optic Tecta."
#
mod <- readRDS(here("output/models/mean_mitodensity_def_OT.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_SX4.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX4.png")
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
#| label: plotmod_potential_OB
#| caption: "Posterior predictive checks for the model of Mitochondrial Potential in Olfactory Bulbs."
#
mod <- readRDS(here("output/models/mean_potential_def_OB.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_SX5.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX5.png")
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
#| label: plotmod_potential_OT
#| caption: "Posterior predictive checks for the model of Mitochondrial Potential in Optic Tecta."
#
mod <- readRDS(here("output/models/mean_potential_def_OT.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_SX6.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX6.png")
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
#| label: plotmod_ros_OB
#| caption: "Posterior predictive checks for the model of ROS Production in Olfactory Bulbs."
#
mod <- readRDS(here("output/models/mean_ros_def_OB.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_SX7.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX7.png")
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
#| label: plotmod_ros_OT
#| caption: "Posterior predictive checks for the model of ROS Production in Optic Tecta."
#
mod <- readRDS(here("output/models/mean_ros_def_OT.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_SX8.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX8.png")
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
#| label: plotmod_dnadamage_OB
#| caption: "Posterior predictive checks for the model of DNA Damage in Olfactory Bulbs."
#
mod <- readRDS(here("output/models/mean_dnadamage_def_OB.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_SX9.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX9.png")
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
#| label: plotmod_dnadamage_OT
#| caption: "Posterior predictive checks for the model of DNA Damage in Optic Tecta."
#
mod <- readRDS(here("output/models/mean_dnadamage_def_OT.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_SX10.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX10.png")
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
#| label: plotmod_peroxidation_OB
#| caption: "Posterior predictive checks for the model of Lipid Peroxidation in Olfactory Bulbs."
#
mod <- readRDS(here("output/models/mean_peroxidation_def_OB.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_SX11.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX11.png")
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
#| label: plotmod_peroxidation_OT
#| caption: "Posterior predictive checks for the model of Lipid Peroxidation in Optic Tecta."
#
mod <- readRDS(here("output/models/mean_peroxidation_def_OT.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_SX12.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX12.png")
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
#| label: plotmod_tD_Chemical
#| caption: "Posterior predictive checks for the model of Detection Latency (t_D) in Chemical trials."
#
mod <- readRDS(here("output/models/t_D_def_Chemical.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_SX1.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX1.png")
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
#| label: plotmod_tD_Visual
#| caption: "Posterior predictive checks for the model of Detection Latency (t_D) in Visual trials."
#
mod <- readRDS(here("output/models/t_D_def_Visual.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_SX2.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX2.png")
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
