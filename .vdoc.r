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
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes, ggh4x, PupillometryR, cowplot, png, grid)
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

```{r, fig-Methods}
#| label: fig-Methods
#| fig.cap: "Panel **A** shows the experimental design of the study. On top, the treatments applied to the eggs. On the bottom the learning task and the brain region extracted together with the physiological analyses performed. In panel **B** the details and measurements of the spatial learning maze."

knitr::include_graphics("./Others/SPAL_METH.svg")

#
#
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
#| label: data_processing
# The result will be the final df with the data for the analysis. To do so, we first estimate the learning slope (choice and errors) for each individual and then merge the data with the mitochondrial data, extracted using the script in extract.R and extraction_finc.R (see R folder). The final df will be saved in (here("output/databases_clean/data_complete.csv") 
refit = FALSE
source(here("R", "data_process.R"))
#
#
#
#| label: sampleSize
# List with the sample sizes from the database (here("output/databases_clean/data_complete.csv") as the sample size per species and group is the same on each task. We used function sample (see func.R) to estimate the sample size per treatment and species.
source(here("R", "func.R"))
#
specie <- c("delicata", "guichenoti")
groups <- c("Red", "Blue")
hormone <- c("CORT", "Control")
temperature <- c("Cold", "Hot")
#
n_list <- list()
#
for(i in 1:length(specie)){
  for(j in 1:length(groups)){
    for(k in 1:length(hormone)){
      for(l in 1:length(temperature)){
        n_list[[paste(specie[i], groups[j], hormone[k], temperature[l], sep = "_")]] <- sample(specie[i], groups[j], hormone[k], temperature[l])
      }
    }
  }
}
#
#
#
#
#
#| label: models
# Fitting the model and extraction of posteriors for both types of task and species using fit_m function (see func.r in R folder). The result everytime the function is used is a df with the posteriors of the model. The functions saves the model automatically in output/models; and when the parameter refit = FALSE then the posteriors are extracted from the model previously written instead of fitting the model again each time.
source(here("R", "func.R"))
# 
# df got from data_process.R
#
fchoice <- (choice_slope__mean + choice_slope__sd ~ cort*temp)
ferrors <- (errors_slope__mean + errors_slope__sd ~ cort*temp)
fROS <- 
#
#
#
# Rename some of the posteriors and make new estimates for the learning rate for the Reversal task doing the same thing we did in the chunk above.
## 1) L. delicata
### Group = red
drr_CORTCold <- deli_rev_red$b_trial_reversal
drr_ControlCold <- (deli_rev_red$'b_trial_reversal:cortControl' + deli_rev_red$b_trial_reversal)
drr_CORTHot <- (deli_rev_red$'b_trial_reversal:tempHot' + deli_rev_red$b_trial_reversal)
drr_ControlHot <- (deli_rev_red$'b_trial_reversal:cortControl:tempHot' + deli_rev_red$b_trial_reversal + deli_rev_red$'b_trial_reversal:cortControl' + deli_rev_red$'b_trial_reversal:tempHot')
### Group = blue
drb_CORTCold <- deli_rev_blue$b_trial_reversal
drb_ControlCold <- (deli_rev_blue$'b_trial_reversal:cortControl' + deli_rev_blue$b_trial_reversal)
drb_CORTHot <- (deli_rev_blue$'b_trial_reversal:tempHot' + deli_rev_blue$b_trial_reversal)
drb_ControlHot <- (deli_rev_blue$'b_trial_reversal:cortControl:tempHot' + deli_rev_blue$b_trial_reversal + deli_rev_blue$'b_trial_reversal:cortControl' + deli_rev_blue$'b_trial_reversal:tempHot')
## 2) L. guichenoti
### Group = red
grr_CORTCold <- guich_rev_red$b_trial_reversal
grr_ControlCold <- (guich_rev_red$'b_trial_reversal:cortControl' + guich_rev_red$b_trial_reversal)
grr_CORTHot <- (guich_rev_red$'b_trial_reversal:tempHot' + guich_rev_red$b_trial_reversal)
grr_ControlHot <- (guich_rev_red$'b_trial_reversal:cortControl:tempHot' + guich_rev_red$b_trial_reversal + guich_rev_red$'b_trial_reversal:cortControl' + guich_rev_red$'b_trial_reversal:tempHot')
### Group = blue
grb_CORTCold <- guich_rev_blue$b_trial_reversal
grb_ControlCold <- (guich_rev_blue$'b_trial_reversal:cortControl' + guich_rev_blue$b_trial_reversal)
grb_CORTHot <- (guich_rev_blue$'b_trial_reversal:tempHot' + guich_rev_blue$b_trial_reversal)
grb_ControlHot <- (guich_rev_blue$'b_trial_reversal:cortControl:tempHot' + guich_rev_blue$b_trial_reversal + guich_rev_blue$'b_trial_reversal:cortControl' + guich_rev_blue$'b_trial_reversal:tempHot')
#
#
#
#| label: tbl-data
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
#| label: fig-learning
#| fig.cap: "Results for L. delicata (A,B) and L. guichenoti (C, D). Panels A and C show the predicted probability of choosing the correct feeder first over trials. The lines represent the mean predicted probability of choosing the correct feeder first on each trial, and the shaded areas indicate the standard deviation of the mean; both were obtained by using the slope and intercept estimates from the posterior distributions. The different colours indicate the different treatments. Panels B and D show the distribution of the estimates of slopes per each treatment. The x-axis represents the slope estimate, and in the y-axis are the density of the estimates. The different colours indicate the different treatments. Points and bars represent the mean and standard deviation of the mean of the estimates, respectively."
#
#
img <- readPNG("./Others/Deli.png")
####### A) Probability of choosing correct figure
## A.1) Get the df with the estimated probability of choosing the correct feeder first over trials (using learning_df from data_process.R)
# Get the df for the intercepts
source(here("R", "func.R"))
post_choice_intercept <- tidy_post_df("choice", "intercept")
# Merge slope and intercept dfs for choice 
learning_df_prob <- merge(learning_df, post_choice_intercept, by = "lizard_id", all = TRUE)
# Create a vector with treatments
treatments <- c("CORT-Cold", "Control-Cold", "CORT-Hot", "Control-Hot")
# Create new matrix and new dfs
treat_matrix <- matrix(NA, nrow = 20, ncol = 40)
fig_df <- data.frame()
# Loop through treatments
for(t in treatments){
  if(t == "CORT-Cold"){
    df <- learning_df_prob %>%
      filter(cort=="CORT", temp=="Cold") %>%
      data.frame()
  } else if(t == "Control-Cold"){
    df <- learning_df_prob %>%
      filter(cort=="Control", temp=="Cold") %>%
      data.frame()
  } else if(t == "CORT-Hot"){
    df <- learning_df_prob %>%
      filter(cort=="CORT", temp=="Hot") %>%
      data.frame()
  } else if(t == "Control-Hot"){
    df <- learning_df_prob %>%
      filter(cort=="Control", temp=="Hot") %>%
      data.frame()
  } else {
    stop("loop wrong")
  }
  # Variables selected
  m <- df$choice_slope__mean
  u <- df$choice_intercept__mean
  # Loop per treatment
  num_individuals <- length(u)
  for(x in 0:40){
    for(j in 1:num_individuals){
      value <- exp(u[j] + m[j] * x) / (1 + exp(u[j] + m[j] * x))
      treat_matrix[j, x] <- value
    }
  }
  treat_df <- as.data.frame(treat_matrix)
  colnames(treat_df) <- paste0("X", 1:40)  # Adjust column names
  treat_df <- gather(treat_df, key = "Trial", value = "Value")  # Reshape data frame
  treat_df$Treatment <- t
  fig_df <- rbind(fig_df, treat_df)
}
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
fig_results_choice <- plot_grid(fig_prob_choice, fig_slope_choice, nrow = 1) +
    annotation_custom(rasterGrob(img), xmin = 0.73, xmax = 0.98, ymin = 0.73, ymax = 0.98)
ggsave("./output/figures/fig_results_choice.png", plot=fig_results_choice, width = 25, height = 18, units = "cm", dpi = 600)
knitr::include_graphics("./output/figures/fig_results.png")
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
fig_errors <- ggplot(errors_df, aes(x = day, y = mean_errors, color = Treatment, fill = Treatment)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1, alpha = 0.075) +
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
fig_results_errors <- plot_grid(fig_errors, fig_slope_errors, nrow = 1) +
    annotation_custom(rasterGrob(img), xmin = 0.73, xmax = 0.98, ymin = 0.73, ymax = 0.98)
ggsave("./output/figures/fig_results_errors.png", plot=fig_results_errors, width = 25, height = 18, units = "cm", dpi = 600)
knitr::include_graphics("./output/figures/fig_results.png")
#
#
#
#| label: ROS-learning
#
plot_ROS_choice <- ggplot(data = clean_df, aes(x = arith.mean_ROS, y = log(choice_slope__mean))) +
  geom_smooth(method = "lm")
print(plot_ROS_choice)
#
plot_ROS_errors <- ggplot(data = clean_df, aes(x = arith.mean_ROS, y = log(errors_slope__mean))) +
  geom_smooth(method = "lm")
print(plot_ROS_errors)
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
# Chunk for calling all the models and making the residuals
# Associative task
## L. delicata
### Red
mod_dar <- readRDS(here("output/models/asso_deli_red.rds"))
resid_dar <- residuals(mod_dar)
### Blue
mod_dab <- readRDS(here("output/models/asso_deli_blue.rds"))
resid_dab <- residuals(mod_dab)
## L. guichenoti
### Red
mod_gar <- readRDS(here("output/models/asso_guich_red.rds"))
resid_gar <- residuals(mod_gar)
### Blue
mod_gab <- readRDS(here("output/models/asso_guich_blue.rds"))
resid_gab <- residuals(mod_gab)
# Reversal task
## L. delicata
### Red
mod_drr <- readRDS(here("output/models/rev_deli_red.rds"))
resid_drr <- residuals(mod_drr)
### Blue
mod_drb <- readRDS(here("output/models/rev_deli_blue.rds"))
resid_drb <- residuals(mod_drb)
## L. guichenoti
### Red
mod_grr <- readRDS(here("output/models/rev_guich_red.rds"))
resid_grr <- residuals(mod_grr)
### Blue
mod_grb <- readRDS(here("output/models/rev_guich_blue.rds"))
resid_grb <- residuals(mod_grb)
#
#
#
#
#
#
#
#
plot(mod_drr)
#
#
#
cat("\\newpage")
#
#
#
#
plot(mod_drb)
#
#
#
cat("\\newpage")
#
#
#
#
#
plot(mod_grr)
#
#
#
cat("\\newpage")
#
#
#
#
plot(mod_grb)
#
#
#
