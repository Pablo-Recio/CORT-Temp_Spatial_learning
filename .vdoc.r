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
# A) Associative task. Model formula: FC_associative ~ Associative_Trial*cort*temp + (1 + Associative_Trial|lizard_id)
## A.1) L. delicata
deli_asso_red <- fit_m("asso", "deli", "red", refit = FALSE)
write.csv(deli_asso_red, file= "./output/Checking/deli_asso_red.csv")
deli_asso_blue <- fit_m("asso", "deli", "blue", refit = FALSE)
write.csv(deli_asso_blue, file= "./output/Checking/deli_asso_blue.csv")
## A.2) L. guichenoti
guich_asso_red <- fit_m("asso", "guich", "red", refit = FALSE)
write.csv(guich_asso_red, file= "./output/Checking/guich_asso_red.csv")
guich_asso_blue <- fit_m("asso", "guich", "blue", refit = FALSE)
write.csv(guich_asso_blue, file= "./output/Checking/guich_asso_blue.csv")
# B) Reversal task. Model formula: FC_reversal ~ trial_reversal*cort*temp + (1 + trial_reversal|lizard_id)
## B.1) L. delicata
deli_rev_red <- fit_m("rev", "deli", "red", refit = FALSE)
write.csv(deli_rev_red, file= "./output/Checking/deli_rev_red.csv")
deli_rev_blue <- fit_m("rev", "deli", "blue", refit = FALSE)
write.csv(deli_rev_blue, file= "./output/Checking/deli_rev_blue.csv")
## B.2) L. guichenoti
guich_rev_red <- fit_m("rev", "guich", "red", refit = FALSE)
write.csv(guich_rev_red, file= "./output/Checking/guich_rev_red.csv")
guich_rev_blue <- fit_m("rev", "guich", "blue", refit = FALSE)
write.csv(guich_rev_blue, file= "./output/Checking/guich_rev_blue.csv")
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
#| label: tbl-bias
#| tbl-cap: "Probability of choosing the correct feeder in the first trial when the correct feeder was blue (Prob Blue) or red (Prob Red) for each species and each treatment. p-values were obtained using pmcmc and test the hypothesis that the probability is >0.33. In bold, those values that are significant (p-value <0.05)"
source(here("R", "func.R"))
# First we estimate the probability of choosing right in the first trial using the intercepts from the posteriors
## 1) L. delicata
### Group = red
probright_drCORTCold <- 1/(1+exp(-deli_asso_red$b_Intercept))
probright_drControlCold <- 1/(1+exp(-(deli_asso_red$b_cortControl + deli_asso_red$b_Intercept)))
probright_drCORTHot <- 1/(1+exp(-(deli_asso_red$b_tempHot + deli_asso_red$b_Intercept)))
probright_drControlHot <- 1/(1+exp(-(deli_asso_red$'b_cortControl:tempHot' + deli_asso_red$b_cortControl + deli_asso_red$b_tempHot + deli_asso_red$b_Intercept)))
### Group = blue
probright_dbCORTCold <- 1/(1+exp(-(deli_asso_blue$b_Intercept)))
probright_dbControlCold <- 1/(1+exp(-(deli_asso_blue$b_cortControl + deli_asso_blue$b_Intercept)))
probright_dbCORTHot <- 1/(1+exp(-(deli_asso_blue$b_tempHot + deli_asso_blue$b_Intercept)))
probright_dbControlHot <- 1/(1+exp(-(deli_asso_blue$'b_cortControl:tempHot' + deli_asso_blue$b_cortControl + deli_asso_blue$b_tempHot + deli_asso_blue$b_Intercept)))
## 2) L. guichenoti
### Group = red
probright_grCORTCold <- 1/(1+exp(-guich_asso_red$b_Intercept))
probright_grControlCold <- 1/(1+exp(-(guich_asso_red$b_cortControl + guich_asso_red$b_Intercept)))
probright_grCORTHot <- 1/(1+exp(-(guich_asso_red$b_tempHot + guich_asso_red$b_Intercept)))
probright_grControlHot <- 1/(1+exp(-(guich_asso_red$'b_cortControl:tempHot' + guich_asso_red$b_cortControl + guich_asso_red$b_tempHot + guich_asso_red$b_Intercept)))
### Group = blue
probright_gbCORTCold <- 1/(1+exp(-guich_asso_blue$b_Intercept))
probright_gbControlCold <- 1/(1+exp(-(guich_asso_blue$b_cortControl + guich_asso_blue$b_Intercept)))
probright_gbCORTHot <- 1/(1+exp(-(guich_asso_blue$b_tempHot + guich_asso_blue$b_Intercept)))
probright_gbControlHot <- 1/(1+exp(-(guich_asso_blue$'b_cortControl:tempHot' + guich_asso_blue$b_cortControl + guich_asso_blue$b_tempHot + guich_asso_blue$b_Intercept)))
# Second we build a df with all the mean probabilities of choosing right in the first trial for each treatment and species
prob_df <- data.frame(
  probright_drCORTCold_t = format_dec(mean(probright_drCORTCold), 3),
  probright_drControlCold_t = format_dec(mean(probright_drControlCold), 3),
  probright_drCORTHot_t = format_dec(mean(probright_drCORTHot), 3),
  probright_drControlHot_t = format_dec(mean(probright_drControlHot), 3),
  probright_dbCORTCold_t = format_dec(mean(probright_dbCORTCold), 3),
  probright_dbControlCold_t = format_dec(mean(probright_dbControlCold), 3),
  probright_dbCORTHot_t = format_dec(mean(probright_dbCORTHot), 3),
  probright_dbControlHot_t = format_dec(mean(probright_dbControlHot), 3),
  probright_grCORTCold_t = format_dec(mean(probright_grCORTCold), 3),
  probright_grControlCold_t = format_dec(mean(probright_grControlCold), 3),
  probright_grCORTHot_t = format_dec(mean(probright_grCORTHot), 3),
  probright_grControlHot_t = format_dec(mean(probright_grControlHot), 3),
  probright_gbCORTCold_t = format_dec(mean(probright_gbCORTCold), 3),
  probright_gbControlCold_t = format_dec(mean(probright_gbControlCold), 3),
  probright_gbCORTHot_t = format_dec(mean(probright_gbCORTHot), 3),
  probright_gbControlHot_t = format_dec(mean(probright_gbControlHot), 3)) %>%
gather(key = "column", value = "value", .) %>%
rename("prob_choice" = "value") %>%
data.frame()
write.csv(prob_df, file= "./output/Checking/prob_df.csv")
# Then we build a df with the statistical test (one-tailed pmcmc) that the probability is >0.33 (the probability of choosing right by chance)
comp_pval <- data.frame(
  probright_drCORTCold_t = format_dec(pmcmc(probright_drCORTCold, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_drControlCold_t = format_dec(pmcmc(probright_drControlCold, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_drCORTHot_t = format_dec(pmcmc(probright_drCORTHot, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_drControlHot_t = format_dec(pmcmc(probright_drControlHot, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_dbCORTCold_t = format_dec(pmcmc(probright_dbCORTCold, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_dbControlCold_t = format_dec(pmcmc(probright_dbControlCold, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_dbCORTHot_t = format_dec(pmcmc(probright_dbCORTHot, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_dbControlHot_t = format_dec(pmcmc(probright_dbControlHot, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_grCORTCold_t = format_dec(pmcmc(probright_grCORTCold, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_grControlCold_t = format_dec(pmcmc(probright_grControlCold, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_grCORTHot_t = format_dec(pmcmc(probright_grCORTHot, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_grControlHot_t = format_dec(pmcmc(probright_grControlHot, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_gbCORTCold_t = format_dec(pmcmc(probright_gbCORTCold, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_gbControlCold_t = format_dec(pmcmc(probright_gbControlCold, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_gbCORTHot_t = format_dec(pmcmc(probright_gbCORTHot, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_gbControlHot_t = format_dec(pmcmc(probright_gbControlHot, null = 0.33, twotail = FALSE, dir=">"), 3)) %>%
gather(key = "column", value = "value", .) %>%
rename("p_value" = "value") %>%
data.frame()
write.csv(comp_pval, file= "./output/Checking/comp_pval.csv")
# Now we merge them and organise them to have the final dt for the table
table_df <- merge(prob_df, comp_pval, by ="column") %>%
  mutate(Specie = gsub(".*_(d|g).*", "\\1",column)) %>%
  mutate(Group = gsub(".*_(dr|db|gr|gb).*", "\\1",column)) %>%
  mutate(Treatment = gsub(".*(CORTCold|ControlCold|CORTHot|ControlHot).*", "\\1",column)) %>% #Last three commands used to get the species, group and treatment from the column names
  mutate(Treatment = factor(Treatment, 
                    levels = c("CORTCold", "ControlCold", "CORTHot", "ControlHot"),
                    labels = c("CORTCold" = "CORT-Cold", "ControlCold" = "Control-Cold", "CORTHot" = "CORT-Hot", "ControlHot" = "Control-Hot"))) %>%
  mutate(Specie = factor(Specie,
                levels = c("d", "g"),
                labels = c("d" = "L. delicata", "g" = "L. guichenoti"))) %>%
  mutate(Group = recode(Group, "dr" = "Red", "db" = "Blue", "gr" = "Red", "gb" = "Blue")) %>%
  select(-column) %>%
  pivot_wider(names_from = Group, values_from = c(prob_choice, p_value)) %>% #To split the df in Red and Blue
  select(Specie, Treatment, prob_choice_Blue, p_value_Blue, prob_choice_Red, p_value_Red) %>%
  arrange(Specie, Treatment) %>% # To arrange the rows the way I wantdata.frame()
data.frame()
write.csv(table_df, file= "./output/Checking/table_df.csv")
#################
colour_table <- flextable(table_df) %>%
  bold(~ `p_value_Blue` < 0.05, ~ `p_value_Blue` + prob_choice_Blue) %>%
  bold(~ `p_value_Red` < 0.05, ~ `p_value_Red` + prob_choice_Red) %>%
  set_table_properties(width = 1) %>%
  align(align="center", part="all") %>% 
  set_header_labels(prob_choice_Blue = "Prob Blue", 
                    p_value_Blue = "p-value Blue", 
                    prob_choice_Red = "Prob Red", 
                    p_value_Red = "p-value Red") %>%
  italic(j = 1, italic = TRUE, part = "body") %>% # To have names of species in italics
  flextable::compose(i = c(2:4,6:8), j = 1, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the first column
  hline(i = c(4), j = c(1:6), part = "body") %>% # To make some horizontal lines
  vline(j = c(2,4), part = "all") %>% # To make some vertical lines on body
autofit()
colour_table
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
cat("\\newpage")
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
