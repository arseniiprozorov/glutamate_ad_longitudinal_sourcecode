library(readxl)   
library(janitor)  
library(jmv)
library(effectsize)
library(interactions)
library(emmeans)

## The Sequential Dynamics of Glutamate and Neurodegeneration in Cognitive Decline  ######
## Arsenii Prozorov 
#ANALYSES PRÉLIMINAIRES :
#Création d’une banque de données

binned_cascade_prediciton <- read_excel("C:/Users/okkam/Desktop/binned_cascade_prediciton.xlsx")
MRS_long <- binned_cascade_prediciton

# Clean the column name
MRS_long <- janitor::clean_names(MRS_long)
names(MRS_long)
sapply(MRS_long,class)



# Force all variables except participant_id to be numeric
MRS_long <- MRS_long %>%
  mutate(across(-participant_id, as.numeric))

# Check your work to ensure they are now "numeric"
sapply(MRS_long, class)






library(dplyr)
library(tidyr)
library(pROC)

# FORCE R TO FORGET THE OLD FUNCTION DEFINITION
if (exists("evaluate_model")) rm(evaluate_model)

# 1. Create a consolidated dataset (Average when both exist, pull single if one is missing)
MRS_combined <- MRS_long %>%
  rowwise() %>% 
  mutate(
    cortical_thickness_base = mean(c(norm_cortical_thickness_dickerson_t_0, norm_cortical_thickness_dickerson_t_pre), na.rm = TRUE),
    sum_hippocampus_base = mean(c(norm_sum_hippocampus_t_0, norm_sum_hippocampus_t_pre), na.rm = TRUE),
    hippocampus_avg_act_base = mean(c(hippocampus_avg_act_t_0, hippocampus_avg_act_t_pre), na.rm = TRUE),
    parietal_sup_l_act_base = mean(c(parietal_sup_l_act_t_0, parietal_sup_l_act_t_pre), na.rm = TRUE),
    glu_precuneus_base = mean(c(glu_precuneus_t_0, glu_precuneus_t_pre), na.rm = TRUE),
    glu_acc_base = mean(c(glu_acc_t_0, glu_acc_t_pre), na.rm = TRUE),
    moca_score_total_base = mean(c(moca_score_total_30_t_0, moca_score_total_30_t_pre), na.rm = TRUE),
    ptau217_base = mean(c(x99453_analyse_plasma_ptau217_t_0, x99453_analyse_plasma_ptau217_t_pre), na.rm = TRUE)
  ) %>%
  ungroup() 

# 2. Define the consolidated predictors and target
target_var <- "decliners"
predictors <- c(
  "cortical_thickness_base", "sum_hippocampus_base", 
  "hippocampus_avg_act_base", "parietal_sup_l_act_base", 
  "glu_precuneus_base", "glu_acc_base", 
  "moca_score_total_base", "ptau217_base"
)

# 3. Redefine the evaluation function with distinct count columns
evaluate_model <- function(vars, data) {
  formula_str <- paste(target_var, "~", paste(vars, collapse = " + "))
  
  # Filter for complete cases on the specifically tested variables
  df_clean <- data %>% 
    select(all_of(c(target_var, vars))) %>% 
    drop_na() %>%
    filter(if_all(all_of(vars), ~ !is.nan(.))) 
  
  n_total <- nrow(df_clean)
  
  if(n_total == 0 || length(unique(df_clean[[target_var]])) < 2) {
    return(NULL)
  }
  
  # Count the true positive events (decliners == 1) in this subset
  n_decliners <- sum(df_clean[[target_var]] == 1)
  
  model <- glm(as.formula(formula_str), data = df_clean, family = binomial)
  pred_probs <- predict(model, type = "response")
  
  roc_obj <- roc(df_clean[[target_var]], pred_probs, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  
  coords_best <- coords(roc_obj, "best", ret = c("specificity", "sensitivity"), 
                        best.method = "youden", transpose = FALSE)
  
  spec <- coords_best$specificity[1]
  sens <- coords_best$sensitivity[1]
  
  return(data.frame(
    Model = paste(vars, collapse = " + "),
    N_Total = n_total,
    N_Decliners = n_decliners, 
    AUC = round(auc_val, 4),
    Sensitivity = round(sens, 4),
    Specificity = round(spec, 4)
  ))
}

# 4. Run the Loops
cat("Evaluating single variables...\n")
results_1var <- bind_rows(lapply(predictors, function(v) evaluate_model(v, MRS_combined))) %>%
  arrange(desc(AUC))

cat("Evaluating 2-variable combinations...\n")
combos_2 <- combn(predictors, 2, simplify = FALSE)
results_2var <- bind_rows(lapply(combos_2, function(v) evaluate_model(v, MRS_combined))) %>%
  arrange(desc(AUC))

cat("Evaluating 3-variable combinations...\n")
combos_3 <- combn(predictors, 3, simplify = FALSE)
results_3var <- bind_rows(lapply(combos_3, function(v) evaluate_model(v, MRS_combined))) %>%
  arrange(desc(AUC))

# 5. Output the Results
cat("\n--- Top Baseline Single Predictors ---\n")
print(head(results_1var, 10))

cat("\n--- Top Baseline Two-Variable Predictors ---\n")
print(head(results_2var, 10))

cat("\n--- Top Baseline Three-Variable Predictors ---\n")
print(head(results_3var, 10))