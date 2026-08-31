

##### Table 1 sensetivity analysis #######
#   pTau217 
mixed_model_ptau217 <- lmer(moca ~ years_from_baseline * plasma_ptau217_z + age_difference + sexe + diagnostic_nick + education + initiale_age +(1 | pscid),  
                            data = MRS_prediction_long)
summary(mixed_model_ptau217)

mixed_model_ptau217 <- lmer(moca ~ years_from_baseline * plasma_ptau217_z * age_difference + sexe + diagnostic_nick + education + initiale_age +(1 | pscid),  
                            data = MRS_prediction_long)
summary(mixed_model_ptau217)

mixed_model_ptau217 <- lmer(moca ~ years_from_baseline * plasma_ptau217_z + years_from_baseline * age_difference + sexe + diagnostic_nick + education + initiale_age +(1 | pscid),  
                            data = MRS_prediction_long)
summary(mixed_model_ptau217)



### Table 2 sensetivity analysis #####
# 1. Plasma p-tau217 (Adjusted)
model_glu_ptau217_cov <- glm(decliner_regression ~ plasma_ptau217_z + sexe + initiale_age, data = MRS_prediction, family = "binomial")
summary(model_glu_ptau217_cov)
roc_ptau217_cov <- roc(model_glu_ptau217_cov$y, fitted(model_glu_ptau217_cov))
auc(roc_ptau217_cov)
coords(roc_ptau217_cov, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_glu_ptau217_cov)

# 2. Precuneus Glutamate (Adjusted)
model_glu_prec_cov <- glm(decliner_regression ~ m_m_precuneus_z + sexe + initiale_age, data = MRS_prediction, family = "binomial")
summary(model_glu_prec_cov)
roc_glu_prec_cov <- roc(model_glu_prec_cov$y, fitted(model_glu_prec_cov))
auc(roc_glu_prec_cov)
coords(roc_glu_prec_cov, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_glu_prec_cov)

# 3. ACC Glutamate (Adjusted)
model_glu_acc_cov <- glm(decliner_regression ~ m_m_acc_z + sexe + initiale_age, data = MRS_prediction, family = "binomial")
summary(model_glu_acc_cov)
roc_glu_acc_cov <- roc(model_glu_acc_cov$y, fitted(model_glu_acc_cov))
auc(roc_glu_acc_cov)
coords(roc_glu_acc_cov, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_glu_acc_cov)

# 4. Cortical Thickness (Adjusted)
model_struc_thick_cov <- glm(decliner_regression ~ cortical_thickness_adsignature_dickson_z + sexe + initiale_age, data = MRS_prediction, family = "binomial")
summary(model_struc_thick_cov)
roc_struc_thick_cov <- roc(model_struc_thick_cov$y, fitted(model_struc_thick_cov))
auc(roc_struc_thick_cov)
coords(roc_struc_thick_cov, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_struc_thick_cov)

# 5. Hippocampal Activation (Adjusted)
model_func_hip_cov <- glm(decliner_regression ~ arsenii_hippocampus_avg_act + sexe + initiale_age, data = MRS_prediction, family = "binomial")
summary(model_func_hip_cov)
roc_func_hip_cov <- roc(model_func_hip_cov$y, fitted(model_func_hip_cov))
auc(roc_func_hip_cov)
coords(roc_func_hip_cov, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_func_hip_cov)
# ACC Glutamate + Plasma p-Tau217 + cortical thickness + covariates 
model_full_step_sig <- glm(decliner_regression ~ plasma_ptau217_z + m_m_acc_z + cortical_thickness_adsignature_dickson_z + sexe + initiale_age, data = MRS_prediction, family = "binomial")
summary(model_full_step_sig)
roc_model_full_step_sig <- roc(model_full_step_sig$y, fitted(model_full_step_sig))
auc(roc_model_full_step_sig)
coords(roc_model_full_step_sig, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_full_step_sig)













############ One leave out cross validation ##############
library(pROC)

# ==============================================================================
# 1. Master LOOCV Function (Computes true out-of-sample metrics)
# ==============================================================================
run_loocv <- function(formula_str, data, model_name) {
  form <- as.formula(formula_str)
  vars <- all.vars(form)
  
  # Extract relevant columns and handle missing data per model
  clean_data <- as.data.frame(data[, vars, drop = FALSE])
  clean_data <- clean_data[complete.cases(clean_data), ]
  
  n <- nrow(clean_data)
  cv_predictions <- numeric(n)
  
  # Leave-One-Out Loop
  for (i in 1:n) {
    train_set <- clean_data[-i, ]
    test_set  <- clean_data[i, ]
    
    fit <- glm(form, data = train_set, family = binomial())
    cv_predictions[i] <- predict(fit, newdata = test_set, type = "response")
  }
  
  # Out-of-sample ROC and optimal threshold (Youden index)
  actual_outcomes <- clean_data[[vars[1]]]
  roc_cv <- pROC::roc(actual_outcomes, cv_predictions, quiet = TRUE)
  coords_cv <- pROC::coords(roc_cv, "best", ret = c("threshold", "specificity", "sensitivity"), best.method = "youden")
  
  if (!is.null(nrow(coords_cv)) && nrow(coords_cv) > 1) {
    coords_cv <- coords_cv[1, ]
  }
  
  predicted_classes <- ifelse(cv_predictions >= coords_cv$threshold, 1, 0)
  acc_cv <- mean(predicted_classes == actual_outcomes)
  
  return(data.frame(
    Model          = model_name,
    N              = n,
    CV_AUC         = round(pROC::auc(roc_cv) * 100, 1),
    CV_Accuracy    = round(acc_cv * 100, 1),
    CV_Sensitivity = round(coords_cv$sensitivity * 100, 1),
    CV_Specificity = round(coords_cv$specificity * 100, 1),
    Optimal_Cutoff = round(coords_cv$threshold, 3)
  ))
}

# ==============================================================================
# 2. Define Your Exact Models
# ==============================================================================
models_to_test <- list(
  # --- Unimodal Models ---
  c("decliner_regression ~ plasma_ptau217_z", 
    "Unimodal: Plasma p-Tau217"),
  c("decliner_regression ~ m_m_precuneus_z", 
    "Unimodal: Precuneus Glutamate"),
  c("decliner_regression ~ m_m_acc_z", 
    "Unimodal: ACC Glutamate"),
  c("decliner_regression ~ cortical_thickness_adsignature_dickson_z", 
    "Unimodal: Cortical Thickness"),
  c("decliner_regression ~ arsenii_hippocampus_avg_act", 
    "Unimodal: Hippocampal Activation"),
  
  # --- Multimodal Model ---
  c("decliner_regression ~ plasma_ptau217_z + m_m_acc_z + cortical_thickness_adsignature_dickson_z", 
    "Multimodal: ACC Glu + p-Tau217 + Cortical Thickness")
)

# ==============================================================================
# 3. Run Validation and Display Summary Table
# ==============================================================================
loocv_results <- do.call(rbind, lapply(models_to_test, function(m) {
  run_loocv(formula_str = m[1], data = MRS_prediction, model_name = m[2])
}))

print(loocv_results, row.names = FALSE)







### Assumptions ####
library(gam)
install.packages("gam")
# 1. Plasma p-tau217
gam_ptau217 <- gam(decliner_regression ~ s(plasma_ptau217_z), 
                   data = MRS_prediction, 
                   family = "binomial")
summary(gam_ptau217)

# 2. Precuneus Glutamate
gam_prec <- gam(decliner_regression ~ s(m_m_precuneus_z), 
                data = MRS_prediction, 
                family = "binomial")
summary(gam_prec)

# 3. ACC Glutamate
gam_acc <- gam(decliner_regression ~ s(m_m_acc_z), 
               data = MRS_prediction, 
               family = "binomial")
summary(gam_acc)

# 4. Cortical Thickness (AD-Signature)
gam_thick <- gam(decliner_regression ~ s(cortical_thickness_adsignature_dickson_z), 
                 data = MRS_prediction, 
                 family = "binomial")
summary(gam_thick)

# 5. Hippocampal Activation
gam_hip <- gam(decliner_regression ~ s(arsenii_hippocampus_avg_act), 
               data = MRS_prediction, 
               family = "binomial")
summary(gam_hip)









names(MRS_prediction)
##### Sensetivty analysis of missing participants ########
# --- 1. Identify Included vs Excluded Participants ---
MRS_prediction$complete_case <- ifelse(
  complete.cases(MRS_prediction[, c("slope_regression_yearly", "m_m_acc_z", 
                                    "m_m_precuneus_z", "plasma_ptau217_z", 
                                    "cortical_thickness_adsignature_dickson_z", 
                                    "arsenii_hippocampus_avg_act", "sexe", 
                                    "diagnostic_nick", "education", "initiale_age", "age_difference")]),
  "Included", "Excluded"
)

# Check sample sizes (should be 48 vs 36)
table(MRS_prediction$complete_case)

# --- 2. Demographic & Baseline Comparisons ---
# Age comparison
cat("\n--- AGE ---\n")
t.test(initiale_age ~ complete_case, data = MRS_prediction)

# Education comparison
cat("\n--- EDUCATION ---\n")
t.test(education ~ complete_case, data = MRS_prediction)

# Sex comparison
cat("\n--- SEX ---\n")
chisq.test(table(MRS_prediction$sexe, MRS_prediction$complete_case))

# Clinical Diagnostic status comparison
cat("\n--- DIAGNOSIS (SCD vs MCI) ---\n")
chisq.test(table(MRS_prediction$diagnostic_nick, MRS_prediction$complete_case))
# --- 3. Cognitive Comparisons ---
# Baseline MoCA comparison
cat("\n--- BASELINE MOCA ---\n")
t.test(initiale_moca_score_total_30 ~ complete_case, data = MRS_prediction)

# Longitudinal cognitive slope comparison
cat("\n--- ANNUALIZED SLOPE ---\n")
t.test(slope_regression_yearly ~ complete_case, data = MRS_prediction)

# Clinically meaningful decliners proportion (Krishnan threshold)
cat("\n--- DECLINER STATUS (KRISHNAN) ---\n")
chisq.test(table(MRS_prediction$decliner_regression, MRS_prediction$complete_case))


# Get Means and SDs
aggregate(cbind(initiale_age, education, initiale_moca_score_total_30, slope_regression_yearly) ~ complete_case, data = MRS_prediction, FUN = function(x) c(mean = mean(x), sd = sd(x)))

# Get Counts for Categorical variables
table(MRS_prediction$sexe, MRS_prediction$complete_case)
table(MRS_prediction$diagnostic_nick, MRS_prediction$complete_case)
table(MRS_prediction$decliner_regression, MRS_prediction$complete_case)







# 1. Overall follow-up duration (Mean, SD, Median, Min, Max)
summary(MRS_prediction$max_years_from_baseline)
sd(MRS_prediction$max_years_from_baseline, na.rm = TRUE)

# 2. Check if follow-up duration differed between Stable and Decliners
t.test(max_years_from_baseline ~ decliner_regression, data = MRS_prediction)

# 3. Test the proportional hazards assumption (Schoenfeld residuals)
cox.zph(surv_acc)
cox.zph(surv_prec)

