
# Moca slope as continous lm

summary(lm(slope_regression_yearly ~ m_m_precuneus , data = MRS_prediction))
summary(lm(slope_regression_yearly ~ m_m_acc, data = MRS_prediction))
summary(lm(slope_regression_yearly ~ m_m_precuneus + I(m_m_precuneus^2) , data = MRS_prediction))
summary(lm(slope_regression_yearly ~ m_m_acc + I(m_m_precuneus^2), data = MRS_prediction))
moca_precuneus_quadratic <- lm(initiale_moca_score_total_30 ~ m_m_precuneus + I(m_m_precuneus^2) , data = MRS_prediction)
AIC(moca_precuneus_quadratic)
moca_precuneus_linear <- lm(initiale_moca_score_total_30 ~ m_m_precuneus, data = MRS_prediction)
AIC(moca_precuneus_linear)
aic_dff <- AIC(moca_precuneus_quadratic) - AIC(moca_precuneus_linear)
aic_dff
summary(lm(initiale_moca_score_total_30 ~ m_m_acc + I(m_m_acc^2) , data = MRS_prediction))
summary(lm(initiale_moca_score_total_30 ~ m_m_acc , data = MRS_prediction))

summary(lm(slope_regression_yearly ~ plasma_ptau217, data = MRS_prediction))
summary(lm(slope_regression_yearly ~ cortical_thickness_adsignature_dickson, data = MRS_prediction))
summary(lm(slope_regression_yearly ~ hipp_mean, data = MRS_prediction))
summary(lm(slope_regression_yearly ~ arsenii_hippocampus_avg_act, data = MRS_prediction))
summary(lm(slope_regression_yearly ~ arsenii_Parietal_Sup_L_act, data = MRS_prediction))

citation("lme4")
citation("lmerTest")







#### Hierarchichal regression ######## (blockwise)
#Block 1: Covariates
#Block 2: Glutamate (Does it predict MoCA on its own?)
#Block 3: Add Activation (Does activation add predictive value above and beyond glutamate?)
#Block 4: Add Thickness.

# 1. Create your clean dataset (stepwise functions still require no missing data)
vars_to_keep <- c("moca", "years_from_baseline", 
                  "m_m_precuneus_z", "m_m_acc_z", "plasma_ptau217_z", 
                  "cortical_thickness_adsignature_dickson_z", "arsenii_hippocampus_avg_act", 
                  "age_difference", "sexe", "diagnostic_nick", "education", "initiale_age", "pscid")

MRS_clean <- MRS_prediction_long[, vars_to_keep]
MRS_clean <- MRS_clean[complete.cases(MRS_clean), ]

# 2. Fit the FULL mixed-effects model
# Note: Place your interacting variables in parentheses multiplied by years_from_baseline
full_mixed_model <- lmer(
  moca ~ years_from_baseline * (m_m_precuneus_z + 
                                  m_m_acc_z + 
                                  plasma_ptau217_z + 
                                  cortical_thickness_adsignature_dickson_z + 
                                  arsenii_hippocampus_avg_act) + 
    age_difference + sexe + diagnostic_nick + education + initiale_age + 
    (1 | pscid), 
  data = MRS_clean
)

# 3. Run the backward stepwise elimination on the lmer model
# reduce.random = FALSE tells it to leave your (1 | pscid) alone and only eliminate fixed effects
step_result <- step(full_mixed_model, reduce.random = FALSE)

# 4. Print the elimination log to see exactly which variables were dropped and in what order
print(step_result)

# 5. Extract the final "winning" model into a new object
final_lmer_model <- get_model(step_result)

# 6. View the summary of your final model
summary(final_lmer_model)

library(lme4)

# 1. Define all variables
vars <- c("moca", "years_from_baseline", "sexe", "diagnostic_nick", "education", 
          "initiale_age", "age_difference", "plasma_ptau217_z", "m_m_acc_z", 
          "m_m_precuneus_z", "arsenii_hippocampus_avg_act", 
          "cortical_thickness_adsignature_dickson_z", "pscid")

# 2. Create the complete dataset
MRS_complete <- MRS_prediction_long[, vars]
MRS_complete <- MRS_complete[complete.cases(MRS_complete), ]

# Block 0: Clinical Baseline
mod_base <- lmer(moca ~ years_from_baseline + sexe + diagnostic_nick + 
                   education + initiale_age + age_difference + 
                   (1 | pscid), 
                 data = MRS_complete, REML = FALSE)

# Block 1: Add Core Pathology (p-tau217)
mod_tau <- lmer(moca ~ years_from_baseline * plasma_ptau217_z + 
                  sexe + diagnostic_nick + education + initiale_age + age_difference + 
                  (1 | pscid), 
                data = MRS_complete, REML = FALSE)

summary(mod_tau)
# Block 2: Add Early Metabolic (Glutamate)
mod_glutamate <- lmer(moca ~ years_from_baseline * plasma_ptau217_z + 
                        years_from_baseline * m_m_acc_z + 
                        years_from_baseline * m_m_precuneus_z + 
                        sexe + diagnostic_nick + education + initiale_age + age_difference + 
                        (1 | pscid), 
                      data = MRS_complete, REML = FALSE)

summary(mod_glutamate)
# Block 3: Add Functional Response (Hippocampal Hyperactivation)
mod_act <- lmer(moca ~ years_from_baseline * plasma_ptau217_z + 
                  years_from_baseline * m_m_acc_z + 
                  years_from_baseline * m_m_precuneus_z + 
                  years_from_baseline * arsenii_hippocampus_avg_act + 
                  sexe + diagnostic_nick + education + initiale_age + age_difference + 
                  (1 | pscid), 
                data = MRS_complete, REML = FALSE)

summary(mod_act)
# Block 4: Add Downstream Neurodegeneration (Cortical Thickness)
mod_thick <- lmer(moca ~ years_from_baseline * plasma_ptau217_z + 
                    years_from_baseline * m_m_acc_z + 
                    years_from_baseline * m_m_precuneus_z + 
                    years_from_baseline * arsenii_hippocampus_avg_act + 
                    years_from_baseline * cortical_thickness_adsignature_dickson_z + 
                    sexe + diagnostic_nick + education + initiale_age + age_difference + 
                    (1 | pscid), 
                  data = MRS_complete, REML = FALSE)
summary(mod_thick)

# Compare the hierarchy
anova(mod_base, mod_tau, mod_glutamate, mod_act, mod_thick)
