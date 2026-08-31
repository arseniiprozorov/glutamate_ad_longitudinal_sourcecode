library(readxl)   
library(janitor)  
library(jmv)
library(effectsize)
library(interactions)
library(emmeans)
library(lme4)
library(lmerTest) 
library(pROC)
library(openxlsx)
library(survival)
library(olsrr)
options(na.action = "na.omit")

## Predicting cognitive change using metabolic, tau, functional and structural predictors  ######
## Arsenii Prozorov 
#ANALYSES PRÉLIMINAIRES :
#Création d’une banque de données
X2026_06_15_dataset_prediction <- read_excel("C:/Users/okkam/Desktop/labo/article 2/Longitudinal_Multimodal_Data_CIMAQ/article_prediction/2026-06-15_dataset_prediction.xlsx")
MRS_prediction <- X2026_06_15_dataset_prediction


# Clean the column name
MRS_prediction <- janitor::clean_names(MRS_prediction)
names(MRS_prediction)
sapply(MRS_prediction,class)
# Convertir  en   factor
MRS_prediction$sexe <- as.factor(MRS_prediction$sexe)
table(MRS_prediction$sexe)
MRS_prediction$diagnostic_nick <- as.factor(MRS_prediction$diagnostic_nick)
table(MRS_prediction$diagnostic_nick)
MRS_prediction$decliners <- as.factor(MRS_prediction$decliners)
table(MRS_prediction$decliners)
# Regrouper les groupes 'HC' et 'SCD' dans un seul groupe 'SCD'
levels(MRS_prediction$diagnostic_nick)[levels(MRS_prediction$diagnostic_nick) == "SCD"] <- "HC"
table(MRS_prediction$diagnostic_nick)


# Creer Hipp mean
MRS_prediction$hipp_mean <- (MRS_prediction$hip_l_nor_icv + MRS_prediction$righ_hip_vol)/2
MRS_prediction$hipp_mean_act <- (MRS_prediction$activation_hippocampus_l + MRS_prediction$activation_hippocampus_r)/2

#  Calculate the temporal gap/age difference between visits
MRS_prediction$age_difference <- abs(MRS_prediction$age_tau - MRS_prediction$initiale_age)

##################  Winzorising ###################
# Define the winsorize function
winsorize_iqr <- function(x, iqr_multiplier) {
  qnt <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- qnt[2] - qnt[1]
  lower <- qnt[1] - iqr_multiplier * iqr
  upper <- qnt[2] + iqr_multiplier * iqr
  x[x < lower] <- lower
  x[x > upper] <- upper
  return(x)}

# Apply winsorization for 1.5 IQR variables
MRS_prediction$memoria_libre_correcte <- winsorize_iqr(MRS_prediction$memoria_libre_correcte, 1.5)
MRS_prediction$face_name_rappel_differe_spectro <- winsorize_iqr(MRS_prediction$face_name_rappel_differe_spectro, 1.5)
## MRS Variables
MRS_prediction$m_m_precuneus <- winsorize_iqr(MRS_prediction$m_m_precuneus, 1.5)
MRS_prediction$m_m_acc <- winsorize_iqr(MRS_prediction$m_m_acc, 1.5)
## Plasma Variables
MRS_prediction$plasma_ptau217 <- winsorize_iqr(MRS_prediction$plasma_ptau217, 1.5)
## Structural MRI Variables
MRS_prediction$hip_l_nor_icv <- winsorize_iqr(MRS_prediction$hip_l_nor_icv, 1.5)
MRS_prediction$cortical_thickness_adsignature_dickson <- winsorize_iqr(MRS_prediction$cortical_thickness_adsignature_dickson, 1.5)
MRS_prediction$hipp_mean <- winsorize_iqr(MRS_prediction$hipp_mean, 1.5)

## fMRI Activation Variables
MRS_prediction$activation_hippocampus_l <- winsorize_iqr(MRS_prediction$activation_hippocampus_l, 1.5)
MRS_prediction$hipp_mean_act <- winsorize_iqr(MRS_prediction$hipp_mean_act, 1.5)
MRS_prediction$arsenii_hippocampus_avg_act <- winsorize_iqr(MRS_prediction$arsenii_hippocampus_avg_act, 1.5)
MRS_prediction$activation_parietal_sup_l <- winsorize_iqr(MRS_prediction$activation_parietal_sup_l, 1.5)
MRS_prediction$arsenii_parietal_sup_l_act <- winsorize_iqr(MRS_prediction$arsenii_parietal_sup_l_act, 1.5)
MRS_prediction$activation_temporal_inf_r <- winsorize_iqr(MRS_prediction$activation_temporal_inf_r, 1.5)



# Descriptives
#sink()
jmv::descriptives(data = MRS_prediction, vars = vars(slope_regression_yearly, m_m_precuneus, m_m_acc, plasma_ptau217,
                                                     cortical_thickness_adsignature_dickson, 
                                                     hipp_mean,hipp_mean_act, arsenii_parietal_sup_l_act),
                  sd = TRUE, iqr = TRUE, skew = TRUE, kurt = TRUE)








#### Calcul de score z
MRS_prediction$m_m_precuneus_z <- scale(MRS_prediction$m_m_precuneus)
MRS_prediction$m_m_acc_z <- scale(MRS_prediction$m_m_acc)
MRS_prediction$plasma_ptau217_z <- scale(MRS_prediction$plasma_ptau217)
MRS_prediction$cortical_thickness_adsignature_dickson_z <- scale(MRS_prediction$cortical_thickness_adsignature_dickson)
MRS_prediction$hipp_mean_z <- scale(MRS_prediction$hipp_mean)
MRS_prediction$activation_hippocampus_l_z <- scale(MRS_prediction$activation_hippocampus_l)
MRS_prediction$hipp_mean_act_z <- scale(MRS_prediction$hipp_mean_act)
MRS_prediction$arsenii_parietal_sup_l_act <- scale(MRS_prediction$arsenii_parietal_sup_l_act)


## Long format 
X2026_06_16_MRS_prediction_long <- read_excel("C:/Users/okkam/Desktop/MRS_prediction_longitudinal_master.xlsx")
MRS_prediction_long <- X2026_06_16_MRS_prediction_long
lapply(MRS_prediction_long,class)
names(MRS_prediction_long)

MRS_prediction_long$arsenii_hippocampus_avg_act <- winsorize_iqr(MRS_prediction_long$arsenii_hippocampus_avg_act, 1.5)
MRS_prediction_long$arsenii_parietal_sup_l_act <- winsorize_iqr(MRS_prediction_long$arsenii_parietal_sup_l_act, 1.5)

# Extract moca slopes for each participant 
raw_individual_models <- lmList(moca ~ years_from_baseline | pscid, data = MRS_prediction_long)
raw_slopes <- coef(raw_individual_models)
raw_slopes$pscid <- rownames(raw_slopes)
raw_slopes$moca_change_3_5_yrs <- raw_slopes$years_from_baseline * 3.5


######################################### Analyses #######################################
names(MRS_prediction)

######### Characterization ###############
#sink("Table 1.txt")
# Table 1: Sociodemographic & Clinical Characteristics
# Table 1.1  Split by Decliners Status 
jmv::descriptives(data = MRS_prediction, vars = vars(initiale_age, education, slope_regression_yearly, age_difference),
                  sd = TRUE, iqr = TRUE, splitBy = decliners, skew = TRUE, kurt = TRUE)

t.test(initiale_age ~ decliners, data = MRS_prediction)

chisq.test(MRS_prediction$sexe, y = MRS_prediction$decliners, correct = TRUE)
table(MRS_prediction$sexe, MRS_prediction$decliners)

t.test(education ~ decliners, data = MRS_prediction)

chisq.test(MRS_prediction$diagnostic_nick, y = MRS_prediction$decliners, correct = TRUE)
table(MRS_prediction$diagnostic_nick, MRS_prediction$decliners)


t.test(slope_regression_yearly ~ decliners, data = MRS_prediction)

t.test(age_difference ~ decliners, data = MRS_prediction)


#  Table 1.2 - Split by Baseline Diagnosis 
jmv::descriptives(data = MRS_prediction, vars = vars(initiale_age, education, slope_regression_yearly, age_difference),
                  sd = TRUE, iqr = TRUE, splitBy = diagnostic_nick, skew = TRUE, kurt = TRUE)

summary(aov(initiale_age ~ diagnostic_nick, data = MRS_prediction))

chisq.test(MRS_prediction$sexe, y = MRS_prediction$diagnostic_nick, correct = TRUE)
table(MRS_prediction$sexe, MRS_prediction$diagnostic_nick)

summary(aov(education ~ diagnostic_nick, data = MRS_prediction))
summary(aov(slope_regression_yearly ~ diagnostic_nick, data = MRS_prediction))
summary(aov(age_difference ~ diagnostic_nick, data = MRS_prediction))



# Table 2: Biomarkers & Neuroimaging Modalities

# Table 2.1 - Split by Decliners Status
jmv::descriptives(data = MRS_prediction, 
                  vars = vars(m_m_acc, m_m_precuneus, plasma_ptau217,
                               hipp_mean, cortical_thickness_adsignature_dickson,
                              arsenii_hippocampus_avg_act_z, arsenii_parietal_sup_l_act),
                  sd = TRUE, iqr = TRUE, splitBy = decliners, skew = TRUE, kurt = TRUE)

# T-tests for Decliners vs Non-Decliners
t.test(m_m_acc ~ decliners, data = MRS_prediction)
t.test(m_m_precuneus ~ decliners, data = MRS_prediction)
t.test(plasma_ptau217 ~ decliners, data = MRS_prediction)
t.test(hipp_mean ~ decliners, data = MRS_prediction)
t.test(cortical_thickness_adsignature_dickson ~ decliners, data = MRS_prediction)
t.test(arsenii_hippocampus_avg_act_z ~ decliners, data = MRS_prediction)
t.test(arsenii_parietal_sup_l_act ~ decliners, data = MRS_prediction)



#  Split by Baseline Diagnosis
#  Descriptive Statistics 
jmv::descriptives(data = MRS_prediction, 
                  vars = vars(m_m_acc, m_m_precuneus, plasma_ptau217,
                              hipp_mean, cortical_thickness_adsignature_dickson,
                              arsenii_hippocampus_avg_act_z, arsenii_parietal_sup_l_act),
                  sd = TRUE, iqr = TRUE, splitBy = diagnostic_nick, skew = TRUE, kurt = TRUE)

# ANOVAs & Post-Hoc Tests for Diagnostic Groups (Harmonized to match descriptives)
anova_acc_diagnostick <- aov(m_m_acc ~ diagnostic_nick, data = MRS_prediction)
summary(anova_acc_diagnostick)
TukeyHSD(anova_acc_diagnostick)

anova_prec_diagnostick <- aov(m_m_precuneus ~ diagnostic_nick, data = MRS_prediction)
summary(anova_prec_diagnostick)
TukeyHSD(anova_prec_diagnostick)

anova_ptau_diagnostick <- aov(plasma_ptau217 ~ diagnostic_nick, data = MRS_prediction)
summary(anova_ptau_diagnostick)
TukeyHSD(anova_ptau_diagnostick)

anova_hipp_vol_diagnostick <- aov(hipp_mean ~ diagnostic_nick, data = MRS_prediction)
summary(anova_hipp_vol_diagnostick)
TukeyHSD(anova_hipp_vol_diagnostick)

summary(aov(cortical_thickness_adsignature_dickson ~ diagnostic_nick, data = MRS_prediction))
summary(aov(arsenii_hippocampus_avg_act_z ~ diagnostic_nick, data = MRS_prediction))
summary(aov(arsenii_parietal_sup_l_act ~ diagnostic_nick, data = MRS_prediction))

#sink()

############################ Objective 1 #######################
names(MRS_prediction_long)
names(MRS_prediction)

citation("lme4")
citation("lmerTest")
############## Mixed-Effects Models ###############
names(MRS_prediction_long)
levels(MRS_prediction_long$diagnostic_nick)

#sink("Table 2.txt")
## Baseline model
mixed_model_moca <- lmer(moca ~ years_from_baseline + sexe + diagnostic_nick + education + initiale_age + (1 | pscid),  
                         data = MRS_prediction_long)
summary(mixed_model_moca)

#   pTau217 
mixed_model_ptau217 <- lmer(moca ~ years_from_baseline * plasma_ptau217_z + age_difference + sexe + diagnostic_nick + education + initiale_age +(1 | pscid),  
                            data = MRS_prediction_long)
summary(mixed_model_ptau217)

#tau217_slopes <- emtrends(mixed_model_ptau217, specs = ~ plasma_ptau217_z, 
#                       var = "years_from_baseline", 
#                       at = list(plasma_ptau217_z = c(-1, 0, 1)))
#summary(tau217_slopes, infer = TRUE)

#anova(mixed_model_ptau217)
#   Glutamate 
mixed_model_precuneus <- lmer(moca ~  m_m_precuneus_z * years_from_baseline + sexe + diagnostic_nick + education + (1 | pscid),  
  data = MRS_prediction_long)
summary(mixed_model_precuneus)

#precuneus_slopes <- emtrends(mixed_model_precuneus, specs = ~ m_m_precuneus_z, 
#                            var = "years_from_baseline", 
#                            at = list(m_m_precuneus_z = c(-1,0,1)))
#summary(precuneus_slopes, infer = TRUE)


#### ACC ##########


mixed_model_acc <- lmer(moca ~ years_from_baseline * m_m_acc_z + sexe + diagnostic_nick + education + initiale_age + (1 | pscid),  
  data = MRS_prediction_long)
summary(mixed_model_acc)

#acc_slopes <- emtrends(mixed_model_acc, specs = ~ m_m_acc_z, 
#                             var = "years_from_baseline", 
#                             at = list(m_m_acc_z = c(-1,0,1)))
#summary(acc_slopes, infer = TRUE)
#anova(mixed_model_acc)





# Structure
mixed_model_thickness <- lmer(
  moca ~ years_from_baseline * cortical_thickness_adsignature_dickson_z + sexe + diagnostic_nick + education + initiale_age + (1 | pscid),  
  data = MRS_prediction_long)
summary(mixed_model_thickness)

#thick_slopes <- emtrends(mixed_model_thickness, specs = ~ cortical_thickness_adsignature_dickson_z, 
#                          var = "years_from_baseline", 
#                          at = list(cortical_thickness_adsignature_dickson_z = c(-1, 0, 1)))
#summary(thick_slopes, infer = TRUE)
#anova(mixed_model_thickness)


mixed_model_hipp_mean <- lmer(moca ~ years_from_baseline * hipp_mean_z  + sexe + diagnostic_nick + education + initiale_age + (1 | pscid),  
  data = MRS_prediction_long)
summary(mixed_model_hipp_mean)

#hipp_vol_slopes <- emtrends(mixed_model_hipp_mean, specs = ~ hipp_mean_z, 
#                         var = "years_from_baseline", 
#                         at = list(hipp_mean_z = c(-1, 1)))
#summary(hipp_vol_slopes, infer = TRUE)
#anova(mixed_model_hipp_mean)




# Activation 
# Hipp
#mixed_model_hipp_mean_act <- lmer(moca ~ years_from_baseline * arsenii_hippocampus_avg_act + + sexe + diagnostic_nick + education + initiale_age + (1 | pscid),  
#  data = MRS_prediction_long)
#summary(mixed_model_hipp_mean_act)

mixed_model_hipp_mean_act <- lmer(moca ~ years_from_baseline * arsenii_hippocampus_avg_act + sexe + diagnostic_nick + education + initiale_age + (1 | pscid),  
                                  data = MRS_prediction_long)
summary(mixed_model_hipp_mean_act)


#hipp_act_slopes <- emtrends(mixed_model_hipp_mean_act, specs = ~ arsenii_hippocampus_avg_act, 
#                        var = "years_from_baseline", 
#                       at = list(hipp_mean_act_z = c(-1, 0, 1)))
#summary(hipp_act_slopes, infer = TRUE)


mixed_model_activation_parietal_l <- lmer(moca ~ years_from_baseline * arsenii_parietal_sup_l_act + + sexe + diagnostic_nick + education +
                                            initiale_age + (1 | pscid),data = MRS_prediction_long)
summary(mixed_model_activation_parietal_l)




########################################### The Multimodal  Model###############################
# Syntax Isaora
#ols_step_backward_p(OLS_mod0,prem=0.1) 
help(package = olsrr)
citation("olsrr")
########### backward stepwise lienar mixed regression ################

### Backward stepwise ###

MRS_step <- na.omit(MRS_prediction[, c("slope_regression_yearly", "m_m_acc_z", 
                                             "m_m_precuneus_z", "plasma_ptau217_z", 
                                             "cortical_thickness_adsignature_dickson_z", 
                                             "arsenii_hippocampus_avg_act", "sexe", 
                                             "diagnostic_nick", "education", "initiale_age", "age_difference")])

#  Full multivariable model & backward selection
full_lm_model <- lm(slope_regression_yearly ~ m_m_acc_z + m_m_precuneus_z + plasma_ptau217_z + 
                      cortical_thickness_adsignature_dickson_z + arsenii_hippocampus_avg_act + 
                      sexe + diagnostic_nick + education + initiale_age + age_difference, 
                    data = MRS_step)

step_result <- ols_step_backward_p(full_lm_model, p_val = 0.25, details = TRUE)


# Output
print(step_result)
summary(step_result$model)


step_aic_result <- ols_step_backward_aic(full_lm_model, details = TRUE)
print(step_aic_result)
########### Post hoc with tertials ########

# Create Tertile Factor Variables (Low, Medium, High)

# Plasma p-Tau217
q_ptau <- quantile(MRS_prediction_long$plasma_ptau217_z, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
MRS_prediction_long$ptau217_tert <- cut(MRS_prediction_long$plasma_ptau217_z, breaks = q_ptau, 
                                        labels = c("Low", "Medium", "High"), include.lowest = TRUE)

# ACC Glutamate
q_acc <- quantile(MRS_prediction_long$m_m_acc_z, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
MRS_prediction_long$acc_tert <- cut(MRS_prediction_long$m_m_acc_z, breaks = q_acc, 
                                    labels = c("Low", "Medium", "High"), include.lowest = TRUE)

# Precuneus Glutamate
q_prec <- quantile(MRS_prediction_long$m_m_precuneus_z, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
MRS_prediction_long$precuneus_tert <- cut(MRS_prediction_long$m_m_precuneus_z, breaks = q_prec, 
                                          labels = c("Low", "Medium", "High"), include.lowest = TRUE)

# Cortical Thickness (AD-Signature)
q_thick <- quantile(MRS_prediction_long$cortical_thickness_adsignature_dickson_z, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
MRS_prediction_long$thickness_tert <- cut(MRS_prediction_long$cortical_thickness_adsignature_dickson_z, breaks = q_thick, 
                                          labels = c("Low", "Medium", "High"), include.lowest = TRUE)

# Hippocampal Activation
q_hip <- quantile(MRS_prediction_long$arsenii_hippocampus_avg_act, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
MRS_prediction_long$hipp_act_tert <- cut(MRS_prediction_long$arsenii_hippocampus_avg_act, breaks = q_hip, 
                                         labels = c("Low", "Medium", "High"), include.lowest = TRUE)

#  Plasma p-Tau217 Tertiles
mod_ptau_tert <- lmer(moca ~ years_from_baseline * ptau217_tert + age_difference + sexe + 
                        diagnostic_nick + education + initiale_age + (1 | pscid), data = MRS_prediction_long)

slopes_ptau <- emtrends(mod_ptau_tert, ~ ptau217_tert, var = "years_from_baseline")
summary(slopes_ptau, infer = TRUE)
pairs(slopes_ptau)

#  Glutamate Tertiles
mod_acc_tert <- lmer(moca ~ years_from_baseline * acc_tert + sexe + diagnostic_nick + 
                       education + initiale_age + (1 | pscid), data = MRS_prediction_long)

slopes_acc <- emtrends(mod_acc_tert, ~ acc_tert, var = "years_from_baseline")
summary(slopes_acc, infer = TRUE)
pairs(slopes_acc)

names(MRS_prediction_long)
#  Precuneus Glutamate Tertiles

mod_prec_tert <- lmer(moca ~ years_from_baseline * precuneus_tert + sexe + diagnostic_nick + 
                        education + (1 | pscid), data = MRS_prediction_long)

slopes_prec <- emtrends(mod_prec_tert, ~ precuneus_tert, var = "years_from_baseline")
summary(slopes_prec, infer = TRUE)
pairs(slopes_prec)


#  AD-Signature Cortical Thickness Tertiles
mod_thick_tert <- lmer(moca ~ years_from_baseline * thickness_tert + sexe + diagnostic_nick + 
                         education + initiale_age + (1 | pscid), data = MRS_prediction_long)

slopes_thick <- emtrends(mod_thick_tert, ~ thickness_tert, var = "years_from_baseline")
summary(slopes_thick, infer = TRUE)
pairs(slopes_thick)


# Hippocampal Activation Tertiles
mod_hip_tert <- lmer(moca ~ years_from_baseline * hipp_act_tert + sexe + diagnostic_nick + 
                       education + initiale_age + (1 | pscid), data = MRS_prediction_long)

# Slopes & Pairwise Tests
slopes_hip <- emtrends(mod_hip_tert, ~ hipp_act_tert, var = "years_from_baseline")
summary(slopes_hip, infer = TRUE)
pairs(slopes_hip)


sink()


################# Objective 2.  Logistic regression ######################
names(MRS_prediction)
## Overal model sig
overall_m_sig <- function(model) {
  chsq <- model$null.deviance - model$deviance
  df <- model$df.null - model$df.residual
  pval <- pchisq(chsq, df, lower.tail = FALSE)
  return(data.frame(Chi_Sq = chsq, df = df, p_value = pval))}
#overall_m_sig(object)


#sink("Table 3.txt")

# plasma_ptau217 
model_glu_ptau217 <- glm(decliner_regression ~ plasma_ptau217_z  , data = MRS_prediction, family = "binomial")
summary(model_glu_ptau217)
roc_ptau217 <- roc(model_glu_ptau217$y, fitted(model_glu_ptau217))
auc(roc_ptau217)
coords(roc_ptau217, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_glu_ptau217)


# Precuneus Glutamate
model_glu_prec <- glm(decliner_regression ~ m_m_precuneus_z , data = MRS_prediction, family = "binomial")
summary(model_glu_prec)
roc_glu_prec <- roc(model_glu_prec$y, fitted(model_glu_prec))
auc(roc_glu_prec)
coords(roc_glu_prec, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_glu_prec)

# ACC Glutamate
model_glu_acc <- glm(decliner_regression ~ m_m_acc_z  , data = MRS_prediction, family = "binomial")
summary(model_glu_acc)
roc_glu_acc <- roc(model_glu_acc$y, fitted(model_glu_acc))
auc(roc_glu_acc)
coords(roc_glu_acc, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_glu_acc)


# Cortical Thickness
model_struc_thick <- glm(decliner_regression ~ cortical_thickness_adsignature_dickson_z  , data = MRS_prediction, family = "binomial")
summary(model_struc_thick)
roc_struc_thick <- roc(model_struc_thick$y, fitted(model_struc_thick))
auc(roc_struc_thick)
coords(roc_struc_thick, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_struc_thick)

# Hippocampal Activation
model_func_hip <- glm(decliner_regression ~ arsenii_hippocampus_avg_act  , data = MRS_prediction, family = "binomial")
summary(model_func_hip)
roc_func_hip <- roc(model_func_hip$y, fitted(model_func_hip))
auc(roc_func_hip)
coords(roc_func_hip, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_func_hip)



# Hippocampal Volume
#model_struc_hip <- glm(decliner_regression ~ hipp_mean_z  , data = MRS_prediction, family = "binomial")
#summary(model_struc_hip)
#roc_struc_hip <- roc(model_struc_hip$y, fitted(model_struc_hip))
#auc(roc_struc_hip)
#coords(roc_struc_hip, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
#overall_m_sig(model_struc_hip)

# Superior Parietal Activation
#model_func_par <- glm(decliner_regression ~ activation_parietal_sup_l_z , data = MRS_prediction, family = "binomial")
#summary(model_func_par)
#roc_func_par <- roc(model_func_par$y, fitted(model_func_par))
#auc(roc_func_par)
#coords(roc_func_par, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
#overall_m_sig(model_func_par)


############ Multimodal model ############

#  ACC Glutamate + Plasma p-Tau217 + cortical thickenss 
model_nocov_step_sig <- glm(decliner_regression ~ plasma_ptau217_z + m_m_acc_z + cortical_thickness_adsignature_dickson_z 
                           , data = MRS_prediction, family = "binomial")
summary(model_nocov_step_sig)

roc_model_nocov_step_sig <- roc(model_nocov_step_sig$y, fitted(model_nocov_step_sig))
auc(roc_model_nocov_step_sig)
coords(roc_model_nocov_step_sig, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_nocov_step_sig)





################# # Objerctive 3 #######################
## Sruvival analysis #####
names(MRS_prediction)
names(MRS_prediction_long)

# Syntax #
# Fit the Cox Proportional Hazards model
#cox_model <- coxph(Surv(time, status) ~ treatment, data = my_data)
# Extract Hazard Ratios (exp(coef)) and 95% Confidence Intervals
#summary(cox_model)
#exp(confint(cox_model))
#Surv(time, status): Defines the survival time and the event indicator (e.g., 1 for event, 0 for censored).
#exp(coef): Represents the Hazard Ratio.If HR = 1, the risk is equal between groups.If HR = 1.5, the event rate is 50% higher at any given moment.

surv_acc <- coxph(Surv(max_years_from_baseline, decliner_regression) ~ m_m_acc_z  + initiale_age, data = MRS_prediction)
summary(surv_acc)
exp(confint(surv_acc))
surv_prec <- coxph(Surv(max_years_from_baseline, decliner_regression) ~ m_m_precuneus_z + initiale_age, data = MRS_prediction)
summary(surv_prec)
exp(confint(surv_prec))

summary(coxph(Surv(max_years_from_baseline, decliner_regression) ~ plasma_ptau217_z + initiale_age, data = MRS_prediction))
#summary(coxph(Surv(max_years_from_baseline, decliner_regression) ~ hipp_mean + initiale_age, data = MRS_prediction))
summary(coxph(Surv(max_years_from_baseline, decliner_regression) ~ cortical_thickness_adsignature_dickson_z + initiale_age, data = MRS_prediction))
summary(coxph(Surv(max_years_from_baseline, decliner_regression) ~ arsenii_hippocampus_avg_act + initiale_age, data = MRS_prediction))
#summary(coxph(Surv(max_years_from_baseline, decliner_regression) ~ activation_parietal_sup_l + initiale_age, data = MRS_prediction))


# full models
cox_full <- coxph(Surv(max_years_from_baseline, decliner_regression) ~ plasma_ptau217_z + m_m_acc_z + cortical_thickness_adsignature_dickson_z 
              + sexe + initiale_age, data = MRS_prediction)
summary(cox_full)
cox_nocov <- coxph(Surv(max_years_from_baseline, decliner_regression) ~ plasma_ptau217_z + m_m_acc_z + cortical_thickness_adsignature_dickson_z 
              , data = MRS_prediction)
summary(cox_nocov)



# split into tertiales 
ptau217_tertiles <- quantile(MRS_prediction$plasma_ptau217_z, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
acc_terials <- quantile(MRS_prediction$m_m_acc_z, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
precuneus_terials <- quantile(MRS_prediction$m_m_precuneus_z, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
thick_tertiles <- quantile(MRS_prediction$cortical_thickness_adsignature_dickson_z, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
hip_act_tertiles <- quantile(MRS_prediction$arsenii_hippocampus_avg_act, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)

# Plasma p-Tau217
MRS_prediction$ptau217_tertiales <- cut(MRS_prediction$plasma_ptau217_z, 
                                        breaks = ptau217_tertiles, labels = c("Low p-Tau217", "Medium p-Tau217", "High p-Tau217"), 
                                        include.lowest = TRUE)

#  acc
MRS_prediction$acc_tertiales <- cut(MRS_prediction$m_m_acc_z, 
                                    breaks = acc_terials, 
                                    labels = c("Low ACC Glu", "Medium ACC Glu", "High ACC Glu"), 
                                    include.lowest = TRUE)

#  precuneus
MRS_prediction$precuneus_tertiales <- cut(MRS_prediction$m_m_precuneus_z, 
                                          breaks = precuneus_terials, 
                                          labels = c("Low Precuneus Glu", "Medium Precuneus Glu", "High Precuneus Glu"), 
                                          include.lowest = TRUE)


# Cortical Thickness (AD-Signature)
MRS_prediction$cortical_thick_tertiales <- cut(MRS_prediction$cortical_thickness_adsignature_dickson_z, 
                                               breaks = thick_tertiles, labels = c("Low Thickness", "Medium Thickness", "High Thickness"), 
                                               include.lowest = TRUE)

# Hippocampal Activation
MRS_prediction$hip_act_tertiales <- cut(MRS_prediction$arsenii_hippocampus_avg_act, 
                                        breaks = hip_act_tertiles, labels = c("Low Hipp Act", "Medium Hipp Act", "High Hipp Act"), 
                                        include.lowest = TRUE)



survfit(Surv(max_years_from_baseline, decliner_regression) ~ acc_tertiales, data = MRS_prediction)

survfit(Surv(max_years_from_baseline, decliner_regression) ~ precuneus_tertiales, data = MRS_prediction)


# Plasma p-Tau217 Tertiles
#fit_ptau217 <- survfit(
#  Surv(max_years_from_baseline, decliner_regression) ~ ptau217_tertiales, 
#  data = MRS_prediction)
#(fit_ptau217)

# Cortical Thickness (AD-Signature) Tertiles
#fit_thick <- survfit(
#  Surv(max_years_from_baseline, decliner_regression) ~ cortical_thick_tertiales, 
#  data = MRS_prediction)
#summary(fit_thick)

# Hippocampal Activation Tertiles
#fit_hip_act <- survfit(Surv(max_years_from_baseline, decliner_regression) ~ hip_act_tertiales, 
#  data = MRS_prediction)
#summary(fit_hip_act)




### Multimodal Keplen Mayer curve
MRS_prediction$multimodal_risk_score <- predict(cox_nocov, newdata = MRS_prediction, type = "lp", na.action = na.pass)

MRS_prediction$multimodal_risk_tertiles <- cut(
  MRS_prediction$multimodal_risk_score,
  breaks = quantile(MRS_prediction$multimodal_risk_score, probs = 0:3/3, na.rm = TRUE),
  labels = c("Low Risk Profile", "Intermediate Risk Profile", "High Risk Profile"),
  include.lowest = TRUE)

# Fit and summarize Kaplan-Meier model
fit_multimodal <- survfit(Surv(max_years_from_baseline, decliner_regression) ~ multimodal_risk_tertiles, data = MRS_prediction)
summary(fit_multimodal)




#sink()

