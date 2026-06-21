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
jmv::descriptives(data = MRS_prediction, vars = vars(m_m_precuneus, m_m_acc, plasma_ptau217,
                                                     cortical_thickness_adsignature_dickson, 
                                                     hipp_mean,hipp_mean_act,activation_parietal_sup_l),
                  sd = TRUE, iqr = TRUE, skew = TRUE, kurt = TRUE, splitBy = converter_to_mci)




#### Calcul de score z
MRS_prediction$m_m_precuneus_z <- scale(MRS_prediction$m_m_precuneus)
MRS_prediction$m_m_acc_z <- scale(MRS_prediction$m_m_acc)
MRS_prediction$plasma_ptau217_z <- scale(MRS_prediction$plasma_ptau217)
MRS_prediction$cortical_thickness_adsignature_dickson_z <- scale(MRS_prediction$cortical_thickness_adsignature_dickson)
MRS_prediction$hipp_mean_z <- scale(MRS_prediction$hipp_mean)
MRS_prediction$activation_hippocampus_l_z <- scale(MRS_prediction$activation_hippocampus_l)
MRS_prediction$hipp_mean_act_z <- scale(MRS_prediction$hipp_mean_act)
MRS_prediction$activation_parietal_sup_l_z <- scale(MRS_prediction$activation_parietal_sup_l)


## Long format 
X2026_06_16_MRS_prediction_long <- read_excel("C:/Users/okkam/Desktop/MRS_prediction_longitudinal_master.xlsx")
MRS_prediction_long <- X2026_06_16_MRS_prediction_long
lapply(MRS_prediction_long,class)

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
#sink("demographic.txt")
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
                              hipp_mean_act, activation_parietal_sup_l),
                  sd = TRUE, iqr = TRUE, splitBy = decliners, skew = TRUE, kurt = TRUE)

# T-tests for Decliners vs Non-Decliners
t.test(m_m_acc ~ decliners, data = MRS_prediction)
t.test(m_m_precuneus ~ decliners, data = MRS_prediction)
t.test(plasma_ptau217 ~ decliners, data = MRS_prediction)
t.test(hipp_mean ~ decliners, data = MRS_prediction)
t.test(cortical_thickness_adsignature_dickson ~ decliners, data = MRS_prediction)
t.test(hipp_mean_act ~ decliners, data = MRS_prediction)
t.test(activation_parietal_sup_l ~ decliners, data = MRS_prediction)



#  Split by Baseline Diagnosis
#  Descriptive Statistics 
jmv::descriptives(data = MRS_prediction, 
                  vars = vars(m_m_acc, m_m_precuneus, plasma_ptau217,
                              hipp_mean, cortical_thickness_adsignature_dickson,
                              hipp_mean_act, activation_parietal_sup_l),
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
summary(aov(hipp_mean_act ~ diagnostic_nick, data = MRS_prediction))
summary(aov(activation_parietal_sup_l ~ diagnostic_nick, data = MRS_prediction))

#sink()

############################ Objective 1 #######################
names(MRS_prediction_long)
names(MRS_prediction)


#sink("objective_1_outputs.txt")
# Moca slope as continous

summary(lm(slope_regression_yearly ~ m_m_precuneus , data = MRS_prediction))
summary(lm(slope_regression_yearly ~ m_m_acc, data = MRS_prediction))

summary(lm(slope_regression_yearly ~ plasma_ptau217, data = MRS_prediction))
summary(lm(slope_regression_yearly ~ cortical_thickness_adsignature_dickson, data = MRS_prediction))
summary(lm(slope_regression_yearly ~ hipp_mean, data = MRS_prediction))
summary(lm(slope_regression_yearly ~ arsenii_hippocampus_avg_act, data = MRS_prediction))
summary(lm(slope_regression_yearly ~ parietal_sup_l_act, data = MRS_prediction))


## same with LMM
names(MRS_prediction_long)

mixed_model_moca <- lmer(moca ~ years_from_baseline + sexe + diagnostic_nick + education + initiale_age + (1 | pscid),  
                         data = MRS_prediction_long)
summary(mixed_model_moca)

# Mixed-Effects Models
#   Glutamate 
mixed_model_precuneus <- lmer(moca ~ years_from_baseline * m_m_precuneus_z + sexe + diagnostic_nick + education + initiale_age + (1 | pscid),  
  data = MRS_prediction_long)
summary(mixed_model_precuneus)

precuneus_slopes <- emtrends(mixed_model_precuneus, specs = ~ m_m_precuneus_z, 
                            var = "years_from_baseline", 
                            at = list(m_m_precuneus_z = c(-1, 0, 1)))
summary(precuneus_slopes, infer = TRUE)
anova(mixed_model_precuneus)



mixed_model_acc <- lmer(moca ~ years_from_baseline * m_m_acc_z + sexe + diagnostic_nick + education + initiale_age + (1 | pscid),  
  data = MRS_prediction_long)
summary(mixed_model_acc)

acc_slopes <- emtrends(mixed_model_acc, specs = ~ m_m_acc_z, 
                             var = "years_from_baseline", 
                             at = list(m_m_acc_z = c(-1, 0, 1)))
summary(acc_slopes, infer = TRUE)
anova(mixed_model_acc)


#   pTau217 
mixed_model_ptau217 <- lmer(moca ~ years_from_baseline * plasma_ptau217_z + age_difference + sexe + diagnostic_nick + education + initiale_age +(1 | pscid),  
  data = MRS_prediction_long)
summary(mixed_model_ptau217)

tau217_slopes <- emtrends(mixed_model_ptau217, specs = ~ plasma_ptau217_z, 
                       var = "years_from_baseline", 
                       at = list(plasma_ptau217_z = c(-1, 0, 1)))
summary(tau217_slopes, infer = TRUE)
anova(mixed_model_ptau217)

# Structure
mixed_model_thickness <- lmer(
  moca ~ years_from_baseline * cortical_thickness_adsignature_dickson_z + sexe + diagnostic_nick + education + initiale_age + (1 | pscid),  
  data = MRS_prediction_long)
summary(mixed_model_thickness)

thick_slopes <- emtrends(mixed_model_thickness, specs = ~ cortical_thickness_adsignature_dickson_z, 
                          var = "years_from_baseline", 
                          at = list(cortical_thickness_adsignature_dickson_z = c(-1, 0, 1)))
summary(thick_slopes, infer = TRUE)
anova(mixed_model_thickness)


mixed_model_hipp_mean <- lmer(moca ~ years_from_baseline * hipp_mean_z  + sexe + diagnostic_nick + education + initiale_age + (1 | pscid),  
  data = MRS_prediction_long)
summary(mixed_model_hipp_mean)

hipp_vol_slopes <- emtrends(mixed_model_hipp_mean, specs = ~ hipp_mean_z, 
                         var = "years_from_baseline", 
                         at = list(hipp_mean_z = c(-1, 0, 1)))
summary(hipp_vol_slopes, infer = TRUE)
anova(mixed_model_hipp_mean)




names(MRS_prediction_long)
# Activation 
# Hipp
#mixed_model_hipp_mean_act <- lmer(moca ~ years_from_baseline * arsenii_hippocampus_avg_act + + sexe + diagnostic_nick + education + initiale_age + (1 | pscid),  
#  data = MRS_prediction_long)
#summary(mixed_model_hipp_mean_act)

mixed_model_hipp_mean_act <- lmer(moca ~ years_from_baseline * hipp_mean_act + sexe + diagnostic_nick + education + initiale_age + (1 | pscid),  
                                  data = MRS_prediction_long)
summary(mixed_model_hipp_mean_act)


hipp_act_slopes <- emtrends(mixed_model_hipp_mean_act, specs = ~ hipp_mean_act_z, 
                        var = "years_from_baseline", 
                        at = list(hipp_mean_act_z = c(-1, 0, 1)))
summary(hipp_act_slopes, infer = TRUE)


mixed_model_activation_parietal_l <- lmer(moca ~ years_from_baseline * activation_parietal_sup_l_z + + sexe + diagnostic_nick + education +
                                            initiale_age + (1 | pscid),data = MRS_prediction_long)
summary(mixed_model_activation_parietal_l)


#sink()


# The Multimodal  Model
mixed_model_multi <- lmer(moca ~ years_from_baseline * (m_m_acc_z + m_m_precuneus_z  + cortical_thickness_adsignature_dickson_z + arsenii_hippocampus_avg_act + plasma_ptau217) + 
    sexe + diagnostic_nick + education + initiale_age + age_difference + (1 | pscid),  
  data = MRS_prediction_long
)
summary(mixed_model_multi)



################# Logistic regression ######################
names(MRS_prediction)
## Overal model sig
overall_m_sig <- function(model) {
  chsq <- model$null.deviance - model$deviance
  df <- model$df.null - model$df.residual
  pval <- pchisq(chsq, df, lower.tail = FALSE)
  return(data.frame(Chi_Sq = chsq, df = df, p_value = pval))}
#overall_m_sig(object)


#sink("glm_models_sensetivity_decliners.txt")


# plasma_ptau217 
model_glu_ptau217 <- glm(decliner_regression ~ plasma_ptau217_z  , data = MRS_prediction, family = "binomial")
summary(model_glu_ptau217)
roc_ptau217 <- roc(model_glu_ptau217$y, fitted(model_glu_ptau217))
auc(roc_ptau217)
coords(roc_ptau217, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_glu_ptau217)


roc_raw_ptau217 <- roc(MRS_prediction$decliner_regression, MRS_prediction$plasma_ptau217)
coords(roc_raw_ptau217, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")

# Precuneus Glutamate
model_glu_prec <- glm(decliner_regression ~ m_m_precuneus_z  , data = MRS_prediction, family = "binomial")
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

roc_raw_acc <- roc(MRS_prediction$decliner_regression, MRS_prediction$m_m_acc)
coords(roc_raw_acc, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")



# Hippocampal Volume
model_struc_hip <- glm(decliner_regression ~ hipp_mean_z  , data = MRS_prediction, family = "binomial")
summary(model_struc_hip)
roc_struc_hip <- roc(model_struc_hip$y, fitted(model_struc_hip))
auc(roc_struc_hip)
coords(roc_struc_hip, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_struc_hip)

# Cortical Thickness
model_struc_thick <- glm(decliner_regression ~ cortical_thickness_adsignature_dickson_z  , data = MRS_prediction, family = "binomial")
summary(model_struc_thick)
roc_struc_thick <- roc(model_struc_thick$y, fitted(model_struc_thick))
auc(roc_struc_thick)
coords(roc_struc_thick, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_struc_thick)

# Hippocampal Activation
model_func_hip <- glm(decliner_regression ~ hipp_mean_act_z  , data = MRS_prediction, family = "binomial")
summary(model_func_hip)
roc_func_hip <- roc(model_func_hip$y, fitted(model_func_hip))
auc(roc_func_hip)
coords(roc_func_hip, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_func_hip)


# Superior Parietal Activation
model_func_par <- glm(decliner_regression ~ activation_parietal_sup_l_z , data = MRS_prediction, family = "binomial")
summary(model_func_par)
roc_func_par <- roc(model_func_par$y, fitted(model_func_par))
auc(roc_func_par)
coords(roc_func_par, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_func_par)



############ 2 Variable Models ##############


# 1. Plasma p-Tau217 and ACC Glutamate
model_ptau217_acc <- glm(decliner_regression ~ plasma_ptau217_z + m_m_acc_z, data = MRS_prediction, family = "binomial")
summary(model_ptau217_acc)
roc_ptau217_acc <- roc(model_ptau217_acc$y, fitted(model_ptau217_acc))
auc(roc_ptau217_acc)
coords(roc_ptau217_acc, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_ptau217_acc)

# 2. Plasma p-Tau217 and Precuneus Glutamate
model_ptau217_prec <- glm(decliner_regression ~ plasma_ptau217_z + m_m_precuneus_z, data = MRS_prediction, family = "binomial")
summary(model_ptau217_prec)
roc_ptau217_prec <- roc(model_ptau217_prec$y, fitted(model_ptau217_prec))
auc(roc_ptau217_prec)
coords(roc_ptau217_prec, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_ptau217_prec)

# 3. Precuneus Glutamate + AD-Signature Cortical Thickness
model_prec_thick <- glm(decliner_regression ~ m_m_precuneus_z + cortical_thickness_adsignature_dickson_z, data = MRS_prediction, family = "binomial")
summary(model_prec_thick)
roc_prec_thick <- roc(model_prec_thick$y, fitted(model_prec_thick))
auc(roc_prec_thick)
coords(roc_prec_thick, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_prec_thick)

# 4. ACC Glutamate + AD-Signature Cortical Thickness
model_acc_thick <- glm(decliner_regression ~ m_m_acc_z + cortical_thickness_adsignature_dickson_z, data = MRS_prediction, family = "binomial")
summary(model_acc_thick)
roc_acc_thick <- roc(model_acc_thick$y, fitted(model_acc_thick))
auc(roc_acc_thick)
coords(roc_acc_thick, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_acc_thick)

# 5. Precuneus Glutamate + Hippocampal Volume
model_prec_hip_vol <- glm(decliner_regression ~ m_m_precuneus_z + hipp_mean_z, data = MRS_prediction, family = "binomial")
summary(model_prec_hip_vol)
roc_prec_hip_vol <- roc(model_prec_hip_vol$y, fitted(model_prec_hip_vol))
auc(roc_prec_hip_vol)
coords(roc_prec_hip_vol, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_prec_hip_vol)

# 6. ACC Glutamate + Hippocampal Volume
model_acc_hip_vol <- glm(decliner_regression ~ m_m_acc_z + hipp_mean_z, data = MRS_prediction, family = "binomial")
summary(model_acc_hip_vol)
roc_acc_hip_vol <- roc(model_acc_hip_vol$y, fitted(model_acc_hip_vol))
auc(roc_acc_hip_vol)
coords(roc_acc_hip_vol, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_acc_hip_vol)

# 7. Precuneus Glutamate + Superior Parietal Activation
model_prec_par_act <- glm(decliner_regression ~ m_m_precuneus_z + activation_parietal_sup_l_z, data = MRS_prediction, family = "binomial")
summary(model_prec_par_act)
roc_prec_par_act <- roc(model_prec_par_act$y, fitted(model_prec_par_act))
auc(roc_prec_par_act)
coords(roc_prec_par_act, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_prec_par_act)

# 8. ACC Glutamate + Superior Parietal Activation
model_acc_par_act <- glm(decliner_regression ~ m_m_acc_z + activation_parietal_sup_l_z, data = MRS_prediction, family = "binomial")
summary(model_acc_par_act)
roc_acc_par_act <- roc(model_acc_par_act$y, fitted(model_acc_par_act))
auc(roc_acc_par_act)
coords(roc_acc_par_act, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_acc_par_act)

# 9. Precuneus Glutamate + Hippocampal Activation
model_prec_hip_act <- glm(decliner_regression ~ m_m_precuneus_z + hipp_mean_act_z, data = MRS_prediction, family = "binomial")
summary(model_prec_hip_act)
roc_prec_hip_act <- roc(model_prec_hip_act$y, fitted(model_prec_hip_act))
auc(roc_prec_hip_act)
coords(roc_prec_hip_act, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_prec_hip_act)

# 10. ACC Glutamate + Hippocampal Activation
model_acc_hip_act <- glm(decliner_regression ~ m_m_acc_z + hipp_mean_act_z, data = MRS_prediction, family = "binomial")
summary(model_acc_hip_act)
roc_acc_hip_act <- roc(model_acc_hip_act$y, fitted(model_acc_hip_act))
auc(roc_acc_hip_act)
coords(roc_acc_hip_act, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_acc_hip_act)

# 11.  Precuneus Glutamate + ACC Glutamate
model_prec_acc <- glm(decliner_regression ~ m_m_precuneus_z + m_m_acc_z, data = MRS_prediction, family = "binomial")
summary(model_prec_acc)
roc_prec_acc <- roc(model_prec_acc$y, fitted(model_prec_acc))
auc(roc_prec_acc)
coords(roc_prec_acc, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_prec_acc)



#sink()






########################### With mixed GLM #########################
library(blme)
library(dplyr)

# Ensure data is cross-sectional
MRS_cross_sectional <- MRS_prediction_long %>%
  distinct(pscid, .keep_all = TRUE)

# Run with a gamma prior on the random effect standard deviation
b_model_ptau <- bglmer(
  decliner_regression ~ plasma_ptau217_z + (1 | pscid), 
  data = MRS_cross_sectional, 
  family = binomial(link = "logit"),
  fixef.prior = normal(sd = 10),                            
  cov.prior = gamma(shape = 2, rate = 0.5, posterior.scale = "sd"), 
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxit = 100000)))

summary(b_model_ptau)




############ One leave out cross validation ##############
library(pROC)

# 1. Define the master LOOCV function
run_loocv <- function(formula_str, data, model_name) {
  form <- as.formula(formula_str)
  vars <- all.vars(form)
  
  # Isolate data for this specific model and drop missing observations
  clean_data <- data[, vars, drop = FALSE]
  clean_data <- clean_data[complete.cases(clean_data), ]
  
  n <- nrow(clean_data)
  cv_predictions <- numeric(n)
  
  # Loop row-by-row (Leave One Out)
  for (i in 1:n) {
    train_set <- clean_data[-i, ]
    test_set  <- clean_data[i, ]
    
    # Train on N-1 subjects
    fit <- glm(form, data = train_set, family = "binomial")
    
    # Predict on the left-out subject
    cv_predictions[i] <- predict(fit, newdata = test_set, type = "response")
  }
  
  # Generate true out-of-sample ROC curves
  actual_outcomes <- clean_data[[vars[1]]]
  roc_cv <- pROC::roc(actual_outcomes, cv_predictions, quiet = TRUE)
  
  # Extract optimal threshold coordinates via Youden Index
  coords_cv <- pROC::coords(roc_cv, "best", ret = c("threshold", "specificity", "sensitivity"), best.method = "youden")
  if (!is.null(nrow(coords_cv)) && nrow(coords_cv) > 1) coords_cv <- coords_cv[1, ]
  
  # Calculate true out-of-sample classification accuracy
  predicted_classes <- ifelse(cv_predictions >= coords_cv$threshold, 1, 0)
  acc_cv <- mean(predicted_classes == actual_outcomes)
  
  # Return metrics as a clean data frame row
  return(data.frame(
    Model  = model_name,
    N_Obs  = n,
    CV_AUC = round(pROC::auc(roc_cv) * 100, 1),
    CV_ACC = round(acc_cv * 100, 1),
    CV_SEN = round(coords_cv$sensitivity * 100, 1),
    CV_SPE = round(coords_cv$specificity * 100, 1)
  ))
}

# 2. Map every single one of your models to its exact formula string
models_to_test <- list(
  # --- Unimodal Models ---
  c("decliner_regression ~ plasma_ptau217_z", "Unimodal: Plasma p-Tau217"),
  c("decliner_regression ~ m_m_precuneus_z", "Unimodal: Precuneus Glutamate"),
  c("decliner_regression ~ m_m_acc_z", "Unimodal: ACC Glutamate"),
  c("decliner_regression ~ hipp_mean_z", "Unimodal: Hippocampal Volume"),
  c("decliner_regression ~ cortical_thickness_adsignature_dickson_z", "Unimodal: Cortical Thickness"),
  c("decliner_regression ~ hipp_mean_act_z", "Unimodal: Hippocampal Activation"),
  c("decliner_regression ~ activation_parietal_sup_l_z", "Unimodal: Superior Parietal Activation"),
  
  # --- Bivariate Models ---
  c("decliner_regression ~ plasma_ptau217_z + m_m_acc_z", "Bivariate: p-Tau217 + ACC Glu"),
  c("decliner_regression ~ plasma_ptau217_z + m_m_precuneus_z", "Bivariate: p-Tau217 + Precuneus Glu"),
  c("decliner_regression ~ m_m_precuneus_z + cortical_thickness_adsignature_dickson_z", "Bivariate: Precuneus Glu + Cortical Thickness"),
  c("decliner_regression ~ m_m_acc_z + cortical_thickness_adsignature_dickson_z", "Bivariate: ACC Glu + Cortical Thickness"),
  c("decliner_regression ~ m_m_precuneus_z + hipp_mean_z", "Bivariate: Precuneus Glu + Hippocampal Volume"),
  c("decliner_regression ~ m_m_acc_z + hipp_mean_z", "Bivariate: ACC Glu + Hippocampal Volume"),
  c("decliner_regression ~ m_m_precuneus_z + activation_parietal_sup_l_z", "Bivariate: Precuneus Glu + Parietal Activation"),
  c("decliner_regression ~ m_m_acc_z + activation_parietal_sup_l_z", "Bivariate: ACC Glu + Parietal Activation"),
  c("decliner_regression ~ m_m_precuneus_z + hipp_mean_act_z", "Bivariate: Precuneus Glu + Hippocampal Activation"),
  c("decliner_regression ~ m_m_acc_z + hipp_mean_act_z", "Bivariate: ACC Glu + Hippocampal Activation"),
  c("decliner_regression ~ m_m_precuneus_z + m_m_acc_z", "Bivariate: Precuneus Glu + ACC Glu")
)

# 3. Execute the cross-validation loop over your dataset
loocv_results_list <- lapply(models_to_test, function(m) {
  run_loocv(formula_str = m[1], data = MRS_prediction, model_name = m[2])
})

# 4. Bind into a clear, unified master results matrix
master_loocv_table <- do.call(rbind, loocv_results_list)
print(master_loocv_table)



################# # Objerctive 3 #######################
## Sruvival analysis #####
library(survival)
names(MRS_prediction)
sapply(MRS_prediction, class)

summary(coxph(Surv(time_to_event, decliner_regression) ~ plasma_ptau217_z + initiale_age, data = MRS_prediction))
summary(coxph(Surv(time_to_event, decliner_regression) ~ m_m_acc_z  + initiale_age, data = MRS_prediction))
summary(coxph(Surv(time_to_event, decliner_regression) ~ m_m_precuneus_z + initiale_age, data = MRS_prediction))
summary(coxph(Surv(time_to_event, decliner_regression) ~ hipp_mean + initiale_age, data = MRS_prediction))
summary(coxph(Surv(time_to_event, decliner_regression) ~ cortical_thickness_adsignature_dickson + initiale_age, data = MRS_prediction))
summary(coxph(Surv(time_to_event, decliner_regression) ~ hipp_mean_act + initiale_age, data = MRS_prediction))
summary(coxph(Surv(time_to_event, decliner_regression) ~ activation_parietal_sup_l + initiale_age, data = MRS_prediction))





summary(coxph(Surv(time_to_event, ever_declined) ~ m_m_acc_z + plasma_ptau217_z + initiale_age, data = MRS_prediction))
summary(coxph(Surv(time_to_event, ever_declined) ~ m_m_acc_z + activation_parietal_sup_l + initiale_age, data = MRS_prediction))
summary(coxph(Surv(time_to_event, ever_declined) ~ m_m_acc_z + hipp_mean_act_z + initiale_age, data = MRS_prediction))
summary(coxph(Surv(time_to_event, ever_declined) ~ m_m_acc_z + cortical_thickness_adsignature_dickson + initiale_age, data = MRS_prediction))
summary(coxph(Surv(time_to_event, ever_declined) ~ m_m_acc_z + hipp_mean + initiale_age, data = MRS_prediction))



# 1. Create the isolated dataset using a completely unique column name
cox_data <- MRS_prediction_long %>%
  group_by(pscid) %>%
  summarize(followup_years = max(years_from_baseline, na.rm = TRUE)) %>%
  left_join(MRS_prediction, by = "pscid")

# 2. Run the Cox models using 'followup_years'
summary(coxph(Surv(followup_years, decliner_regression) ~ plasma_ptau217_z + initiale_age, data = cox_data))

summary(coxph(Surv(followup_years, decliner_regression) ~ m_m_acc_z + initiale_age, data = cox_data))

summary(coxph(Surv(followup_years, decliner_regression) ~ m_m_precuneus_z + initiale_age, data = cox_data))

summary(coxph(Surv(followup_years, decliner_regression) ~ hipp_mean_z + initiale_age, data = cox_data))

summary(coxph(Surv(followup_years, decliner_regression) ~ cortical_thickness_adsignature_dickson_z + initiale_age, data = cox_data))

summary(coxph(Surv(followup_years, decliner_regression) ~ hipp_mean_act_z + initiale_age, data = cox_data))

summary(coxph(Surv(followup_years, decliner_regression) ~ activation_parietal_sup_l_z + initiale_age, data = cox_data))




# 1. Overall Time to 50% Decline for the whole cohort
survfit(Surv(followup_years, decliner_regression) ~ 1, data = cox_data)

# 2. Time to 50% Decline split by your ACC Glutamate groups
survfit(Surv(followup_years, decliner_regression) ~ acc_strata, data = cox_data)

survfit(Surv(followup_years, decliner_regression) ~ prec_strata, data = cox_data)