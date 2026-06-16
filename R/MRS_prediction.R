library(readxl)   
library(janitor)  
library(jmv)
library(effectsize)
library(interactions)
library(emmeans)

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
MRS_prediction$activation_hippocampus_r <- winsorize_iqr(MRS_prediction$hipp_mean_act, 1.5)
MRS_prediction$activation_parietal_sup_l <- winsorize_iqr(MRS_prediction$activation_parietal_sup_l, 1.5)
MRS_prediction$activation_temporal_inf_r <- winsorize_iqr(MRS_prediction$activation_temporal_inf_r, 1.5)



# Descriptives
#sink()
jmv::descriptives(data = MRS_prediction, vars = vars(m_m_precuneus, m_m_acc, plasma_ptau217,
                                                     cortical_thickness_adsignature_dickson, 
                                                     hipp_mean,hipp_mean_act,activation_parietal_sup_l),
                  sd = TRUE, iqr = TRUE, skew = TRUE, kurt = TRUE)


######### Analyses ######
## Characterization ##
# Table 1: Sociodemographic & Clinical Characteristics
# Table 1.1  Split by Decliners Status 
jmv::descriptives(data = MRS_prediction, vars = vars(initiale_age, education, slope_moca_raw, age_difference),
                  sd = TRUE, iqr = TRUE, splitBy = decliners, skew = TRUE, kurt = TRUE)

t.test(initiale_age ~ decliners, data = MRS_prediction)

chisq.test(MRS_prediction$sexe, y = MRS_prediction$decliners, correct = TRUE)
table(MRS_prediction$sexe, MRS_prediction$decliners)

t.test(education ~ decliners, data = MRS_prediction)

chisq.test(MRS_prediction$diagnostic_nick, y = MRS_prediction$decliners, correct = TRUE)
table(MRS_prediction$diagnostic_nick, MRS_prediction$decliners)

t.test(slope_moca_raw ~ decliners, data = MRS_prediction)
t.test(age_difference ~ decliners, data = MRS_prediction)


#  Table 1.2 - Split by Baseline Diagnosis 
jmv::descriptives(data = MRS_prediction, vars = vars(initiale_age, education, slope_moca_raw, age_difference),
                  sd = TRUE, iqr = TRUE, splitBy = diagnostic_nick, skew = TRUE, kurt = TRUE)

summary(aov(initiale_age ~ diagnostic_nick, data = MRS_prediction))

chisq.test(MRS_prediction$sexe, y = MRS_prediction$diagnostic_nick, correct = TRUE)
table(MRS_prediction$sexe, MRS_prediction$diagnostic_nick)

summary(aov(education ~ diagnostic_nick, data = MRS_prediction))
summary(aov(slope_moca_raw ~ diagnostic_nick, data = MRS_prediction))
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


########## Objective 1 #############
names(MRS_prediction)
# Moca slope as continous
summary(lm(slope_moca_raw ~ m_m_precuneus, data = MRS_prediction))
summary(lm(slope_moca_raw ~ m_m_acc, data = MRS_prediction))
summary(lm(slope_moca_raw ~ plasma_ptau217, data = MRS_prediction))
summary(lm(slope_moca_raw ~ cortical_thickness_adsignature_dickson, data = MRS_prediction))
summary(lm(slope_moca_raw ~ hipp_mean, data = MRS_prediction))
summary(lm(slope_moca_raw ~ hipp_mean_act, data = MRS_prediction))
summary(lm(slope_moca_raw ~ activation_parietal_sup_l, data = MRS_prediction))


## same with LMM
library(lme4)
library(lmerTest) 

X2026_06_15_MRS_prediction_long <- read_excel("C:/Users/okkam/Desktop/labo/article 2/Longitudinal_Multimodal_Data_CIMAQ/article_prediction/2026-06-15_MRS_prediction_long.xlsx")
MRS_prediction_long <- X2026_06_15_MRS_prediction_long

names(MRS_prediction_long)
# Mixed-Effects Models

#   Glutamate 
mixed_model_precuneus <- lmer(moca ~ years_from_baseline * m_m_precuneus + (1 | pscid),  
  data = MRS_prediction_long)
summary(mixed_model_precuneus)

mixed_model_acc <- lmer(moca ~ years_from_baseline * m_m_acc + (1 | pscid),  
  data = MRS_prediction_long)
summary(mixed_model_acc)

#   pTau217 
mixed_model_ptau217 <- lmer(moca ~ years_from_baseline * plasma_ptau217 + (1 | pscid),  
  data = MRS_prediction_long)
summary(mixed_model_ptau217)

# 4. Cortical Thickness AD Signature Model
mixed_model_thickness <- lmer(
  moca ~ age_centered * cortical_thickness_adsignature_dickson_centered + (age_centered | pscid),  
  data = MRS_long_mixed
)
summary(mixed_model_thickness)

# 5. Normalized Left Hippocampus Volume Model
mixed_model_hip_l_vol <- lmer(
  moca ~ age_centered * hip_l_nor_icv_centered + (age_centered | pscid),  
  data = MRS_long_mixed
)
summary(mixed_model_hip_l_vol)

# 6. Left Hippocampus Activation Model
mixed_model_activation_hip_l <- lmer(
  moca ~ age_centered * activation_hippocampus_l_centered + (age_centered | pscid),  
  data = MRS_long_mixed
)
summary(mixed_model_activation_hip_l)

# 7. Right Hippocampus Activation Model
mixed_model_activation_hip_r <- lmer(
  moca ~ age_centered * activation_hippocampus_r_centered + (age_centered | pscid),  
  data = MRS_long_mixed
)
summary(mixed_model_activation_hip_r)

# 8. Left Superior Parietal Activation Model
mixed_model_activation_parietal_l <- lmer(
  moca ~ age_centered * activation_parietal_sup_l_centered + (age_centered | pscid),  
  data = MRS_long_mixed
)
summary(mixed_model_activation_parietal_l)




## Logistic regression
library(pROC)

# ==========================================
# 1 Variable Models
# ==========================================

# plasma_ptau217 
model_glu_ptau217 <- glm(decliners ~ plasma_ptau217, data = MRS_prediction, family = "binomial")
summary(model_glu_ptau217)
roc_ptau217 <- roc(model_glu_ptau217$y, fitted(model_glu_ptau217))
auc(roc_ptau217)
coords(roc_ptau217, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")

# Precuneus Glutamate
model_glu_prec <- glm(decliners ~ m_m_precuneus, data = MRS_prediction, family = "binomial")
summary(model_glu_prec)
roc_glu_prec <- roc(model_glu_prec$y, fitted(model_glu_prec))
auc(roc_glu_prec)
coords(roc_glu_prec, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")

# ACC Glutamate
model_glu_acc <- glm(decliners ~ m_m_acc, data = MRS_prediction, family = "binomial")
summary(model_glu_acc)
roc_glu_acc <- roc(model_glu_acc$y, fitted(model_glu_acc))
auc(roc_glu_acc)
coords(roc_glu_acc, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")

# Hippocampal Volume
model_struc_hip <- glm(decliners ~ hip_l_nor_icv, data = MRS_prediction, family = "binomial")
summary(model_struc_hip)
roc_struc_hip <- roc(model_struc_hip$y, fitted(model_struc_hip))
auc(roc_struc_hip)
coords(roc_struc_hip, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")

# Cortical Thickness
model_struc_thick <- glm(decliners ~ cortical_thickness_adsignature_dickson, data = MRS_prediction, family = "binomial")
summary(model_struc_thick)
roc_struc_thick <- roc(model_struc_thick$y, fitted(model_struc_thick))
auc(roc_struc_thick)
coords(roc_struc_thick, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")

# Hippocampal Activation
model_func_hip <- glm(decliners ~ activation_hippocampus_l, data = MRS_prediction, family = "binomial")
summary(model_func_hip)
roc_func_hip <- roc(model_func_hip$y, fitted(model_func_hip))
auc(roc_func_hip)
coords(roc_func_hip, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")

# Superior Parietal Activation
model_func_par <- glm(decliners ~ activation_parietal_sup_l, data = MRS_prediction, family = "binomial")
summary(model_func_par)
roc_func_par <- roc(model_func_par$y, fitted(model_func_par))
auc(roc_func_par)
coords(roc_func_par, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")


# ==========================================
# 2 Variable Models
# ==========================================

# plasma_ptau217 and glu acc
model_glu_ptau217_acc <- glm(decliners ~ plasma_ptau217 + m_m_acc, data = MRS_prediction, family = "binomial")
summary(model_glu_ptau217_acc)
roc_ptau217_acc <- roc(model_glu_ptau217_acc$y, fitted(model_glu_ptau217_acc))
auc(roc_ptau217_acc)
coords(roc_ptau217_acc, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")

# plasma_ptau217 and glu prec
model_glu_ptau217_prec <- glm(decliners ~ plasma_ptau217 + m_m_precuneus, data = MRS_prediction, family = "binomial")
summary(model_glu_ptau217_prec)
roc_ptau217_prec <- roc(model_glu_ptau217_prec$y, fitted(model_glu_ptau217_prec))
auc(roc_ptau217_prec)
coords(roc_ptau217_prec, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")

# Glu acc + parietal activation 
model_metab_func <- glm(decliners ~ m_m_acc + activation_parietal_sup_l, data = MRS_prediction, family = "binomial")
summary(model_metab_func)
roc_metab_func <- roc(model_metab_func$y, fitted(model_metab_func))
auc(roc_metab_func)
coords(roc_metab_func, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")

# Glu precuneus + Hipp volume 
model_metab_struc_hip <- glm(decliners ~ m_m_precuneus + hip_l_nor_icv, data = MRS_prediction, family = "binomial")
summary(model_metab_struc_hip)
roc_metab_struc_hip <- roc(model_metab_struc_hip$y, fitted(model_metab_struc_hip))
auc(roc_metab_struc_hip)
coords(roc_metab_struc_hip, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")

# Glu precuneus + cortical thickness 
model_metab_struc_thick <- glm(decliners ~ m_m_precuneus + cortical_thickness_adsignature_dickson, data = MRS_prediction, family = "binomial")
summary(model_metab_struc_thick)
roc_metab_struc_thick <- roc(model_metab_struc_thick$y, fitted(model_metab_struc_thick))
auc(roc_metab_struc_thick)
coords(roc_metab_struc_thick, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")


# ==========================================
# 3 Variable Models
# ==========================================

# 3 variables with tau (ptau217 + glu acc + parietal)
model_tau_glu_acc_parietal <- glm(decliners ~ plasma_ptau217 + m_m_acc + activation_parietal_sup_l, data = MRS_prediction, family = "binomial")
summary(model_tau_glu_acc_parietal)
roc_tau_glu_acc_parietal <- roc(model_tau_glu_acc_parietal$y, fitted(model_tau_glu_acc_parietal))
auc(roc_tau_glu_acc_parietal)
coords(roc_tau_glu_acc_parietal, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")

# 3 variables with tau (ptau217 + hipp vol + parietal)
model_tau_hipp_vol <- glm(decliners ~ plasma_ptau217 + hip_l_nor_icv + activation_parietal_sup_l, data = MRS_prediction, family = "binomial")
summary(model_tau_hipp_vol)
roc_tau_hipp_vol <- roc(model_tau_hipp_vol$y, fitted(model_tau_hipp_vol))
auc(roc_tau_hipp_vol)
coords(roc_tau_hipp_vol, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")












