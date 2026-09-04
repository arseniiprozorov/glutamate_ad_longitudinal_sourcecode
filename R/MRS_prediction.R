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

# Regrouper les groupes 'HC' et 'SCD' dans un seul groupe 'SCD'
levels(MRS_prediction$diagnostic_nick)[levels(MRS_prediction$diagnostic_nick) == "SCD"] <- "HC"
table(MRS_prediction$diagnostic_nick)

# Creer Hipp mean
MRS_prediction$hipp_mean <- (MRS_prediction$hip_l_nor_icv + MRS_prediction$righ_hip_vol)/2
#MRS_prediction$hipp_mean_act <- (MRS_prediction$activation_hippocampus_l + MRS_prediction$activation_hippocampus_r)/2

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
#MRS_prediction$activation_hippocampus_l <- winsorize_iqr(MRS_prediction$activation_hippocampus_l, 1.5)
#MRS_prediction$hipp_mean_act <- winsorize_iqr(MRS_prediction$hipp_mean_act, 1.5)
MRS_prediction$arsenii_hippocampus_avg_act <- winsorize_iqr(MRS_prediction$arsenii_hippocampus_avg_act, 1.5)
#MRS_prediction$activation_parietal_sup_l <- winsorize_iqr(MRS_prediction$activation_parietal_sup_l, 1.5)
MRS_prediction$arsenii_parietal_sup_l_act <- winsorize_iqr(MRS_prediction$arsenii_parietal_sup_l_act, 1.5)
MRS_prediction$activation_temporal_inf_r <- winsorize_iqr(MRS_prediction$activation_temporal_inf_r, 1.5)



# Descriptives
#sink()
jmv::descriptives(data = MRS_prediction, vars = vars(slope_regression_yearly, sustained_decliner, ever_declined, decliner_regression,
                                                     m_m_precuneus, m_m_acc, plasma_ptau217,
                                                     cortical_thickness_adsignature_dickson, hipp_mean, 
                                                     arsenii_hippocampus_avg_act, arsenii_parietal_sup_l_act),
                  sd = TRUE, iqr = TRUE, skew = TRUE, kurt = TRUE)

jmv::descriptives(data = MRS_prediction, vars = vars(slope_regression_yearly, m_m_precuneus, m_m_acc, plasma_ptau217,
                                                     cortical_thickness_adsignature_dickson, 
                                                     hipp_mean,arsenii_hippocampus_avg_act, arsenii_parietal_sup_l_act),
                  , splitBy =  sustained_decliner , sd = TRUE, iqr = TRUE, skew = TRUE, kurt = TRUE)





#### Calcul de score z
MRS_prediction$m_m_precuneus_z <- scale(MRS_prediction$m_m_precuneus)
MRS_prediction$m_m_acc_z <- scale(MRS_prediction$m_m_acc)
MRS_prediction$plasma_ptau217_z <- scale(MRS_prediction$plasma_ptau217)
MRS_prediction$cortical_thickness_adsignature_dickson_z <- scale(MRS_prediction$cortical_thickness_adsignature_dickson)
MRS_prediction$hipp_mean_z <- scale(MRS_prediction$hipp_mean)
#MRS_prediction$activation_hippocampus_l_z <- scale(MRS_prediction$activation_hippocampus_l)
MRS_prediction$arsenii_hippocampus_avg_act <- scale(MRS_prediction$arsenii_hippocampus_avg_act)
MRS_prediction$arsenii_parietal_sup_l_act <- scale(MRS_prediction$arsenii_parietal_sup_l_act)


## Long format 
X2026_06_16_MRS_prediction_long <- read_excel("C:/Users/okkam/Desktop/MRS_prediction_longitudinal_master.xlsx")
MRS_prediction_long <- X2026_06_16_MRS_prediction_long
lapply(MRS_prediction_long,class)
names(MRS_prediction_long)

MRS_prediction_long$arsenii_hippocampus_avg_act <- winsorize_iqr(MRS_prediction_long$arsenii_hippocampus_avg_act, 1.5)
MRS_prediction_long$arsenii_parietal_sup_l_act <- winsorize_iqr(MRS_prediction_long$arsenii_parietal_sup_l_act, 1.5)


# Convert to factor longitudinal
MRS_prediction_long$diagnostic_nick <- as.factor(MRS_prediction_long$diagnostic_nick)
MRS_prediction_long$sexe <- as.factor(MRS_prediction_long$sexe)

# Extract moca slopes for each participant 
raw_individual_models <- lmList(moca ~ years_from_baseline | pscid, data = MRS_prediction_long)
raw_slopes <- coef(raw_individual_models)
raw_slopes$pscid <- rownames(raw_slopes)
raw_slopes$moca_change_3_5_yrs <- raw_slopes$years_from_baseline * 3.5


######################################### Analyses #######################################
names(MRS_prediction)

######### Characterization #############


## change all decliners to sustained decliners !!!!


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
sapply(MRS_prediction_long,class)

#sink("Table 2.txt")
## Baseline model
mixed_model_moca <- lmer(moca ~ years_from_baseline + sexe + diagnostic_nick + education + initiale_age + (1 | pscid),  
                         REML = FALSE, data = MRS_prediction_long)
summary(mixed_model_moca)

#   pTau217 
mixed_model_ptau217 <- lmer(moca ~ years_from_baseline * plasma_ptau217_z + age_difference + sexe + diagnostic_nick + education + initiale_age +(1 | pscid),  
                            , REML = FALSE, data = MRS_prediction_long)
summary(mixed_model_ptau217)

#tau217_slopes <- emtrends(mixed_model_ptau217, specs = ~ plasma_ptau217_z, 
#                       var = "years_from_baseline", 
#                       at = list(plasma_ptau217_z = c(-1, 0, 1)))
#summary(tau217_slopes, infer = TRUE)

#anova(mixed_model_ptau217)

#   Glutamate 
mixed_model_precuneus <- lmer(moca ~  m_m_precuneus_z * years_from_baseline + sexe + diagnostic_nick + education + (1 | pscid),  
  ,REML = FALSE, data = MRS_prediction_long)
summary(mixed_model_precuneus)

#precuneus_slopes <- emtrends(mixed_model_precuneus, specs = ~ m_m_precuneus_z, 
#                            var = "years_from_baseline", 
#                            at = list(m_m_precuneus_z = c(-1,0,1)))
#summary(precuneus_slopes, infer = TRUE)


#### ACC ##########


mixed_model_acc <- lmer(moca ~ years_from_baseline * m_m_acc_z + sexe + diagnostic_nick + education + initiale_age + (1 | pscid),  
                        REML = FALSE, data = MRS_prediction_long)
summary(mixed_model_acc)

#acc_slopes <- emtrends(mixed_model_acc, specs = ~ m_m_acc_z, 
#                             var = "years_from_baseline", 
#                             at = list(m_m_acc_z = c(-1,0,1)))
#summary(acc_slopes, infer = TRUE)
#anova(mixed_model_acc)





# Structure
mixed_model_thickness <- lmer(
  moca ~ years_from_baseline * cortical_thickness_adsignature_dickson_z + sexe + diagnostic_nick + education + initiale_age + (1 | pscid),  
  REML = FALSE, data = MRS_prediction_long)
summary(mixed_model_thickness)

#thick_slopes <- emtrends(mixed_model_thickness, specs = ~ cortical_thickness_adsignature_dickson_z, 
#                          var = "years_from_baseline", 
#                          at = list(cortical_thickness_adsignature_dickson_z = c(-1, 0, 1)))
#summary(thick_slopes, infer = TRUE)
#anova(mixed_model_thickness)


mixed_model_hipp_mean <- lmer(moca ~ years_from_baseline * hipp_mean_z  + sexe + diagnostic_nick + education + initiale_age + (1 | pscid),  
                              REML = FALSE, data = MRS_prediction_long)
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
                                  REML = FALSE, data = MRS_prediction_long)
summary(mixed_model_hipp_mean_act)


#hipp_act_slopes <- emtrends(mixed_model_hipp_mean_act, specs = ~ arsenii_hippocampus_avg_act, 
#                        var = "years_from_baseline", 
#                       at = list(hipp_mean_act_z = c(-1, 0, 1)))
#summary(hipp_act_slopes, infer = TRUE)


mixed_model_activation_parietal_l <- lmer(moca ~ years_from_baseline * arsenii_parietal_sup_l_act + + sexe + diagnostic_nick + education +
                                            initiale_age + (1 | pscid),REML = FALSE, data = MRS_prediction_long)
summary(mixed_model_activation_parietal_l)





########### Post hoc with tertials ########

# Create Tertile Factor Variables (Low, Medium, High)
# Plasma p-Tau217
q_ptau <- quantile(MRS_prediction_long$plasma_ptau217_z, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
MRS_prediction_long$ptau217_tert <- cut(MRS_prediction_long$plasma_ptau217_z, breaks = q_ptau, 
                                        labels = c("Low", "Medium", "High"), include.lowest = TRUE)

q_acc <- quantile(MRS_prediction_long$m_m_acc_z, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
MRS_prediction_long$acc_tert <- cut(MRS_prediction_long$m_m_acc_z, breaks = q_acc, 
                                    labels = c("Low", "Medium", "High"), include.lowest = TRUE)

q_prec <- quantile(MRS_prediction_long$m_m_precuneus_z, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
MRS_prediction_long$precuneus_tert <- cut(MRS_prediction_long$m_m_precuneus_z, breaks = q_prec, 
                                          labels = c("Low", "Medium", "High"), include.lowest = TRUE)

q_thick <- quantile(MRS_prediction_long$cortical_thickness_adsignature_dickson_z, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
MRS_prediction_long$thickness_tert <- cut(MRS_prediction_long$cortical_thickness_adsignature_dickson_z, breaks = q_thick, 
                                          labels = c("Low", "Medium", "High"), include.lowest = TRUE)

q_hip <- quantile(MRS_prediction_long$arsenii_hippocampus_avg_act, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
MRS_prediction_long$hipp_act_tert <- cut(MRS_prediction_long$arsenii_hippocampus_avg_act, breaks = q_hip, 
                                         labels = c("Low", "Medium", "High"), include.lowest = TRUE)

# Run models
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



########################################### The Multimodal  Model###############################
# Syntax Isaora
#ols_step_backward_p(OLS_mod0,prem=0.1) 
help(package = olsrr)
citation("olsrr")
########### backward stepwise lienar mixed regression ################

### Arsenii Backward stepwise ###

MRS_step <- na.omit(MRS_prediction[, c("slope_regression_yearly","sustained_decliner",  "m_m_acc_z", 
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

# With AIC
#step_aic_result <- ols_step_backward_aic(full_lm_model, details = TRUE)
#print(step_aic_result)
#summary(step_aic_result)



############### Habib theory driven ###########
names(MRS_prediction_long)
MRS_core_long <- na.omit(MRS_prediction_long[, c("pscid", "moca", "years_from_baseline", "m_m_acc_z", 
                                                 "m_m_precuneus_z", "plasma_ptau217_z", 
                                                 "cortical_thickness_adsignature_dickson_z", 
                                                 "arsenii_hippocampus_avg_act", "sexe", 
                                                 "diagnostic_nick", "education", "initiale_age", "age_difference")])
# Step 1: null model 
base_model <- lmer(
  moca ~ years_from_baseline * (plasma_ptau217_z + cortical_thickness_adsignature_dickson_z) + 
    age_difference + sexe + diagnostic_nick + education + initiale_age + 
    (1 | pscid),data = MRS_core_long, REML = FALSE)
summary(base_model)


# Step 2a ACC Glutamate Addition
base_model_acc <- lmer(
  moca ~ years_from_baseline * (plasma_ptau217_z + cortical_thickness_adsignature_dickson_z + m_m_acc_z) + 
    age_difference + sexe + diagnostic_nick + education + initiale_age + 
    (1 | pscid),data = MRS_core_long, REML = FALSE)
summary(base_model_acc)
anova(base_model, base_model_acc)


# Step 2b Precuneus Glutamate Addition 
base_model_prec <- lmer(
  moca ~ years_from_baseline * (plasma_ptau217_z + cortical_thickness_adsignature_dickson_z + m_m_precuneus_z) + 
    age_difference + sexe + diagnostic_nick + education + initiale_age + 
    (1 | pscid),data = MRS_core_long, REML = FALSE)
summary(base_model_prec)
anova(base_model, base_model_prec)


# Step 2c Hippocampal fMRI Activation Addition
base_model_act <- lmer(moca ~ years_from_baseline * (plasma_ptau217_z + cortical_thickness_adsignature_dickson_z + arsenii_hippocampus_avg_act) + 
    age_difference + sexe + diagnostic_nick + education + initiale_age + 
    (1 | pscid),data = MRS_core_long, REML = FALSE)
summary(base_model_act)
anova(base_model, base_model_act)

##  step 3  Full model
full_lmm_model <- lmer(moca ~ years_from_baseline * (plasma_ptau217_z + cortical_thickness_adsignature_dickson_z + 
      m_m_acc_z + m_m_precuneus_z + arsenii_hippocampus_avg_act) + 
    age_difference + sexe + diagnostic_nick + education + initiale_age + 
    (1 | pscid),data = MRS_core_long, REML = FALSE)
summary(full_lmm_model)
anova(base_model, full_lmm_model)

## Step 4 backward elimination
step_lmm_res <- lmerTest::step(full_lmm_model, 
  reduce.fixed = TRUE, reduce.random = FALSE, alpha.remove = 0.05)
print(step_lmm_res)
final_lmm_model <- get_model(step_lmm_res)
summary(final_lmm_model)


#sink()


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


# Arsenii: decliners defined based on slope

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




# Habib: decliners defined based on flipping to delciner and staying decliner


# plasma_ptau217 
model_glu_ptau217 <- glm(sustained_decliner ~ plasma_ptau217_z + age_difference  , data = MRS_prediction, family = "binomial")
summary(model_glu_ptau217)
overall_m_sig(model_glu_ptau217)

roc_ptau217 <- roc(model_glu_ptau217$y, fitted(model_glu_ptau217))
auc(roc_ptau217)

coords(roc_ptau217, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")


# Precuneus Glutamate
model_glu_prec <- glm(sustained_decliner ~ m_m_precuneus_z , data = MRS_prediction, family = "binomial")
summary(model_glu_prec)
overall_m_sig(model_glu_prec)

roc_glu_prec <- roc(model_glu_prec$y, fitted(model_glu_prec))
auc(roc_glu_prec)

coords(roc_glu_prec, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")

# ACC Glutamate
model_glu_acc <- glm(sustained_decliner ~ m_m_acc_z  , data = MRS_prediction, family = "binomial")
summary(model_glu_acc)
overall_m_sig(model_glu_acc)

roc_glu_acc <- roc(model_glu_acc$y, fitted(model_glu_acc))
auc(roc_glu_acc)

coords(roc_glu_acc, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")



# Cortical Thickness
model_struc_thick <- glm(sustained_decliner ~ cortical_thickness_adsignature_dickson_z  , data = MRS_prediction, family = "binomial")
summary(model_struc_thick)
overall_m_sig(model_struc_thick)

roc_struc_thick <- roc(model_struc_thick$y, fitted(model_struc_thick))
auc(roc_struc_thick)
coords(roc_struc_thick, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")


# Hippocampal Activation
model_func_hip <- glm(sustained_decliner ~ arsenii_hippocampus_avg_act  , data = MRS_prediction, family = "binomial")
summary(model_func_hip)
overall_m_sig(model_func_hip)

roc_func_hip <- roc(model_func_hip$y, fitted(model_func_hip))
auc(roc_func_hip)

coords(roc_func_hip, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")




############ Multimodal model ############

# Arsenii multimodal.  ACC Glutamate + Plasma p-Tau217 + cortical thickenss 
model_nocov_step_sig <- glm(decliner_regression ~ plasma_ptau217_z  + cortical_thickness_adsignature_dickson_z + m_m_acc_z
                           , data = MRS_prediction, family = "binomial")
summary(model_nocov_step_sig)

roc_model_nocov_step_sig <- roc(model_nocov_step_sig$y, fitted(model_nocov_step_sig))
auc(roc_model_nocov_step_sig)
coords(roc_model_nocov_step_sig, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
overall_m_sig(model_nocov_step_sig)





# Habib multimodal. 
# Step 1. Null  model
model_glm_null <- glm(sustained_decliner ~ plasma_ptau217_z + cortical_thickness_adsignature_dickson_z,
                          data = MRS_step, family = "binomial")
summary(model_glm_null)
overall_m_sig(model_glm_null)

roc_null <- roc(model_glm_null$y, fitted(model_glm_null))
auc(roc_null)

ci_roc_null  <- ci.auc(roc_null, method = "bootstrap", boot.n = 1000)
print(ci_roc_null)

coords(roc_null, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
ci.coords(roc_null, x = 0.1916574, ret = c("sensitivity", "specificity"), 
          method = "bootstrap", boot.n = 1000)


# Step 2a. Null + ACC 
model_glm_null_acc <- glm(sustained_decliner ~ plasma_ptau217_z + cortical_thickness_adsignature_dickson_z + m_m_acc_z,
                          data = MRS_step, family = "binomial")
summary(model_glm_null_acc)
overall_m_sig(model_glm_null_acc)

roc_null_acc <- roc(model_glm_null_acc$y, fitted(model_glm_null_acc))
auc(roc_null_acc)

ci_roc_null_acc  <- ci.auc(roc_null_acc, method = "bootstrap", boot.n = 1000)
print(ci_roc_null_acc)

anova(model_glm_null, model_glm_null_acc, test = "Chisq")
roc.test(roc_null, roc_null_acc, method = "bootstrap", paired = TRUE, boot.n = 1000)

coords(roc_null_acc, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
ci.coords(roc_null_acc, x = 0.5374003, ret = c("sensitivity", "specificity"), 
          method = "bootstrap", boot.n = 1000)



# Step 2b. Null +  precuneus 
model_glm_null_prec <- glm(sustained_decliner ~ plasma_ptau217_z + cortical_thickness_adsignature_dickson_z + m_m_precuneus_z,
  data = MRS_step, family = "binomial")
summary(model_glm_null_prec)
overall_m_sig(model_glm_null_prec)

roc_null_prec <- roc(model_glm_null_prec$y, fitted(model_glm_null_prec))
auc(roc_null_prec)


ci_roc_null_prec  <- ci.auc(roc_null_prec, method = "bootstrap", boot.n = 1000)
print(ci_roc_null_prec)

anova(model_glm_null, model_glm_null_prec, test = "Chisq")
roc.test(roc_null, roc_null_prec, method = "bootstrap", paired = TRUE, boot.n = 1000)


coords(roc_null_prec, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
ci.coords(roc_null_prec, x = 0.235148, ret = c("sensitivity", "specificity"), 
          method = "bootstrap", boot.n = 1000)


# Step 2c. Null +  activation hipp 
model_glm_null_act <- glm(sustained_decliner ~ plasma_ptau217_z + cortical_thickness_adsignature_dickson_z + arsenii_hippocampus_avg_act,
                     data = MRS_step, family = "binomial")
summary(model_glm_null_act)
overall_m_sig(model_glm_null_act)

roc_null_act <- roc(model_glm_null_act$y, fitted(model_glm_null_act))
auc(roc_null_act)

ci_roc_null_act  <- ci.auc(roc_null_act, method = "bootstrap", boot.n = 1000)
print(ci_roc_null_act)

anova(model_glm_null, model_glm_null_act, test = "Chisq")
roc.test(roc_null, roc_null_act, method = "bootstrap", paired = TRUE, boot.n = 1000)

coords(roc_null_act, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
ci.coords(roc_null_act, x = 0.174571, ret = c("sensitivity", "specificity"), 
          method = "bootstrap", boot.n = 1000)



# Step 3. Full model
model_glm_full <- glm(sustained_decliner ~ plasma_ptau217_z + cortical_thickness_adsignature_dickson_z + m_m_acc_z + m_m_precuneus_z + arsenii_hippocampus_avg_act,
  data = MRS_step, 
  family = "binomial")
summary(model_glm_full)
overall_m_sig(model_glm_full)


roc_full <- roc(model_glm_full$y, fitted(model_glm_full))
auc(roc_full)

ci_roc_full  <- ci.auc(roc_full, method = "bootstrap", boot.n = 1000)
print(ci_roc_full)

anova(model_glm_null, model_glm_full, test = "Chisq")
roc.test(roc_null, roc_full, method = "bootstrap", paired = TRUE, boot.n = 1000)

coords(roc_full, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")
ci.coords(roc_full, x = 0.3600494, ret = c("sensitivity", "specificity"), 
          method = "bootstrap", boot.n = 1000)







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

surv_acc <- coxph(Surv(time_to_sustained_decline , sustained_decliner) ~ m_m_acc_z  + initiale_age, data = MRS_prediction)
summary(surv_acc)
exp(confint(surv_acc))

surv_prec <- coxph(Surv(time_to_sustained_decline, sustained_decliner) ~ m_m_precuneus_z + initiale_age, data = MRS_prediction)
summary(surv_prec)
exp(confint(surv_prec))
summary(coxph(Surv(time_to_sustained_decline, sustained_decliner) ~ plasma_ptau217_z + initiale_age, data = MRS_prediction))
summary(coxph(Surv(time_to_sustained_decline, sustained_decliner) ~ cortical_thickness_adsignature_dickson_z + initiale_age, data = MRS_prediction))
summary(coxph(Surv(time_to_sustained_decline, sustained_decliner) ~ arsenii_hippocampus_avg_act + initiale_age, data = MRS_prediction))

#summary(coxph(Surv(max_years_from_baseline, decliner_regression) ~ activation_parietal_sup_l + initiale_age, data = MRS_prediction))
#summary(coxph(Surv(max_years_from_baseline, decliner_regression) ~ hipp_mean + initiale_age, data = MRS_prediction))


#### Multimodal Cox
# Null model
cox_null <- coxph(Surv(time_to_sustained_decline, sustained_decliner) ~ plasma_ptau217_z  + cortical_thickness_adsignature_dickson_z, data = MRS_prediction)
summary(cox_null)

# + individual biomarkers
cox_null_acc <- coxph(Surv(time_to_sustained_decline, sustained_decliner) ~ plasma_ptau217_z  + cortical_thickness_adsignature_dickson_z + 
                        m_m_acc_z, data = MRS_prediction)
summary(cox_null_acc)

cox_null_prec <- coxph(Surv(time_to_sustained_decline, sustained_decliner) ~ plasma_ptau217_z  + cortical_thickness_adsignature_dickson_z  +
                    m_m_precuneus_z , data = MRS_prediction)
summary(cox_null_prec)

cox_null_act <- coxph(Surv(time_to_sustained_decline, sustained_decliner) ~ plasma_ptau217_z  + cortical_thickness_adsignature_dickson_z +
                     arsenii_hippocampus_avg_act, data = MRS_prediction)
summary(cox_null)

# Full model 
cox_full <- coxph(Surv(time_to_sustained_decline, sustained_decliner) ~ plasma_ptau217_z  + cortical_thickness_adsignature_dickson_z + m_m_acc_z +
                     m_m_precuneus_z + arsenii_hippocampus_avg_act, data = MRS_prediction)
summary(cox_full)







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



survfit(Surv(time_to_sustained_decline, sustained_decliner) ~ acc_tertiales, data = MRS_prediction)

survfit(Surv(time_to_sustained_decline, sustained_decliner) ~ precuneus_tertiales, data = MRS_prediction)


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





### Multimodal Keplen Mayer curve Arsenii
MRS_prediction$multimodal_risk_score <- predict(cox_nocov, newdata = MRS_prediction, type = "lp", na.action = na.pass)

MRS_prediction$multimodal_risk_tertiles <- cut(
  MRS_prediction$multimodal_risk_score,
  breaks = quantile(MRS_prediction$multimodal_risk_score, probs = 0:3/3, na.rm = TRUE),
  labels = c("Low Risk Profile", "Intermediate Risk Profile", "High Risk Profile"),
  include.lowest = TRUE)

# Fit and summarize Kaplan-Meier model
fit_multimodal <- survfit(Surv(max_years_from_baseline, decliner_regression) ~ multimodal_risk_tertiles, data = MRS_prediction)
summary(fit_multimodal)






library(survival)

# 1. Align time variable (event time for decliners, last visit for censored)
MRS_prediction$surv_time <- ifelse(
  MRS_prediction$sustained_decliner == 1,
  MRS_prediction$time_to_sustained_decline,
  MRS_prediction$max_years_from_baseline
)
MRS_prediction$surv_time <- pmax(MRS_prediction$surv_time, 0.1)

# 2. Fit Multimodal Model (cox_full) and generate risk scores
cox_full <- coxph(
  Surv(surv_time, sustained_decliner) ~ plasma_ptau217_z + cortical_thickness_adsignature_dickson_z + 
    m_m_acc_z + m_m_precuneus_z + arsenii_hippocampus_avg_act, 
  data = MRS_prediction
)

MRS_prediction$risk_score <- predict(cox_full, newdata = MRS_prediction, type = "lp")

MRS_prediction$risk_tertile <- cut(
  MRS_prediction$risk_score,
  breaks = quantile(MRS_prediction$risk_score, probs = 0:3/3, na.rm = TRUE),
  labels = c("Low Risk", "Intermediate Risk", "High Risk"),
  include.lowest = TRUE
)

# 3. Median Time-to-Decline: Multimodal Risk Strata
fit_multimodal_strata <- survfit(Surv(surv_time, sustained_decliner) ~ risk_tertile, data = MRS_prediction)
print(fit_multimodal_strata)

# 4. Median Time-to-Decline: Significant Unimodal Stratum (ACC Glutamate)
acc_tertiles <- quantile(MRS_prediction$m_m_acc_z, probs = 0:3/3, na.rm = TRUE)
MRS_prediction$acc_tertile_cat <- cut(
  MRS_prediction$m_m_acc_z,
  breaks = acc_tertiles,
  labels = c("Low ACC Glu", "Medium ACC Glu", "High ACC Glu"),
  include.lowest = TRUE
)

fit_acc_strata <- survfit(Surv(surv_time, sustained_decliner) ~ acc_tertile_cat, data = MRS_prediction)
print(fit_acc_strata)

# 5. Median Time-to-Decline: Progressor Subset Only (Observed vs. Midpoint)
decliners <- subset(MRS_prediction, sustained_decliner == 1)
decliners$midpoint_time <- (pmax(0, decliners$time_to_sustained_decline - 2) + decliners$time_to_sustained_decline) / 2

cat("\n--- Progressors Only (n =", nrow(decliners), ") ---\n")
cat("Median First Documented Visit (R):", median(decliners$time_to_sustained_decline, na.rm = TRUE), "years\n")
cat("IQR First Documented Visit (R):", quantile(decliners$time_to_sustained_decline, probs = c(0.25, 0.75), na.rm = TRUE), "years\n")
cat("Median Interval Midpoint Onset:", median(decliners$midpoint_time, na.rm = TRUE), "years\n")
cat("IQR Interval Midpoint Onset:", quantile(decliners$midpoint_time, probs = c(0.25, 0.75), na.rm = TRUE), "years\n")




