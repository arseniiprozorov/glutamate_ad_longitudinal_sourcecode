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

X2026_05_02_glu_moca_struc_func <- read_excel("C:/Users/okkam/Desktop/labo/article 2/Longitudinal_Multimodal_Data_CIMAQ/2026-05-02_glu_moca_struc_func.xlsx")
MRS_long <- X2026_05_02_glu_moca_struc_func

# Clean the column name
MRS_long <- janitor::clean_names(MRS_long)
names(MRS_long)
sapply(MRS_long,class)
# Convertir  en   factor
MRS_long$sexe <- as.factor(MRS_long$sexe)
table(MRS_long$sexe)
MRS_long$diagnostic_nick <- as.factor(MRS_long$diagnostic_nick)
table(MRS_long$diagnostic_nick)
MRS_long$cluster_moca_raw <- as.factor(MRS_long$cluster_moca_raw)
table(MRS_long$cluster_moca_raw)
MRS_long$decliners <- as.factor(MRS_long$decliners)
table(MRS_long$decliners)
# Regrouper les groupes 'HC' et 'SCD' dans un seul groupe 'SCD'
levels(MRS_long$diagnostic_nick)[levels(MRS_long$diagnostic_nick) == "SCD"] <- "HC"
table(MRS_long$diagnostic_nick)

# Create traj glu variable (either increasing or decreasing)
#  precuneus
MRS_long$traj_glu_prec <- ifelse(MRS_long$percent_change_glu_prec > 0, "Increasing",
                                 ifelse(MRS_long$percent_change_glu_prec < 0, "Decreasing", NA))

# 2.ACC
MRS_long$traj_glu_acc <- ifelse(MRS_long$percent_change_glu_acc > 0, "Increasing",
                                ifelse(MRS_long$percent_change_glu_acc < 0, "Decreasing", NA))

# 3. Convert both to factors 
MRS_long$traj_glu_prec <- as.factor(MRS_long$traj_glu_prec)
MRS_long$traj_glu_acc <- as.factor(MRS_long$traj_glu_acc)


#Descriptives
jmv::descriptives(data = MRS_long, vars = vars(sexe,   diagnostic_nick,education,  t1_age,   t2_age, 
                                               age_difference, t1_glu_prec, t2_glu_prec, hipp_difference,
                                               percent_change_glu_prec,t1_glu_acc, t2_glu_acc, percent_change_glu_acc, t1_hipp_e_tiv,
                                               slope_moca_raw, intercept_moca_raw, cluster_moca_raw, change_over_3_5,decliners,  
                                               memor_rappel_libre_nombre_reponses_correctes_t1,memor_rappel_libre_nombre_reponses_correctes_t2,
                                               visage_visages_score_rappel_differe_9_t1, visage_visages_score_rappel_differe_9_t2),
                  sd = TRUE, iqr = TRUE, skew = TRUE, kurt = TRUE)

jmv::descriptives(data = MRS_long, vars = vars(t1_cortical_thickness_dickson, t2_cortical_thickness_dickson, thickness_difference),
                  sd = TRUE, iqr = TRUE, skew = TRUE, kurt = TRUE)


jmv::descriptives(data = MRS_long, vars = vars(hippocampus_l_average_t1, 
                                               hippocampus_l_average_t2, 
                                               parietal_sup_l_average_t1, 
                                               parietal_sup_l_average_t2),
                  sd = TRUE, iqr = TRUE, skew = TRUE, kurt = TRUE)


# Completness of dataset 


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
##  MRS Variables
MRS_long$t1_glu_prec <- winsorize_iqr(MRS_long$t1_glu_prec, 1.5)
MRS_long$t2_glu_prec <- winsorize_iqr(MRS_long$t2_glu_prec, 1.5)
MRS_long$percent_change_glu_prec <- winsorize_iqr(MRS_long$percent_change_glu_prec, 1.5)
MRS_long$t1_glu_acc <- winsorize_iqr(MRS_long$t1_glu_acc, 1.5)
MRS_long$t2_glu_acc <- winsorize_iqr(MRS_long$t2_glu_acc, 1.5)
MRS_long$percent_change_glu_acc <- winsorize_iqr(MRS_long$percent_change_glu_acc, 1.5)
##  Structural MRI Variables
MRS_long$t1_hipp_e_tiv <- winsorize_iqr(MRS_long$t1_hipp_e_tiv, 1.5)
MRS_long$t2_hipp_e_tiv <- winsorize_iqr(MRS_long$t2_hipp_e_tiv, 1.5)
MRS_long$hipp_difference <- winsorize_iqr(MRS_long$hipp_difference, 1.5)
MRS_long$t1_cortical_thickness_dickson <- winsorize_iqr(MRS_long$t1_cortical_thickness_dickson, 1.5)
MRS_long$t2_cortical_thickness_dickson <- winsorize_iqr(MRS_long$t2_cortical_thickness_dickson, 1.5)
MRS_long$thickness_difference <- winsorize_iqr(MRS_long$thickness_difference, 1.5)
## fMRI Activation Variables
MRS_long$parietal_sup_l_average_t1 <- winsorize_iqr(MRS_long$parietal_sup_l_average_t1, 1.5)
MRS_long$parietal_sup_l_average_t2 <- winsorize_iqr(MRS_long$parietal_sup_l_average_t2, 1.5)
MRS_long$parietal_sup_l_diff <- winsorize_iqr(MRS_long$parietal_sup_l_diff, 1.5)
MRS_long$hippocampus_l_average_t1 <- winsorize_iqr(MRS_long$hippocampus_l_average_t1, 1.5)
MRS_long$hippocampus_l_average_t2 <- winsorize_iqr(MRS_long$hippocampus_l_average_t2, 1.5)
MRS_long$hipp_diff_activation <- winsorize_iqr(MRS_long$hipp_diff_activation, 1.5)
MRS_long$precuneus_t1_avg_act <- winsorize_iqr(MRS_long$precuneus_t1_avg_act, 1.5)
MRS_long$precuneus_t2_avg_act <- winsorize_iqr(MRS_long$precuneus_t2_avg_act, 1.5)
MRS_long$precuneus_diff <- winsorize_iqr(MRS_long$precuneus_diff, 1.5)
MRS_long$acc_t1_avg_activation <- winsorize_iqr(MRS_long$acc_t1_avg_activation, 1.5)
MRS_long$acc_t2_avg_activation <- winsorize_iqr(MRS_long$acc_t2_avg_activation, 1.5)
MRS_long$acc_diff <- winsorize_iqr(MRS_long$acc_diff, 1.5)
MRS_long$hipp_l_diff_activation <- ((MRS_long$hippocampus_l_average_t2 - MRS_long$hippocampus_l_average_t1)/MRS_long$hippocampus_l_average_t1)/MRS_long$age_difference
MRS_long$hipp_l_diff_activation <- winsorize_iqr(MRS_long$hipp_l_diff_activation, 1.5)



names(MRS_long)
### Characteristiques cliniques
#Table 1
#Table 1.1 - decliners
jmv::descriptives(data = MRS_long, vars = vars(t1_age,education,slope_moca_raw,age_difference),
                  sd = TRUE, iqr = TRUE, splitBy = decliners, skew = TRUE, kurt = TRUE)

t.test(t1_age ~ decliners, MRS_long)

chisq.test(MRS_long$sexe, y = MRS_long$decliners, correct = TRUE,
           simulate.p.value = FALSE, B = 2000)
table(MRS_long$sexe, MRS_long$decliners)

t.test(education ~ decliners, MRS_long)

chisq.test(MRS_long$diagnostic_nick, y = MRS_long$decliners, correct = TRUE,
           simulate.p.value = FALSE, B = 2000)
table(MRS_long$diagnostic_nick, MRS_long$decliners)

t.test(slope_moca_raw ~ decliners, MRS_long)

t.test(age_difference ~ decliners, MRS_long)


#Table 1.2 - diagnostick
jmv::descriptives(data = MRS_long, vars = vars(t1_age,education,slope_moca_raw,age_difference),
                  sd = TRUE, iqr = TRUE, splitBy = diagnostic_nick, skew = TRUE, kurt = TRUE)

summary(aov(t1_age ~ diagnostic_nick, MRS_long))

chisq.test(MRS_long$sexe, y = MRS_long$diagnostic_nick, correct = TRUE,
           simulate.p.value = FALSE, B = 2000)
table(MRS_long$sexe, MRS_long$diagnostic_nick)

summary(aov(education ~ diagnostic_nick, MRS_long))

summary(aov(slope_moca_raw ~ diagnostic_nick, MRS_long))

summary(aov(age_difference ~ diagnostic_nick, MRS_long))


TukeyHSD()

### table 2
jmv::descriptives(data = MRS_long, vars = vars(t1_glu_acc,t1_glu_prec,t1_hipp_e_tiv,t1_cortical_thickness_dickson,hippocampus_l_average_t1,
                                               hippocampus_avg_activation_t1, parietal_sup_l_average_t1),
                  sd = TRUE, iqr = TRUE, splitBy = decliners, skew = TRUE, kurt = TRUE)
t.test(t1_glu_acc ~ decliners, MRS_long)
t.test(t1_glu_prec ~ decliners, MRS_long)

t.test(t1_hipp_e_tiv ~ decliners, MRS_long)
t.test(t1_cortical_thickness_dickson ~ decliners, MRS_long)
t.test(hippocampus_l_average_t1 ~ decliners, MRS_long)
t.test(hippocampus_avg_activation_t1 ~ decliners, MRS_long)
t.test(parietal_sup_l_average_t1 ~ decliners, MRS_long)


## Table 2.2 neuroimaging and diagnostick
jmv::descriptives(data = MRS_long, vars = vars(t1_glu_acc,t1_glu_prec,t1_hipp_e_tiv,t1_cortical_thickness_dickson,hippocampus_l_average_t1,
                                               hippocampus_avg_activation_t1, parietal_sup_l_average_t1),
                  sd = TRUE, iqr = TRUE, splitBy = diagnostic_nick, skew = TRUE, kurt = TRUE)

anova_acc_diagnostick <- aov(t1_glu_acc ~ diagnostic_nick, MRS_long)
summary(anova_acc_diagnostick)
TukeyHSD(anova_acc_diagnostick)

summary(aov(t1_glu_prec ~ diagnostic_nick, MRS_long))

anova_hipp_struc_diagnostick <- aov(t1_hipp_e_tiv ~ diagnostic_nick, MRS_long)
summary(anova_hipp_struc_diagnostick)
TukeyHSD(anova_hipp_struc_diagnostick)

summary(aov(t1_cortical_thickness_dickson ~ diagnostic_nick, MRS_long))
summary(aov(hippocampus_l_average_t1 ~ diagnostic_nick, MRS_long))
summary(aov(hippocampus_avg_activation_t1 ~ diagnostic_nick, MRS_long))
summary(aov(parietal_sup_l_average_t1 ~ diagnostic_nick, MRS_long))



### Group analyses with lmm
library(lmerTest)
names(MRS_long)
table(MRS_long$diagnostic_nick)

MRS_long <- as.data.frame(MRS_long)
# Regrouper les groupes 'HC' et 'SCD' dans un seul groupe 'HC'
levels(MRS_long$diagnostic_nick)[levels(MRS_long$diagnostic_nick) == "SCD"] <- "HC"
table(MRS_long$diagnostic_nick)
# Regrouper les groupes 'HC' et 'SCD' dans un seul groupe 'SCD'
levels(MRS_long$diagnostic_nick)[levels(MRS_long$diagnostic_nick) %in% c("MCI", "SCD+")] <- "SCD+ and MCI"
table(MRS_long$diagnostic_nick)


## variables to convert 
MRS_prec_wide <- MRS_long[, c("pscid", "diagnostic_nick","decliners", "t1_age", "t1_glu_prec", "t2_glu_prec","age_difference")]
MRS_acc_wide <- MRS_long[, c("pscid", "diagnostic_nick","decliners", "t1_age", "t1_glu_acc", "t2_glu_acc","age_difference")]

# Convert to long format
MRS_prec_long <- reshape(MRS_prec_wide, direction = "long",varying = c("t1_glu_prec", "t2_glu_prec"), 
                         v.names = "glutamate_prec", timevar = "time",  times = c("T1", "T2"),               
                         idvar = "pscid")  
MRS_prec_long$time <- as.factor(MRS_prec_long$time)
MRS_prec_long$years_elapsed <- ifelse(MRS_prec_long$time == "T1", 0, MRS_prec_long$age_difference)

MRS_acc_long <- reshape(MRS_acc_wide, direction = "long",varying = c("t1_glu_acc", "t2_glu_acc"), 
                         v.names = "glutamate_acc", timevar = "time",  times = c("T1", "T2"),               
                         idvar = "pscid")   
MRS_acc_long$time <- as.factor(MRS_acc_long$time)
MRS_acc_long$years_elapsed <- ifelse(MRS_acc_long$time == "T1", 0, MRS_acc_long$age_difference)


# Prec Glu
lmm_prec <- lmer(glutamate_prec ~ years_elapsed * diagnostic_nick + decliners + t1_age + (1 | pscid), 
                   data = MRS_prec_long)
summary(lmm_prec)

prec_stats <- emtrends(lmm_prec, ~ diagnostic_nick, var = "years_elapsed")
summary(prec_stats, infer = TRUE)

# ACC Glu
lmm_acc <- lmer(glutamate_acc ~ years_elapsed * diagnostic_nick + decliners + t1_age + (1 | pscid), 
                 data = MRS_acc_long)
summary(lmm_acc)

acc_stats <- emtrends(lmm_acc, ~ diagnostic_nick, var = "years_elapsed")
summary(acc_stats, infer = TRUE)


names(MRS_long)
### Regression
## Structure
summary(lm(t2_hipp_e_tiv ~ t1_glu_prec, data = MRS_long))
summary(lm(t2_hipp_e_tiv ~ t1_glu_acc, data = MRS_long))
summary(lm(t2_cortical_thickness_dickson ~  t1_glu_prec, data = MRS_long))
summary(lm(t2_cortical_thickness_dickson  ~  t1_glu_acc, data = MRS_long))

## Activation
names(MRS_long)


summary(lm(hipp_l_diff_activation ~ t1_glu_prec, data = MRS_long))
summary(lm(hipp_l_diff_activation ~ t1_glu_acc, data = MRS_long))
summary(lm(parietal_sup_l_diff ~ t1_glu_prec, data = MRS_long))
summary(lm(parietal_sup_l_diff ~ t1_glu_acc, data = MRS_long))


## Cognition
summary(lm(slope_moca_scaled ~ t1_glu_prec, data = MRS_long))
summary(lm(slope_moca_scaled ~ t1_glu_acc, data = MRS_long))
summary(lm(slope_moca_scaled ~ t1_glu_prec + I(t1_glu_prec^2), data = MRS_long))
summary(lm(slope_moca_scaled ~ t1_glu_acc + I(t1_glu_acc^2), data = MRS_long))



## Between themselves
summary(lm(slope_moca_raw ~ hipp_l_diff_activation, data = MRS_long))
summary(lm(slope_moca_raw ~ hipp_l_diff_activation, data = MRS_long))
summary(lm(t2_hipp_e_tiv ~ hipp_l_diff_activation, data = MRS_long))
summary(lm(t2_cortical_thickness_dickson ~ hipp_l_diff_activation, data = MRS_long))

## Mediation
library(mediation)
library(lavaan)
# Center 
MRS_long$t1_glu_acc_c <- as.numeric(scale(MRS_long$t1_glu_acc, scale = FALSE))

sem_linear_model_acc <- 't2_hipp_e_tiv ~ a1 * t1_glu_acc_c
  slope_moca_scaled ~ cp1 * t1_glu_acc_c + b1 * t2_hipp_e_tiv
ind_eff := a1* b1'

fit_sem_lin_acc <- sem(sem_linear_model_acc,data = MRS_long, missing = "fiml",
                                 fixed.x   = FALSE, estimator = "MLR")

summary(fit_sem_lin_acc, standardized = TRUE, ci = TRUE)
parameterEstimates(fit_sem_lin_acc, standardized = TRUE, ci = TRUE)


## Sensirivity bootstrap 
fit_sem_lin_acc_boot <- sem(sem_linear_model_acc,data = MRS_long,missing = "fiml",
  fixed.x = FALSE,estimator = "ML",se = "bootstrap",bootstrap = 5000)



## precuneus
MRS_long$t1_glu_prec_c <- as.numeric(scale(MRS_long$t1_glu_prec, scale = FALSE))

sem_linear_model_prec <- 't2_hipp_e_tiv ~ a1 * t1_glu_prec_c
slope_moca_scaled ~ cp1 * t1_glu_prec_c + b1  * t2_hipp_e_tiv
  ind_structural := a1 * b1'

fit_sem_lin_prec <- sem(sem_linear_model_prec,data      = MRS_long,
  missing   = "fiml",fixed.x   = FALSE,estimator = "MLR")

summary(fit_sem_lin_prec, standardized = TRUE, ci = TRUE)
parameterEstimates(fit_sem_lin_prec, standardized = TRUE, ci = TRUE)



## Precuneus quadratic 
MRS_long$t1_glu_prec_c2 <- MRS_long$t1_glu_prec_c^2

sem_quadratic_model_prec <- 't2_hipp_e_tiv ~ a1 * t1_glu_prec_c
  slope_moca_scaled ~ cp1 * t1_glu_prec_c +cp2 * t1_glu_prec_c2 + b1  * t2_hipp_e_tiv
  ind_structural := a1 * b1'

fit_sem_quad_prec <- sem(sem_quadratic_model_prec,data      = MRS_long,
  missing   = "fiml",fixed.x   = FALSE,estimator = "MLR")

summary(fit_sem_quad_prec, standardized = TRUE, ci = TRUE)
parameterEstimates(fit_sem_quad_prec, standardized = TRUE, ci = TRUE)




############ Logistic regressions
library(pROC)

# Precuneus GLutamate
model_glu_prec <- glm(decliners ~ t1_glu_prec, data = MRS_long, family = "binomial")
summary(model_glu_prec)
roc_glu_prec <- roc(model_glu_prec$y, fitted(model_glu_prec))
auc(roc_glu_prec)
coords(roc_glu_prec, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")

# ACC Glutamate
model_glu_acc <- glm(decliners ~ t1_glu_acc, data = MRS_long, family = "binomial")
summary(model_glu_acc)
roc_glu_acc <- roc(model_glu_acc$y, fitted(model_glu_acc))
auc(roc_glu_acc)
coords(roc_glu_acc, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")

# Hippocampal Volume
model_struc_hip <- glm(decliners ~ t1_hipp_e_tiv, data = MRS_long, family = "binomial")
summary(model_struc_hip)
roc_struc_hip <- roc(model_struc_hip$y, fitted(model_struc_hip))
auc(roc_struc_hip)
coords(roc_struc_hip, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")

# Cortical Thickness
model_struc_thick <- glm(decliners ~ t1_cortical_thickness_dickson, data = MRS_long, family = "binomial")
summary(model_struc_thick)
roc_struc_thick <- roc(model_struc_thick$y, fitted(model_struc_thick))
auc(roc_struc_thick)
coords(roc_struc_thick, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")

# Hippocampal Activation
model_func_hip <- glm(decliners ~ hippocampus_l_average_t1, data = MRS_long, family = "binomial")
summary(model_func_hip)
roc_func_hip <- roc(model_func_hip$y, fitted(model_func_hip))
auc(roc_func_hip)
coords(roc_func_hip, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")

# Superior Parietal Activation
model_func_par <- glm(decliners ~ parietal_sup_l_average_t1, data = MRS_long, family = "binomial")
summary(model_func_par)
roc_func_par <- roc(model_func_par$y, fitted(model_func_par))
auc(roc_func_par)
coords(roc_func_par, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")


## Combined Models (2 variables)

#  Glu acc + parietal activaiton 
model_metab_func <- glm(decliners ~ t1_glu_acc + parietal_sup_l_average_t1, data = MRS_long, family = "binomial")
summary(model_metab_func)
roc_metab_func <- roc(model_metab_func$y, fitted(model_metab_func))
auc(roc_metab_func)
coords(roc_metab_func, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")

# Glu precuneus  + Hipp volume 
model_metab_struc_hip <- glm(decliners ~ t1_glu_prec + t1_hipp_e_tiv, data = MRS_long, family = "binomial")
summary(model_metab_struc_hip)
roc_metab_struc_hip <- roc(model_metab_struc_hip$y, fitted(model_metab_struc_hip))
auc(roc_metab_struc_hip)
coords(roc_metab_struc_hip, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")

# Glu precuneus + cortical thickness 
model_metab_struc_thick <- glm(decliners ~ t1_glu_prec + t1_cortical_thickness_dickson, data = MRS_long, family = "binomial")
summary(model_metab_struc_thick)
roc_metab_struc_thick <- roc(model_metab_struc_thick$y, fitted(model_metab_struc_thick))
auc(roc_metab_struc_thick)
coords(roc_metab_struc_thick, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")




## Combined models (3 variables )

model_3var_1 <- glm(decliners ~ t1_glu_prec + t1_hipp_e_tiv + parietal_sup_l_average_t1, data = MRS_long, family = "binomial")
summary(model_3var_1)
roc_3var_1 <- roc(model_3var_1$y, fitted(model_3var_1))
auc(roc_3var_1)
coords(roc_3var_1, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")


model_3var_2 <- glm(decliners ~ t1_glu_prec + t1_glu_acc + parietal_sup_l_average_t1, data = MRS_long, family = "binomial")
summary(model_3var_2)
roc_3var_2 <- roc(model_3var_2$y, fitted(model_3var_2))
auc(roc_3var_2)
coords(roc_3var_2, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")

model_3var_3 <- glm(decliners ~ t1_hipp_e_tiv + hippocampus_l_average_t1 + parietal_sup_l_average_t1, data = MRS_long, family = "binomial")
summary(model_3var_3)
roc_3var_3 <- roc(model_3var_3$y, fitted(model_3var_3))
auc(roc_3var_3)
coords(roc_3var_3, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="youden")




### Survival analysis
library(survival)
library(survminer)

MRS_long$decliners_numeric <- as.numeric(as.character(MRS_long$decliners))

# Cox model
cox_metab_func <- coxph(Surv(age_change_moca, decliners_numeric) ~ t1_glu_acc + parietal_sup_l_average_t1, data = MRS_long)

# 3. View the results
summary(cox_metab_func)



# Cox model
cox_metab_func <- coxph(Surv(age_change_moca, decliners_numeric) ~ t1_glu_prec, data = MRS_long)

# 3. View the results
summary(cox_metab_func)








