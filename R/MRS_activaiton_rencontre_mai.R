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



### Regression
names(MRS_long)
MRS_long$hipp_l_diff_activation <- ((MRS_long$hippocampus_l_average_t2 - MRS_long$hippocampus_l_average_t1)/MRS_long$hippocampus_l_average_t1)/MRS_long$age_difference
MRS_long$hipp_l_diff_activation <- winsorize_iqr(MRS_long$hipp_l_diff_activation, 1.5)


summary(lm(hipp_diff_activation ~ t1_glu_acc, data = MRS_long))
summary(lm(hipp_l_diff_activation ~ t1_glu_prec + t1_age, data = MRS_long))
summary(lm(hipp_l_diff_activation ~ t1_glu_acc, data = MRS_long))

summary(lm(hipp_l_diff_activation ~ t1_glu_prec, data = MRS_long))
summary(lm(hipp_l_diff_activation ~ t1_glu_acc, data = MRS_long))
summary(lm(parietal_sup_l_diff ~ t1_glu_prec, data = MRS_long))
summary(lm(parietal_sup_l_diff ~ t1_glu_acc, data = MRS_long))



summary(lm(hipp_l_diff_activation ~ hipp_difference + I(hipp_difference^2), data = MRS_long))
summary(lm(hipp_l_diff_activation ~ thickness_difference + I(thickness_difference^2), data = MRS_long))
summary(lm(hipp_difference ~  hipp_l_diff_activation  + I(hipp_l_diff_activation^2), data = MRS_long))
summary(lm(thickness_difference  ~  hipp_l_diff_activation + I(hipp_l_diff_activation^2), data = MRS_long))






