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

X2026_04_04_glu_moca_struc <- read_excel("C:/Users/okkam/Desktop/labo/article 2/Longitudinal_Multimodal_Data_CIMAQ/2026-04-04_glu_moca_struc.xlsx")
MRS_long <- X2026_04_04_glu_moca_struc

# Clean the column name
MRS_long <- janitor::clean_names(MRS_long)
names(MRS_long)
sapply(MRS_long,class)
# Convertir  en   factor
MRS_long$sexe <- as.factor(MRS_long$sexe)
table(MRS_long$sexe)
MRS_long$diagnostic_nick <- as.factor(MRS_long$diagnostic_nick)
levels(MRS_long$diagnostic_nick)
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

MRS_long$t1_glu_prec <- winsorize_iqr(MRS_long$t1_glu_prec, 1.5)
MRS_long$t2_glu_prec <- winsorize_iqr(MRS_long$t2_glu_prec, 1.5)
MRS_long$percent_change_glu_prec <- winsorize_iqr(MRS_long$percent_change_glu_prec, 1.5)
MRS_long$t1_glu_acc <- winsorize_iqr(MRS_long$t1_glu_acc, 1.5)
MRS_long$t2_glu_acc <- winsorize_iqr(MRS_long$t2_glu_acc, 1.5)
MRS_long$percent_change_glu_acc <- winsorize_iqr(MRS_long$percent_change_glu_acc, 1.5)
MRS_long$t1_hipp_e_tiv <- winsorize_iqr(MRS_long$t1_hipp_e_tiv, 1.5)
MRS_long$t2_hipp_e_tiv  <- winsorize_iqr(MRS_long$t2_hipp_e_tiv , 1.5)
MRS_long$hipp_difference  <- winsorize_iqr(MRS_long$hipp_difference , 1.5)
MRS_long$t2_hipp_e_tiv  <- winsorize_iqr(MRS_long$t2_hipp_e_tiv , 1.5)
MRS_long$t1_cortical_thickness_dickson    <- winsorize_iqr(MRS_long$t1_cortical_thickness_dickson, 1.5)
MRS_long$thickness_difference  <- winsorize_iqr(MRS_long$thickness_difference   , 1.5)
############################## ANALYSES PRINCIPALES ###########################
names(MRS_long)
sapply(MRS_long,class)
levels(MRS_long$diagnostic_nick)
levels(MRS_long$decliners)
levels(MRS_long$cluster_moca_raw)
table(MRS_long$traj_glu_prec)
table(MRS_long$traj_glu_acc)




# 1. ANOVA for ACC Percentage Change
anova_acc <- aov(percent_change_glu_acc ~ decliners, data = MRS_long)
summary(anova_acc)

# 2. ANOVA for Precuneus Percentage Change
anova_prec <- aov(percent_change_glu_prec ~ decliners, data = MRS_long)
summary(anova_prec)



############## Characterising Glu ################
jmv::descriptives(data = MRS_long, vars = vars(sexe, diagnostic_nick, education, t1_age, 
                                               t1_glu_prec, t2_glu_prec, t1_hipp_e_tiv,hipp_difference,
                                               percent_change_glu_prec, t1_glu_acc, t2_glu_acc, percent_change_glu_acc, 
                                               slope_moca_raw, intercept_moca_raw, cluster_moca_raw),
                  splitBy = "traj_glu_prec",  sd = TRUE, iqr = TRUE, skew = TRUE, kurt = TRUE)


jmv::descriptives(data = MRS_long,  vars = vars(sexe, diagnostic_nick, education, t1_age, 
                                                t1_glu_prec, t2_glu_prec, t1_hipp_e_tiv,hipp_difference,
                                                percent_change_glu_prec, t1_glu_acc, t2_glu_acc, percent_change_glu_acc, 
                                                slope_moca_raw, intercept_moca_raw, cluster_moca_raw),
                  splitBy = "traj_glu_acc", sd = TRUE, iqr = TRUE, skew = TRUE, kurt = TRUE)


## Chi qarre ###
help("chisq.test")
## glu between regions
chisq.test(MRS_long$traj_glu_prec, y = MRS_long$traj_glu_acc, correct = TRUE,
           simulate.p.value = FALSE, B = 2000)
table(MRS_long$traj_glu_prec, MRS_long$traj_glu_acc)

## glu between decliners in both regions
chisq.test(MRS_long$traj_glu_prec, y = MRS_long$decliners, correct = TRUE,
           simulate.p.value = FALSE, B = 2000)
table(MRS_long$traj_glu_prec, MRS_long$decliners)

chisq_acc_decliners <- chisq.test(MRS_long$traj_glu_acc, y = MRS_long$decliners, correct = TRUE,
           simulate.p.value = FALSE, B = 2000)
chisq_acc_decliners
table(MRS_long$traj_glu_acc, MRS_long$decliners)
chisq_acc_decliners$stdres
addmargins(table(MRS_long$traj_glu_acc, MRS_long$decliners))
fisher.test(MRS_long$traj_glu_acc, MRS_long$decliners)

## glu between clinical groups  in both regions
chisq.test(MRS_long$traj_glu_prec, y = MRS_long$diagnostic_nick, correct = TRUE,
           simulate.p.value = FALSE, B = 2000)
table(MRS_long$traj_glu_prec, MRS_long$diagnostic_nick)

chisq.test(MRS_long$traj_glu_acc, y = MRS_long$diagnostic_nick, correct = TRUE,
           simulate.p.value = FALSE, B = 2000)
table(MRS_long$traj_glu_acc, MRS_long$diagnostic_nick)

chisq.test(MRS_long$traj_glu_prec, y = MRS_long$cluster_moca_raw, correct = TRUE,
           simulate.p.value = FALSE, B = 2000)
table(MRS_long$traj_glu_prec, MRS_long$cluster_moca_raw)

chisq.test(MRS_long$traj_glu_acc, y = MRS_long$cluster_moca_raw, correct = TRUE,
           simulate.p.value = FALSE, B = 2000)
table(MRS_long$traj_glu_acc, MRS_long$diagnostic_nick)




## Glu and structure
sapply(MRS_long,class)
jmv::descriptives(data = MRS_long, vars = vars( t1_hipp_e_tiv ,hipp_difference,t1_cortical_thickness_dickson,thickness_difference), 
                  sd = TRUE, iqr = TRUE, skew = TRUE, kurt = TRUE, hist = TRUE)

## test with glutamate trajectories
## Hipp
t.test(t1_hipp_e_tiv ~ traj_glu_prec, data = MRS_long,
       alternative = c("two.sided"),
       mu = 0, var.equal = FALSE,
       conf.level = 0.95)

t.test(t1_hipp_e_tiv ~ traj_glu_acc, data = MRS_long,
       alternative = c("two.sided"),
       mu = 0, var.equal = FALSE,
       conf.level = 0.95)

t.test(hipp_difference ~ traj_glu_prec, data = MRS_long,
       alternative = c("two.sided"),
       mu = 0, var.equal = FALSE,
       conf.level = 0.95)

t.test(hipp_difference ~ traj_glu_acc, data = MRS_long,
       alternative = c("two.sided"),
       mu = 0, var.equal = FALSE,
       conf.level = 0.95)


## Thickness
t.test(t1_cortical_thickness_dickson ~ traj_glu_prec, data = MRS_long,
       alternative = c("two.sided"),
       mu = 0, var.equal = FALSE,
       conf.level = 0.95)

t.test(t1_cortical_thickness_dickson  ~ traj_glu_acc, data = MRS_long,
       alternative = c("two.sided"),
       mu = 0, var.equal = FALSE,
       conf.level = 0.95)

t.test(thickness_difference  ~ traj_glu_prec, data = MRS_long,
       alternative = c("two.sided"),
       mu = 0, var.equal = FALSE,
       conf.level = 0.95)

t.test(thickness_difference  ~ traj_glu_acc, data = MRS_long,
       alternative = c("two.sided"),
       mu = 0, var.equal = FALSE,
       conf.level = 0.95)


## GLM ###
library(pROC)
# Standard Logistic Regression - extremely high std error
logistic_decliners_t1acc <- glm(decliners ~ t1_glu_acc, 
                 data = MRS_long, 
                 family = binomial)
predicted_decliners_t1acc <- predict(logistic_decliners_t1acc, type = "response")
summary(logistic_decliners_t1acc)
exp(coef(logistic_decliners_t1acc)) 
exp(confint(logistic_decliners_t1acc))
roc_acc <- roc(logistic_decliners_t1acc$y, predicted_decliners_t1acc)
auc(roc_acc)
ci.auc(roc_acc)

# Precuneus
logistic_decliners_t1prec <- glm(decliners ~ t1_glu_prec, 
            data = MRS_long, 
            family = binomial)
predicted_decliners_t1prec <- predict(logistic_decliners_t1prec, type = "response")
summary(logistic_decliners_t1prec)
exp(coef(logistic_decliners_t1prec)) 
exp(confint(logistic_decliners_t1prec))
roc_prec <- roc(logistic_decliners_t1prec$y, predicted_decliners_t1prec)
auc(roc_prec)
ci.auc(roc_prec)

# Combined 
logistic_decliners_t1glu <- glm(decliners ~ t1_glu_prec +t1_glu_acc,
            data = MRS_long, 
            family = binomial)
predicted_decliners_t1glu <- predict(logistic_decliners_t1glu, type = "response")
summary(logistic_decliners_t1glu)
exp(coef(logistic_decliners_t1glu)) 
exp(confint(logistic_decliners_t1glu))










############## Decliners and moca slope ###################
jmv::descriptives(data = MRS_long, vars = vars(sexe,   diagnostic_nick,education,  t1_age,   t2_age, 
                                               age_difference, t1_glu_prec, t2_glu_prec, t1_hipp_e_tiv, hipp_difference,
                                               percent_change_glu_prec,t1_glu_acc, t2_glu_acc, percent_change_glu_acc, 
                                               slope_moca_raw, intercept_moca_raw, cluster_moca_raw, change_over_3_5,decliners,  
                                               memor_rappel_libre_nombre_reponses_correctes_t1,memor_rappel_libre_nombre_reponses_correctes_t2,
                                               visage_visages_score_rappel_differe_9_t1, visage_visages_score_rappel_differe_9_t2),
                  , splitBy = "decliners", sd = TRUE, iqr = TRUE)



# Regression
lin_moca_t1glu_acc <- lm(slope_moca_raw ~ t1_glu_acc + t1_age + sexe, data = MRS_long)
summary(lin_moca_t1glu_acc)
lin_moca_t1glu_prec <- lm(slope_moca_raw ~ t1_glu_prec + t1_age + sexe, data = MRS_long)
summary(lin_moca_t1glu_prec)
# The model with the linear term AND the quadratic term
lin_moca_t1glu_prec_quad <- lm(slope_moca_raw ~ t1_glu_prec + I(t1_glu_prec^2) + t1_age + sexe, data = MRS_long)
summary(lin_moca_t1glu_prec_quad)
lin_moca_t1glu_acc_quad <- lm(slope_moca_raw ~ t1_glu_acc + I(t1_glu_acc^2) + t1_age + sexe, data = MRS_long)
summary(lin_moca_t1glu_acc_quad)


lin_moca_t2glu_acc <- lm(slope_moca_raw ~ t2_glu_acc + t1_age + sexe, data = MRS_long)
summary(lin_moca_t2glu_acc)
lin_moca_t2glu_prec <- lm(slope_moca_raw ~ t2_glu_prec + t1_age + sexe, data = MRS_long)
summary(lin_moca_t2glu_prec)


## regression structure
lin_t2hipp_t1glu_acc <- lm(t2_hipp_e_tiv ~ t1_glu_acc + t1_age + sexe, data = MRS_long)
summary(lin_t2hipp_t1glu_acc)
lin_t2hipp_t1glu_acc <- lm(t2_cortical_thickness_dickson ~ t1_glu_acc + t1_age + sexe, data = MRS_long)
summary(lin_t2hipp_t1glu_acc)
lin_t2thick_t1glu_prec <- lm(t2_hipp_e_tiv ~ t1_glu_prec + t1_age + sexe, data = MRS_long)
summary(lin_t2thick_t1glu_prec)
lin_t2thick_t1glu_prec <- lm(t2_cortical_thickness_dickson ~ t1_glu_prec + t1_age + sexe, data = MRS_long)
summary(lin_t2thick_t1glu_prec)

lin_t2hipp_t1glu_acc <- lm(t2_hipp_e_tiv ~ t1_glu_acc, data = MRS_long)
summary(lin_t2hipp_t1glu_acc)
lin_t2thick_t1glu_prec <- lm(t2_cortical_thickness_dickson ~ t1_glu_prec, data = MRS_long)
summary(lin_t2thick_t1glu_prec)





# --- MEAN-CENTERING PREDICTORS ---
MRS_long$t1_glu_acc_c <- scale(MRS_long$t1_glu_acc, center = TRUE, scale = FALSE)
MRS_long$t1_glu_prec_c <- scale(MRS_long$t1_glu_prec, center = TRUE, scale = FALSE)

# 1. T2 Hippocampus ~ T1 ACC Glutamate
quad_t2hipp_t1acc <- lm(t2_hipp_e_tiv ~ t1_glu_acc_c + I(t1_glu_acc_c^2) + t1_age + sexe, data = MRS_long)
summary(quad_t2hipp_t1acc)

# 2. T2 Cortical Thickness ~ T1 ACC Glutamate
quad_t2thick_t1acc <- lm(t2_cortical_thickness_dickson ~ t1_glu_acc_c + I(t1_glu_acc_c^2) + t1_age + sexe, data = MRS_long)
summary(quad_t2thick_t1acc)

# 3. T2 Hippocampus ~ T1 Precuneus Glutamate
quad_t2hipp_t1prec <- lm(t2_hipp_e_tiv ~ t1_glu_prec_c + I(t1_glu_prec_c^2) + t1_age + sexe, data = MRS_long)
summary(quad_t2hipp_t1prec)

# 4. T2 Cortical Thickness ~ T1 Precuneus Glutamate
quad_t2thick_t1prec <- lm(t2_cortical_thickness_dickson ~ t1_glu_prec_c + I(t1_glu_prec_c^2) + t1_age + sexe, data = MRS_long)
summary(quad_t2thick_t1prec)



# --- MEAN-CENTERING PREDICTORS ---
MRS_long$t2_hipp_c <- scale(MRS_long$t2_hipp_e_tiv, center = TRUE, scale = FALSE)
MRS_long$t2_thick_c <- scale(MRS_long$t2_cortical_thickness_dickson, center = TRUE, scale = FALSE)

# 1. T1 ACC Glutamate ~ T2 Hippocampus
flip_t1acc_t2hipp <- lm(t1_glu_acc ~ t2_hipp_c + I(t2_hipp_c^2), data = MRS_long)
summary(flip_t1acc_t2hipp)

# 2. T1 ACC Glutamate ~ T2 Cortical Thickness
flip_t1acc_t2thick <- lm(t1_glu_acc ~ t2_thick_c + I(t2_thick_c^2), data = MRS_long)
summary(flip_t1acc_t2thick)

# 3. T1 Precuneus Glutamate ~ T2 Hippocampus
flip_t1prec_t2hipp <- lm(t1_glu_prec ~ t2_hipp_c + I(t2_hipp_c^2), data = MRS_long)
summary(flip_t1prec_t2hipp)

# 4. T1 Precuneus Glutamate ~ T2 Cortical Thickness
flip_t1prec_t2thick <- lm(t1_glu_prec ~ t2_thick_c + I(t2_thick_c^2), data = MRS_long)
summary(flip_t1prec_t2thick)



############# test with decliners and sexe  #######
t.test(hipp_difference ~ decliners, data = MRS_long,
       alternative = c("two.sided"),
       mu = 0, var.equal = FALSE,
       conf.level = 0.95)

t.test(t1_hipp_e_tiv ~ decliners, data = MRS_long,
       alternative = c("two.sided"),
       mu = 0, var.equal = FALSE,
       conf.level = 0.95)

t.test(hipp_difference ~ sexe, data = MRS_long,
       alternative = c("two.sided"),
       mu = 0, var.equal = FALSE,
       conf.level = 0.95)


## test ANOVAS   #####
help("aov")
help('TukeyHSD')

summary(aov(hipp_difference ~ diagnostic_nick, data = MRS_long))
summary(aov(hipp_difference ~ cluster_moca_raw , data = MRS_long))
t1_hipp_diagnostic_nick <- aov(t1_hipp_e_tiv ~ diagnostic_nick, data = MRS_long)
TukeyHSD(t1_hipp_diagnostic_nick)
t1_hipp_cluster <- aov(t1_hipp_e_tiv ~ cluster_moca_raw , data = MRS_long)
TukeyHSD(t1_hipp_cluster)












