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



### Composite scores

## Glutamate 
MRS_long$z_glu_prec_t1 <- as.numeric(scale(MRS_long$t1_glu_prec)) 
MRS_long$z_glu_acc_t1  <- as.numeric(scale(MRS_long$t1_glu_acc))
MRS_long$z_glu_prec_t2 <- as.numeric(scale(MRS_long$t2_glu_prec))
MRS_long$z_glu_acc_t2  <- as.numeric(scale(MRS_long$t2_glu_acc))
MRS_long$comp_glu_t1 <- (MRS_long$z_glu_prec_t1 + MRS_long$z_glu_acc_t1) / 2
MRS_long$comp_glu_t2 <- (MRS_long$z_glu_prec_t2 + MRS_long$z_glu_acc_t2) / 2
MRS_long$comp_glu_change <- MRS_long$comp_glu_t2 - MRS_long$comp_glu_t1

# Structure
MRS_long$z_thick_t1 <- as.numeric(scale(MRS_long$t1_cortical_thickness_dickson))
MRS_long$z_hipp_t1  <- as.numeric(scale(MRS_long$t1_hipp_e_tiv))
MRS_long$z_thick_t2 <- as.numeric(scale(MRS_long$t2_cortical_thickness_dickson))
MRS_long$z_hipp_t2  <- as.numeric(scale(MRS_long$t2_hipp_e_tiv))
MRS_long$comp_struct_t1 <- (MRS_long$z_thick_t1 + MRS_long$z_hipp_t1) / 2
MRS_long$comp_struct_t2 <- (MRS_long$z_thick_t2 + MRS_long$z_hipp_t2) / 2
MRS_long$comp_struct_change <- MRS_long$comp_struct_t2 - MRS_long$comp_struct_t1

# Activation
MRS_long$z_act_hipp_t1 <- as.numeric(scale(MRS_long$hippocampus_avg_activation_t1))
MRS_long$z_act_par_t1  <- as.numeric(scale(MRS_long$parietal_sup_l_average_t1))
MRS_long$z_act_hipp_t2 <- as.numeric(scale(MRS_long$hippocampus_avg_activation_t2))
MRS_long$z_act_par_t2  <- as.numeric(scale(MRS_long$parietal_sup_l_average_t2))
MRS_long$comp_act_t1 <- (MRS_long$z_act_hipp_t1 + MRS_long$z_act_par_t1) / 2
MRS_long$comp_act_t2 <- (MRS_long$z_act_hipp_t2 + MRS_long$z_act_par_t2) / 2
MRS_long$comp_act_change <- MRS_long$comp_act_t2 - MRS_long$comp_act_t1

# Cognition
MRS_long$z_moca_t1 <- as.numeric(scale(MRS_long$intercept_moca_raw))
MRS_long$z_moca_slope <- as.numeric(scale(MRS_long$slope_moca_raw))

##########//////######## Reproducing article 1 ###########################
names(MRS_long)
#### Centering and Squaring for MRS_long Polynomial Models
MRS_long$t1_hipp_e_tiv_c <- scale(MRS_long$t1_hipp_e_tiv, center = TRUE, scale = FALSE)
MRS_long$t1_hipp_e_tiv_sq <- MRS_long$t1_hipp_e_tiv_c^2
MRS_long$t1_cortical_thickness_dickson_c <- scale(MRS_long$t1_cortical_thickness_dickson, center = TRUE, scale = FALSE)
MRS_long$t1_cortical_thickness_dickson_sq <- MRS_long$t1_cortical_thickness_dickson_c^2
MRS_long$parietal_sup_l_average_t1_c <- scale(MRS_long$parietal_sup_l_average_t1, center = TRUE, scale = FALSE)
MRS_long$parietal_sup_l_average_t1_sq <- MRS_long$parietal_sup_l_average_t1_c^2
MRS_long$hippocampus_l_average_t1_c <- scale(MRS_long$hippocampus_l_average_t1, center = TRUE, scale = FALSE)
MRS_long$hippocampus_l_average_t1_sq <- MRS_long$hippocampus_l_average_t1_c^2
MRS_long$hippocampus_avg_activation_t1_c <- scale(MRS_long$hippocampus_avg_activation_t1, center = TRUE, scale = FALSE)
MRS_long$hippocampus_avg_activation_t1_sq <- MRS_long$hippocampus_avg_activation_t1_c^2
MRS_long$t1_glu_acc_c <- scale(MRS_long$t1_glu_acc, center = TRUE, scale = FALSE)
MRS_long$t1_glu_acc_sq <- MRS_long$t1_glu_acc_c^2
MRS_long$t1_glu_prec_c <- scale(MRS_long$t1_glu_prec, center = TRUE, scale = FALSE)
MRS_long$t1_glu_prec_sq <- MRS_long$t1_glu_prec_c^2
MRS_long$slope_moca_raw_c <- scale(MRS_long$slope_moca_raw, center = TRUE, scale = FALSE)
MRS_long$slope_moca_raw_sq <- MRS_long$slope_moca_raw_c^2


######################### Polynomial analyses  ####################################

#### structure and memory ######


# Model 1 t1_glu_acc ~ hipp mean
model_lin_acc_hipp_mean <- lm(t1_glu_acc ~ t1_hipp_e_tiv_c, data = MRS_long)
summary(model_lin_acc_hipp_mean)
AIC(model_lin_acc_hipp_mean)
model_quad_acc_hipp_mean <- lm(t1_glu_acc ~ t1_hipp_e_tiv_c + t1_hipp_e_tiv_sq, data = MRS_long)
summary(model_quad_acc_hipp_mean)
AIC(model_quad_acc_hipp_mean)


# Model 2 t1_glu_prec ~ hipp mean 
model_lin_prec_hipp_mean <- lm(t1_glu_prec ~ t1_hipp_e_tiv_c, data = MRS_long)
summary(model_lin_prec_hipp_mean)
AIC(model_lin_prec_hipp_mean)
model_quad_prec_hipp_mean <- lm(t1_glu_prec ~ t1_hipp_e_tiv_c + t1_hipp_e_tiv_sq, data = MRS_long)
summary(model_quad_prec_hipp_mean)
AIC(model_quad_prec_hipp_mean)


# Model 3 t1_glu_acc ~ thickness
model_lin_acc_thick <- lm(t1_glu_acc ~ t1_cortical_thickness_dickson_c, data = MRS_long)
summary(model_lin_acc_thick)
AIC(model_lin_acc_thick)
model_quad_acc_thick <- lm(t1_glu_acc ~ t1_cortical_thickness_dickson_c + t1_cortical_thickness_dickson_sq, data = MRS_long)
summary(model_quad_acc_thick)
AIC(model_quad_acc_thick)


# Model 4 t1_glu_prec ~ thickness
model_lin_prec_thick <- lm(t1_glu_prec ~ t1_cortical_thickness_dickson_c, data = MRS_long)
summary(model_lin_prec_thick)
AIC(model_lin_prec_thick)
model_quad_prec_thick <- lm(t1_glu_prec ~ t1_cortical_thickness_dickson_c + t1_cortical_thickness_dickson_sq, data = MRS_long)
summary(model_quad_prec_thick)
AIC(model_quad_prec_thick)


#### Activaiton ######

## Model 5 ACC ~ activaiton parietal
model_lin_acc_sup_act_rev <- lm(t1_glu_acc ~ parietal_sup_l_average_t1, data = MRS_long)
summary(model_lin_acc_sup_act_rev)
AIC(model_lin_acc_sup_act_rev)

model_quad_acc_sup_act_rev <- lm(t1_glu_acc ~ parietal_sup_l_average_t1 + parietal_sup_l_average_t1_sq, data = MRS_long)
summary(model_quad_acc_sup_act_rev)
AIC(model_quad_acc_sup_act_rev)

## Model 6 precuneus ~ activaiton parietal
model_lin_prec_sup_act_rev <- lm(t1_glu_prec ~ parietal_sup_l_average_t1, data = MRS_long)
summary(model_lin_prec_sup_act_rev)
AIC(model_lin_prec_sup_act_rev)

model_quad_prec_sup_act_rev <- lm(t1_glu_prec ~ parietal_sup_l_average_t1 + parietal_sup_l_average_t1_sq, data = MRS_long)
summary(model_quad_prec_sup_act_rev)
AIC(model_quad_prec_sup_act_rev)


# Model 7 acc ~ hipp activation
model_lin_acc_hipp_act <- lm(t1_glu_acc ~ hippocampus_l_average_t1_c, data = MRS_long)
summary(model_lin_acc_hipp_act)
AIC(model_lin_acc_hipp_act)
model_quad_acc_hipp_act <- lm(t1_glu_acc ~ hippocampus_l_average_t1_c + hippocampus_l_average_t1_sq, data = MRS_long)
summary(model_quad_acc_hipp_act)
AIC(model_quad_acc_hipp_act)

# Model 8 precuneus  ~ hipp activation
model_lin_prec_hipp_act <- lm(t1_glu_prec ~ hippocampus_l_average_t1_c, data = MRS_long)
summary(model_lin_prec_hipp_act)
AIC(model_lin_prec_hipp_act)
model_quad_prec_hipp_act <- lm(t1_glu_prec ~ hippocampus_l_average_t1_c + hippocampus_l_average_t1_sq, data = MRS_long)
summary(model_quad_prec_hipp_act)
AIC(model_quad_prec_hipp_act)

# Model 7 (Mean) acc ~ hipp activation mean
model_lin_acc_hipp_act_mean <- lm(t1_glu_acc ~ hippocampus_avg_activation_t1_c, data = MRS_long)
summary(model_lin_acc_hipp_act_mean)
AIC(model_lin_acc_hipp_act_mean)
model_quad_acc_hipp_act_mean <- lm(t1_glu_acc ~ hippocampus_avg_activation_t1_c + hippocampus_avg_activation_t1_sq, data = MRS_long)
summary(model_quad_acc_hipp_act_mean)
AIC(model_quad_acc_hipp_act_mean)

# Model 8 (Mean) precuneus  ~ hipp activation mean
model_lin_prec_hipp_act_mean <- lm(t1_glu_prec ~ hippocampus_avg_activation_t1_c, data = MRS_long)
summary(model_lin_prec_hipp_act_mean)
AIC(model_lin_prec_hipp_act_mean)
model_quad_prec_hipp_act_mean <- lm(t1_glu_prec ~ hippocampus_avg_activation_t1_c + hippocampus_avg_activation_t1_sq, data = MRS_long)
summary(model_quad_prec_hipp_act_mean)
AIC(model_quad_prec_hipp_act_mean)

### Memory ###

# Model 9 memoria ~ ACC Glu
model_lin_acc_memor <- lm(slope_moca_raw ~ t1_glu_acc_c, data = MRS_long)
summary(model_lin_acc_memor)
AIC(model_lin_acc_memor)
model_quad_acc_memor <- lm(slope_moca_raw ~ t1_glu_acc_c + t1_glu_acc_sq, data = MRS_long)
summary(model_quad_acc_memor)
AIC(model_quad_acc_memor)

# Model 10 memoria ~ Precuneus Glu
model_lin_prec_memor <- lm(slope_moca_raw ~ t1_glu_prec_c, data = MRS_long)
summary(model_lin_prec_memor)
AIC(model_lin_prec_memor)
model_quad_prec_memor <- lm(slope_moca_raw ~ t1_glu_prec_c + t1_glu_prec_sq, data = MRS_long)
summary(model_quad_prec_memor)
AIC(model_quad_prec_memor)



#################### Mediation ###############
library(mediation)
MRS_t1<- na.omit(MRS_long[c("comp_struct_t1", "comp_glu_t1", "comp_act_t1")])
## Reproducing Article 1 # Mean Structure -> Mean Glu -> Mean Activation
# M ~ X
model_m_comb <- lm(comp_glu_t1 ~ comp_struct_t1, data = MRS_t1)
summary(model_m_comb)

# Y ~ X + M
model_y_comb <- lm(comp_act_t1 ~ comp_struct_t1 + comp_glu_t1, data = MRS_t1)
summary(model_y_comb)

# Mediate
med_comb <- mediate(model_m_comb, model_y_comb, 
                    treat = "comp_struct_t1", mediator = "comp_glu_t1", 
                    boot = TRUE, sims = 5000)
summary(med_comb)
plot(med_comb)


#########################
names(MRS_long)

### Activation
summary(aov(hipp_diff_activation ~ diagnostic_nick, data = MRS_long))
summary(aov(hipp_diff_activation ~ cluster_moca_raw , data = MRS_long))
summary(aov(hipp_diff_activation ~ decliners , data = MRS_long))
summary(aov(hipp_diff_activation ~ sexe , data = MRS_long))
summary(aov(hipp_diff_activation ~ traj_glu_prec , data = MRS_long))
summary(aov(hipp_diff_activation ~ traj_glu_acc , data = MRS_long))


summary(aov(parietal_sup_l_diff ~ diagnostic_nick, data = MRS_long))
summary(aov(parietal_sup_l_diff ~ cluster_moca_raw , data = MRS_long))
summary(aov(parietal_sup_l_diff ~ decliners , data = MRS_long))
summary(aov(parietal_sup_l_diff ~ sexe , data = MRS_long))
summary(aov(hipp_diff_activation ~ traj_glu_prec , data = MRS_long))
summary(aov(hipp_diff_activation ~ traj_glu_acc , data = MRS_long))

summary(ttest)
summary(aov(parietal_sup_l_diff ~ cluster_moca_raw , data = MRS_long))
summary(aov(parietal_sup_l_diff ~ decliners , data = MRS_long))
summary(aov(parietal_sup_l_diff ~ sexe , data = MRS_long))

\