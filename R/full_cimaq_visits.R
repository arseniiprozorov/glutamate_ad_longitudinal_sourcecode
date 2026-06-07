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

X01_multimodal_longitudinal <- read_excel("C:/Users/okkam/Desktop/labo/article 2/Longitudinal_Multimodal_Data_CIMAQ/article longitudinal/long format/01_multimodal_longitudinal.xlsx")
MRS_long <- X01_multimodal_longitudinal


# Clean the column name
MRS_long <- janitor::clean_names(MRS_long)
names(MRS_long)
sapply(MRS_long,class)
# Convertir  en   factor
MRS_long$sex <- as.factor(MRS_long$sex)
table(MRS_long$sex)
MRS_long$diagnostic_nick <- as.factor(MRS_long$diagnostic_nick)
table(MRS_long$diagnostic_nick)

# Regrouper les groupes 'HC' et 'SCD' dans un seul groupe 'SCD'
levels(MRS_long$diagnostic_nick)[levels(MRS_long$diagnostic_nick) == "SCD"] <- "HC"
levels(MRS_long$diagnostic_nick)[levels(MRS_long$diagnostic_nick) == "HC"] <- "CU"
table(MRS_long$diagnostic_nick)



# Descriptives
#sink("descriptives_output_all_vars.txt")
jmv::descriptives(data = MRS_long, vars = vars(education, age, 
    moca_score_total_30, slope_score_moca, intercept_score_moca, x99453_analyse_plasma_ptau217, 
    glu_precuneus, glu_acc, norm_left_hippocampus, norm_sum_hippocampus, 
    cortical_thickness_dickerson, norm_cortical_thickness_dickerson, hippocampus_l_act, 
    hippocampus_r_act, hippocampus_avg_act, parietal_sup_l_act, parietal_sup_r_act, 
    parietal_sup_avg_act, precuneus_l_act, precuneus_avg_act, acc_avg_act),
  splitBy = vars(visit),sd = TRUE, iqr = TRUE, skew = TRUE, kurt = TRUE)
#sink()





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
MRS_long$glu_precuneus <- winsorize_iqr(MRS_long$glu_precuneus, 1.5)
MRS_long$glu_acc <- winsorize_iqr(MRS_long$glu_acc, 1.5)

##  Structural MRI Variables
MRS_long$norm_left_hippocampus <- winsorize_iqr(MRS_long$norm_left_hippocampus, 1.5)
MRS_long$norm_sum_hippocampus <- winsorize_iqr(MRS_long$norm_sum_hippocampus, 1.5)
MRS_long$cortical_thickness_dickerson <- winsorize_iqr(MRS_long$cortical_thickness_dickerson, 1.5)
MRS_long$norm_cortical_thickness_dickerson <- winsorize_iqr(MRS_long$norm_cortical_thickness_dickerson, 1.5)

## fMRI Activation Variables
MRS_long$hippocampus_l_act <- winsorize_iqr(MRS_long$hippocampus_l_act, 1.5)
MRS_long$hippocampus_r_act <- winsorize_iqr(MRS_long$hippocampus_r_act, 1.5)
MRS_long$hippocampus_avg_act <- winsorize_iqr(MRS_long$hippocampus_avg_act, 1.5)
MRS_long$parietal_sup_l_act <- winsorize_iqr(MRS_long$parietal_sup_l_act, 1.5)
MRS_long$parietal_sup_r_act <- winsorize_iqr(MRS_long$parietal_sup_r_act, 1.5)
MRS_long$parietal_sup_avg_act <- winsorize_iqr(MRS_long$parietal_sup_avg_act, 1.5)
MRS_long$precuneus_l_act <- winsorize_iqr(MRS_long$precuneus_l_act, 1.5)
MRS_long$precuneus_avg_act <- winsorize_iqr(MRS_long$precuneus_avg_act, 1.5)
MRS_long$acc_avg_act <- winsorize_iqr(MRS_long$acc_avg_act, 1.5)






###################### Analyses ####################
names(MRS_long)
#LMM# 
library(lmerTest)
model_glu_prec <- lmer( glu_precuneus ~ visit * diagnostic_nick + age + sex + education + (1 | participant_id), 
  data = MRS_long)
summary(model_glu_prec)

model_glu_acc <- lmer( glu_acc ~ visit * diagnostic_nick + age + sex + education + (1 | participant_id), 
                        data = MRS_long)
summary(model_glu_acc)


# Convert "t02" string into the number 2
MRS_long$visit_num <- as.numeric(gsub("t", "", MRS_long$visit))





