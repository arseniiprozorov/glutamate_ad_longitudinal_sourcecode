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






library(dplyr)
library(tidyr)

# ---------------------------------------------------------
# Step 1: Group the relative timeline and average adjacent visits
# ---------------------------------------------------------
binned_cascade <- MRS_relative_discrete %>%
  filter(relative_visit >= -2 & relative_visit <= 2) %>%
  mutate(
    sequence_window = case_when(
      relative_visit %in% c(-2, -1) ~ "T_Pre",
      relative_visit == 0           ~ "T_0",
      relative_visit %in% c(1, 2)   ~ "T_Post"
    )
  ) %>%
  group_by(participant_id, sequence_window) %>%
  # Average the values for participants with multiple scans in a window
  summarise(across(c(norm_cortical_thickness_dickerson,
                     norm_sum_hippocampus,
                     hippocampus_avg_act,
                     parietal_sup_l_act,
                     glu_precuneus,
                     glu_acc,
                     moca_score_total_30),
                   ~ mean(., na.rm = TRUE)), .groups = "drop") %>%
  # Clean up NaN values created by averaging NAs
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA, .))) %>%
  
  # ---------------------------------------------------------
# Step 2: Reshape to wide format for direct linear regressions
# ---------------------------------------------------------
pivot_wider(
  names_from = sequence_window,
  values_from = c(norm_cortical_thickness_dickerson,
                  norm_sum_hippocampus,
                  hippocampus_avg_act,
                  parietal_sup_l_act,
                  glu_precuneus,
                  glu_acc,
                  moca_score_total_30)
)

# ---------------------------------------------------------
# Step 3: Run the Core Sequence Models
# ---------------------------------------------------------

print("=== ACT I: UPSTREAM (T_Pre Structure -> T_0 Metabolism & Activation) ===")

# Does historical cortical thickness dictate the index baseline?
print("--- T_Pre Cortical Thickness -> T_0 Glutamate ---")
summary(lm(glu_precuneus_T_0 ~ norm_cortical_thickness_dickerson_T_Pre, data = binned_cascade))
print("--- T_Pre Cortical Thickness -> T_0 Parietal Activation ---")
summary(lm(parietal_sup_l_act_T_0 ~ norm_cortical_thickness_dickerson_T_Pre, data = binned_cascade))

# Does historical hippocampal volume dictate the index baseline?
print("--- T_Pre Hippocampal Volume -> T_0 Glutamate ---")
summary(lm(glu_precuneus_T_0 ~ norm_sum_hippocampus_T_Pre, data = binned_cascade))
print("--- T_Pre Hippocampal Volume -> T_0 Hippocampal Activation ---")
summary(lm(hippocampus_avg_act_T_0 ~ norm_sum_hippocampus_T_Pre, data = binned_cascade))


print("=== ACT II: DOWNSTREAM (T_0 Metabolism -> T_Post Structure, Activation & Cognition) ===")

# Does the index glutamate surge protect future structure?
print("--- T_0 Glutamate -> T_Post Cortical Thickness ---")
summary(lm(norm_cortical_thickness_dickerson_T_Post ~ glu_precuneus_T_0, data = binned_cascade))
print("--- T_0 Glutamate -> T_Post Hippocampal Volume ---")
summary(lm(norm_sum_hippocampus_T_Post ~ glu_precuneus_T_0, data = binned_cascade))

# Does the index glutamate surge predict future functional network states?
print("--- T_0 Glutamate -> T_Post Parietal Activation ---")
summary(lm(parietal_sup_l_act_T_Post ~ glu_precuneus_T_0, data = binned_cascade))
print("--- T_0 Glutamate -> T_Post Hippocampal Activation ---")
summary(lm(hippocampus_avg_act_T_Post ~ glu_precuneus_T_0, data = binned_cascade))

# Does the index glutamate surge preserve global cognition?
print("--- T_0 Glutamate -> T_Post MoCA Total ---")
summary(lm(moca_score_total_30_T_Post ~ glu_precuneus_T_0, data = binned_cascade))




# ---------------------------------------------------------
# Step 3: Run the Core Sequence Models (With Autoregressive Controls)
# ---------------------------------------------------------

print("=== ACT I: UPSTREAM (T_Pre Structure -> T_0 Metabolism & Activation) ===")
# Note: Act I is already an unconfounded historical prediction model, as T_Pre occurs 
# chronologically before T_0, and historical glutamate data is not present to control for.

print("--- T_Pre Cortical Thickness -> T_0 Glutamate ---")
summary(lm(glu_precuneus_T_0 ~ norm_cortical_thickness_dickerson_T_Pre, data = binned_cascade))

print("--- T_Pre Cortical Thickness -> T_0 Parietal Activation ---")
summary(lm(parietal_sup_l_act_T_0 ~ norm_cortical_thickness_dickerson_T_Pre, data = binned_cascade))

print("--- T_Pre Hippocampal Volume -> T_0 Glutamate ---")
summary(lm(glu_precuneus_T_0 ~ norm_sum_hippocampus_T_Pre, data = binned_cascade))

print("--- T_Pre Hippocampal Volume -> T_0 Hippocampal Activation ---")
summary(lm(hippocampus_avg_act_T_0 ~ norm_sum_hippocampus_T_Pre, data = binned_cascade))


print("=== ACT II: DOWNSTREAM WITH HISTORICAL CONTROLS (Controlling for T_Pre) ===")
# These models test if T0 Glutamate predicts change relative to the absolute pre-disease baseline

print("--- T_0 Glutamate -> T_Post Cortical Thickness [Ctrl: T_Pre Thickness] ---")
summary(lm(norm_cortical_thickness_dickerson_T_Post ~ glu_precuneus_T_0 + norm_cortical_thickness_dickerson_T_Pre, data = binned_cascade))

print("--- T_0 Glutamate -> T_Post Hippocampal Volume [Ctrl: T_Pre Volume] ---")
summary(lm(norm_sum_hippocampus_T_Post ~ glu_precuneus_T_0 + norm_sum_hippocampus_T_Pre, data = binned_cascade))

print("--- T_0 Glutamate -> T_Post Parietal Activation [Ctrl: T_Pre Activation] ---")
summary(lm(parietal_sup_l_act_T_Post ~ glu_precuneus_T_0 + parietal_sup_l_act_T_Pre, data = binned_cascade))

print("--- T_0 Glutamate -> T_Post Hippocampal Activation [Ctrl: T_Pre Activation] ---")
summary(lm(hippocampus_avg_act_T_Post ~ glu_precuneus_T_0 + hippocampus_avg_act_T_Pre, data = binned_cascade))

print("--- T_0 Glutamate -> T_Post MoCA Total [Ctrl: T_Pre MoCA] ---")
summary(lm(moca_score_total_30_T_Post ~ glu_precuneus_T_0 + moca_score_total_30_T_Pre, data = binned_cascade))


print("=== ACT II: DOWNSTREAM WITH CONCURRENT CONTROLS (Controlling for T_0) ===")
# These models test if T0 Glutamate predicts subsequent change relative to the exact moment of the scan

print("--- T_0 Glutamate -> T_Post Cortical Thickness [Ctrl: T_0 Thickness] ---")
summary(lm(norm_cortical_thickness_dickerson_T_Post ~ glu_precuneus_T_0 + norm_cortical_thickness_dickerson_T_0, data = binned_cascade))

print("--- T_0 Glutamate -> T_Post Parietal Activation [Ctrl: T_0 Activation] ---")
summary(lm(parietal_sup_l_act_T_Post ~ glu_precuneus_T_0 + parietal_sup_l_act_T_0, data = binned_cascade))

print("--- T_0 Glutamate -> T_Post MoCA Total [Ctrl: T_0 MoCA] ---")
summary(lm(moca_score_total_30_T_Post ~ glu_precuneus_T_0 + moca_score_total_30_T_0, data = binned_cascade))





