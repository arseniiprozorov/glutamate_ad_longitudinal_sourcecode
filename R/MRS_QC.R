library(readxl)
library(janitor)  
library(ggplot2)


Glu_T1_T2_QC_R_Precuneus <- read_excel("C:/Users/okkam/Desktop/labo/article 2/Longitudinal_Multimodal_Data_CIMAQ/Glu/QC/Glu_T1_T2_QC_R_Precuneus.xlsx")
Glu_T1_T2_QC_R_ACC <- read_excel("C:/Users/okkam/Desktop/labo/article 2/Longitudinal_Multimodal_Data_CIMAQ/Glu/QC/Glu_T1_T2_QC_R_ACC.xlsx")

MRS_QC_Precuneus <- Glu_T1_T2_QC_R_Precuneus
MRS_QC_ACC <- Glu_T1_T2_QC_R_ACC


# Clean the column name
MRS_QC_Precuneus <- janitor::clean_names(MRS_QC_Precuneus)
MRS_QC_ACC <- janitor::clean_names(MRS_QC_ACC)
names(MRS_QC_Precuneus)
names(MRS_QC_ACC)


### Rapport signal/bruit (vérifie que le SNR global est à peu près constant au cours du temps, et en particulier le SNR par métabolite majeur (NAA, Cr)
## Create databases
SNR_constance_Precuneus <- MRS_QC_Precuneus[, c("participant_id", "visite", "glu_snr", "naa_snr", "cr_snr","p_ch_snr")]
jmv::descriptives(data = SNR_constance_Precuneus, vars = vars(glu_snr, naa_snr, cr_snr, p_ch_snr), sd = TRUE,range = TRUE,skew = TRUE,kurt = TRUE)
jmv::descriptives(data = SNR_constance_Precuneus, vars = vars(glu_snr, naa_snr, cr_snr, p_ch_snr),splitBy = "visite", sd = TRUE,range = TRUE,skew = TRUE,kurt = TRUE)
glu_snr_precuneus_t1 <- SNR_constance_Precuneus$glu_snr[SNR_constance_Precuneus$visite == 1]
glu_snr_precuneus_t2 <- SNR_constance_Precuneus$glu_snr[SNR_constance_Precuneus$visite == 2]
naa_snr_precuneus_t1 <- SNR_constance_Precuneus$naa_snr[SNR_constance_Precuneus$visite == 1]
naa_snr_precuneus_t2 <- SNR_constance_Precuneus$naa_snr[SNR_constance_Precuneus$visite == 2]
cr_snr_precuneus_t1 <- SNR_constance_Precuneus$cr_snr[SNR_constance_Precuneus$visite == 1]
cr_snr_precuneus_t2 <- SNR_constance_Precuneus$cr_snr[SNR_constance_Precuneus$visite == 2]
pcho_snr_precuneus_t1 <- SNR_constance_Precuneus$p_ch_snr[SNR_constance_Precuneus$visite == 1]
pcho_snr_precuneus_t2 <- SNR_constance_Precuneus$p_ch_snr[SNR_constance_Precuneus$visite == 2]



SNR_constance_ACC <- MRS_QC_ACC[, c("participant_id", "visite", "glu_snr", "naa_snr", "cr_snr","p_ch_snr")]
jmv::descriptives(data = SNR_constance_ACC, vars = vars(glu_snr, naa_snr, cr_snr, p_ch_snr), sd = TRUE,range = TRUE,skew = TRUE,kurt = TRUE)
jmv::descriptives(data = SNR_constance_ACC, vars = vars(glu_snr, naa_snr, cr_snr, p_ch_snr),splitBy = "visite", sd = TRUE,range = TRUE,skew = TRUE,kurt = TRUE)
glu_snr_acc_t1 <- SNR_constance_ACC$glu_snr[SNR_constance_ACC$visite == 1]
glu_snr_acc_t2 <- SNR_constance_ACC$glu_snr[SNR_constance_ACC$visite == 2]
naa_snr_acc_t1 <- SNR_constance_ACC$naa_snr[SNR_constance_ACC$visite == 1]
naa_snr_acc_t2 <- SNR_constance_ACC$naa_snr[SNR_constance_ACC$visite == 2]
cr_snr_acc_t1 <- SNR_constance_ACC$cr_snr[SNR_constance_ACC$visite == 1]
cr_snr_acc_t2 <- SNR_constance_ACC$cr_snr[SNR_constance_ACC$visite == 2]
pcho_snr_acc_t1 <- SNR_constance_ACC$p_ch_snr[SNR_constance_ACC$visite == 1]
pcho_snr_acc_t2 <- SNR_constance_ACC$p_ch_snr[SNR_constance_ACC$visite == 2]


#  Run the Paired T-test
# Precuneus
t.test(glu_snr_precuneus_t1, glu_snr_precuneus_t2, paired = TRUE)
t.test(naa_snr_precuneus_t1, naa_snr_precuneus_t2, paired = TRUE)
t.test(cr_snr_precuneus_t1, cr_snr_precuneus_t2, paired = TRUE)
t.test(pcho_snr_precuneus_t1, pcho_snr_precuneus_t2, paired = TRUE)
# ACC
t.test(glu_snr_acc_t1, glu_snr_acc_t2, paired = TRUE)
t.test(naa_snr_acc_t1, naa_snr_acc_t2, paired = TRUE)
t.test(cr_snr_acc_t1, cr_snr_acc_t2, paired = TRUE)
t.test(pcho_snr_acc_t1, pcho_snr_acc_t2, paired = TRUE)



#### Stabilité des ratio NAA/Cr, Cho/Cr, Glu/Cr

### PRECUNEUS
## Create databases
Ratio_stabilite_Precuneus <- MRS_QC_Precuneus[, c("participant_id", "visite", "glu_cr_p_cr", "naa_cr_p_cr", "gpc_cr_p_cr", "p_ch_cr_p_cr")]

# Calculate Total Choline Ratio (GPC + PCh)
Ratio_stabilite_Precuneus$t_cho_cr_p_cr <- Ratio_stabilite_Precuneus$gpc_cr_p_cr + Ratio_stabilite_Precuneus$p_ch_cr_p_cr

# Run Descriptives
jmv::descriptives(data = Ratio_stabilite_Precuneus, vars = vars(glu_cr_p_cr, naa_cr_p_cr, t_cho_cr_p_cr), sd = TRUE, range = TRUE, skew = TRUE, kurt = TRUE)
jmv::descriptives(data = Ratio_stabilite_Precuneus, vars = vars(glu_cr_p_cr, naa_cr_p_cr, t_cho_cr_p_cr), splitBy = "visite", sd = TRUE, range = TRUE, skew = TRUE, kurt = TRUE)

# Extract vectors for T-Tests
glu_ratio_precuneus_t1 <- Ratio_stabilite_Precuneus$glu_cr_p_cr[Ratio_stabilite_Precuneus$visite == 1]
glu_ratio_precuneus_t2 <- Ratio_stabilite_Precuneus$glu_cr_p_cr[Ratio_stabilite_Precuneus$visite == 2]
naa_ratio_precuneus_t1 <- Ratio_stabilite_Precuneus$naa_cr_p_cr[Ratio_stabilite_Precuneus$visite == 1]
naa_ratio_precuneus_t2 <- Ratio_stabilite_Precuneus$naa_cr_p_cr[Ratio_stabilite_Precuneus$visite == 2]
tcho_ratio_precuneus_t1 <- Ratio_stabilite_Precuneus$t_cho_cr_p_cr[Ratio_stabilite_Precuneus$visite == 1]
tcho_ratio_precuneus_t2 <- Ratio_stabilite_Precuneus$t_cho_cr_p_cr[Ratio_stabilite_Precuneus$visite == 2]


### ACC
## Create databases
Ratio_stabilite_ACC <- MRS_QC_ACC[, c("participant_id", "visite", "glu_cr_p_cr", "naa_cr_p_cr", "gpc_cr_p_cr", "p_ch_cr_p_cr")]

# Calculate Total Choline Ratio (GPC + PCh)
Ratio_stabilite_ACC$t_cho_cr_p_cr <- Ratio_stabilite_ACC$gpc_cr_p_cr + Ratio_stabilite_ACC$p_ch_cr_p_cr

# Run Descriptives
jmv::descriptives(data = Ratio_stabilite_ACC, vars = vars(glu_cr_p_cr, naa_cr_p_cr, t_cho_cr_p_cr), sd = TRUE, range = TRUE, skew = TRUE, kurt = TRUE)
jmv::descriptives(data = Ratio_stabilite_ACC, vars = vars(glu_cr_p_cr, naa_cr_p_cr, t_cho_cr_p_cr), splitBy = "visite", sd = TRUE, range = TRUE, skew = TRUE, kurt = TRUE)

# Extract vectors for T-Tests
glu_ratio_acc_t1 <- Ratio_stabilite_ACC$glu_cr_p_cr[Ratio_stabilite_ACC$visite == 1]
glu_ratio_acc_t2 <- Ratio_stabilite_ACC$glu_cr_p_cr[Ratio_stabilite_ACC$visite == 2]
naa_ratio_acc_t1 <- Ratio_stabilite_ACC$naa_cr_p_cr[Ratio_stabilite_ACC$visite == 1]
naa_ratio_acc_t2 <- Ratio_stabilite_ACC$naa_cr_p_cr[Ratio_stabilite_ACC$visite == 2]
tcho_ratio_acc_t1 <- Ratio_stabilite_ACC$t_cho_cr_p_cr[Ratio_stabilite_ACC$visite == 1]
tcho_ratio_acc_t2 <- Ratio_stabilite_ACC$t_cho_cr_p_cr[Ratio_stabilite_ACC$visite == 2]


# Run the Paired T-test
# Precuneus
t.test(glu_ratio_precuneus_t1, glu_ratio_precuneus_t2, paired = TRUE)
t.test(naa_ratio_precuneus_t1, naa_ratio_precuneus_t2, paired = TRUE)
t.test(tcho_ratio_precuneus_t1, tcho_ratio_precuneus_t2, paired = TRUE)

# ACC
t.test(glu_ratio_acc_t1, glu_ratio_acc_t2, paired = TRUE)
t.test(naa_ratio_acc_t1, naa_ratio_acc_t2, paired = TRUE)
t.test(tcho_ratio_acc_t1, tcho_ratio_acc_t2, paired = TRUE)




#### Stabilité des Concentrations Absolues (mM)

### PRECUNEUS
## Create databases
# Note: We include cr_m_m here to check if the reference itself is stable
mM_stabilite_Precuneus <- MRS_QC_Precuneus[, c("participant_id", "visite", "glu_m_m", "naa_m_m", "cr_m_m", "gpc_m_m", "p_ch_m_m")]

# Calculate Total Choline mM (GPC + PCh)
mM_stabilite_Precuneus$t_cho_m_m <- mM_stabilite_Precuneus$gpc_m_m + mM_stabilite_Precuneus$p_ch_m_m

# Run Descriptives
jmv::descriptives(data = mM_stabilite_Precuneus, vars = vars(glu_m_m, naa_m_m, cr_m_m, t_cho_m_m), sd = TRUE, range = TRUE, skew = TRUE, kurt = TRUE)
jmv::descriptives(data = mM_stabilite_Precuneus, vars = vars(glu_m_m, naa_m_m, cr_m_m, t_cho_m_m), splitBy = "visite", sd = TRUE, range = TRUE, skew = TRUE, kurt = TRUE)

# Extract vectors for T-Tests
glu_mm_precuneus_t1 <- mM_stabilite_Precuneus$glu_m_m[mM_stabilite_Precuneus$visite == 1]
glu_mm_precuneus_t2 <- mM_stabilite_Precuneus$glu_m_m[mM_stabilite_Precuneus$visite == 2]

naa_mm_precuneus_t1 <- mM_stabilite_Precuneus$naa_m_m[mM_stabilite_Precuneus$visite == 1]
naa_mm_precuneus_t2 <- mM_stabilite_Precuneus$naa_m_m[mM_stabilite_Precuneus$visite == 2]

cr_mm_precuneus_t1  <- mM_stabilite_Precuneus$cr_m_m[mM_stabilite_Precuneus$visite == 1]
cr_mm_precuneus_t2  <- mM_stabilite_Precuneus$cr_m_m[mM_stabilite_Precuneus$visite == 2]

tcho_mm_precuneus_t1 <- mM_stabilite_Precuneus$t_cho_m_m[mM_stabilite_Precuneus$visite == 1]
tcho_mm_precuneus_t2 <- mM_stabilite_Precuneus$t_cho_m_m[mM_stabilite_Precuneus$visite == 2]


### ACC
## Create databases
mM_stabilite_ACC <- MRS_QC_ACC[, c("participant_id", "visite", "glu_m_m", "naa_m_m", "cr_m_m", "gpc_m_m", "p_ch_m_m")]

# Calculate Total Choline mM (GPC + PCh)
mM_stabilite_ACC$t_cho_m_m <- mM_stabilite_ACC$gpc_m_m + mM_stabilite_ACC$p_ch_m_m

# Run Descriptives
jmv::descriptives(data = mM_stabilite_ACC, vars = vars(glu_m_m, naa_m_m, cr_m_m, t_cho_m_m), sd = TRUE, range = TRUE, skew = TRUE, kurt = TRUE)
jmv::descriptives(data = mM_stabilite_ACC, vars = vars(glu_m_m, naa_m_m, cr_m_m, t_cho_m_m), splitBy = "visite", sd = TRUE, range = TRUE, skew = TRUE, kurt = TRUE)

# Extract vectors for T-Tests
glu_mm_acc_t1 <- mM_stabilite_ACC$glu_m_m[mM_stabilite_ACC$visite == 1]
glu_mm_acc_t2 <- mM_stabilite_ACC$glu_m_m[mM_stabilite_ACC$visite == 2]

naa_mm_acc_t1 <- mM_stabilite_ACC$naa_m_m[mM_stabilite_ACC$visite == 1]
naa_mm_acc_t2 <- mM_stabilite_ACC$naa_m_m[mM_stabilite_ACC$visite == 2]

cr_mm_acc_t1  <- mM_stabilite_ACC$cr_m_m[mM_stabilite_ACC$visite == 1]
cr_mm_acc_t2  <- mM_stabilite_ACC$cr_m_m[mM_stabilite_ACC$visite == 2]

tcho_mm_acc_t1 <- mM_stabilite_ACC$t_cho_m_m[mM_stabilite_ACC$visite == 1]
tcho_mm_acc_t2 <- mM_stabilite_ACC$t_cho_m_m[mM_stabilite_ACC$visite == 2]


# Run the Paired T-test
# Precuneus
t.test(glu_mm_precuneus_t1, glu_mm_precuneus_t2, paired = TRUE)
t.test(naa_mm_precuneus_t1, naa_mm_precuneus_t2, paired = TRUE)
t.test(cr_mm_precuneus_t1, cr_mm_precuneus_t2, paired = TRUE) # Important Check!
t.test(tcho_mm_precuneus_t1, tcho_mm_precuneus_t2, paired = TRUE)

# ACC
t.test(glu_mm_acc_t1, glu_mm_acc_t2, paired = TRUE)
t.test(naa_mm_acc_t1, naa_mm_acc_t2, paired = TRUE)
t.test(cr_mm_acc_t1, cr_mm_acc_t2, paired = TRUE) # Important Check!
t.test(tcho_mm_acc_t1, tcho_mm_acc_t2, paired = TRUE)



#### Coefficient de Variation (CV) par métabolite par sujet

# 1. Define Function to calculate Mean Intra-Subject CV (%)
# Formula: Mean of [ (Difference / sqrt(2)) / Mean ] * 100
CV_function <- function(t1, t2) {
  diffs <- abs(t1 - t2)
  means <- (t1 + t2) / 2
  individual_cvs <- ( (diffs / sqrt(2)) / means ) * 100
  return(mean(individual_cvs))}
# Precuneus mM
cv_prec_glu_mm <- CV_function(glu_mm_precuneus_t1, glu_mm_precuneus_t2)
cv_prec_naa_mm <- CV_function(naa_mm_precuneus_t1, naa_mm_precuneus_t2)
cv_prec_cho_mm <- CV_function(tcho_mm_precuneus_t1, tcho_mm_precuneus_t2)

summary(cv_prec_glu_mm)
summary(cv_prec_naa_mm)
summary(cv_prec_cho_mm)

# ACC mM
cv_acc_glu_mm <- CV_function(glu_mm_acc_t1, glu_mm_acc_t2)
cv_acc_naa_mm <- CV_function(naa_mm_acc_t1, naa_mm_acc_t2)
cv_acc_cho_mm <- mean( (abs(tcho_mm_acc_t1 - tcho_mm_acc_t2) / sqrt(2)) / ((tcho_mm_acc_t1 + tcho_mm_acc_t2) / 2) * 100, na.rm = TRUE )

summary(cv_acc_glu_mm)
summary(cv_acc_naa_mm)
summary(cv_acc_cho_mm)




# Si possible calcul un Z-score intra-sujet qui donne une idée sur la variabilité longitudinale
#  Define Function for Longitudinal Z-Score
# Formula: (My_Change - Average_Change) / SD_of_Changes
longitudinal_z_function <- function(t1, t2) {
  diffs <- t2 - t1
  mean_diff <- mean(diffs, na.rm = TRUE)
  sd_diff   <- sd(diffs, na.rm = TRUE)
  z_scores <- (diffs - mean_diff) / sd_diff
  return(z_scores)}


#  Calculate Z-Scores for Precuneus (mM)

z_prec_glu_mm <- longitudinal_z_function(glu_mm_precuneus_t1, glu_mm_precuneus_t2)
z_prec_naa_mm <- longitudinal_z_function(naa_mm_precuneus_t1, naa_mm_precuneus_t2)
z_prec_cho_mm <- longitudinal_z_function(tcho_mm_precuneus_t1, tcho_mm_precuneus_t2)

summary(z_prec_glu_mm)
summary(z_prec_naa_mm)
summary(z_prec_cho_mm)


#  Calculate Z-Scores for ACC (mM)

z_acc_glu_mm <- longitudinal_z_function(glu_mm_acc_t1, glu_mm_acc_t2)
z_acc_naa_mm <- longitudinal_z_function(naa_mm_acc_t1, naa_mm_acc_t2)
z_acc_cho_mm <- longitudinal_z_function(tcho_mm_acc_t1, tcho_mm_acc_t2)

summary(z_acc_glu_mm)
summary(z_acc_naa_mm)
summary(z_acc_cho_mm)



# 1. Create a data frame for plotting
# We combine both regions into one table
plot_data <- data.frame(
  Region = c(rep("Precuneus", length(z_prec_glu_mm)), 
             rep("ACC", length(z_acc_glu_mm))),
  Z_Score = c(z_prec_glu_mm, z_acc_glu_mm)
)

# 2. Plot
ggplot(plot_data, aes(x = Region, y = Z_Score, fill = Region)) +
  # Add the boxplot
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  # Add individual points (jittered slightly so they don't overlap)
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  # Add Z-score threshold lines (-2 and +2) for reference
  geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", color = "red") +
  geom_hline(yintercept = c(-3, 3), linetype = "dotted", color = "darkred") +
  # Formatting
  theme_minimal() +
  labs(title = "Longitudinal Variability of Glutamate (mM)",
       subtitle = "Z-Scores of Change (T2 - T1)",
       y = "Z-Score (Standard Deviations from Mean Change)",
       x = "Brain Region") +
  theme(legend.position = "none")
