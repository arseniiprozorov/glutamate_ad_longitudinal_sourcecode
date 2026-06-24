########## Figures ############
library(dplyr)
library(ggplot2)
library(emmeans)
library(dplyr)
library(ggplot2)
library(emmeans)


library(emmeans)
library(ggplot2)
library(dplyr)



library(ggplot2)

# 1. Convert the binary 0/1 decliner variable into a clean factor for the legend
MRS_prediction$decliner_factor <- factor(MRS_prediction$decliner_regression, 
                                         levels = c(0, 1), 
                                         labels = c("Stable", "Declined"))

# 2. Build the plot
ggplot(MRS_prediction, aes(x = m_m_precuneus, y = initiale_moca_score_total_30)) +
  
  # Map color to diagnosis (MCI/SCD+) and shape to longitudinal outcome
  geom_point(aes(color = diagnostic_nick, shape = decliner_factor), 
             alpha = 0.8, size = 3) +
  
  # The overarching quadratic biological curve (kept black to anchor the colored points)
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), 
              color = "black", fill = "gray80", linewidth = 1.2) +
  
  theme_minimal(base_size = 14) +
  
  # Custom poster-ready colors (colorblind-friendly) and shapes (16=circle, 17=triangle)
  scale_color_manual(values = c("MCI" = "#E69F00", "SCD+" = "#56B4E9")) +
  scale_shape_manual(values = c("Stable" = 16, "Declined" = 17)) + 
  
  labs(
    title = "Baseline Precuneus Glutamate vs. Initial MoCA",
    subtitle = "",
    x = "Precuneus Glutamate ",
    y = "Initial MoCA Score ",
    color = "Baseline Diagnosis",
    shape = "Longitudinal Outcome"
  ) +
  
  coord_cartesian(ylim = c(15, 30))


# ==============================================================================
# 1. MODIFIED MASTER VISUALIZATION ENGINE (With Ribbons & 24-30 scale)
# ==============================================================================
generate_trajectory_plot <- function(model_obj, var_name, title_text, legend_text) {
  
  # Build spec formula
  preds_formula <- as.formula(paste("~ years_from_baseline *", var_name))
  
  # Generate line coordinates (emmeans automatically includes CLs)
  at_grid <- list(years_from_baseline = seq(0, 6, by = 0.1))
  at_grid[[var_name]] <- c(-1, 0, 1)
  
  preds_df <- emmeans(model_obj, specs = preds_formula, at = at_grid) %>% as.data.frame()
  
  # Map grouping factor
  preds_df$Biomarker_Level <- factor(preds_df[[var_name]], 
                                     levels = c(-1, 0, 1), 
                                     labels = c("Low (-1 SD)", "Mean (0 SD)", "High (+1 SD)"))
  
  # Color/Fill Palette
  pal_colors <- c("Low (-1 SD)" = "#D55E00", 
                  "Mean (0 SD)" = "#737373", 
                  "High (+1 SD)" = "#3B5998")
  
  # Build ggplot
  p <- ggplot() +
    # Uniform gray background spaghetti data
    geom_line(data = MRS_prediction_long, 
              aes(x = years_from_baseline, y = moca, group = pscid), 
              color = "grey80", alpha = 0.3, linewidth = 0.4) +
    
    # ADDED: Error Shadow Ribbons (placed before lines so they appear behind)
    geom_ribbon(data = preds_df, 
                aes(x = years_from_baseline, ymin = lower.CL, ymax = upper.CL, fill = Biomarker_Level), 
                alpha = 0.15) +
    
    # Solid prediction paths
    geom_line(data = preds_df, 
              aes(x = years_from_baseline, y = emmean, group = Biomarker_Level, color = Biomarker_Level), 
              linewidth = 2) + 
    
    scale_color_manual(values = pal_colors, name = legend_text) +
    scale_fill_manual(values = pal_colors, name = legend_text) + # Required for the ribbon fill
    
    scale_x_continuous(limits = c(0, 6), breaks = c(0, 2, 4, 6)) +
    # FIXED SCALE 24-30
    coord_cartesian(ylim = c(24, 28)) + 
    labs(
      title = title_text,
      subtitle = "",
      x = "Years from Baseline",
      y = "MoCA Total Score"
    ) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), legend.position = "right")
  
  return(p)
}

# ==============================================================================
# 2. RUN THE ENGINE
# ==============================================================================
# All calls remain the same as your previous script
plot_precuneus <- generate_trajectory_plot(mixed_model_precuneus, "m_m_precuneus_z", "Precuneus Glutamate", "Baseline Precuneus Glu")
plot_acc <- generate_trajectory_plot(mixed_model_acc, "m_m_acc_z", "ACC Glutamate", "Baseline ACC Glutamate")
plot_ptau217 <- generate_trajectory_plot(mixed_model_ptau217, "plasma_ptau217_z", "Plasma p-Tau217", "Baseline Plasma p-Tau")
plot_thickness <- generate_trajectory_plot(mixed_model_thickness, "cortical_thickness_adsignature_dickson_z", "Cortical Thickness", "Baseline Thickness")
plot_hipp_vol <- generate_trajectory_plot(mixed_model_hipp_mean, "hipp_mean_z", "Hippocampal Volume", "Baseline Hipp Volume")
plot_hipp_act <- generate_trajectory_plot(mixed_model_hipp_mean_act, "hipp_mean_act_z", "Hippocampal Activation", "Baseline Hipp Activation")
plot_parietal_act <- generate_trajectory_plot(mixed_model_activation_parietal_l, "activation_parietal_sup_l_z", "Left Parietal Activation", "Baseline Parietal Act")

# ==============================================================================
# 3. DISPLAY PLOTS
# ==============================================================================
print(plot_ptau217)
print(plot_acc)
print(plot_precuneus)
print(plot_thickness)
print(plot_hipp_vol)
print(plot_hipp_act)
print(plot_parietal_act)















### ROC curves 

### ROC curves for Unimodal Significant Predictors

library(pROC)

# 1. Initialize the plot with the strongest individual predictor (ACC Glutamate)
plot(roc_glu_acc, 
     col = "#1f77b4", # Professional Blue
     lwd = 3, 
     main = "ROC Curves: Significant Predictors of Cognitive Decline",
     xlab = "1 - Specificity", 
     ylab = "Sensitivity",
     legacy.axes = TRUE) # Ensures standard 0 to 1 (1-Specificity) axis mapping

# 2. Add the remaining two significant unadjusted models
plot(roc_glu_prec, add = TRUE, col = "#2ca02c", lwd = 3)     # Professional Green
plot(roc_struc_thick, add = TRUE, col = "#ff7f0e", lwd = 3)  # Professional Orange

# 3. Add the diagonal identity/reference line (chance level)
abline(a = 0, b = 1, lty = 2, col = "darkgrey")

# 4. Add a clean legend matching your exact AUC values
legend("bottomright", 
       legend = c(
         paste0("ACC Glutamate (AUC = ", round(auc(roc_glu_acc), 3), ")"),
         paste0("Precuneus Glutamate (AUC = ", round(auc(roc_glu_prec), 3), ")"),
         paste0("Cortical Thickness (AUC = ", round(auc(roc_struc_thick), 3), ")")
       ),
       col = c("#1f77b4", "#2ca02c", "#ff7f0e"), 
       lwd = 3, 
       cex = 0.9, 
       bty = "n") # Removes the bounding border box






## Combined models 2 variables 

library(pROC)

# OPTIONAL: Run this line if you want all 3 plots side-by-side in one layout window
# par(mfrow = c(1, 3))


### =========================================================================
### PLOT 1: Fluid Pathology & Regional Metabolic Interactions
### =========================================================================

# Initialize with the top-performing hybrid model (p-Tau217 + ACC Glu)
plot(roc_ptau217_acc, 
     col = "#1f77b4", lwd = 3, legacy.axes = TRUE,
     main = "Pathology & Metabolism Models",
     xlab = "1 - Specificity", ylab = "Sensitivity")

# Add secondary combinations
plot(roc_ptau217_prec, add = TRUE, col = "#aec7e8", lwd = 3) # Lighter Blue
plot(roc_prec_acc, add = TRUE, col = "#9467bd", lwd = 3)     # Purple

# Baseline chance line
abline(a = 0, b = 1, lty = 2, col = "darkgrey")

# Legend
legend("bottomright", 
       legend = c(
         paste0("p-Tau217 + ACC Glu (AUC = ", round(auc(roc_ptau217_acc), 3), ")"),
         paste0("p-Tau217 + Precuneus Glu (AUC = ", round(auc(roc_ptau217_prec), 3), ")"),
         paste0("Precuneus Glu + ACC Glu (AUC = ", round(auc(roc_prec_acc), 3), ")")
       ),
       col = c("#1f77b4", "#aec7e8", "#9467bd"), lwd = 3, cex = 0.8, bty = "n")


### =========================================================================
### PLOT 2: Metabolic-Structural Models (Glutamate + Atrophy)
### =========================================================================

# Initialize with the strongest structural combo (Precuneus Glu + Cortical Thickness)
plot(roc_prec_thick, 
     col = "#2ca02c", lwd = 3, legacy.axes = TRUE,
     main = "Metabolic-Structural Models",
     xlab = "1 - Specificity", ylab = "Sensitivity")

# Add other structural combinations
plot(roc_prec_hip_vol, add = TRUE, col = "#98df8a", lwd = 3) # Lighter Green
plot(roc_acc_thick, add = TRUE, col = "#ff7f0e", lwd = 3)    # Orange
plot(roc_acc_hip_vol, add = TRUE, col = "#ffbb78", lwd = 3)  # Lighter Orange

# Baseline chance line
abline(a = 0, b = 1, lty = 2, col = "darkgrey")

# Legend
legend("bottomright", 
       legend = c(
         paste0("Precuneus Glu + Thick (AUC = ", round(auc(roc_prec_thick), 3), ")"),
         paste0("Precuneus Glu + Hip Vol (AUC = ", round(auc(roc_prec_hip_vol), 3), ")"),
         paste0("ACC Glu + Thick (AUC = ", round(auc(roc_acc_thick), 3), ")"),
         paste0("ACC Glu + Hip Vol (AUC = ", round(auc(roc_acc_hip_vol), 3), ")")
       ),
       col = c("#2ca02c", "#98df8a", "#ff7f0e", "#ffbb78"), lwd = 3, cex = 0.8, bty = "n")


### =========================================================================
### PLOT 3: Metabolic-Functional Models (Glutamate + fMRI Activation)
### =========================================================================

# Initialize with the strongest functional combo (Precuneus Glu + Parietal Activation)
plot(roc_prec_par_act, 
     col = "#d62728", lwd = 3, legacy.axes = TRUE,
     main = "Metabolic-Functional Models",
     xlab = "1 - Specificity", ylab = "Sensitivity")

# Add remaining functional combinations
plot(roc_acc_par_act, add = TRUE, col = "#ff9896", lwd = 3)  # Lighter Red/Pink
plot(roc_prec_hip_act, add = TRUE, col = "#e377c2", lwd = 3) # Magenta/Pink
plot(roc_acc_hip_act, add = TRUE, col = "#f7b6d2", lwd = 3)  # Lighter Pink

# Baseline chance line
abline(a = 0, b = 1, lty = 2, col = "darkgrey")

# Legend
legend("bottomright", 
       legend = c(
         paste0("Precuneus Glu + Parietal Act (AUC = ", round(auc(roc_prec_par_act), 3), ")"),
         paste0("ACC Glu + Parietal Act (AUC = ", round(auc(roc_acc_par_act), 3), ")"),
         paste0("Precuneus Glu + Hip Act (AUC = ", round(auc(roc_prec_hip_act), 3), ")"),
         paste0("ACC Glu + Hip Act (AUC = ", round(auc(roc_acc_hip_act), 3), ")")
       ),
       col = c("#d62728", "#ff9896", "#e377c2", "#f7b6d2"), lwd = 3, cex = 0.75, bty = "n")


# Reset layout window to default single view when done
# par(mfrow = c(1, 1))








names(MRS_prediction)

library(ggplot2)
library(dplyr)

raw_mm_column <- "m_m_acc" 

# 1. Prepare your data
plot_data <- MRS_prediction %>%
  mutate(
    Clinical_Status = factor(decliner_regression, levels = c(0, 1), labels = c("Stable", "Decliner"))
  )

# =========================================================================
# GRAPH 1: Logistic Regression Probability Curve (mM version)
# =========================================================================
ggplot(plot_data, aes(x = .data[[raw_mm_column]], y = decliner_regression)) +
  geom_point(aes(color = Clinical_Status), 
             position = position_jitter(height = 0.03, width = 0), 
             size = 2.5, alpha = 0.7) +
  stat_smooth(method = "glm", method.args = list(family = "binomial"), 
              se = TRUE, color = "#2C3E50", fill = "gray80", alpha = 0.4) +
  scale_color_manual(values = c("Stable" = "#56B4E9", "Decliner" = "#FDB863")) +
  labs(
    title = "Predicted Probability of Cognitive Decline by ACC Glutamate",
    subtitle = "Absolute metabolic concentrations drive progressive longitudinal risk",
    x = "Baseline ACC Glutamate Concentration (mM)",
    y = "Predicted Probability of Longitudinal Regression",
    color = "Clinical Outcome"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, color = "gray40", hjust = 0.5),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 20, unit = "pt") # Prevents clipping
  )

# =========================================================================
# GRAPH 2: Clinical Separation Boxplot (mM version)
# =========================================================================
ggplot(plot_data, aes(x = Clinical_Status, y = .data[[raw_mm_column]], fill = Clinical_Status)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.4, color = "gray30") +
  geom_jitter(aes(color = Clinical_Status), width = 0.15, size = 2.5, alpha = 0.8) +
  scale_fill_manual(values = c("Stable" = "#56B4E9", "Decliner" = "#FDB863")) +
  scale_color_manual(values = c("Stable" = "#56B4E9", "Decliner" = "#FDB863")) +
  labs(
    title = "Baseline ACC Glutamate Separation Profiles",
    subtitle = "High Specificity (96.7%): Intact metabolic signatures strictly rule out decline",
    x = "Longitudinal Clinical Outcome (8-Year Window)",
    y = "Baseline ACC Glutamate Concentration (mM)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, color = "gray40", hjust = 0.5)
  )





# 1. Define the exact column name as a string
raw_tau_column <- "plasma_ptau217"

# =========================================================================
# GRAPH 1: Logistic Regression Probability Curve (p-Tau217 version)
# =========================================================================
ggplot(plot_data, aes(x = .data[[raw_tau_column]], y = decliner_regression)) +
  geom_point(aes(color = Clinical_Status), 
             position = position_jitter(height = 0.03, width = 0), 
             size = 2.5, alpha = 0.7) +
  stat_smooth(method = "glm", method.args = list(family = "binomial"), 
              se = TRUE, color = "#2C3E50", fill = "gray80", alpha = 0.4) +
  scale_color_manual(values = c("Stable" = "#56B4E9", "Decliner" = "#FDB863")) +
  labs(
    title = "Predicted Probability of Cognitive Decline by Plasma p-Tau217",
    subtitle = "Positive relationship: Systemic pathology trends with progressive risk",
    x = "Baseline Plasma p-Tau217 Concentration (pg/mL)",
    y = "Predicted Probability of Longitudinal Regression",
    color = "Clinical Outcome"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, color = "gray40", hjust = 0.5),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 20, unit = "pt")
  )

# =========================================================================
# GRAPH 2: Clinical Separation Boxplot (p-Tau217 version)
# =========================================================================
ggplot(plot_data, aes(x = Clinical_Status, y = .data[[raw_tau_column]], fill = Clinical_Status)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.4, color = "gray30") +
  geom_jitter(aes(color = Clinical_Status), width = 0.15, size = 2.5, alpha = 0.8) +
  scale_fill_manual(values = c("Stable" = "#56B4E9", "Decliner" = "#FDB863")) +
  scale_color_manual(values = c("Stable" = "#56B4E9", "Decliner" = "#FDB863")) +
  labs(
    title = "Baseline Plasma p-Tau217 Separation Profiles",
    subtitle = "High Specificity (92.9%): Elevated systemic pathology strictly rules in decline",
    x = "Longitudinal Clinical Outcome (8-Year Window)",
    y = "Baseline Plasma p-Tau217 Concentration (pg/mL)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, color = "gray40", hjust = 0.5)
  )








library(survival)
library(survminer)
library(dplyr)
library(survival)
library(survminer)
library(dplyr)

# 1. Categorize continuous ACC glutamate into tertiles inside your cox_data
cox_data <- cox_data %>%
  mutate(
    acc_tertile = ntile(m_m_acc_z, 3),
    acc_strata = case_when(
      acc_tertile == 1 ~ "Low ACC Glu",
      acc_tertile == 2 ~ "Medium ACC Glu",
      acc_tertile == 3 ~ "High ACC Glu"
    ),
    acc_strata = factor(acc_strata, levels = c("Low ACC Glu", "Medium ACC Glu", "High ACC Glu"))
  )

# 2. Fit the survival curve matching your Cox model variables
fit_acc <- survfit(Surv(followup_years, decliner_regression) ~ acc_strata, data = cox_data)

# 3. Plot the ACC graph
ggsurvplot(
  fit_acc,
  data = cox_data,
  conf.int = TRUE,          
  conf.int.alpha = 0.15,      # Slightly lower alpha prevents overlapping cloud confusion
  censor = TRUE,            
  censor.shape = 43,        
  censor.size = 4.5,
  palette = c("#FDB863", "#999999", "#56B4E9"), 
  title = "Cognitive Stability by ACC Glutamate Levels",
  xlab = "Years from Baseline",
  ylab = "Probability of Stability",
  legend.title = "Metabolic Strata",
  legend = "top",
  ggtheme = theme_minimal() 
)










# 1. Categorize continuous Precuneus glutamate into tertiles inside your cox_data
cox_data <- cox_data %>%
  mutate(
    prec_tertile = ntile(m_m_precuneus_z, 3),
    prec_strata = case_when(
      prec_tertile == 1 ~ "Low Precuneus Glu",
      prec_tertile == 2 ~ "Medium Precuneus Glu",
      prec_tertile == 3 ~ "High Precuneus Glu"
    ),
    prec_strata = factor(prec_strata, levels = c("Low Precuneus Glu", "Medium Precuneus Glu", "High Precuneus Glu"))
  )

# 2. Fit the survival curve matching your Cox model variables
fit_prec <- survfit(Surv(followup_years, decliner_regression) ~ prec_strata, data = cox_data)

# 3. Plot the Precuneus graph
ggsurvplot(
  fit_prec,
  data = cox_data,
  conf.int = TRUE,          
  conf.int.alpha = 0.15,     
  censor = TRUE,            
  censor.shape = 43,        
  censor.size = 4.5,
  palette = c("#FDB863", "#999999", "#56B4E9"), 
  title = "Cognitive Stability by Precuneus Glutamate Levels",
  xlab = "Years from Baseline",
  ylab = "Probability of Stability",
  legend.title = "Metabolic Strata",
  legend = "top",
  ggtheme = theme_minimal() 
)