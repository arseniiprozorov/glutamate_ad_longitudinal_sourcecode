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




# 1. Convert the binary 0/1 decliner variable into a clean factor for the legend
MRS_prediction$decliner_factor <- factor(MRS_prediction$decliner_regression, 
                                         levels = c(0, 1), 
                                         labels = c("Stable", "Declined"))

# 2. Build the plot
ggplot(MRS_prediction, aes(x = m_m_precuneus, y = initiale_moca_score_total_30)) +
  
  # The overarching global quadratic curve (kept as a baseline reference)
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), 
              color = "black", fill = "gray80", linewidth = 1.2, alpha = 0.5) +
  
  # NEW: Subgroup LINEAR curves mapping color to diagnosis and line style to decline
  # Changed formula to y ~ x to force straight, stable trendlines
  stat_smooth(aes(color = diagnostic_nick, linetype = decliner_factor),
              method = "lm", formula = y ~ x, 
              se = FALSE, linewidth = 1) +
  
  # Map color to diagnosis and shape to longitudinal outcome
  geom_point(aes(color = diagnostic_nick, shape = decliner_factor), 
             alpha = 0.8, size = 3) +
  
  theme_minimal(base_size = 14) +
  
  # FIXED: Added the third group ("HC") to explicitly color and label them
  scale_color_manual(values = c("MCI" = "#E69F00", "SCD+" = "#56B4E9", "HC" = "gray50")) +
  scale_shape_manual(values = c("Stable" = 16, "Declined" = 17)) + 
  
  labs(
    title = "Baseline Precuneus Glutamate vs. Initial MoCA",
    subtitle = "Subgroup Linear Trends vs. Global Quadratic Trend",
    x = "Precuneus Glutamate ",
    y = "Initial MoCA Score ",
    color = "Baseline Diagnosis",
    shape = "Outcome",
    linetype = "Outcome"
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


formula(mixed_model_hipp_mean_act)
PATH_FIG.precuneus = "Slopes_precuneus.tiff"
if(file.exists(PATH_FIG.precuneus)){file.remove(PATH_FIG.precuneus)}

##########
graphics.off()

plot_precuneus <- interactions::interact_plot(
  mixed_model_precuneus,               # <--- Your lmer model
  pred        = years_from_baseline,   # <--- Time on the X-axis
  modx        = m_m_precuneus_z,       # <--- Glutamate creates the different lines
  modx.values = c(-1, 0, 1),           
  modx.labels = c("Low Glu (-1 SD)", "Average Glu (0 SD)", "High Glu (+1 SD)"),
  legend.main = "Precuneus Glutamate",
  plot.points = TRUE,                  # Plots the raw longitudinal data points
  
  # Engine Visual Style
  interval    = TRUE,                
  int.alpha   = 0.15,                
  line.thickness = 2,                
  point.alpha = 0.3,                 
  colors      = c("#D55E00", "#737373", "#3B5998") 
)

plot_precuneus <- plot_precuneus + 
  theme_minimal(base_size = 12) +
  labs(
    x = "Years from Baseline", 
    y = "MoCA Total Score"
  ) +
  theme(panel.grid.minor = element_blank(), legend.position = "right")

############
ggsave(PATH_FIG.precuneus, plot = plot_precuneus, width = 13, height = 7, units = "cm", dpi=300)
print(plot_precuneus)










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







############################### Survival plots #######################


library(survival)
library(survminer)

# 1. Lock in the factor levels so Low is plotted first and High is plotted last
MRS_prediction$acc_tertiales <- factor(MRS_prediction$acc_tertiales, 
                                       levels = c("Low ACC Glu", "Medium ACC Glu", "High ACC Glu"))

MRS_prediction$precuneus_tertiales <- factor(MRS_prediction$precuneus_tertiales, 
                                             levels = c("Low Precuneus Glu", "Medium Precuneus Glu", "High Precuneus Glu"))


# ==============================================================================
# ACC SURVIVAL PLOT
# ==============================================================================

# Fit the survival curve using the calculated tertiles
fit_acc <- survfit(Surv(max_years_from_baseline, decliner_regression) ~ acc_tertiales, data = MRS_prediction)

# Plot the ACC graph
ggsurvplot(
  fit_acc,
  data = MRS_prediction,
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


# ==============================================================================
# PRECUNEUS SURVIVAL PLOT
# ==============================================================================

# Fit the survival curve using the calculated tertiles
fit_prec <- survfit(Surv(max_years_from_baseline, decliner_regression) ~ precuneus_tertiales, data = MRS_prediction)

# Plot the Precuneus graph
ggsurvplot(
  fit_prec,
  data = MRS_prediction,
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
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  install.packages("ragg")
  
  
  
  ########################### Article figures ####################
  ########################### Article Figure 1 ####################
  library(ggplot2)
  library(emmeans)
  library(dplyr)
  library(patchwork)
  library(ragg)
  
  # ==============================================================================
  # 1. Helper Function: Plot Adjusted Trajectories for Tertiles
  # ==============================================================================
  
  plot_tertile_traj <- function(model, df, tert_var, plot_title, y_limits = c(24, 28)) {
    
    # a) Generate model-adjusted marginal predictions (Year 0 to 6)
    pred_grid <- emmip(
      model, 
      as.formula(paste(tert_var, "~ years_from_baseline")), 
      at = list(years_from_baseline = seq(0, 6, by = 0.5)),
      plotit = FALSE, 
      CIs = TRUE
    )
    
    # Ensure clean factor levels
    pred_grid$Tertile <- factor(pred_grid[[tert_var]], levels = c("Low", "Medium", "High"))
    
    # b) Prepare raw data for background spaghetti plots
    valid_rows <- !is.na(df[[tert_var]]) & !is.na(df$moca)
    df_plot <- df[valid_rows, ]
    df_plot$Tertile <- factor(df_plot[[tert_var]], levels = c("Low", "Medium", "High"))
    
    # c) High-Contrast Palette (Red = At Risk, Blue = Medium, Green = Protective)
    tert_colors <- c("Low" = "#E64B35", "Medium" = "#4DBBD5", "High" = "#00A087")
    
    # d) Build ggplot
    p <- ggplot() +
      geom_line(
        data = df_plot,
        aes(x = years_from_baseline, y = moca, group = pscid, color = Tertile),
        alpha = 0.15, linewidth = 0.4
      ) +
      geom_ribbon(
        data = pred_grid,
        aes(x = years_from_baseline, ymin = LCL, ymax = UCL, fill = Tertile),
        alpha = 0.20, color = NA
      ) +
      geom_line(
        data = pred_grid,
        aes(x = years_from_baseline, y = yvar, color = Tertile),
        linewidth = 1.2
      ) +
      scale_color_manual(values = tert_colors) +
      scale_fill_manual(values = tert_colors) +
      scale_x_continuous(breaks = seq(0, 6, 1), limits = c(0, 6)) +
      scale_y_continuous(breaks = seq(min(y_limits), max(y_limits), by = 1)) +
      coord_cartesian(ylim = y_limits) + 
      labs(
        title = plot_title,
        x = "Years from Baseline",
        y = "MoCA Score",
        color = "Baseline Tertile",
        fill  = "Baseline Tertile"
      ) +
      theme_classic(base_size = 12, base_family = "Arial") +
      theme(
        plot.title       = element_text(face = "bold", size = 13),
        axis.title       = element_text(face = "bold"),
        legend.position  = "right",
        legend.title     = element_text(face = "bold", size = 11),
        legend.text      = element_text(size = 10),
        plot.margin      = margin(10, 10, 10, 10)
      )
    
    return(p)
  }
  
  # ==============================================================================
  # 2. Build Individual Panels in Requested Order:
  # Order: Glutamate (A, B) -> fMRI (C) -> sMRI (D) -> Plasma p-tau217 (E)
  # ==============================================================================
  p_acc   <- plot_tertile_traj(mod_acc_tert,   MRS_prediction_long, "acc_tert",       "A. ACC Glutamate")
  p_prec  <- plot_tertile_traj(mod_prec_tert,  MRS_prediction_long, "precuneus_tert", "B. Precuneus Glutamate")
  p_hip   <- plot_tertile_traj(mod_hip_tert,   MRS_prediction_long, "hipp_act_tert",  "C. Hippocampal Activation")
  p_thick <- plot_tertile_traj(mod_thick_tert, MRS_prediction_long, "thickness_tert", "D. Cortical Thickness")
  p_ptau  <- plot_tertile_traj(mod_ptau_tert,  MRS_prediction_long, "ptau217_tert",   "E. Plasma p-Tau217")
  
  # ==============================================================================
  # 3. Patchwork Arrangement & Save
  # ==============================================================================
  combo_traj <- (p_acc | p_prec) / 
    (p_hip | p_thick) / 
    (p_ptau | plot_spacer()) + 
    plot_layout(guides = "collect") & 
    theme(
      legend.position = "right",
      legend.box.margin = margin(10, 8, 0, 0)
    )
  
  # Set target directory
  out_dir <- "C:/Users/okkam/Desktop/labo/article 2/A&D"
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Base R tiff export (1200 DPI, LZW compression)
  tiff(
    filename    = file.path(out_dir, "fig1_longitudinal_tertile_trajectories_1200dpi.tiff"),
    width       = 10.0, 
    height      = 11.5, 
    units       = "in", 
    res         = 1200, 
    compression = "lzw"
  )
  print(combo_traj)
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  library(ggplot2)
  library(emmeans)
  library(patchwork)
  
  # ==============================================================================
  # 1. Helper Function: Plot Adjusted Trajectories for Tertiles
  # ==============================================================================
  
  plot_tertile_traj <- function(model, df, tert_var, plot_title, y_limits = c(23, 29)) {
    
    # a) Generate model-adjusted marginal predictions (Year 0 to 6)
    pred_grid <- emmip(
      model, 
      as.formula(paste(tert_var, "~ years_from_baseline")), 
      at = list(years_from_baseline = seq(0, 6, by = 0.5)),
      plotit = FALSE, 
      CIs = TRUE
    )
    
    # Ensure clean factor levels
    pred_grid$Tertile <- factor(pred_grid[[tert_var]], levels = c("Low", "Medium", "High"))
    
    # b) Prepare raw data for background spaghetti plots
    valid_rows <- !is.na(df[[tert_var]]) & !is.na(df$moca)
    df_plot <- df[valid_rows, ]
    df_plot$Tertile <- factor(df_plot[[tert_var]], levels = c("Low", "Medium", "High"))
    
    # c) High-Contrast Palette (Red = At Risk, Blue = Medium, Green = Protective)
    tert_colors <- c("Low" = "#E64B35", "Medium" = "#4DBBD5", "High" = "#00A087")
    
    # d) Build ggplot
    p <- ggplot() +
      geom_line(
        data = df_plot,
        aes(x = years_from_baseline, y = moca, group = pscid, color = Tertile),
        alpha = 0.15, linewidth = 0.4
      ) +
      geom_ribbon(
        data = pred_grid,
        aes(x = years_from_baseline, ymin = LCL, ymax = UCL, fill = Tertile),
        alpha = 0.20, color = NA
      ) +
      geom_line(
        data = pred_grid,
        aes(x = years_from_baseline, y = yvar, color = Tertile),
        linewidth = 1.2
      ) +
      scale_color_manual(values = tert_colors) +
      scale_fill_manual(values = tert_colors) +
      scale_x_continuous(breaks = seq(0, 6, 1), limits = c(0, 6)) +
      scale_y_continuous(breaks = seq(min(y_limits), max(y_limits), by = 1)) +
      coord_cartesian(ylim = y_limits) + 
      labs(
        title = plot_title,
        x = "Years from Baseline",
        y = "MoCA Score",
        color = "Baseline Tertile",
        fill  = "Baseline Tertile"
      ) +
      theme_classic(base_size = 12, base_family = "Arial") +
      theme(
        plot.title       = element_text(face = "bold", size = 13),
        axis.title       = element_text(face = "bold"),
        legend.position  = "right",
        legend.title     = element_text(face = "bold", size = 11),
        legend.text      = element_text(size = 10),
        plot.margin      = margin(10, 10, 10, 10)
      )
    
    return(p)
  }
  
  # ==============================================================================
  # 2. Build Individual Panels (Titles Matched to Grid Order)
  # ==============================================================================
  
  p_ptau  <- plot_tertile_traj(mod_ptau_tert,  MRS_prediction_long, "ptau217_tert",   "A. Plasma p-Tau217")
  p_acc   <- plot_tertile_traj(mod_acc_tert,   MRS_prediction_long, "acc_tert",       "B. ACC Glutamate")
  p_prec  <- plot_tertile_traj(mod_prec_tert,  MRS_prediction_long, "precuneus_tert", "C. Precuneus Glutamate")
  p_thick <- plot_tertile_traj(mod_thick_tert, MRS_prediction_long, "thickness_tert", "D. Cortical Thickness")
  p_hip   <- plot_tertile_traj(mod_hip_tert,   MRS_prediction_long, "hipp_act_tert",  "E. Hippocampal Activation")
  
  # ==============================================================================
  # 3. Patchwork Arrangement & Save
  # Layout:
  # [ A. p-Tau217 ]   [ B. ACC Glu    ]
  # [ C. Prec Glu ]   [ D. Thickness  ]
  # [ E. Hipp Act ]   [  (spacer)     ]
  # ==============================================================================
  
  combo_traj <- (p_ptau | p_acc) / 
    (p_prec | p_thick) / 
    (p_hip  | plot_spacer()) + 
    plot_layout(guides = "collect") & 
    theme(
      legend.position = "right",
      legend.box.margin = margin(10, 8, 0, 0)
    )
  
  # Set target directory
  out_dir <- "C:/Users/okkam/Desktop/labo/article 2/A&D"
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Base R tiff export (1200 DPI, LZW compression)
  tiff(
    filename    = file.path(out_dir, "fig1_longitudinal_tertile_trajectories_1200dpi.tiff"),
    width       = 10.0, 
    height      = 11.5, 
    units       = "in",
    res         = 1200, 
    compression = "lzw"
  )
  print(combo_traj)
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  
  #### Figure 2 ######
  #### Figure 2 ######
  library(pROC)
  library(ggplot2)
  library(patchwork)
  
  # ==============================================================================
  # 0. Helper: Extract and Format ROC Dataframe with Explicit (0,0) and (1,1)
  # ==============================================================================
  get_roc_df <- function(roc_obj, model_name) {
    df <- data.frame(
      fpr = 1 - roc_obj$specificities,
      tpr = roc_obj$sensitivities,
      model = model_name
    )
    # Ensure strict (0,0) and (1,1) boundaries
    df <- rbind(data.frame(fpr = 0, tpr = 0, model = model_name),
                df,
                data.frame(fpr = 1, tpr = 1, model = model_name))
    # Sort strictly by FPR then TPR for clean geom_path drawing
    df <- df[order(df$fpr, df$tpr), ]
    rownames(df) <- NULL
    return(df)
  }
  
  # 1. Build tidy dataframes for Panel A (Unimodal)
  df_acc   <- get_roc_df(roc_glu_acc, "ACC Glutamate")
  df_prec  <- get_roc_df(roc_glu_prec, "Precuneus Glutamate")
  df_thick <- get_roc_df(roc_struc_thick, "Cortical Thickness")
  df_ptau  <- get_roc_df(roc_ptau217, "Plasma p-tau217")
  df_hip   <- get_roc_df(roc_func_hip, "Hippocampal Act")
  
  df_unimodal <- rbind(df_acc, df_prec, df_thick, df_ptau, df_hip)
  df_unimodal$model <- factor(df_unimodal$model, levels = c(
    "ACC Glutamate", "Precuneus Glutamate", "Cortical Thickness", "Plasma p-tau217", "Hippocampal Act"
  ))
  
  # 2. Build tidy dataframe for Panel B (Multimodal)
  df_multi <- get_roc_df(roc_model_nocov_step_sig, "Multimodal")
  
  # ==============================================================================
  # 1. Panel A: Unimodal Plot
  # ==============================================================================
  
  unimodal_colors <- c(
    "ACC Glutamate"       = "#00A087", # Protective Green/Teal
    "Precuneus Glutamate" = "#4DBBD5", # Medium Blue/Cyan
    "Cortical Thickness"  = "#3C5488", # Deep Navy
    "Plasma p-tau217"     = "#F39B7F", # Coral Accent
    "Hippocampal Act"     = "#7E6148"  # Muted Slate
  )
  
  youden_uni <- data.frame(
    model = c("ACC Glutamate", "Precuneus Glutamate", "Cortical Thickness", "Plasma p-tau217", "Hippocampal Act"),
    fpr   = c(1 - 0.967, 1 - 0.787, 1 - 0.790, 1 - 0.929, 1 - 0.418),
    tpr   = c(0.333,     0.591,     0.591,     0.333,     0.810)
  )
  youden_uni$model <- factor(youden_uni$model, levels = levels(df_unimodal$model))
  
  p_roc_uni <- ggplot() +
    # Chance diagonal
    geom_segment(
      aes(x = 0, xend = 1, y = 0, yend = 1),
      color = "grey65", linetype = "dashed", linewidth = 0.7
    ) +
    # ROC Curves via geom_path
    geom_path(
      data = df_unimodal,
      aes(x = fpr, y = tpr, color = model),
      linewidth = 1.0
    ) +
    # Optimal Operating Points (Youden)
    geom_point(
      data = youden_uni,
      aes(x = fpr, y = tpr, color = model),
      size = 2.8, shape = 19, show.legend = FALSE
    ) +
    scale_color_manual(values = unimodal_colors) +
    scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
    scale_y_continuous(expand = c(0.01, 0.01), breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
    labs(
      title = "A. Unimodal Models",
      x = "1 - Specificity",
      y = "Sensitivity",
      color = NULL
    ) +
    theme_classic(base_size = 12, base_family = "Arial") +
    theme(
      plot.title        = element_text(face = "bold", size = 13),
      axis.title        = element_text(face = "bold", size = 11),
      axis.text         = element_text(color = "black", size = 10),
      legend.position   = c(0.62, 0.22),
      legend.text       = element_text(size = 9, face = "bold"),
      legend.background = element_rect(fill = alpha("white", 0.85), color = "grey85", linewidth = 0.5),
      legend.key.height = unit(0.42, "cm"),
      plot.margin       = margin(10, 15, 10, 10)
    )
  
  # ==============================================================================
  # 2. Panel B: Multimodal Plot
  # ==============================================================================
  
  youden_multi <- data.frame(
    fpr = 1 - 0.900,
    tpr = 0.643
  )
  
  p_roc_multi <- ggplot() +
    # Chance diagonal
    geom_segment(
      aes(x = 0, xend = 1, y = 0, yend = 1),
      color = "grey65", linetype = "dashed", linewidth = 0.7
    ) +
    # Multimodal Curve
    geom_path(
      data = df_multi,
      aes(x = fpr, y = tpr),
      color = "#E64B35", linewidth = 1.3
    ) +
    # Youden Point
    geom_point(
      data = youden_multi,
      aes(x = fpr, y = tpr),
      color = "#E64B35", size = 3.4, shape = 19
    ) +
    scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
    scale_y_continuous(expand = c(0.01, 0.01), breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
    labs(
      title = "B. Multimodal Model",
      x = "1 - Specificity",
      y = "Sensitivity"
    ) +
    theme_classic(base_size = 12, base_family = "Arial") +
    theme(
      plot.title      = element_text(face = "bold", size = 13),
      axis.title      = element_text(face = "bold", size = 11),
      axis.text       = element_text(color = "black", size = 10),
      legend.position = "none",
      plot.margin     = margin(10, 10, 10, 15)
    )
  
  # ==============================================================================
  # 3. Patchwork Arrangement & TIFF Export
  # ==============================================================================
  
  combo_roc <- p_roc_uni | p_roc_multi
  
  out_dir <- "C:/Users/okkam/Desktop/labo/article 2/A&D"
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  tiff(
    filename    = file.path(out_dir, "fig2_roc_curves_unimodal_multimodal_1200dpi.tiff"),
    width       = 10.5,
    height      = 5.5,
    units       = "in",
    res         = 1200,
    compression = "lzw"
  )
  print(combo_roc)
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  ############ Figure 3 #########
  library(survival)
  library(survminer)
  library(ggplot2)
  library(patchwork)
  
  # ==============================================================================
  # 1. Survival Fits & Explicit Factor/Color Alignment
  # Level 1 = Red (#E64B35, High Risk / Low Glu)
  # Level 2 = Cyan (#4DBBD5, Intermediate)
  # Level 3 = Green (#00A087, Low Risk / High Glu)
  # ==============================================================================
  shared_labels <- c(
    "Low Glutamate / High Multimodal Risk",
    "Intermediate",
    "High Glutamate / Low Multimodal Risk"
  )
  
  # Unimodal factor coding
  MRS_prediction$acc_tert_clean <- factor(
    MRS_prediction$acc_tertiales,
    labels = shared_labels
  )
  
  MRS_prediction$prec_tert_clean <- factor(
    MRS_prediction$precuneus_tertiales,
    labels = shared_labels
  )
  
  # Multimodal Risk Tertiles
  MRS_prediction$multimodal_risk <- predict(cox_nocov, newdata = MRS_prediction, type = "lp")
  risk_tertiles <- quantile(MRS_prediction$multimodal_risk, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  
  # Cut and align multimodal factor levels to match Red -> Blue -> Green
  MRS_prediction$multimodal_risk_clean <- cut(
    MRS_prediction$multimodal_risk,
    breaks = risk_tertiles,
    labels = c(shared_labels[3], shared_labels[2], shared_labels[1]),
    include.lowest = TRUE
  )
  
  MRS_prediction$multimodal_risk_clean <- factor(
    MRS_prediction$multimodal_risk_clean,
    levels = shared_labels
  )
  
  # Fit survival models
  fit_acc_km   <- survfit(Surv(max_years_from_baseline, decliner_regression) ~ acc_tert_clean, data = MRS_prediction)
  fit_prec_km  <- survfit(Surv(max_years_from_baseline, decliner_regression) ~ prec_tert_clean, data = MRS_prediction)
  fit_multi_km <- survfit(Surv(max_years_from_baseline, decliner_regression) ~ multimodal_risk_clean, data = MRS_prediction)
  
  # ==============================================================================
  # 2. General KM Panel Function
  # ==============================================================================
  plot_km_panel <- function(fit, df, title_text) {
    p <- ggsurvplot(
      fit,
      data = df,
      size = 1.1,
      palette = c("#E64B35", "#4DBBD5", "#00A087"), # Red = High Risk, Cyan = Intermediate, Green = Low Risk
      censor.size = 2.2,
      censor.shape = 3,
      pval = FALSE,
      conf.int = FALSE,
      xlab = "Follow-up Duration (Years)",
      ylab = "Cognitive Stability Probability",
      title = title_text,
      legend.title = "Stratum / Prognostic Profile",
      legend.labs = shared_labels,
      ggtheme = theme_classic(base_size = 12, base_family = "Arial") +
        theme(
          plot.title   = element_text(face = "bold", size = 13),
          axis.title   = element_text(face = "bold", size = 11),
          axis.text    = element_text(color = "black", size = 10),
          legend.title = element_text(face = "bold", size = 10),
          legend.text  = element_text(size = 9.5),
          plot.margin  = margin(10, 10, 10, 10)
        )
    )
    return(p$plot)
  }
  
  # ==============================================================================
  # 3. Build Panels
  # ==============================================================================
  p_acc  <- plot_km_panel(fit_acc_km,   MRS_prediction, "A. ACC Glutamate")
  p_prec <- plot_km_panel(fit_prec_km,  MRS_prediction, "B. Precuneus Glutamate")
  p_mult <- plot_km_panel(fit_multi_km, MRS_prediction, "C. Multimodal Risk Profile")
  
  # ==============================================================================
  # 4. Patchwork Layout with Shared Bottom Legend & TIFF Export
  # ==============================================================================
  combo_km <- (p_acc | p_prec | p_mult) +
    plot_layout(guides = "collect") &
    theme(
      legend.position   = "bottom",
      legend.box.margin = margin(10, 0, 0, 0)
    )
  
  # Display in Plots viewer
  print(combo_km)
  
  # Save TIFF (1200 DPI, LZW Compression)
  out_dir <- "C:/Users/okkam/Desktop/labo/article 2/A&D"
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  tiff(
    filename    = file.path(out_dir, "fig3_kaplan_meier_survival_curves_1200dpi.tiff"),
    width       = 14.5,
    height      = 5.8,
    units       = "in",
    res         = 1200,
    compression = "lzw"
  )
  print(combo_km)
  dev.off()