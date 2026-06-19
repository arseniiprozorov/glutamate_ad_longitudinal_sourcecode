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
      subtitle = "Trajectories with 95% Confidence Intervals",
      x = "Years from Baseline",
      y = "MoCA Total Score (24-30)"
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

























library(emmeans)
library(ggplot2)
library(dplyr)

# ==============================================================================
# 1. DEFINE THE MASTER VISUALIZATION ENGINE
# ==============================================================================
generate_trajectory_plot <- function(model_obj, var_name, title_text, legend_text, 
                                     y1 = 23.4, y2 = 22.7, y3 = 22.0) {
  
  # A. Dynamic list allocation for emtrends/emmeans slicing
  at_list <- list()
  at_list[[var_name]] <- c(-1, 0, 1)
  
  # Build spec formulas dynamically from strings
  spec_formula <- as.formula(paste("~", var_name))
  preds_formula <- as.formula(paste("~ years_from_baseline *", var_name))
  
  # B. Extract trends/slopes with math expression parsing
  slopes <- emtrends(model_obj, specs = spec_formula, 
                     var = "years_from_baseline", at = at_list)
  slopes_summary <- summary(slopes, infer = TRUE)
  
  slopes_summary <- slopes_summary %>%
    mutate(
      type = c("Low (-1 SD)", "Mean (0 SD)", "High (+1 SD)"),
      label = paste0("'", type, ":' ~ beta == ", round(years_from_baseline.trend, 2), 
                     " ~~~~~ p == ", round(p.value, 3))
    )
  
  # C. Generate line coordinates across a continuous grid (0 to 6 years)
  at_grid <- list(years_from_baseline = seq(0, 6, by = 0.1))
  at_grid[[var_name]] <- c(-1, 0, 1)
  
  preds_df <- emmeans(model_obj, specs = preds_formula, at = at_grid) %>% as.data.frame()
  
  # Map grouping factor
  preds_df$Biomarker_Level <- factor(preds_df[[var_name]], 
                                     levels = c(-1, 0, 1), 
                                     labels = c("Low (-1 SD)", "Mean (0 SD)", "High (+1 SD)"))
  
  # D. Build ggplot following exact user visual layout parameters
  p <- ggplot() +
    # Uniform gray background spaghetti data
    geom_line(data = MRS_prediction_long, 
              aes(x = years_from_baseline, y = moca, group = pscid), 
              color = "grey80", alpha = 0.3, linewidth = 0.4) +
    geom_point(data = MRS_prediction_long, 
               aes(x = years_from_baseline, y = moca), 
               color = "grey75", alpha = 0.4, size = 1) +
    
    # Solid, non-graded model prediction paths
    geom_line(data = preds_df, 
              aes(x = years_from_baseline, y = emmean, group = Biomarker_Level, color = Biomarker_Level), 
              linewidth = 2) + 
    
    # Exact publication hex color mappings
    scale_color_manual(values = c("Low (-1 SD)" = "#D55E00",  
                                  "Mean (0 SD)" = "#737373",  
                                  "High (+1 SD)" = "#3B5998"), 
                       name = legend_text) +
    
    # Mathematical expression text overlays
    annotate("text", x = 0.2, y = y1, label = slopes_summary$label[1], parse = TRUE, hjust = 0, color = "#D55E00", size = 4.5) +
    annotate("text", x = 0.2, y = y2, label = slopes_summary$label[2], parse = TRUE, hjust = 0, color = "#737373", size = 4.5) +
    annotate("text", x = 0.2, y = y3, label = slopes_summary$label[3], parse = TRUE, hjust = 0, color = "#3B5998", size = 4.5) +
    
    scale_x_continuous(limits = c(0, 6), breaks = c(0, 2, 4, 6)) +
    coord_cartesian(ylim = c(21.5, 30)) + 
    labs(
      title = title_text,
      subtitle = "Model-derived trajectories evaluated at three distinct baseline levels",
      x = "Years from Baseline",
      y = "MoCA Total Score (0-30)"
    ) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), legend.position = "right")
  
  return(p)
}


# ==============================================================================
# 2. RUN THE ENGINE FOR ALL UNIMODAL PLOTS
# ==============================================================================

# Note: Adjust y1, y2, y3 coordinates as needed per plot if prediction lines 
# slide upward or downward and overlap your text labels.

# --- TRACK 1: METABOLIC (Glutamate) ---
plot_precuneus <- generate_trajectory_plot(
  model_obj   = mixed_model_precuneus, 
  var_name    = "m_m_precuneus_z", 
  title_text  = "Continuous Cognitive Trajectories: Precuneus Glutamate", 
  legend_text = "Baseline Precuneus Glu",
  y1 = 23.4, y2 = 22.7, y3 = 22.0
)

plot_acc <- generate_trajectory_plot(
  model_obj   = mixed_model_acc, 
  var_name    = "m_m_acc_z", 
  title_text  = "Continuous Cognitive Trajectories: ACC Glutamate", 
  legend_text = "Baseline ACC Glutamate",
  y1 = 23.4, y2 = 22.7, y3 = 22.0
)


# --- TRACK 2: FLUID (Pathology) ---
plot_ptau217 <- generate_trajectory_plot(
  model_obj   = mixed_model_ptau217, 
  var_name    = "plasma_ptau217_z", 
  title_text  = "Continuous Cognitive Trajectories: Plasma p-Tau217", 
  legend_text = "Baseline Plasma p-Tau",
  y1 = 23.4, y2 = 22.7, y3 = 22.0
)


# --- TRACK 3: STRUCTURAL (Atrophy) ---
plot_thickness <- generate_trajectory_plot(
  model_obj   = mixed_model_thickness, 
  var_name    = "cortical_thickness_adsignature_dickson_z", 
  title_text  = "Continuous Cognitive Trajectories: Cortical Thickness", 
  legend_text = "Baseline Thickness",
  y1 = 23.4, y2 = 22.7, y3 = 22.0
)

plot_hipp_vol <- generate_trajectory_plot(
  model_obj   = mixed_model_hipp_mean, 
  var_name    = "hipp_mean_z", 
  title_text  = "Continuous Cognitive Trajectories: Hippocampal Volume", 
  legend_text = "Baseline Hipp Volume",
  y1 = 23.4, y2 = 22.7, y3 = 22.0
)


# --- TRACK 4: FUNCTIONAL (fMRI Activation) ---
plot_hipp_act <- generate_trajectory_plot(
  model_obj   = mixed_model_hipp_mean_act, 
  var_name    = "hipp_mean_act_z", 
  title_text  = "Continuous Cognitive Trajectories: Hippocampal Activation", 
  legend_text = "Baseline Hipp Activation",
  y1 = 23.4, y2 = 22.7, y3 = 22.0
)

plot_parietal_act <- generate_trajectory_plot(
  model_obj   = mixed_model_activation_parietal_l, 
  var_name    = "activation_parietal_sup_l_z", 
  title_text  = "Continuous Cognitive Trajectories: Left Parietal Activation", 
  legend_text = "Baseline Parietal Act",
  y1 = 23.4, y2 = 22.7, y3 = 22.0
)


# ==============================================================================
# 3. DISPLAY OR SAVE CHOSEN TARGET PLOTS
# ==============================================================================
# To view a plot instantly in RStudio, just call its object:
print(plot_ptau217)
print(plot_acc)
print(plot_precuneus)
print(plot_thickness)
print(plot_hipp_vol)
print(plot_hipp_act)
print(plot_parietal_act)





# Example for your Precuneus model
model_sum <- summary(mixed_model_precuneus)
interaction_coef <- round(model_sum$coefficients["years_from_baseline:m_m_precuneus_z", "Estimate"], 3)
interaction_p    <- round(model_sum$coefficients["years_from_baseline:m_m_precuneus_z", "Pr(>|t|)"], 4)

cat("Precuneus Glutamate Interaction Slope:", interaction_coef, "\n")
cat("Precuneus Glutamate Interaction p-value:", interaction_p, "\n")



library(emmeans)

# --- 1. Hippocampal Activation ---
# This calculates the slope of MoCA decline at -1 SD, Mean, and +1 SD
hipp_slopes <- emtrends(mixed_model_hipp_mean_act, 
                        specs = ~ hipp_mean_act_z, 
                        var = "years_from_baseline", 
                        at = list(hipp_mean_act_z = c(-1, 0, 1)))

# This table shows the p-value for the slope of each group
summary(hipp_slopes, infer = TRUE)

# --- 2. Parietal Activation ---
parietal_slopes <- emtrends(mixed_model_activation_parietal_l, 
                            specs = ~ activation_parietal_sup_l_z, 
                            var = "years_from_baseline", 
                            at = list(activation_parietal_sup_l_z = c(-1, 0, 1)))

summary(parietal_slopes, infer = TRUE)






library(emmeans)

# --- 1. Hippocampal Activation ---
# This calculates the slope of MoCA decline at -1 SD, Mean, and +1 SD
hipp_slopes <- emtrends(mixed_model_hipp_mean_act, 
                        specs = ~ hipp_mean_act_z, 
                        var = "years_from_baseline", 
                        at = list(hipp_mean_act_z = c(-1, 0, 1)))

# This table shows the p-value for the slope of each group
summary(hipp_slopes, infer = TRUE)

# --- 2. Parietal Activation ---
parietal_slopes <- emtrends(mixed_model_activation_parietal_l, 
                            specs = ~ activation_parietal_sup_l_z, 
                            var = "years_from_baseline", 
                            at = list(activation_parietal_sup_l_z = c(-1, 0, 1)))

summary(parietal_slopes, infer = TRUE)

#### ACC ####
library(dplyr)
library(ggplot2)
library(emmeans)

# 1. Extract slopes from the continuous ACC model
acc_slopes <- emtrends(mixed_model_acc, ~ m_m_acc_z, 
                       var = "years_from_baseline", 
                       at = list(m_m_acc_z = c(-1, 0, 1)))
acc_slopes_summary <- summary(acc_slopes, infer = TRUE)

# 2. Wrap text in single quotes for clean math parsing
acc_slopes_summary <- acc_slopes_summary %>%
  mutate(
    type = c("Low (-1 SD)", "Mean (0 SD)", "High (+1 SD)"),
    label = paste0("'", type, ":' ~ beta == ", round(years_from_baseline.trend, 2), 
                   " ~~~~~ p == ", round(p.value, 3))
  )

# 3. Generate model prediction lines and define discrete factor levels
acc_preds <- emmeans(mixed_model_acc, ~ years_from_baseline * m_m_acc_z,
                     at = list(years_from_baseline = seq(0, 6, by = 0.1), 
                               m_m_acc_z = c(-1, 0, 1))) %>% 
  as.data.frame()

acc_preds$Glutamate <- factor(acc_preds$m_m_acc_z, 
                              levels = c(-1, 0, 1), 
                              labels = c("Low (-1 SD)", "Mean (0 SD)", "High (+1 SD)"))

# 4. SIMPLIFIED DISCRETE PLOT FOR ACC
ggplot() +
  # Raw background data set to a neutral, uniform gray to eliminate the gradient
  geom_line(data = MRS_prediction_long, 
            aes(x = years_from_baseline, y = moca, group = pscid), 
            color = "grey80", alpha = 0.3, linewidth = 0.4) +
  geom_point(data = MRS_prediction_long, 
             aes(x = years_from_baseline, y = moca), 
             color = "grey75", alpha = 0.4, size = 1) +
  
  # Model prediction lines colored dynamically by their discrete factor level
  geom_line(data = acc_preds, 
            aes(x = years_from_baseline, y = emmean, group = Glutamate, color = Glutamate), 
            linewidth = 2) + 
  
  # Define exactly three solid, non-graded colors manually
  scale_color_manual(values = c("Low (-1 SD)" = "#D55E00",   # Orange/Red
                                "Mean (0 SD)" = "#737373",  # Neutral Dark Gray
                                "High (+1 SD)" = "#3B5998"), # Blue
                     name = "Baseline Glutamate") +
  
  # Text annotations matching the exact line colors
  annotate("text", x = 0.2, y = 23.4, label = acc_slopes_summary$label[1], parse = TRUE, hjust = 0, color = "#D55E00", size = 4.5) +
  annotate("text", x = 0.2, y = 22.7, label = acc_slopes_summary$label[2], parse = TRUE, hjust = 0, color = "#737373", size = 4.5) +
  annotate("text", x = 0.2, y = 22.0, label = acc_slopes_summary$label[3], parse = TRUE, hjust = 0, color = "#3B5998", size = 4.5) +
  
  scale_x_continuous(limits = c(0, 6), breaks = c(0, 2, 4, 6)) +
  coord_cartesian(ylim = c(21.5, 30)) + 
  
  labs(
    title = "Continuous Cognitive Trajectories: ACC Glutamate",
    subtitle = "Model-derived trajectories evaluated at three distinct baseline levels",
    x = "Years from Baseline",
    y = "MoCA Total Score (0-30)"
  ) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), legend.position = "right")




########### Tau ############
library(dplyr)
library(ggplot2)
library(emmeans)

# 1. Extract slopes from the continuous p-Tau217 model
ptau_slopes <- emtrends(mixed_model_ptau217, ~ plasma_ptau217_z, 
                        var = "years_from_baseline", 
                        at = list(plasma_ptau217_z = c(-1, 0, 1)))
ptau_slopes_summary <- summary(ptau_slopes, infer = TRUE)

# 2. Wrap text in single quotes for clean math parsing
ptau_slopes_summary <- ptau_slopes_summary %>%
  mutate(
    type = c("Low (-1 SD)", "Mean (0 SD)", "High (+1 SD)"),
    label = paste0("'", type, ":' ~ beta == ", round(years_from_baseline.trend, 2), 
                   " ~~~~~ p == ", round(p.value, 3))
  )

# 3. Generate model prediction lines and define discrete factor levels
ptau_preds <- emmeans(mixed_model_ptau217, ~ years_from_baseline * plasma_ptau217_z,
                      at = list(years_from_baseline = seq(0, 6, by = 0.1), 
                                plasma_ptau217_z = c(-1, 0, 1))) %>% 
  as.data.frame()

ptau_preds$pTau217 <- factor(ptau_preds$plasma_ptau217_z, 
                             levels = c(-1, 0, 1), 
                             labels = c("Low (-1 SD)", "Mean (0 SD)", "High (+1 SD)"))

# 4. PLOT FOR PLASMA p-TAU217
ggplot() +
  # Raw background data set to neutral gray
  geom_line(data = MRS_prediction_long, 
            aes(x = years_from_baseline, y = moca, group = pscid), 
            color = "grey80", alpha = 0.3, linewidth = 0.4) +
  geom_point(data = MRS_prediction_long, 
             aes(x = years_from_baseline, y = moca), 
             color = "grey75", alpha = 0.4, size = 1) +
  
  # Model prediction lines
  geom_line(data = ptau_preds, 
            aes(x = years_from_baseline, y = emmean, group = pTau217, color = pTau217), 
            linewidth = 2) + 
  
  # Color Mapping: Low p-Tau is blue (protective), High p-Tau is orange (risk)
  scale_color_manual(values = c("Low (-1 SD)" = "#3B5998",   # Blue
                                "Mean (0 SD)" = "#737373",  # Neutral Dark Gray
                                "High (+1 SD)" = "#D55E00"), # Orange/Red
                     name = "Baseline p-Tau217") +
  
  # Text annotations matching the respective line color assignments
  annotate("text", x = 0.2, y = 23.4, label = ptau_slopes_summary$label[1], parse = TRUE, hjust = 0, color = "#3B5998", size = 4.5) +
  annotate("text", x = 0.2, y = 22.7, label = ptau_slopes_summary$label[2], parse = TRUE, hjust = 0, color = "#737373", size = 4.5) +
  annotate("text", x = 0.2, y = 22.0, label = ptau_slopes_summary$label[3], parse = TRUE, hjust = 0, color = "#D55E00", size = 4.5) +
  
  scale_x_continuous(limits = c(0, 6), breaks = c(0, 2, 4, 6)) +
  coord_cartesian(ylim = c(21.5, 30)) + 
  
  labs(
    title = "Continuous Cognitive Trajectories: Plasma p-Tau217",
    subtitle = "Model-derived trajectories evaluated at three distinct baseline levels",
    x = "Years from Baseline",
    y = "MoCA Total Score (0-30)"
  ) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), legend.position = "right")



### Structure

# 1. Extract slopes from the continuous Cortical Thickness model
thick_slopes <- emtrends(mixed_model_thickness, ~ cortical_thickness_adsignature_dickson_z, 
                         var = "years_from_baseline", 
                         at = list(cortical_thickness_adsignature_dickson_z = c(-1, 0, 1)))
thick_slopes_summary <- summary(thick_slopes, infer = TRUE)

# 2. Wrap text in single quotes for clean math parsing
thick_slopes_summary <- thick_slopes_summary %>%
  mutate(
    type = c("Low (-1 SD)", "Mean (0 SD)", "High (+1 SD)"),
    label = paste0("'", type, ":' ~ beta == ", round(years_from_baseline.trend, 2), 
                   " ~~~~~ p == ", round(p.value, 3))
  )

# 3. Generate model prediction lines and define discrete factor levels
thick_preds <- emmeans(mixed_model_thickness, ~ years_from_baseline * cortical_thickness_adsignature_dickson_z,
                       at = list(years_from_baseline = seq(0, 6, by = 0.1), 
                                 cortical_thickness_adsignature_dickson_z = c(-1, 0, 1))) %>% 
  as.data.frame()

thick_preds$Thickness <- factor(thick_preds$cortical_thickness_adsignature_dickson_z, 
                                levels = c(-1, 0, 1), 
                                labels = c("Low (-1 SD)", "Mean (0 SD)", "High (+1 SD)"))

# 4. PLOT FOR CORTICAL THICKNESS
ggplot() +
  # Raw background data set to neutral gray
  geom_line(data = MRS_prediction_long, 
            aes(x = years_from_baseline, y = moca, group = pscid), 
            color = "grey80", alpha = 0.3, linewidth = 0.4) +
  geom_point(data = MRS_prediction_long, 
             aes(x = years_from_baseline, y = moca), 
             color = "grey75", alpha = 0.4, size = 1) +
  
  # Model prediction lines
  geom_line(data = thick_preds, 
            aes(x = years_from_baseline, y = emmean, group = Thickness, color = Thickness), 
            linewidth = 2) + 
  
  # Color Mapping: Low thickness is orange (at-risk), High thickness is blue (protective)
  scale_color_manual(values = c("Low (-1 SD)" = "#D55E00",   # Orange/Red
                                "Mean (0 SD)" = "#737373",  # Neutral Dark Gray
                                "High (+1 SD)" = "#3B5998"), # Blue
                     name = "Baseline Thickness") +
  
  # Text annotations matching the respective line color assignments
  annotate("text", x = 0.2, y = 23.4, label = thick_slopes_summary$label[1], parse = TRUE, hjust = 0, color = "#D55E00", size = 4.5) +
  annotate("text", x = 0.2, y = 22.7, label = thick_slopes_summary$label[2], parse = TRUE, hjust = 0, color = "#737373", size = 4.5) +
  annotate("text", x = 0.2, y = 22.0, label = thick_slopes_summary$label[3], parse = TRUE, hjust = 0, color = "#3B5998", size = 4.5) +
  
  scale_x_continuous(limits = c(0, 6), breaks = c(0, 2, 4, 6)) +
  coord_cartesian(ylim = c(21.5, 30)) + 
  
  labs(
    title = "Continuous Cognitive Trajectories: AD-Signature Cortical Thickness",
    subtitle = "Model-derived trajectories evaluated at three distinct baseline levels",
    x = "Years from Baseline",
    y = "MoCA Total Score (0-30)"
  ) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), legend.position = "right")


  











  
  library(tidyr); library(dplyr); library(lme4); library(lmerTest); library(emmeans); library(ggplot2)

# --- DATA PREP ---
MRS_decline_long <- MRS_long %>%
  mutate(decliners = factor(decliners, levels = c("0", "1"))) %>%
  select(pscid, decliners, age_difference, t1_glu_acc, t2_glu_acc, t1_glu_prec, t2_glu_prec) %>%
  pivot_longer(cols = matches("glu_"), names_to = c("time", "region"), names_pattern = "(t[12])_glu_(.*)") %>%
  mutate(years_elapsed = ifelse(time == "t1", 0, age_difference))

# --- PRECUNEUS MODEL ---
mod_prec_dec <- lmer(value ~ decliners * years_elapsed + (1 | pscid), 
                     data = filter(MRS_decline_long, region == "prec"))
stats_prec_dec <- summary(emtrends(mod_prec_dec, ~ decliners, var = "years_elapsed"), infer = TRUE)
lab_prec_dec <- stats_prec_dec %>% mutate(label = paste0("beta == ", round(years_elapsed.trend, 2), "~~p == ", round(p.value, 3)))

# --- ACC MODEL ---
mod_acc_dec <- lmer(value ~ decliners * years_elapsed + (1 | pscid), 
                    data = filter(MRS_decline_long, region == "acc"))
stats_acc_dec <- summary(emtrends(mod_acc_dec, ~ decliners, var = "years_elapsed"), infer = TRUE)
lab_acc_dec <- stats_acc_dec %>% mutate(label = paste0("beta == ", round(years_elapsed.trend, 2), "~~p == ", round(p.value, 3)))




ggplot(filter(MRS_decline_long, region == "prec"), aes(x = years_elapsed, y = value, color = decliners)) +
  geom_line(aes(group = pscid), alpha = 0.2) + geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", linewidth = 1.5, se = TRUE, aes(fill = decliners)) +
  facet_wrap(~ decliners) +
  geom_text(data = lab_prec_dec, aes(x = 1.5, y = Inf, label = label), vjust = 2, parse = TRUE, inherit.aes = FALSE) +
  scale_x_continuous(limits = c(0, 3)) +
  scale_color_manual(values = c("0" = "#3B5998", "1" = "#D55E00")) +
  scale_fill_manual(values = c("0" = "#3B5998", "1" = "#D55E00")) +
  labs(title = "Precuneus Glutamate: Decliners (1) vs. Non-Decliners (0)", x = "Years from Baseline", y = "Glutamate (mM)") +
  theme_minimal() + theme(legend.position = "none")







ggplot(filter(MRS_decline_long, region == "acc"), aes(x = years_elapsed, y = value, color = decliners)) +
  # Individual patient trajectories
  geom_line(aes(group = pscid), alpha = 0.2, linewidth = 0.4) + 
  geom_point(alpha = 0.5) +
  # Group-level annualized slope with 95% CI
  geom_smooth(method = "lm", linewidth = 1.5, se = TRUE, aes(fill = decliners)) +
  # Separate into Decliners (1) and Non-Decliners (0)
  facet_wrap(~ decliners) +
  # Add Stats from the ACC model (Beta = change per year)
  geom_text(data = lab_acc_dec, aes(x = 1.5, y = Inf, label = label), 
            vjust = 2, color = "black", parse = TRUE, inherit.aes = FALSE) +
  # Visual formatting
  scale_x_continuous(limits = c(0, 3), breaks = c(0, 1, 2, 3)) +
  scale_color_manual(values = c("0" = "#3B5998", "1" = "#D55E00")) +
  scale_fill_manual(values = c("0" = "#3B5998", "1" = "#D55E00")) +
  labs(
    title = "ACC Glutamate: Decliners (1) vs. Non-Decliners (0)",
    subtitle = "Annualized trajectories; beta represents change in mM per year",
    x = "Years from Baseline (Capped at 3 Years)",
    y = "Glutamate (mM)"
  ) +
  theme_minimal() + 
  theme(legend.position = "none", panel.grid.minor = element_blank())













######## Regressions 

library(ggplot2)
library(tidyr)
library(dplyr)

# 1. Prepare the data for faceting by region
plot_data_t1 <- MRS_long %>%
  select(pscid, diagnostic_nick, slope_moca_raw, t1_glu_prec, t1_glu_acc) %>%
  pivot_longer(cols = c(t1_glu_prec, t1_glu_acc), 
               names_to = "region", 
               values_to = "t1_glutamate") %>%
  mutate(region = ifelse(region == "t1_glu_prec", "Precuneus", "ACC"),
         diagnostic_nick = factor(diagnostic_nick, levels = c("HC", "SCD+", "MCI")))

# 2. Create the regression plot
ggplot(plot_data_t1, aes(x = t1_glutamate, y = slope_moca_raw)) +
  # Add points colored by clinical group
  geom_point(aes(color = diagnostic_nick), alpha = 0.6, size = 2.5) +
  # Add the linear regression line for the whole cohort
  geom_smooth(method = "lm", color = "black", linewidth = 1.2, se = TRUE) +
  # Facet by region (Precuneus vs. ACC)
  facet_wrap(~ region, scales = "free_x") +
  # Styling to match the paper
  scale_color_manual(values = c("HC" = "#3B5998", "SCD+" = "#E69F00", "MCI" = "#D55E00")) +
  labs(
    title = "Baseline Glutamate (T1) vs. Cognitive Decline (MoCA Slope)",
    subtitle = "",
    x = "Baseline Glutamate Concentration (mM)",
    y = "MoCA Slope (Change in Score per Year)",
    color = "Diagnosis"
  ) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())





# 1. Prepare the data for faceting by region
plot_data_t1 <- MRS_long %>%
  select(pscid, diagnostic_nick, slope_moca_raw, t1_glu_prec, t1_glu_acc) %>%
  pivot_longer(cols = c(t1_glu_prec, t1_glu_acc), 
               names_to = "region", 
               values_to = "t1_glutamate") %>%
  mutate(region = ifelse(region == "t1_glu_prec", "Precuneus", "ACC"),
         diagnostic_nick = factor(diagnostic_nick, levels = c("HC", "SCD+", "MCI")))

# 2. Create the regression plot
ggplot(plot_data_t1, aes(x = t1_glutamate, y = slope_moca_raw)) +
  # Add points colored by clinical group
  geom_point(aes(color = diagnostic_nick), alpha = 0.6, size = 2.5) +
  # Add the linear regression line for the whole cohort
  geom_smooth(method = "lm", color = "black", linewidth = 1.2, se = TRUE) +
  # Facet by region (Precuneus vs. ACC)
  facet_wrap(~ region, scales = "free_x") +
  # Styling to match the paper
  scale_color_manual(values = c("HC" = "#3B5998", "SCD+" = "#E69F00", "MCI" = "#D55E00")) +
  labs(
    title = "Baseline Glutamate (T1) vs. Cognitive Decline (MoCA Slope)",
    subtitle = "",
    x = "Baseline Glutamate Concentration (mM)",
    y = "MoCA Slope (Change in Score per Year)",
    color = "Diagnosis"
  ) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())





# 1. Prepare the data for faceting by region
plot_data_t2 <- MRS_long %>%
  select(pscid, diagnostic_nick, slope_moca_raw, t2_glu_prec, t2_glu_acc) %>%
  pivot_longer(cols = c(t2_glu_prec, t2_glu_acc), 
               names_to = "region", 
               values_to = "t2_glutamate") %>%
  mutate(region = ifelse(region == "t2_glu_prec", "Precuneus", "ACC"),
         diagnostic_nick = factor(diagnostic_nick, levels = c("HC", "SCD+", "MCI")))

# 2. Create the regression plot
ggplot(plot_data_t2, aes(x = t2_glutamate, y = slope_moca_raw)) +
  # Add points colored by clinical group
  geom_point(aes(color = diagnostic_nick), alpha = 0.6, size = 2.5) +
  # Add the linear regression line for the whole cohort
  geom_smooth(method = "lm", color = "black", linewidth = 1.2, se = TRUE) +
  # Facet by region (Precuneus vs. ACC)
  facet_wrap(~ region, scales = "free_x") +
  # Styling to match the paper
  scale_color_manual(values = c("HC" = "#3B5998", "SCD+" = "#E69F00", "MCI" = "#D55E00")) +
  labs(
    title = "Glutamate (T2) vs. Cognitive Decline (MoCA Slope)",
    subtitle = "",
    x = "T2 Glutamate Concentration (mM)",
    y = "MoCA Slope (Change in Score per Year)",
    color = "Diagnosis"
  ) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())








#### Regressions #####



library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Prepare data (ensure groups are ordered)
plot_data_t1_groups <- MRS_long %>%
  select(pscid, diagnostic_nick, slope_moca_raw, t1_glu_prec, t1_glu_acc) %>%
  pivot_longer(cols = c(t1_glu_prec, t1_glu_acc), 
               names_to = "region", 
               values_to = "t1_glutamate") %>%
  mutate(region = ifelse(region == "t1_glu_prec", "Precuneus", "ACC"),
         diagnostic_nick = factor(diagnostic_nick, levels = c("HC", "SCD+", "MCI")))

# 2. Create group-wise regression plot
ggplot(plot_data_t1_groups, aes(x = t1_glutamate, y = slope_moca_raw, color = diagnostic_nick, fill = diagnostic_nick)) +
  geom_point(alpha = 0.6, size = 2.5) +
  # Add separate regression lines for each diagnostic group
  geom_smooth(method = "lm", alpha = 0.2, linewidth = 1.2) +
  # Facet by region
  facet_wrap(~ region, scales = "free_x") +
  # Styling
  scale_color_manual(values = c("HC" = "#3B5998", "SCD+" = "#E69F00", "MCI" = "#D55E00")) +
  scale_fill_manual(values = c("HC" = "#3B5998", "SCD+" = "#E69F00", "MCI" = "#D55E00")) +
  labs(
    title = "Group-Wise Associations: T1 Glutamate vs. MoCA Slope",
    subtitle = "Testing if the baseline glutamate-cognition link differs by stage",
    x = "Baseline Glutamate Concentration (mM)",
    y = "MoCA Slope (Change in Score per Year)",
    color = "Diagnosis",
    fill = "Diagnosis"
  ) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())






### Structual ####
library(ggplot2)
library(dplyr)
library(patchwork) # This package is magic for combining separate plots

# Ensure your diagnostic group is ordered properly
MRS_long <- MRS_long %>%
  mutate(diagnostic_nick = factor(diagnostic_nick, levels = c("HC", "SCD+", "MCI")))

# =========================================================================
# 1. HELPER FUNCTION: Builds a plot with specific labels and calculates stats
# =========================================================================
create_scatter <- function(data, x_var, y_var, x_lab, y_lab) {
  
  # A. Calculate the overall model stats for the whole sample
  mod <- lm(data[[y_var]] ~ data[[x_var]], data = data)
  beta_val <- round(summary(mod)$coefficients[2, "Estimate"], 3)
  p_val <- summary(mod)$coefficients[2, "Pr(>|t|)"]
  
  # B. Format the p-value nicely
  p_text <- ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", round(p_val, 3)))
  
  # C. Create the label string (using the unicode character for beta: β)
  stat_label <- paste0("\u03B2 = ", beta_val, " | ", p_text)
  
  # D. Build the plot
  ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(aes(color = diagnostic_nick), alpha = 0.6, size = 2.5) +
    geom_smooth(method = "lm", color = "black", fill = "gray70", alpha = 0.3, linewidth = 1.2) +
    scale_color_manual(values = c("HC" = "#3B5998", "SCD+" = "#E69F00", "MCI" = "#D55E00")) +
    
    # Specific X and Y labels applied here!
    labs(x = x_lab, y = y_lab) + 
    
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      plot.margin = margin(15, 15, 15, 15)
    ) +
    # Stamp the stats dynamically in the top-left corner (-Inf, Inf)
    annotate("text", x = -Inf, y = Inf, label = stat_label, 
             hjust = -0.1, vjust = 1.5, size = 4.5, fontface = "bold", color = "black")
}

# =========================================================================
# 2. GENERATE THE 4 INDIVIDUAL PLOTS
# =========================================================================
p1 <- create_scatter(MRS_long, "t1_glu_acc", "t2_hipp_e_tiv", 
                     "Baseline ACC Glutamate", "T2 Hippocampal Volume")

p2 <- create_scatter(MRS_long, "t1_glu_acc", "t2_cortical_thickness_dickson", 
                     "Baseline ACC Glutamate", "T2 Cortical Thickness")

p3 <- create_scatter(MRS_long, "t1_glu_prec", "t2_hipp_e_tiv", 
                     "Baseline Precuneus Glutamate", "T2 Hippocampal Volume")

p4 <- create_scatter(MRS_long, "t1_glu_prec", "t2_cortical_thickness_dickson", 
                     "Baseline Precuneus Glutamate", "T2 Cortical Thickness")

# =========================================================================
# 3. COMBINE THEM INTO A 2x2 GRID USING PATCHWORK
# =========================================================================
# The (p1 | p2) / (p3 | p4) syntax tells patchwork to put 1 and 2 on top, 3 and 4 on bottom
final_plot <- (p1 | p2) / (p3 | p4) + 
  plot_layout(guides = 'collect') + # This merges the 4 legends into one single legend on the side
  plot_annotation(
    title = "Associations Between Baseline Glutamate and Follow-up Brain Structure",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

# View the final result
final_plot







######## Regressions (Quadratic Fit)

library(ggplot2)
library(tidyr)
library(dplyr)

# 1. Prepare the data for faceting by region
plot_data_t1 <- MRS_long %>%
  select(pscid, diagnostic_nick, slope_moca_raw, t1_glu_prec, t1_glu_acc) %>%
  pivot_longer(cols = c(t1_glu_prec, t1_glu_acc), 
               names_to = "region", 
               values_to = "t1_glutamate") %>%
  mutate(region = ifelse(region == "t1_glu_prec", "Precuneus", "ACC"),
         diagnostic_nick = factor(diagnostic_nick, levels = c("HC", "SCD+", "MCI")))

# 2. Create the quadratic regression plot
ggplot(plot_data_t1, aes(x = t1_glutamate, y = slope_moca_raw)) +
  # Add points colored by clinical group
  geom_point(aes(color = diagnostic_nick), alpha = 0.6, size = 2.5) +
  # Add the QUADRATIC regression line for the whole cohort
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "black", linewidth = 1.2, se = TRUE) +
  # Facet by region (Precuneus vs. ACC)
  facet_wrap(~ region, scales = "free_x") +
  # Styling to match the paper
  scale_color_manual(values = c("HC" = "#3B5998", "SCD+" = "#E69F00", "MCI" = "#D55E00")) +
  labs(
    title = "Baseline Glutamate (T1) vs. Cognitive Decline (MoCA Slope)",
    subtitle = "Quadratic Fit",
    x = "Baseline Glutamate Concentration (mM)",
    y = "MoCA Slope (Change in Score per Year)",
    color = "Diagnosis"
  ) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())












##### Activation #########








library(ggplot2)
library(tidyr)
library(dplyr)

# 1. Prepare the data for faceting by region
plot_data_t1 <- MRS_long %>%
  # Swap out slope_moca_raw for hipp_l_diff_activation
  select(pscid, diagnostic_nick, hipp_l_diff_activation, t1_glu_prec, t1_glu_acc) %>%
  pivot_longer(cols = c(t1_glu_prec, t1_glu_acc), 
               names_to = "region", 
               values_to = "t1_glutamate") %>%
  mutate(region = ifelse(region == "t1_glu_prec", "Precuneus", "ACC"))

# 2. Create the quadratic regression plot
ggplot(plot_data_t1, aes(x = t1_glutamate, y = hipp_l_diff_activation)) +
  # Add points colored by clinical group
  geom_point(aes(color = diagnostic_nick), alpha = 0.6, size = 2.5) +
  
  # Add the QUADRATIC regression line matching your exact lm() formula
  geom_smooth(method = "lm", formula = y ~ x, color = "black", linewidth = 1.2, se = TRUE) +
  
  # Facet by region (Precuneus vs. ACC)
  facet_wrap(~ region, scales = "free_x") +
  
  # Styling (Updated to account for your combined SCD+ and MCI group)
  scale_color_manual(values = c("HC" = "#3B5998", 
                                "SCD+" = "#E69F00", 
                                "MCI" = "#D55E00", 
                                "SCD+ and MCI" = "#D55E00")) +
  labs(
    title = "",
    subtitle = "",
    x = "Baseline Glutamate Concentration (mM)",
    y = "Annualized Change in Left Hippocampal Activation",
    color = "Diagnosis"
  ) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())














### ROC curves 

library(pROC)

# Set up the plotting area
# We start by plotting the strongest individual predictor (Cortical Thickness)
plot(roc_struc_thick, 
     col = "#1f77b4", # Professional Blue
     lwd = 3, 
     main = "ROC Curves: Significant Predictors of Cognitive Decline",
     xlab = "1 - Specificity", 
     ylab = "Sensitivity",
     legacy.axes = TRUE) # This ensures the x-axis is 0 to 1 (1-Specificity)

# Add the other three significant models
plot(roc_struc_hip, add = TRUE, col = "#ff7f0e", lwd = 3)    # Orange
plot(roc_glu_prec, add = TRUE, col = "#2ca02c", lwd = 3)     # Green
plot(roc_glu_acc, add = TRUE, col = "#d62728", lwd = 3)      # Red

# Add a reference line (chance level)
abline(a = 0, b = 1, lty = 2, col = "darkgrey")

# Add a legend with AUC values for clarity
legend("bottomright", 
       legend = c(
         paste0("Cortical Thickness (AUC = ", round(auc(roc_struc_thick), 3), ")"),
         paste0("Hippocampal Volume (AUC = ", round(auc(roc_struc_hip), 3), ")"),
         paste0("Precuneus Glutamate (AUC = ", round(auc(roc_glu_prec), 3), ")"),
         paste0("ACC Glutamate (AUC = ", round(auc(roc_glu_acc), 3), ")")
       ),
       col = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"), 
       lwd = 3, 
       cex = 0.8, 
       bty = "n") # Removes the box around the legend




## Combined models 2 variables 

# Set up the plotting area with your top-performing model
# Metab (ACC) + Func (Parietal)
plot(roc_metab_func, 
     col = "#1f77b4", # Professional Blue
     lwd = 3, 
     main = "Multimodal ROC Curves: Top 2-Variable Predictors",
     xlab = "1 - Specificity", 
     ylab = "Sensitivity",
     legacy.axes = TRUE) # Standard 0 to 1 axis

# Add the other two winning models
plot(roc_metab_struc_hip, add = TRUE, col = "#ff7f0e", lwd = 3)    # Orange
plot(roc_metab_struc_thick, add = TRUE, col = "#2ca02c", lwd = 3)  # Green

# Add a reference line (chance level)
abline(a = 0, b = 1, lty = 2, col = "darkgrey")

# Add a professional legend
legend("bottomright", 
       legend = c(
         paste0("Glu ACC + Parietal L Activaiton (AUC = ", round(auc(roc_metab_func), 3), ")"),
         paste0("Glu Precuneus + Hip Vol (AUC = ", round(auc(roc_metab_struc_hip), 3), ")"),
         paste0("Glu Precuneus + Thickness (AUC = ", round(auc(roc_metab_struc_thick), 3), ")")
       ),
       col = c("#1f77b4", "#ff7f0e", "#2ca02c"), 
       lwd = 3, 
       cex = 0.75, 
       bty = "n")




## 3 variables 
library(pROC)

# Set up the plotting area with the top-performing 3-variable model
# Model 1: Metab + Struc + Func (AUC = 0.784)
plot(roc_3var_1, 
     col = "#1f77b4", # Professional Blue
     lwd = 3, 
     main = "Exploratory Multimodal ROC Curves: 3-Variable Predictors",
     xlab = "1 - Specificity", 
     ylab = "Sensitivity",
     legacy.axes = TRUE) # Standard 0 to 1 axis

# Add the other two winning models
# Model 2: Metab(x2) + Func (AUC = 0.771)
plot(roc_3var_2, add = TRUE, col = "#ff7f0e", lwd = 3)    # Orange

# Model 3: Struc + Func(x2) (AUC = 0.767)
plot(roc_3var_3, add = TRUE, col = "#2ca02c", lwd = 3)    # Green

# Add a reference line (chance level)
abline(a = 0, b = 1, lty = 2, col = "darkgrey")

# Add a professional legend
legend("bottomright", 
       legend = c(
         paste0("Glu Prec + Hip Vol + Parietal Act (AUC = ", round(auc(roc_3var_1), 3), ")"),
         paste0("Glu Prec + Glu ACC + Parietal Act (AUC = ", round(auc(roc_3var_2), 3), ")"),
         paste0("Hip Vol + Hip Act + Parietal Act (AUC = ", round(auc(roc_3var_3), 3), ")")
       ),
       col = c("#1f77b4", "#ff7f0e", "#2ca02c"), 
       lwd = 3, 
       cex = 0.75, 
       bty = "n") # Removes the box around the legend









# Install these if you haven't already: 
# install.packages(c("survival", "survminer", "dplyr"))

library(survival)
library(survminer)
library(dplyr)

# 1. Categorize the continuous precuneus glutamate variable into tertiles
# This creates the distinct "Low", "Medium", and "High" groups for the plot
MRS_long <- MRS_long %>%
  mutate(
    glu_tertile = ntile(t1_glu_prec, 3),
    glu_strata = case_when(
      glu_tertile == 1 ~ "Low Prec Glu",
      glu_tertile == 2 ~ "Medium Prec Glu",
      glu_tertile == 3 ~ "High Prec Glu"
    ),
    # Convert to a factor to ensure the legend displays in the correct order
    glu_strata = factor(glu_strata, levels = c("Low Prec Glu", "Medium Prec Glu", "High Prec Glu"))
  )

# 2. Fit the survival curve specifically for plotting
# Note: We use survfit here to generate the survival probabilities for the plot
fit_plot <- survfit(Surv(age_change_moca, decliners_numeric) ~ glu_strata, data = MRS_long)

# 3. Generate the graph using ggsurvplot
ggsurvplot(
  fit_plot,
  data = MRS_long,
  conf.int = TRUE,          # Adds the shaded confidence intervals
  conf.int.alpha = 0.3,     # Adjusts the transparency of the shading
  censor = TRUE,            # Adds the '+' censor marks for stability
  censor.shape = 43,        # Sets censor shape to the standard '+'
  censor.size = 4.5,
  # Approximating the orange, grey, and blue color palette from your image
  palette = c("#FDB863", "#999999", "#56B4E9"), 
  title = "",
  xlab = "Years",
  ylab = "Probability of Stability",
  legend.title = "",
  legend = "top",
  ggtheme = theme_minimal() # Provides the clean, white grid background
)














library(survival)
library(survminer)
library(dplyr)

# 1. Categorize the continuous ACC glutamate variable into tertiles
MRS_long <- MRS_long %>%
  mutate(
    glu_tertile = ntile(t1_glu_acc, 3),
    glu_strata = case_when(
      glu_tertile == 1 ~ "Low ACC Glu",
      glu_tertile == 2 ~ "Medium ACC Glu",
      glu_tertile == 3 ~ "High ACC Glu"
    ),
    # CORRECTED: The levels now perfectly match the strings from case_when above
    glu_strata = factor(glu_strata, levels = c("Low ACC Glu", "Medium ACC Glu", "High ACC Glu"))
  )

# 2. Fit the survival curve specifically for plotting
fit_plot <- survfit(Surv(age_change_moca, decliners_numeric) ~ glu_strata, data = MRS_long)

# 3. Generate the graph using ggsurvplot
ggsurvplot(
  fit_plot,
  data = MRS_long,
  conf.int = TRUE,          
  conf.int.alpha = 0.3,     
  censor = TRUE,            
  censor.shape = 43,        
  censor.size = 4.5,
  palette = c("#FDB863", "#999999", "#56B4E9"), 
  # UPDATED TITLE: Reflects the shift to the Anterior Cingulate
  title = "",
  xlab = "Years",
  ylab = "Probability of Stability",
  legend.title = "Strata",
  legend = "top",
  ggtheme = theme_minimal() 
)