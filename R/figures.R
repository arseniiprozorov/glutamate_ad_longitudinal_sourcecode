########## Figures ############
library(vcd)

# 1. Prepare data (Percentages of Decliners within each Trajectory)
acc_prop <- as.data.frame(prop.table(table(MRS_long$traj_glu_acc, MRS_long$decliners), margin = 1) * 100)
colnames(acc_prop) <- c("Trajectory", "Outcome", "Percentage")

# 2. Plot
library(ggplot2)
ggplot(acc_prop, aes(x = Trajectory, y = Percentage, fill = Outcome)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("0" = "#999999", "1" = "#d7191c"), 
                    labels = c("Stable", "Declining")) +
  theme_minimal() +
  labs(title = "Clinical Outcome by ACC Glutamate Trajectory",
       subtitle = "p = 0.0069 (Chi-Square)",
       y = "Percentage of Group (%)",
       x = "ACC Glutamate Trajectory (2 Years)") +
  geom_text(aes(label = paste0(round(Percentage), "%")), 
            position = position_dodge(width = 0.9), vjust = -0.5)








#### Figure 2 




library(tidyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)

# 1. Reshape ACC data and annualize the time axis
# We filter age_difference here if you want to strictly exclude data past 3 years, 
# or just cap the plot later to keep the model power.
MRS_acc_long <- MRS_long %>%
  mutate(diagnostic_nick = factor(diagnostic_nick, levels = c("HC", "SCD+", "MCI"))) %>%
  select(pscid, diagnostic_nick, sexe, education, 
         t1_age, t2_age, age_difference,
         t1_glu_acc, t2_glu_acc) %>%
  pivot_longer(
    cols = c(t1_glu_acc, t2_glu_acc),
    names_to = "visit",
    values_to = "glu_acc"
  ) %>%
  mutate(years_elapsed = ifelse(visit == "t1_glu_acc", 0, age_difference))

# 2. Run the Annualized LMM for ACC
# Replicating Model 1: neurotransmitter ~ diagnostic group x years
model_acc_annual <- lmer(glu_acc ~ diagnostic_nick * years_elapsed + (1 | pscid), 
                         data = MRS_acc_long)

# 3. Extract annualized slopes (Beta = change in mmol/L per year)
acc_stats <- emtrends(model_acc_annual, ~ diagnostic_nick, var = "years_elapsed")
acc_summary <- summary(acc_stats, infer = TRUE)

# 4. Create labels for the plot
acc_labels <- acc_summary %>%
  mutate(label = paste0("beta == ", round(years_elapsed.trend, 2), 
                        "~~p == ", round(p.value, 3)))

# 5. Create the faceted trajectory plot capped at 3 years
ggplot(MRS_acc_long, aes(x = years_elapsed, y = glu_acc, color = diagnostic_nick)) +
  geom_line(aes(group = pscid), alpha = 0.2, linewidth = 0.4) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", linewidth = 1.5, se = TRUE, aes(fill = diagnostic_nick)) +
  facet_wrap(~ diagnostic_nick) +
  
  # Add Stats (Annualized Beta and P-value)
  geom_text(data = acc_labels, aes(x = 1.5, y = Inf, label = label), 
            vjust = 2, color = "black", parse = TRUE, inherit.aes = FALSE) +
  
  # Cap the X-axis at 3 years for visual clarity
  scale_x_continuous(limits = c(0, 3), breaks = c(0, 1, 2, 3)) +
  
  scale_color_manual(values = c("HC" = "#3B5998", "SCD+" = "#E69F00", "MCI" = "#D55E00")) +
  scale_fill_manual(values = c("HC" = "#3B5998", "SCD+" = "#E69F00", "MCI" = "#D55E00")) +
  labs(
    title = "Annualized Glutamate Trajectories: ACC",
    subtitle = "Faceted by clinical stage; Slopes (Beta) represent change per year",
    x = "Years from Baseline ",
    y = "ACC Glutamate (mM)"
  ) +
  theme_minimal() +
  theme(legend.position = "none", panel.grid.minor = element_blank())





library(tidyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)

# 1. Reshape Precuneus data and annualize the time axis
MRS_prec_long <- MRS_long %>%
  mutate(diagnostic_nick = factor(diagnostic_nick, levels = c("HC", "SCD+", "MCI"))) %>%
  select(pscid, diagnostic_nick, sexe, education, 
         t1_age, t2_age, age_difference,
         t1_glu_prec, t2_glu_prec) %>%
  pivot_longer(
    cols = c(t1_glu_prec, t2_glu_prec),
    names_to = "visit",
    values_to = "glu_prec"
  ) %>%
  mutate(years_elapsed = ifelse(visit == "t1_glu_prec", 0, age_difference))

# 2. Run the Annualized LMM for Precuneus
# Replicating Model 1 from the paper: neurotransmitter ~ diagnosis * years
model_prec_annual <- lmer(glu_prec ~ diagnostic_nick * years_elapsed + (1 | pscid), 
                          data = MRS_prec_long)

# 3. Extract annualized slopes (Beta = change in mmol/L per year)
prec_stats <- emtrends(model_prec_annual, ~ diagnostic_nick, var = "years_elapsed")
prec_summary <- summary(prec_stats, infer = TRUE)

# 4. Create labels for the plot
prec_labels <- prec_summary %>%
  mutate(label = paste0("beta == ", round(years_elapsed.trend, 2), 
                        "~~p == ", round(p.value, 3)))

# 5. Create the faceted trajectory plot capped at 3 years
ggplot(MRS_prec_long, aes(x = years_elapsed, y = glu_prec, color = diagnostic_nick)) +
  geom_line(aes(group = pscid), alpha = 0.2, linewidth = 0.4) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", linewidth = 1.5, se = TRUE, aes(fill = diagnostic_nick)) +
  facet_wrap(~ diagnostic_nick) +
  
  # Add Stats (Annualized Beta and P-value)
  geom_text(data = prec_labels, aes(x = 1.5, y = Inf, label = label), 
            vjust = 2, color = "black", parse = TRUE, inherit.aes = FALSE) +
  
  # Cap the X-axis at 3 years for visual clarity
  scale_x_continuous(limits = c(0, 3), breaks = c(0, 1, 2, 3)) +
  
  scale_color_manual(values = c("HC" = "#3B5998", "SCD+" = "#E69F00", "MCI" = "#D55E00")) +
  scale_fill_manual(values = c("HC" = "#3B5998", "SCD+" = "#E69F00", "MCI" = "#D55E00")) +
  labs(
    title = "Annualized Glutamate Trajectories: Precuneus",
    subtitle = "Faceted by clinical stage; Slopes (Beta) represent change per year",
    x = "Years from Baseline",
    y = "Precuneus Glutamate (mM)"
  ) +
  theme_minimal() 
  
  
  
  






  
  
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