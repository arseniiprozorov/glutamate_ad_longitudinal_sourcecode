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

