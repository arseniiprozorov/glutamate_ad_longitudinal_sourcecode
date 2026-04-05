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

jmv::descriptives(data = MRS_long, vars = vars(sexe,   diagnostic_nick,education,  t1_age,   t2_age, 
                                               age_difference, t1_glu_prec, t2_glu_prec, t1_hipp_e_tiv, hipp_difference,
                                               percent_change_glu_prec,t1_glu_acc, t2_glu_acc, percent_change_glu_acc, 
                                               slope_moca_raw, intercept_moca_raw, cluster_moca_raw, change_over_3_5,decliners,  
                                               memor_rappel_libre_nombre_reponses_correctes_t1,memor_rappel_libre_nombre_reponses_correctes_t2,
                                               visage_visages_score_rappel_differe_9_t1, visage_visages_score_rappel_differe_9_t2),
                  , splitBy = "decliners", sd = TRUE, iqr = TRUE)






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

############################## ANALYSES PRINCIPALES ###########################
names(MRS_long)
sapply(MRS_long,class)






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




############## Decliners and moca slope ###################




################## Charctarize hipp ############################
sapply(MRS_long,class)
jmv::descriptives(data = MRS_long, vars = vars( t1_hipp_e_tiv,t2_hipp_e_tiv ,hipp_difference), 
                  sd = TRUE, iqr = TRUE, skew = TRUE, kurt = TRUE, hist = TRUE)


## test with glutamate trajectories
t.test(hipp_difference ~ traj_glu_prec, data = MRS_long,
       alternative = c("two.sided"),
       mu = 0, var.equal = FALSE,
       conf.level = 0.95)

t.test(hipp_difference ~ traj_glu_acc, data = MRS_long,
       alternative = c("two.sided"),
       mu = 0, var.equal = FALSE,
       conf.level = 0.95)

t.test(t1_hipp_e_tiv ~ traj_glu_prec, data = MRS_long,
       alternative = c("two.sided"),
       mu = 0, var.equal = FALSE,
       conf.level = 0.95)

t.test(t1_hipp_e_tiv ~ traj_glu_acc, data = MRS_long,
       alternative = c("two.sided"),
       mu = 0, var.equal = FALSE,
       conf.level = 0.95)

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


## test correlations ##
library(psych)
library(corrplot)
help('cor.test')
cor_vars <- c("t1_glu_prec", "t2_glu_prec", "percent_change_glu_prec", 
  "t1_glu_acc", "t2_glu_acc", "percent_change_glu_acc",
  "slope_moca_raw", "intercept_moca_raw", 
  "t1_hipp_e_tiv", "t2_hipp_e_tiv", "hipp_difference", 
  "education", "t1_age")

cor_results <- corr.test(MRS_long[, cor_vars], use = "pairwise", 
                         method = "pearson",  adjust = "none")
r_matrix <- cor_results$r
p_matrix <- cor_results$p
r_matrix
p_matrix

corrplot(r_matrix, method = "color",type = "upper",order = "hclust", 
         hclust.method = "ward.D2", p.mat = p_matrix, sig.level = 0.05,
         insig = "pch",     addCoef.col = "red",
         number.cex = 0.7)














