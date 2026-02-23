## Qoala-T: Estimations of MRI Qoala-T using BrainTime model

# Code to reproduce step 4 of our Qoala-T Tool
# Copyright (C) 2017-2019 Lara Wierenga - Leiden University, Brain and Development Research Center
# 
# This package contains data and R code for use of the Qoala-T tool based on the BrainTime model:
#   
#   title: Qoala-T: A supervised-learning tool for quality control of FreeSurfer segmented MRI data
# author:
#  - name:   Klapwijk, E.T., van de Kamp, F., Meulen, M., Peters, S. and Wierenga, L.M.
# https://doi.org/10.1016/j.neuroimage.2019.01.014
#
# If you have any question or suggestion, dont hesitate to get in touch:
# https://github.com/Qoala-T/QC/issues

# 1. LOAD PACKAGES (DMwR successfully removed!)
packages <- c("caret", "corrplot", "gbm", "plyr", "randomForest", "e1071",
              "pROC", "dplyr", "pbkrtest", "car", "doParallel", "ROSE", "repmis")

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
lapply(packages, library, character.only = TRUE)

# -----------------------------------------------------------------
# 2. SET FOLDERS & DATASET NAME
# -----------------------------------------------------------------
inputFolder <- "/Volumes/TOSHIBA_EXT/02-Raw_data-anat/longitudinal_freesurfer_149/reconall_output/"
outputFolder <- paste0(inputFolder, "Output_Qoala_T/") 
ifelse(dir.exists(outputFolder), FALSE, dir.create(outputFolder))

setwd(inputFolder)
dataset_name <- "CrossSectional_Glu_FS7" 

# -----------------------------------------------------------------
# 3. LOAD YOUR DATASET
# -----------------------------------------------------------------
stats2Table <- read.csv(paste("FreeSurfer_Output_", dataset_name, ".csv", sep=""), header=T, row.names=1)
test_data <- stats2Table

# -----------------------------------------------------------------
# 4. LOAD QOALA-T MODEL FROM GITHUB 
# -----------------------------------------------------------------
githubURL <- "https://github.com/Qoala-T/QC/blob/master/Qoala_T_model.Rdata?raw=true"
rf.tune <- get(load(url(githubURL)))

# -----------------------------------------------------------------
# 5. MATCH COLUMNS & PREDICT
# -----------------------------------------------------------------
dataset_colnames <- names(rf.tune$trainingData)[-ncol(rf.tune$trainingData)]
testing <- test_data[,dataset_colnames]
testing <- testing[complete.cases(testing),]

rf.pred <-  predict(rf.tune,testing)
rf.probs <- predict(rf.tune,testing,type="prob") 

# -----------------------------------------------------------------
# 6. SAVE OUTPUT CSV
# -----------------------------------------------------------------
Qoala_T_predictions <- data.frame(matrix(ncol = 4, nrow = nrow(rf.probs)))                              
colnames(Qoala_T_predictions) = c('participant_id','Scan_QoalaT', 'Recommendation', 'manual_QC_adviced') 

Qoala_T_predictions$participant_id <- row.names(rf.probs)
Qoala_T_predictions$Scan_QoalaT <- rf.probs$Include*100 
Qoala_T_predictions$Recommendation <- rf.pred
Qoala_T_predictions$manual_QC_adviced <- ifelse(Qoala_T_predictions$Scan_QoalaT<70&Qoala_T_predictions$Scan_QoalaT>30,"yes","no")
Qoala_T_predictions <- Qoala_T_predictions[order(Qoala_T_predictions$Scan_QoalaT, Qoala_T_predictions$participant_id),]

csv_Qoala_T_predictions = paste0(outputFolder, 'Qoala_T_predictions_model_based_', dataset_name, '.csv')
write.csv(Qoala_T_predictions, file = csv_Qoala_T_predictions, row.names=F)
print("CSV Saved!")

# -----------------------------------------------------------------
# 7. PLOT RESULTS 
# -----------------------------------------------------------------
excl_rate <- table(Qoala_T_predictions$Recommendation)
fill_colour <- rev(c("#88A825","#CF4A30"))
font_size <- 12
text_col <- "Black"

p <- ggplot(Qoala_T_predictions, aes(x=Scan_QoalaT,y=1,col=Recommendation)) +  
  annotate("rect", xmin=30, xmax=70, ymin=1.12, ymax=.88, alpha=0.2, fill="#777777") +
  geom_jitter(alpha=.8,height=.1,size=6) +
  ggtitle(paste("Qoala-T estimation model based ",dataset_name,"\nMean Qoala-T Score = ",round(mean(Qoala_T_predictions$Scan_QoalaT),1),sep="")) + 
  annotate("text", x=20, y=1.15, label=paste("Excluded = ",as.character(round(excl_rate[1]))," scans",sep="")) + 
  annotate("text", x=80, y=1.15, label=paste("Included = ",as.character(round(excl_rate[2]))," scans",sep="")) + 
  scale_colour_manual(values=fill_colour) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.text.x = element_text (size = font_size,color=text_col), axis.text.y = element_blank(),
        axis.title.x = element_text (size = font_size,color=text_col), axis.title.y = element_blank(), 
        axis.ticks=element_blank(), plot.title=element_text (size =16,color=text_col,hjust=.5))
print(p) 

filename<- paste0(outputFolder, "Figure_Rating_model_based_", dataset_name, ".pdf")
dev.copy(pdf,filename,width=30/2.54, height=20/2.54)
dev.off()
print("Plot Saved! Quality Control Complete.")