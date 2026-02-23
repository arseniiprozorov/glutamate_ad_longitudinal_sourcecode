##################
# STATS2TABLE.R #
##################
# Written by Olga Veth - s2067579 - University of Leiden
# Created on 30-09-2019
# Edited by Lara M Wierenga on 22-02-21
# Most Recent update: 22-02-21
# Version 4.0

# Edited for FS7
##################
# MERGED WITH QOALA-T AND ROBUST BINDING FIXES

library(dplyr)

datasetDir <- "/Volumes/TOSHIBA_EXT/02-Raw_data-anat/longitudinal_freesurfer_149/reconall_output/" # Note the trailing slash
setwd(datasetDir) 
dataset_name <- "CrossSectional_Glu_FS7" 

readAseg <- function(){
  aseg_file <- data.frame(read.table(paste("./stats/aseg.stats", sep=""), row.names=1))[,c(3,4)] 
  asegTable <- t(data.frame(aseg_file[,1], row.names = aseg_file[,2])) 
  return (asegTable)
}

readMetaAseg <- function(){
  aseg_meta <- readLines("./stats/aseg.stats", n=35)[15:34] # FS7 update
  meta1 <- gsub("# ", "", aseg_meta)
  meta <- t(data.frame(strsplit(meta1, ",")))[,c(2,4)]
  metaTable <- t(data.frame(meta[,2]))
  colnames(metaTable) <- meta[,1]
  return(metaTable)
}

readAparc <- function(value){
  sides <- c("lh", "rh")
  ifelse((value == "area"), pos <- 1, pos <- 2)
  for (x in 1:length(sides)){
    areaThickness <- as.data.frame(read.table(paste("./stats/", sides[x], ".aparc.stats", sep=""), row.names=1))[, c(2,4)]
    rowValues <- rownames(areaThickness)
    meta <- readLines(paste("./stats/", sides[x], ".aparc.stats", sep=""))[c(20, 21)] 
    meta1 <- gsub("# ", "", meta)
    meta2 <- t(data.frame(strsplit(meta1, ",")))[, c(2,4)]
    meta3 <- data.frame(meta2[pos,2])
    value2 <- gsub(" ", "", meta2[pos,1])
    colnames(meta3) <- paste(sides[x], "_", value2, "_" , value, sep="")
    extra <- t(matrix(areaThickness[,pos]))
    colnames(extra) <- paste(sides[x], "_", rowValues, "_", value, sep="")
    ifelse(x==1, aparcTable <- cbind(extra, meta3), aparcTable <- cbind(aparcTable, extra, meta3))
  }
  return(aparcTable)
}

readFiles <- function(){
  asegTable <- readAseg()
  metaTable <- readMetaAseg()
  areaAparc <- readAparc("area")
  thickAparc <- readAparc("thickness")
  subjectTable <- cbind(asegTable, metaTable, areaAparc, thickAparc) 
  subjectTable <- data.frame(subjectTable)
  return (subjectTable)
}

preprocTable <- function(subjectTable){
  # FIX: Kept SurfaceHoles for Qoala-T
  removeCols <- c("*.WM-hypointensities$","*.WM.hypointensities$", "*pole*", "*bankssts*", 
                  "VentricleChoroidVol", "*CerebralWhiteMatterVol", 
                  "SegVolFile.mri.aseg.mgz.", "*CorticalWhiteMatterVol")
  remove <- grep(paste(removeCols, collapse="|"), colnames(subjectTable))
  if(length(remove) > 0) subjectTable <- subjectTable[, -remove]
  
  colnames(subjectTable) <- gsub("^X\\.", "", colnames(subjectTable))
  colnames(subjectTable) <- gsub("_\\.", "_", colnames(subjectTable))
  colnames(subjectTable) <- gsub("-", ".", colnames(subjectTable))
  colnames(subjectTable) <- gsub(" ", "", colnames(subjectTable))
  
  colnames(subjectTable)[which(colnames(subjectTable) == "eTIV")] <- "EstimatedTotalIntraCranialVol"
  colnames(subjectTable)[which(colnames(subjectTable) %in% c("rd.Ventricle", "th.Ventricle", "5th.Ventricle"))] <- c("X4th.Ventricle", "X3rd.Ventricle", "X5th.Ventricle") 
  
  return(subjectTable)
}

main <- function(){
  subjects <- c()
  first <- T
  subjectDirs <- unique(list.dirs('.', recursive=FALSE)) 
  
  # Ensure we only process cross-sectional (ignore anything with .long or .base)
  subjectDirs <- subjectDirs[!grepl("\\.long\\.", subjectDirs) & !grepl("\\.base", subjectDirs, ignore.case=TRUE)]
  
  # FIX: Loop starts at 1, not 0
  for (x in 1:length(subjectDirs)){
    setwd(paste0(datasetDir, substring(subjectDirs[x], 3), "/"))
    
    if (file.exists("./stats/aseg.stats")){
      subjectTable <- readFiles()
      subjectTable <- preprocTable(subjectTable)
      
      if (first == T){
        stats2Table <- subjectTable
        subjects <- c(subjects, substring(subjectDirs[x], 3))
        first <- F
      } else {
        # FIX: Robust row binding
        stats2Table <- dplyr::bind_rows(stats2Table, subjectTable)
        subjects <- c(subjects, substring(subjectDirs[x], 3))
      } 
    }
  }
  
  stats2Table <- data.frame(stats2Table)
  rownames(stats2Table) <- subjects
  
  # FS7 Updates
  colnames(stats2Table)[which(names(stats2Table) == "Left.Thalamus")] <- "Left.Thalamus.Proper"
  colnames(stats2Table)[which(names(stats2Table) == "Right.Thalamus")] <- "Right.Thalamus.Proper"
  
  # These variables might be missing entirely if it was a missing column, so we check if they exist first
  if("BrainSegVolNotVent" %in% colnames(stats2Table)) {
    stats2Table$BrainSegVolNotVentSurf <- stats2Table$BrainSegVolNotVent
  }
  if("SupraTentorialVolNotVent" %in% colnames(stats2Table)) {
    stats2Table$SupraTentorialVolNotVentVox <- stats2Table$SupraTentorialVolNotVent
  }
  
  setwd(datasetDir)
  write.csv(stats2Table, paste0("FreeSurfer_Output_", dataset_name,".csv"))
  print("Extraction complete!")
}

main()