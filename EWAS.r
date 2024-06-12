##---------------------------------------------------------------------#
##
## Title: Run EWAS
##
## Purpose of script: Run linear regression within each cell type
##                    
##                    This scripted is adapted from one co-authored
##                    by EJH and EMW for the MRC schizophrenia project
##
##                   
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# to run: sbatch..... 

# NeuN+ models

# pathology: lm(DNAm ~ QCmetrics$Group + QCmetrics$Pathology + QCmetrics$Age + QCmetrics$Sex + QCmetrics$Batch


# group*age: DNAm ~ QCmetrics$Group + QCmetrics$Age*QCmetrics$Group + QCmetrics$Age + QCmetrics$Sex + QCmetrics$Batch



# NeuN- models- need to include individual as random effect

#lmer(DNAm ~ QCmetrics$Group + QCmetrics$Pathology + QCmetrics$Age + QCmetrics$Sex + QCmetrics$Batch (1|pheno$Individual_ID), REML = FALSE)


#lmer(DNAm ~ QCmetrics$Group + QCmetrics$Age*QCmetrics$Group + QCmetrics$Age + QCmetrics$Sex + QCmetrics$Batch (1|pheno$Individual_ID), REML = FALSE)


# include null models of each with anova to test that pathology/interaction term is improving the model fit



#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(lme4)
library(lmerTest)
library(dplyr)


#----------------------------------------------------------------------#
# DEFINE ANALYSIS FUNCTION
#----------------------------------------------------------------------#

# if cell type = neun+ source functions without random terms
# if cell type = neun- source functions w/ random terms

runEWAS<-function(row,QCmetrics){
  
  nullLM<-lm(row ~ QCmetrics$Group + QCmetrics$Age + QCmetrics$Sex + QCmetrics$Batch)
  
  interactionLM<-lm(row ~ QCmetrics$Group + QCmetrics$Age*QCmetrics$Group + QCmetrics$Age + QCmetrics$Sex + QCmetrics$Batch)
  
  
  # extract case control main effect and cell proportion effect
  return(c(summary(interactionLM)$coefficients["QCmetrics$GroupWT",c(1,2,4)],
           summary(interactionLM)$coefficients["QCmetrics$Age",c(1,2,4)],
           summary(interactionLM)$coefficients["QCmetrics$SexM",c(1,2,4)],
           summary(interactionLM)$coefficients["QCmetrics$GroupWT:QCmetrics$Age",c(1,2,4)],
           
           
           # extract cell specific case control effect
           summary(nullLM)$coefficients["QCmetrics$GroupWT",c(1,2,4)],
           summary(nullLM)$coefficients["QCmetrics$Age",c(1,2,4)],
           summary(nullLM)$coefficients["QCmetrics$SexM",c(1,2,4)]))
  
}



#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(doParallel)

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

args<-commandArgs(trailingOnly = TRUE)
projDir <- args[1]
#projDir <- "/lustre/projects/Research_Project-191406/cellSortedEWAS"
cellType <- args[2]
#cellType <- "NEUNpos"


normData<-file.path(projDir, "2_normalised/normalisedData.rdat")

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(projDir)
load(normData)


QCmetrics$Group <- as.factor(QCmetrics$Group)
QCmetrics$Sex <- as.factor(QCmetrics$Sex)
QCmetrics$Batch <- as.factor(QCmetrics$Batch)

print(paste0("running EWAS on ", cellType, " cell type..."))
## subset to cell type samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell_Type == cellType),]


# subset beta matrix to cell type specific samples
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]

# take top 100 rows for debugging
#betasSub <- celltypeNormbeta[1:100,]
#row <- colMeans(celltypeNormbeta)

#----------------------------------------------------------------------#
# INTITATE PARALLEL ENV
#----------------------------------------------------------------------#

nCores<-detectCores()
cl <- makeCluster(nCores-1)
registerDoParallel(cl)
clusterExport(cl, list("runEWAS"))

outtab<-matrix(data = parRapply(cl, celltypeNormbeta, runEWAS, QCmetrics), ncol = 9, byrow = TRUE)
#outtab<-matrix(data = parRapply(cl, betasSub, runEWAS, QCmetrics), ncol = 21, byrow = TRUE)


rownames(outtab)<-rownames(celltypeNormbeta)
#rownames(outtab)<-rownames(betasSub)
colnames(outtab)<-c("GroupWT_coeff", "GroupWT_SE", "GroupWT_P", 
                    "Age_coeff", "Age_SE", "Age_P",
                    "SexM_coeff", "SexM_SE", "SexM_P",
                    "GroupWT:QCmetrics$Age_coeff", "GroupWT:QCmetrics$SE", "GroupWT:QCmetrics$Age_P",
                    
                    "nullGroupWT_coeff", "nullGroupWT_SE", "nullGroupWT_P", 
                    "nullAge_coeff", "nullAge_SE", "nullAge_P",
                    "nullSexM_coeff", "nullSexM_SE", "nullSexM_P") 


save(outtab, file = "3_analysis/results/EWASout.rdat")
