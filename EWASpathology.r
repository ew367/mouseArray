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

# within cell type models

# null model to check gorup term improves the model

# excluded age (Jon's advice)



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
  
  nullLM<-lm(row ~ QCmetrics$Sex + QCmetrics$Batch)
  
  fullLM<-lm(row ~ QCmetrics$Pathology + QCmetrics$Sex + QCmetrics$Batch)
  
  
  # extract case control main effect and sex effect
  return(c(summary(fullLM)$coefficients["QCmetrics$Pathology",c(1,2,4)],
           summary(fullLM)$coefficients["QCmetrics$SexM",c(1,2,4)],
           
           # extract sex effect
           summary(nullLM)$coefficients["QCmetrics$SexM",c(1,2,4)],
           
           # run anova
           anova(fullLM, nullLM)[2,"Pr(>F)"]))
  
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

outtab<-matrix(data = parRapply(cl, celltypeNormbeta, runEWAS, QCmetrics), ncol = 10, byrow = TRUE)


rownames(outtab)<-rownames(celltypeNormbeta)
#rownames(outtab)<-rownames(betasSub)
colnames(outtab)<-c("GroupWT_coeff", "GroupWT_SE", "GroupWT_P", 
                    "SexM_coeff", "SexM_SE", "SexM_P",
                
                    "nullSexM_coeff", "nullSexM_SE", "nullSexM_P",
                    
                    "anovoP") 

filePath <- paste0("3_analysis/results/", cellType, "EWASpathOut.rdat")
save(outtab, file = filePath)



#----------------------------------------------------------------------#
# Explore results
#----------------------------------------------------------------------#


bonfP <- 0.05/nrow(outtab)

# for each col, how many resultsa are lower than bonfP

countSig <- function(x){
  sig <- sum(x < bonfP)
  return(sig)
}

apply(outtab[,1:10], 2, countSig) 
