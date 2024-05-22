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
  
  modelLM<-lm(row ~ QCmetrics$Phenotype + QCmetrics$Cell.Proportions + QCmetrics$CCDNAmAge + QCmetrics$Sex + QCmetrics$Tissue.Centre)
  nullCT<-lm(row ~ QCmetrics$Phenotype +  QCmetrics$CCDNAmAge + QCmetrics$Sex + QCmetrics$Tissue.Centre)
  
  # extract case control main effect and cell proportion effect
  return(c(summary(modelLM)$coefficients["QCmetrics$PhenotypeSchizophrenia",c(1,2,4)],
           summary(modelLM)$coefficients["QCmetrics$Cell.Proportions",c(1,2,4)],
           
           # extract cell specific case control effect
           summary(nullCT)$coefficients["QCmetrics$PhenotypeSchizophrenia",c(1,2,4)]))
}



#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(doParallel)

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

args<-commandArgs(trailingOnly = TRUE)
#dataDir <- args[1]
#cellType <- args[2]
cellType <- "NEUNpos"


normData<-file.path(dataDir, "2_normalised/normalisedData.rdat")

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(dataDir)
load(normData)
elisa <- read.csv("0_metadata/AB42ELISAconcentrations.csv", stringsAsFactors = F) # pathology data


# add in pathology data
elisa$Individual_ID <- as.character(elisa$Sample_ID)
elisa$Pathology <- elisa$Mean_N3

QCmetrics <- left_join(QCmetrics, elisa %>% dplyr::select(Individual_ID, Pathology))
QCmetrics$Group <- as.factor(QCmetrics$Group)

print(paste0("running EWAS on ", cellType, " cell type..."))
## subset to cell type samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell_Type == cellType),]


# subset beta matrix to cell type specific samples
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]

# take top 100 rows for debugging
betasSub <- celltypeNormbeta[1:100,]
meanBetas <- colMeans(celltypeNormbeta)

#----------------------------------------------------------------------#
# INTITATE PARALLEL ENV
#----------------------------------------------------------------------#

nCores<-detectCores()
cl <- makeCluster(nCores-1)
registerDoParallel(cl)
clusterExport(cl, list("runEWAS"))

outtab<-matrix(data = parRapply(cl, celltypeNormbeta, runEWAS, QCmetrics), ncol = 9, byrow = TRUE)
#outtab<-matrix(data = parRapply(cl, betasSub, runEWAS, QCmetrics), ncol = 9, byrow = TRUE)


rownames(outtab)<-rownames(celltypeNormbeta)
#rownames(outtab)<-rownames(betasSub)
colnames(outtab)<-c("SCZ_coeff", "SCZ_SE", "SCZ_P", 
                    paste0(cellType,"_coeff"), paste0(cellType,"_SE"), paste0(cellType,"_P"),
                    "nullCT_SCZ_coeff", "nullCT_SCZ_SE", "nullCT_SCZ_P") 


save(outtab, file = file.path(paste0(resPath, cellType,"LM.rdata")))
