##---------------------------------------------------------------------#
##
## Title: Normalise the data for downstream analyses
##
## Purpose of script: Filter the data using the QC metrics and cell type
##                    checks
##
##                    Remove flagged probes
##
##                    Normalise the data within celltype and save
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#


#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

args<-commandArgs(trailingOnly = TRUE)
dataDir <- args[1]

setwd(dataDir)

print("loading packages...")
library(data.table)
library(wateRmelon)
library(ENmix)


source("config.r")

#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#

# get raw betas for all samples
load(file = file.path(QCDir, "mraw.rdat"))
rawbetas <- getB(mraw)

# load detP for removing failed probes
load(file = file.path(QCDir, "detP.rdat"))

# load manifest for removing sex probes and mfg_flagged probes
man <- fread(manifest, skip=7, fill=TRUE, data.table=F)

# load QCmetrics and QC summary to remove failed samples
if(ctCheck){
  QCSum <- read.csv(file.path(QCDir, "passQCStatusStage3AllSamples.csv"), stringsAsFactors = F)
  QCSum <- na.omit(QCSum)
  passQC <- QCSum$Basename[QCSum$passQCS3]
  load(file = file.path(QCDir, "QCmetricsPostCellTypeChecks.rdat"))
} else {
  QCSum <- read.csv(file.path(QCDir, "passQCStatusStage1AllSamples.csv"), stringsAsFactors = F)
  QCSum <- na.omit(QCSum)
  passQC <- QCSum$Basename[QCSum$passQCS1]
  load(file = file.path(QCDir, "QCmetrics.rdat"))
}


QCmetrics <- QCmetrics[QCmetrics$Basename %in% passQC,]
     

#----------------------------------------------------------------------#
# REMOVE FAILED PROBES AND SAMPLES
#----------------------------------------------------------------------#

print("filtering SNPs...")
betas<-rawbetas[-grep("rs", rownames(rawbetas)),]

print("filtering flagged probes...")
flagged.probes<-man$IlmnID[man$MFG_Change_Flagged == TRUE]
betas <- betas[!row.names(betas) %in% flagged.probes,]

print("filtering sex and MT probes...")
#auto.probes<-man$IlmnID[man$CHR != "X" & man$CHR != "Y" & man$CHR != "MT"]
#betas<-betas[row.names(betas) %in% auto.probes,]

print("filtering detP failed probes...")
failedProbes <- rownames(detP)[((rowSums(detP > pFiltProbeThresh)/ncol(detP)) * 100) > pFiltSampleThresh]
betas<-betas[!row.names(betas) %in% failedProbes,]


print("filtering samples failed any stage of QC...")
betas <- betas[,passQC]

# filter mraw object to retain intens info for normalisation
mrawPass <- mraw[row.names(betas), colnames(betas)]



#----------------------------------------------------------------------#
# NORMALISE (WITHIN CELL TYPE FOR CELL SORTED DATA)
#----------------------------------------------------------------------#

#### This needs to be checked still!!
## also make sure that cols/rows are in same order before saving

cellTypes<-unique(QCmetrics$Cell_Type)

if(ctCheck){
  print("normalising within cell type...")
  
  celltypeNormbeta<-matrix(NA, nrow = nrow(assays(mrawPass)$Meth), ncol = ncol(assays(mrawPass)$Meth))
  rownames(celltypeNormbeta)<-rownames(betas)
  colnames(celltypeNormbeta)<-colnames(betas)
  for(each in cellTypes){
    print(each)
    index<-which(QCmetrics$Cell_Type == each)
    if(length(index) > 2){
      celltypeNormbeta[,index]<-as.matrix(adjustedDasen(mns = assays(mrawPass)$Meth[,index], uns = assays(mrawPass)$Unmeth[,index], onetwo = mrawPass@elementMetadata$Infinium_Design_Type, chr = mrawPass@elementMetadata$chr))
    }
  }
  
  save(celltypeNormbeta, QCmetrics,  file = file.path(normDir, "normalisedData.rdat"))
  print("normalised QC object saved")
} else{
  print("normalising bulk tissue...")
  normBeta <- adjustedDasen(mns = assays(mrawPass)$Meth, uns = assays(mrawPass)$Unmeth, onetwo = mrawPass@elementMetadata$Infinium_Design_Type, chr = mrawPass@elementMetadata$chr, cores=1)
  save(normbeta, QCmetrics, file = file.path(normDir, "normalisedData.rdat"))
}



