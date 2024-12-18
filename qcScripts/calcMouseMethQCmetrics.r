##---------------------------------------------------------------------#
##
## Title: Calculate QC metrics for mouse methylation array data
##
## Purpose of script: Calculate standard QC metrics to be used to
##                    filter samples prior to normalisation
##
##                    The QC.rmd script uses these metrics to output
##                    a quality control html report
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# parameters and relative paths etc are loaded from the config.r 
# file in the dataect folder


#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

args<-commandArgs(trailingOnly = TRUE)
dataDir <- args[1]
refDir <- args[2]

print("loading packages...")

library(ENmix)
library(SummarizedExperiment)
library(dplyr)
library(stringr)
library(data.table)


# define directories

setwd(dataDir)

manifest <- file.path(refDir, "MouseMethylation-12v1-0_A2.csv")
pheno <- "0_metadata/sampleSheet.csv"
idatPath <- "1_raw"
normDir <- "2_normalised"
QCDir <- "2_normalised/QC"
configFile <- file.path(dataDir, "config.r") 


# load config.r
if (file.exists(configFile)) {
  source(configFile)
} else {
  stop("config.r file does not exist.")
}

# load sampleSheet
if (!file.exists(pheno)) {
  stop("Phenotype file not found at ", pheno)
}


# load manifest
man <- fread(manifest, skip=7, fill=TRUE, data.table=F)

# read in samplesheet
sampleSheet <- read.csv(pheno, stringsAsFactors = F)

# Exclude empty wells if defined in sampleSheet
if(exists("empty")){
  sampleSheet <- sampleSheet[!sampleSheet$Basename %in% empty,]
}

for(i in projVar){
  if(is.character(sampleSheet[,i])){
    sampleSheet[,i] <- trimws(sampleSheet[,i])
  }
}
#----------------------------------------------------------------------#
# LOAD IDATS TO RGSET
#----------------------------------------------------------------------#

if(file.exists(file = file.path(QCDir, "rgSet.rdat"))){
  print("Loading rgSet")
  load(file = file.path(QCDir, "rgSet.rdat"))
  rgSet <- rgSet[,sampleSheet$Basename]
  print(paste0(ncol(rgSet), " samples loaded"))
} else{
  rgSet <- readidat(path = idatPath, manifestfile=manifest ,recursive = TRUE)
  save(rgSet, file=file.path(QCDir, "rgSet.rdat"))
  print("rgSet created and saved")
}



#----------------------------------------------------------------------#
# Calculate intensities
#----------------------------------------------------------------------#

if(file.exists(file = file.path(QCDir, "mraw.rdat"))){
  print("Loading mraw object")
  load(file = file.path(QCDir, "mraw.rdat"))
} else{
  mraw <- getmeth(rgSet)
  save(mraw, file=file.path(QCDir, "mraw.rdat"))
  print("mraw object created and saved")
}


m_intensities <- assays(mraw)$Meth
u_intensities <- assays(mraw)$Unmeth

M.median <- apply(m_intensities, 2, median)
U.median <- apply(u_intensities, 2, median)

M.mean <- apply(m_intensities, 2, mean)
U.mean <- apply(u_intensities, 2, mean)


M <- as.data.frame(M.median)
M$M.mean <- M.mean
M$Basename <- rownames(M)

U <-as.data.frame(U.median)
U$U.mean <- U.mean
U$Basename <- rownames(U)


# make QC metrics object
QCmetrics <- left_join(sampleSheet, M, by = "Basename")
QCmetrics <- left_join(QCmetrics, U, by = "Basename")

QCmetrics$IntensityPass <- ifelse(QCmetrics$M.median > 2000 & QCmetrics$U.median > 2000, TRUE, FALSE)


#----------------------------------------------------------------------#
# BISULPHITE CONVERSION
#----------------------------------------------------------------------#

if(file.exists(file = file.path(QCDir, "bsCon.rdat"))){
  print("Loading bsCon object")
  load(file = file.path(QCDir, "bsCon.rdat"))
} else{

ctrls <- metadata(rgSet)$ictrl

# subset to just mouse specific probes
bs.type1 <- ctrls$Address[ctrls$Type == "BISULFITE CONVERSION I"][11:15]
bs.type2 <- ctrls$Address[ctrls$Type == "BISULFITE CONVERSION II"][4:6]

bs.green.type1 <- assays(rgSet)$Green[bs.type1,]
bs.green.type2 <- assays(rgSet)$Green[bs.type2,]
bs.red.type1 <- assays(rgSet)$Red[bs.type1,]
bs.red.type2 <- assays(rgSet)$Red[bs.type2,]

red.med <- apply(bs.red.type1, 1, median)
green.med <- apply(bs.green.type1, 1, median)

redGreater<-rowMeans(bs.red.type1) > rowMeans(bs.green.type1)

BScon1<-rbind(
  bs.red.type1[which(redGreater),] / ( bs.red.type1[which(redGreater),] + bs.green.type1[which(redGreater),] ), 
  bs.green.type1[which(!redGreater),] / ( bs.red.type1[which(!redGreater),] + bs.green.type1[which(!redGreater),] )
)

BScon2 <- bs.red.type2 / ( bs.red.type2 + bs.green.type2 )

BSconAll<-rbind(BScon1, BScon2)
BScon.med<-apply(BSconAll, 2, median)
BSconAll<-rbind(BSconAll, BScon.med)*100

BSconAll <- BSconAll[,QCmetrics$Basename]
BScon.med <- BScon.med[QCmetrics$Basename]*100

save(BSconAll, file=file.path(QCDir, "bsCon.rdat"))
print("bsCon object created and saved")

}


QCmetrics$BsCon <- BSconAll["BScon.med",]
QCmetrics$BsConPass <- ifelse(QCmetrics$BsCon > bsConThresh, TRUE, FALSE)
  


#----------------------------------------------------------------------#
# P FILTER
#----------------------------------------------------------------------#

if(file.exists(file = file.path(QCDir, "detP.rdat"))){
  print("Loading detP object")
  load(file = file.path(QCDir, "detP.rdat"))
} else{
  detP <- calcdetP(rgSet)
  save(detP, file=file.path(QCDir, "detP.rdat"))
  print("detP object created and saved")
}


# check if any samples have > 1 percent of probes with a detection p value of > pFiltThresh
pfiltdf <- data.frame(matrix(ncol = 2, nrow = nrow(QCmetrics)))
colnames(pfiltdf) <- c("Basename", "PercProbesFail")

for(i in 1:ncol(detP)){
  pfiltdf$Basename[i] <- colnames(detP)[i]
  pfiltdf$PercProbesFail[i] <- sum(detP[,i] > pFiltProbeThresh)/nrow(detP)*100
}

pfiltdf$PfiltPass <- ifelse(pfiltdf$PercProbes < 1, TRUE, FALSE)

QCmetrics <- left_join(QCmetrics, pfiltdf, by = "Basename")


# check if any probes fail in more than pFiltSampleThresh of samples
# samples with low intensity/bscon values or fail pfilt sample check above are excluded first

goodsamps <- QCmetrics$IntensityPass & QCmetrics$BsCon > bsConThresh & QCmetrics$PfiltPass
detP <- detP[,goodsamps]
failedProbes <- rownames(detP)[((rowSums(detP > pFiltProbeThresh)/ncol(detP)) * 100) > pFiltSampleThresh]
save(failedProbes, file = file.path(QCDir, "failedProbes.rdat"))


#----------------------------------------------------------------------#
# SEX CHECK (ENMIX)
#----------------------------------------------------------------------#

if(file.exists(file = file.path(QCDir, "sexPred.rdat"))){
  print("Loading sexPred object")
  load(file = file.path(QCDir, "sexPred.rdat"))
} else{
  sexPred <- predSex(rgSet)
  colnames(sexPred) <- c("Basename", "PredSex")
  save(sexPred, file=file.path(QCDir, "sexPred.rdat"))
  print("sexPred object created and saved")
}


QCmetrics <- left_join(QCmetrics, sexPred, by = "Basename")

if(sexCheck){
  QCmetrics$sexPass <- ifelse(QCmetrics$PredSex == QCmetrics$Sex, TRUE, FALSE)
}


#----------------------------------------------------------------------#
# PCA OF BETAS
#----------------------------------------------------------------------#


if(file.exists(file = file.path(QCDir, "PCAbetas.rdat"))){
  print("Loading PCA betas object")
  load(file = file.path(QCDir, "PCAbetas.rdat"))
  
} else {
# filter to autosomal probes and passed intens check samples only
  auto.probes<-man$IlmnID[man$CHR != "X" & man$CHR != "Y" & man$CHR != "MT"]
  rawbetas <- getB(mraw)
  rawbetas<-rawbetas[row.names(rawbetas) %in% auto.probes, QCmetrics$IntensityPass]

#run PCA  
  pca <- prcomp(t(na.omit(rawbetas)))
  save(pca, file = file.path(QCDir, "PCAbetas.rdat"))
  
}
  
  betas.scores = pca$x
  colnames(betas.scores) = paste(colnames(betas.scores), '_betas', sep='')
  
  betas.pca<-pca$sdev^2/sum(pca$sdev^2)
  betas.scores<-betas.scores[match(QCmetrics$Basename, QCmetrics$Basename[QCmetrics$IntensityPass]),]
  rownames(betas.scores)<-QCmetrics$Basename	
# only save PCs which explain > 1% of the variance
  QCmetrics<-cbind(QCmetrics,betas.scores[,which(betas.pca > 0.01)])




#----------------------------------------------------------------------#
# FLAG SAMPLES THAT FAIL THE FIRST QC STAGE
#----------------------------------------------------------------------#

# only exclude samples where reported sex != predicted Sex if sexCheck is turned on in config.r
if(sexCheck){
       QCmetrics$PassQC1 <- QCmetrics$IntensityPass & QCmetrics$PfiltPass & QCmetrics$BsConPass & QCmetrics$sexPass
} else {
       QCmetrics$PassQC1 <- QCmetrics$IntensityPass & QCmetrics$PfiltPass & QCmetrics$BsConPass
}

# Create QC pass/fail object
keepCols <- c("Basename", "Individual_ID", projVar,
              "IntensityPass", "PfiltPass", "BsConPass", "sexPass",
              "PassQC1")

QCSum<-QCmetrics[, colnames(QCmetrics) %in% keepCols]



#----------------------------------------------------------------------#
# SAVE AND CLOSE
#----------------------------------------------------------------------#

save(QCmetrics, file=file.path(QCDir, "QCmetrics.rdat"))
write.csv(QCSum, file.path(QCDir, "passQCStatusStage1AllSamples.csv"), row.names = F)

print("QC objects created and saved")



