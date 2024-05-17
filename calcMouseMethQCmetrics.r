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
# file in the project folder


# MAY BE GOOD TO INCLUDE SESAME PROBE MASKING HERE TOO???? - in normalisation script

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#
print("loading packages...")

library(ENmix)
library(SummarizedExperiment)
library(dplyr)
library(plotrix)
library(factoextra)
library(stringr)


source("config.r")


#----------------------------------------------------------------------#
# LOAD IDATS TO RGSET
#----------------------------------------------------------------------#

sampleSheet <- read.csv(pheno, stringsAsFactors = F)
sampleSheet$Cell_Type <- trimws(sampleSheet$Cell_Type)

if(file.exists(file = file.path(QCDir, "rgSet.rdat"))){
  print("Loading rgSet")
  load(file = file.path(QCDir, "rgSet.rdat"))
} else{
  rgSet <- readidat(path = idatPath ,manifestfile=manifest ,recursive = TRUE)
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
failedProbes <- rownames(detP)[((rowSums(detP > pFiltProbeThresh)/ncol(detP)) * 100) > pFiltSampleThresh]


#----------------------------------------------------------------------#
# BISULPHITE CONVERSION
#----------------------------------------------------------------------#

if(file.exists(file = file.path(QCDir, "bsCon.rdat"))){
  print("Loading bsCon object")
  load(file = file.path(QCDir, "bsCon.rdat"))
} else{

ctrls <- metadata(rgSet)$ictrl

# subset to just mouse specific probes
bs.type1 <- ctrls$Address[ctrls$Type == "BISULFITE CONVERSION I"][11:20]
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


QCmetrics$BsCon <- BScon.med
QCmetrics$BsConPass <- ifelse(QCmetrics$BsCon > bsConThresh, TRUE, FALSE)
  


  

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


QCmetrics <- left_join(QCmetrics, sexPred)
QCmetrics$sexPass <- ifelse(QCmetrics$PredSex == QCmetrics$Sex, TRUE, FALSE)


#----------------------------------------------------------------------#
# REMOVE SAMPLES/PROBES THAT FAIL THE FIRST QC STAGE
#----------------------------------------------------------------------#

sampleSheet <- sampleSheet[sampleSheet$Basename %in% QCmetrics$Basename[QCmetrics$IntensityPass & QCmetrics$PfiltPass & QCmetrics$BsConPass],]

rgSet <- rgSet[ ,QCmetrics$Basename[QCmetrics$IntensityPass & QCmetrics$PfiltPass & QCmetrics$BsConPass]]

#remove failed probes
rgSet <- rgSet[!rgSet@elementMetadata$Name %in% failedProbes, ]

QCmetrics$PassQC1 <- QCmetrics$IntensityPass & QCmetrics$PfiltPass & QCmetrics$BsConPass


#----------------------------------------------------------------------#
# SAVE AND CLOSE
#----------------------------------------------------------------------#

save(QCmetrics, file=file.path(QCDir, "QCmetrics.rdat"))
