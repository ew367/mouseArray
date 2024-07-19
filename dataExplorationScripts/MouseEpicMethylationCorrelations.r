##---------------------------------------------------------------------#
##
## Title: Human - Mouse DMP methylation
##
## Purpose of script: correlate the results from the EWAS with otehr studies
##
##
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# This script currently compares AGE ewas DMPs to BDR dmps
# the BDR data it uses is not normalised

#----------------------------------------------------------------------#
# DEFINE PARAMETERS ETC.
#----------------------------------------------------------------------#

library(dplyr)
library(ggplot2)
library(reshape2)

source("config.r")


cellType <- "NEUNpos"
hsCellType <- "NeuN+"
#model <- "path"
model <- ""
var <- "nullAge"
#var <- "nullSexM"
#colVar <- "Sex"

varP <- paste0(var, "_P")


# mouse EWAS res files
resOut <- paste0(cellType, "EWAS", model, "out.rdat")
resFile <- file.path("3_analysis/results/", resOut)


# mouse methylation data
normData <- file.path(normDir, "normalisedData.rdat")


# comparison dataset files
epicBetas <- "otherDatasets/rawbetas.rdat"
epicQCsum <- "/lustre/projects/Research_Project-MRC190311/DNAm/sortedBDREPICv1/2_gds/QCmetrics/passQCStatusStage3AllSamples.csv" 


#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#

# load epic dataset to compare EWAS DMPS to
epicQC <- read.csv(epicQCsum)
load(epicBetas)

# subset to passed samples  only
rawbetas <- rawbetas[, epicQC$passQCS3]
epicQC <- epicQC[epicQC$passQCS3,]
epicQC$Cell.type <- as.character(epicQC$Cell.type)

# subset to celltype betas
epicCellQC <- as.character(epicQC[which(epicQC$Cell.type == hsCellType), "Basename"])
epicBetas <- as.data.frame(rawbetas[, epicCellQC])
epicBetas$HumanMean <- rowMeans(epicBetas) # calc mean

epicBetas$epicProbeID <- row.names(epicBetas)



# load manifest
epicMatched <- read.csv("0_metadata/epicMatchedManifest.csv")


# load mouse EWAS results
load(resFile)
outtab <- as.data.frame(outtab)


# load mouse methylation data
load(normData)
celltypeNormbeta <- as.data.frame(celltypeNormbeta)


# subset to cell type of interest
QCmetrics<-QCmetrics[which(QCmetrics$Cell_Type == cellType),]
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]

# calcualte mean
celltypeNormbeta$MouseMean <- rowMeans(celltypeNormbeta)
celltypeNormbeta$Probe_ID <- row.names(celltypeNormbeta)



#----------------------------------------------------------------------#
# ANALYSE
#----------------------------------------------------------------------#

# compare methylation levels at matched probes

species <- as.data.frame(rbind(cbind(colnames(epicBetas), "human"), cbind(colnames(celltypeNormbeta), "mouse")))
colnames(species) <- c("variable", "species")

resDf <- left_join(epicMatched %>% dplyr::select(Probe_ID, epicProbeID), celltypeNormbeta) # add mouse methylation
resDf <-  left_join(resDf, epicBetas)# add human methlyation
resDf <- resDf[complete.cases(resDf),] # keep only sites have both for
resDf <- melt(resDf) # melt
resDf <- left_join(resDf, species) #add species


#ggplot(plotdf, aes(x=Probe_ID, y=value, group=species))+
 # geom_point()




#----------------------------------------------------------------------#
# CORRELATION PLOTS
#----------------------------------------------------------------------#

mouseDf <- left_join(epicMatched %>% dplyr::select(Probe_ID, epicProbeID), celltypeNormbeta %>% dplyr::select(Probe_ID, MouseMean))
#mouseDf <- melt(mouseDf)
#colnames(mouseDf)[3] <- "mouseBasenames"

humanDf <- left_join(epicMatched %>% dplyr::select(Probe_ID, epicProbeID), epicBetas%>% dplyr::select(epicProbeID, HumanMean))
#humanDf <- melt(humanDf)
#colnames(humanDf)[3] <- "humanBasenames"

corDf <- left_join(mouseDf, humanDf)
corDf <- corDf[complete.cases(corDf),]

# all matched probes
ggplot(corDf, aes(x=MouseMean, y=HumanMean))+
  geom_point()+
  ggtitle(paste0("All matched probes - ", cellType, " R=", signif(cor.test(corDf$MouseMean, corDf$HumanMean)$estimate, 3)))


# plot just sig for term of interest

dmps <- row.names(outtab %>% filter(!!as.symbol(varP) < 0.05/nrow(outtab)))

sigCorDf <- corDf[corDf$Probe_ID %in% dmps,]

ggplot(sigCorDf, aes(x=MouseMean, y=HumanMean))+
  geom_point()+
  ggtitle(paste0("Age DMP matched probes - ", cellType, " R=", signif(cor.test(sigCorDf$MouseMean, sigCorDf$HumanMean)$estimate, 3)))











