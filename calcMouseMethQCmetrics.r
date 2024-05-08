##---------------------------------------------------------------------#
##
## Title: Calculate QC metrics for mouse methylation array data
##
## Purpose of script: 
##
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#



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

# if Rgset exists load - else create
#rgSet <- readidat(path = idatPath ,manifestfile=manifest ,recursive = TRUE)
#save(rgSet, file=file.path(normDir, "rgSet.rdat"))
load(file.path(normDir, "rgSet.rdat"))


#----------------------------------------------------------------------#
# Calculate intensities
#----------------------------------------------------------------------#


mraw <- getmeth(rgSet)

betas<-getB(mraw)

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

detP <- calcdetP(rgSet) 

# check if any samples have > 1 percent of probes with a detection p value of > pFiltThresh
pfiltdf <- data.frame(matrix(ncol = 2, nrow = nrow(QCmetrics)))
colnames(pfiltdf) <- c("Basename", "PercProbesFail")

for(i in 1:ncol(detP)){
  pfiltdf$Basename[i] <- colnames(detP)[i]
  pfiltdf$PercProbesFail[i] <- sum(detP[,i] > pFiltProbeThresh)/nrow(detP)*100
}

pfiltdf$PfiltPass <- ifelse(pfiltdf$PercProbes < 1, TRUE, FALSE)

QCmetrics <- left_join(QCmetrics, pfiltdf, by = "Basename")


# do this just prior to normalising?
# check if any probes fail in more than pFiltSampleThresh of samples
#failedProbes <- rownames(detP)[((rowSums(detP > pFiltProbeThresh)/ncol(detP)) * 100) > pFiltSampleThresh]

#remove failed probes
#rgSet <- rgSet[!rgSet@elementMetadata$Name %in% failedProbes, ]



#----------------------------------------------------------------------#
# BISULPHITE CONVERSION
#----------------------------------------------------------------------#

ctrls <- metadata(rgSet)$ictrl
ctrls <- ctrls[ctrls$Address %in% rownames(rgSet),]
ctrl_r <- assays(rgSet)$Red[ctrls$Address,]
ctrl_g <- assays(rgSet)$Green[ctrls$Address,]


# from here: https://emea.support.illumina.com/bulletins/2021/07/infinium-mouse-methylation-beadchip-genomestudio-controls-interp.html
green1 <- c("BS1-396C_MUS", "BS1-396U_MUS", "BS1-140C_MUS", "BS1-140U_MUS")
red1 <- c("BS1-409C_MUS", "BS1-409U_MUS", "BS1-318C_MUS", "BS1-318U_MUS", "BS1-317C_MUS", "BS1-317U_MUS")
red2 <- c("BS2-330_MUS", "BS2-505_MUS", "BS2-649_MUS")

#green
cc=ctrls[(ctrls$ExtendedType %in% green1),]
I_green=colMedians(ctrl_g[cc$Address,])

#red
cc=ctrls[(ctrls$ExtendedType %in% red1),]
I_red=colMedians(ctrl_r[cc$Address,])

cc=ctrls[(ctrls$ExtendedType %in% red2),]
II_red=colMedians(ctrl_r[cc$Address,])


#BSI.betas <- I_red/(I_red + I_green)
#Bisulphite <- BSI.betas*100
#hist(Bisulphite, xlab = "Median % BS conversion", main = "Bisulphite Converstion Statistics")
#png("Bisulphite_conversion.png", width = 800, height = 600)
#hist(Bisulphite, xlab = "Median % BS conversion", main = "Bisulphite Converstion Statistics")
#dev.off()


# green
cc <- as.data.frame(ctrls[(ctrls$ExtendedType %in% green1),])
I_green <- t(ctrl_g[cc$Address,])
colnames(I_green) <- c("high", "low", "background", "background")

# calc SE

#for each col in I_green, calc SEx1.96

#thresh <- calcThresh(I_green)

#test <- c()

#for(i in 1:ncol(ctrlType)){
 # threshVal <- 1.96*std.error(ctrlType[,i])
  #meanVal <- mean(ctrlType[,i])
  #print(threshVal)
  #print(meanVal)
  
  #for(j in 1:nrow(ctrlType)){
   # if(ctrlType[j,i] < meanVal+threshVal & ctrlType[j,i] > meanVal-threshVal) {
    #  print("TRUE")
    #}
 # }
#}

#for(i in 1:nrow(I_green)){
 # for(j in 1:ncol(I_green)){
    
  #  if(I_green[i,j] > )
  #}
#}
  

I_green <- reshape2::melt(I_green)
colnames(I_green) <- c("Basename", "ExpectedIntensity", "Intensity")
I_green$ExpectedIntensity <- as.factor(I_green$ExpectedIntensity)
I_green$channel <- "green_I"


#ggplot(I_green, aes(ExpectedIntensity, Intensity)) +
 # geom_boxplot()+
  #ggtitle("Green_I")


# red I
cc=ctrls[(ctrls$ExtendedType %in% red1),]
I_red=t(ctrl_r[cc$Address,])
colnames(I_red) <- c("high", "medium", "low", "background", "background", "background")
I_red <- reshape2::melt(I_red)
colnames(I_red) <- c("Basename", "ExpectedIntensity", "Intensity")
I_red$ExpectedIntensity <- as.factor(I_red$ExpectedIntensity)
I_red$channel <- "red_I"

#ggplot(I_red, aes(ExpectedIntensity, Intensity)) + 
 # geom_boxplot() +
 # ggtitle("Red_I")


# red II
cc=ctrls[(ctrls$ExtendedType %in% red2),]
II_red=t(ctrl_r[cc$Address,])
colnames(II_red) <- c("high", "high", "high")
II_red <- reshape2::melt(II_red)
colnames(II_red) <- c("Basename", "ExpectedIntensity", "Intensity")
II_red$ExpectedIntensity <- as.factor(II_red$ExpectedIntensity)
II_red$channel <- "red_II"

#ggplot(II_red, aes(ExpectedIntensity, Intensity)) + 
 # geom_boxplot() +
  #ggtitle("Red_II")

#combine together
plotdf <- rbind(I_green, I_red, II_red)
colnames(plotdf)[4] <- "ProbeType"
plotdf$ExpectedIntensity <- factor(plotdf$ExpectedIntensity, levels = c("high", "medium", "low" , "background"))
ggplot(plotdf, aes(ExpectedIntensity, Intensity, fill = ProbeType))+
  geom_boxplot()+
  ggtitle("Bisulphite Conversion Control Probes")+
  scale_fill_manual(values=c("#00CC66", "#CC3300", "#FF00CC"))
#ggsave("Bisulphite Conversion Control Probes.png")

# plot after removing failed samples

# use this to subset plotdf to only passed samples
passedOnly <- QCmetrics$Basename[QCmetrics$IntensityPass & QCmetrics$PfiltPass]

# if ctrol probes are > +/- 1.96×SE from mean then flag??
# do this on each probe type individually??

plotdfPass <- plotdf[plotdf$Basename %in% passedOnly,]

ggplot(plotdfPass, aes(ExpectedIntensity, Intensity, fill = ProbeType))+
  geom_boxplot()+
  ggtitle("Bisulphite Conversion Control Probes passed Samples Only")+
  scale_fill_manual(values=c("#00CC66", "#CC3300", "#FF00CC"))


findoutlier <- function(x) {
  return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
}

# for each col in each probeType, find if outlier on any col


#----------------------------------------------------------------------#
# REMOVE FAILED SAMPLES HERE??
#----------------------------------------------------------------------#


#QCmetrics$Basename[!(QCmetrics$IntensityPass & QCmetrics$PfiltPass & QCmetrics$sexPass)]





#----------------------------------------------------------------------#
# SEX CHECK (ENMIX)
#----------------------------------------------------------------------#

sex <- predSex(rgSet)
colnames(sex) <- c("Basename", "PredSex")

QCmetrics <- left_join(QCmetrics, sex)
QCmetrics$sexPass <- ifelse(QCmetrics$PredSex == QCmetrics$Sex, TRUE, FALSE)



#----------------------------------------------------------------------#
# CLUSTERING
#----------------------------------------------------------------------#

sigma <- apply(betas, 1, sd)

# replace NAN with 0.5 in 16 cases
#test <-betas[complete.cases(betas),]

pca.res <- prcomp(t(betas))

#----------------------------------------------------------------------#
# Normalise within cell type (move to another script)
#----------------------------------------------------------------------#


#----------------------------------------------------------------------#
# SAVE AND EXIT
#----------------------------------------------------------------------#

# QCMetrics
# RgSet

# mraw
# detP
# pc.res
# sigma


# bscon stuff
