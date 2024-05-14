##---------------------------------------------------------------------#
##
## Title: explore BScon metrics for mouse methylation array data
##
## Purpose of script: The bisulphite conversion statistic cannot be calculated the ##                   same way as the DNAm QC pipeline (brainFans repo). This is 
##                    because: a) the libraries are not set up to work with the 
##                    mouse array and b) the computed values appear to be lower.
##
##                    This script explores different options for calculating a 
##                    bscon statistic that can be included in the standard QC 
##                    pipline
##                    
##
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# First section adapted from code for epicV2 (from Alice via slack 13/05/2024)
# This does not seem to work very well as despite the data looking good 
# on other QC metrics, almost all samples would appear to fail at 
# threshold of 80%, which seems unlikely as the intensity values look good.

# THE enmix package outputs some bisulphite QC plots etc but the QCinfo function 
# is broken
# https://rdrr.io/bioc/ENmix/man/QCinfo.html

# The SeSame package does not have bscon output for mouse array
# https://www.bioconductor.org/packages/release/bioc/vignettes/sesame/inst/doc/nonhuman.html
# NB: IT DOES MASK possibly 'dodgy' probes

# the minfi package appears to now have some external support for mouse array data
# https://github.com/chiaraherzog/IlluminaMouseMethylationanno.12.v1.mm10
# but the files won't download from github


#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

# run using R version 4.2.1 (2022-06-23) -- "Funny-Looking Kid"

library(ENmix)
library(SummarizedExperiment)
library(sesame)
remotes::install_github("chiaraherzog/IlluminaMouseMethylationmanifest")
library(IlluminaMouseMethylationmanifest)
remotes::install_github("https://github.com/chiaraherzog/IlluminaMouseMethylationanno.12.v1.mm10")
library(IlluminaMouseMethylationanno.12.v1.mm10)

source("config.r")

#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#

sampleSheet <- read.csv(pheno, stringsAsFactors = F)

# if Rgset exists load - else create
#rgSet <- readidat(path = idatPath ,manifestfile=manifest ,recursive = TRUE)
#save(rgSet, file=file.path(normDir, "rgSet.rdat"))
load(file.path(normDir, "rgSet.rdat"))

load("QCmetrics.rdat")


# reovme failed samples
passedOnly <- QCmetrics$Basename[QCmetrics$IntensityPass & QCmetrics$PfiltPass]
rgSet <- rgSet[, passedOnly]

#mraw <- getmeth(rgSet)
#betas<-getB(mraw)
load("betas.rdat")

#----------------------------------------------------------------------#
# USING CODE DESIGNED FOR EPIC V2
#----------------------------------------------------------------------#

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

par(mfrow=c(1,2))
boxplot(t(bs.red.type1), col=c('black','red')[factor(red.med>10000)], xlab='Type 1 bisulfite probes', ylab='Intensity', main='Red channel')
abline(h=10000, lty=3)
boxplot(t(bs.green.type1), col=c('black','green')[factor(green.med>10000)], xlab='Type 1 bisulfite probes', ylab='Intensity', main='Green channel')
abline(h=10000, lty=3)

redGreater<-rowMeans(bs.red.type1) > rowMeans(bs.green.type1)

BScon1<-rbind(
  bs.red.type1[which(redGreater),] / ( bs.red.type1[which(redGreater),] + bs.green.type1[which(redGreater),] ), 
  bs.green.type1[which(!redGreater),] / ( bs.red.type1[which(!redGreater),] + bs.green.type1[which(!redGreater),] )
)

BScon2 <- bs.red.type2 / ( bs.red.type2 + bs.green.type2 )

BSconAll<-rbind(BScon1, BScon2)
BScon.med<-apply(BSconAll, 2, median)
hist(BScon.med)
BSconAll<-rbind(BSconAll, BScon.med)*100


boxplot(t(BSconAll), ylab = "% conversion", axes = FALSE)
axis(2)
axis(1, c("Type I", "Type II", "Median"), at = c(length(bs.type1)/2,length(bs.type1)+(length(bs.type2)/2)+1, nrow(BSconAll)))
abline(v = length(bs.type1)+0.5, lty = 2)
abline(v = nrow(BSconAll)-0.5, lty = 2)
box()
#*

boxplot(t(BSconAll))
abline(h=80, lty=3)



#----------------------------------------------------------------------#
# ENmix package
#----------------------------------------------------------------------#


QCinfo(rgSet) # default parameters 

# output (- not working, no recent updates)
# NA  samples with percentage of low quality CpG value greater than  0.05  or bisulfite intensity less than  NaN
#NA  CpGs with percentage of low quality value greater than  0.05
#Ploting qc_sample.jpg ...Error in plot.window(...) : need finite 'ylim' values
#In addition: Warning messages:
 # 1: In min(x) : no non-missing arguments to min; returning Inf
#2: In max(x) : no non-missing arguments to max; returning -Inf


#----------------------------------------------------------------------#
# seSame package
#----------------------------------------------------------------------#

# prep codes "TQCDPB" recommended: https://github.com/zwdzwd/sesame/blob/devel/vignettes/nonhuman.Rmd
# T	inferStrain	Set strain-specific mask (mouse)
# Q	qualityMask	Mask probes of poor design
# C	inferInfiniumIChannel	Infer channel for Infinium-I probes
# D	dyeBiasNL	Dye bias correction (non-linear)
# P	pOOBAH	Detection p-value masking using oob
# P	pOOBAH	Detection p-value masking using oob


betas = openSesame(idatPath, prep="TQCDPB", func = getBetas) # this just outputs a betas matrix
# 296070    204 - this is larger than the ENmix matrix..
qcs = openSesame(idatPath, prep="TQCDPB", func=sesameQC_calcStats)
head(do.call(rbind, lapply(qcs, as.data.frame))) # combine to df for easy viewing
colnames(do.call(rbind, lapply(qcs, as.data.frame))) # view what output has been generated - no bscon stats


#all_sdfs = mclapply(searchIDATprefixes(idatPath), readIDATpair, mc.cores=16)
load("allSamplesSDFsSesame.rdat")

# this doesn't appear to have worked as expected as lots of the slots are filled with NA's


#bisConversionControl(sset, use.median = FALSE) # doesn't work... stops as not epic or 450k.

# code from: https://www.bioconductor.org/packages/devel/bioc/vignettes/sesame/inst/doc/sesame.html

#bisConversionControl <- function(sset, use.median=FALSE) {
  
# stopifnot(sset@platform %in% c('EPIC','HM450'))
#  extC <- sesameDataGet(paste0(sset@platform, '.probeInfo'))$typeI.extC
# extT <- sesameDataGet(paste0(sset@platform, '.probeInfo'))$typeI.extT
#  prbs <- rownames(oobG(sset))
# extC <- intersect(prbs, extC)
# extT <- intersect(prbs, extT)
#  if (use.median) {
#   median(oobG(sset)[extC,], na.rm=TRUE) /
#      median(oobG(sset)[extT,], na.rm=TRUE)
# } else {
#  mean(oobG(sset)[extC,], na.rm=TRUE) /
#   mean(oobG(sset)[extT,], na.rm=TRUE)
# }
#}



#----------------------------------------------------------------------#
# MINFI
#----------------------------------------------------------------------#
# can create RGset but can't access github packages
minfiRGset <- read.metharray.exp(base = idatPath, targets = sampleSheet, force = TRUE)
minfiRGset@annotation <-  c(array = "IlluminaMouseMethylation", annotation = "12.v1.mm10")



#----------------------------------------------------------------------#
# plot distributions of the control Probes
#----------------------------------------------------------------------#

# just use BS I and BS II


# get the control probe data from the rgset
ctrls <- metadata(rgSet)$ictrl

# get the Adresses of the BScon probes
bs.type1 <- ctrls$Address[ctrls$Type == "BISULFITE CONVERSION I"][11:20]
bs.type2 <- ctrls$Address[ctrls$Type == "BISULFITE CONVERSION II"][4:6]

# use the adresses to get the green/red values
bs.green.type1 <- assays(rgSet)$Green[bs.type1,]
bs.green.type2 <- assays(rgSet)$Green[bs.type2,]
bs.red.type1 <- assays(rgSet)$Red[bs.type1,]
bs.red.type2 <- assays(rgSet)$Red[bs.type2,]

# calculate the beta values for every sample
BSI.betas <- bs.red.type1/(bs.red.type1 + bs.green.type1)*100 # 10 probes
BSII.betas <- bs.red.type2/(bs.red.type2 + bs.green.type2)*100 # 3 probes

# HSA is human only and MUS is mouse only
type1CU <- as.data.frame(ctrls[ctrls$Type == "BISULFITE CONVERSION I",])[11:20,] 
#no C/U letters for type2 control probes
type1CU$letter <- substr(type1CU$ExtendedType, 8, 8)

# plot
densityPlot(BSI.betas, main = "Bisulphite I", sampGroups = type1CU$letter)
densityPlot(BSII.betas, main = "Bisulphite II")




