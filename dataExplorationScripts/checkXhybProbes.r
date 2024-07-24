##---------------------------------------------------------------------#
##
## Title: Look at cross hybridisign probes and sig Sex diff DMPs
##
## Purpose of script: 
##                    
##                 
##
##                   
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# there are 235 probes that cross hybridise on the epic array and
# are included on the mouse array


#----------------------------------------------------------------------#
# LOAD PACKAGES AND IMPORT DATA
#----------------------------------------------------------------------#

# load EWAS results and extract sex DMPs
cellType <- "NEUNneg"
load(paste0("3_analysis/results/", cellType, "EWASout.rdat"))
outtab <- as.data.frame(outtab)
sexDMPs <- row.names(outtab %>% filter(nullSexM_P < 0.05/nrow(outtab)))



# load epic Matched Manifest
v1 <- read.csv("0_metadata/epicMatchedManifest.csv", stringsAsFactors = F)

# load file containing cross hybridising probes from Alice and subset to those on
# mouse array
xhyb <- read.csv("/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/Bulk/3_analysis/linearRegression/Sex/crossHybr/BLAT_cross_perc90.csv", stringsAsFactors = F)
xhyb <- xhyb[xhyb$qName %in% v1$epicProbeID,]


#----------------------------------------------------------------------#
# 
#----------------------------------------------------------------------#

se

xhybSexDMPs <- 
