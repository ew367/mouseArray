##---------------------------------------------------------------------#
##
## Title: Human - Mouse CETYGO overlap
##
## Purpose of script: use human brain ref panal to deconvolute 
##                    cell sorted mouse data
##
##
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# https://github.com/ejh243/CETYGO/wiki/Deconvolution-of-brain-cell-types

# according to https://www.illumina.com/content/dam/illumina/gcs/assembled-assets/marketing-literature/infinium-mouse-methylation-array-data-sheet-370-2020-002/infinium-mouse-methylation-array-data-sheet-370-2020-002.pdf
# There is a "Human MethylationEPIC liftover"


# https://zwdzwd.github.io/InfiniumAnnotation#mouse this contains a file called "Mouse array design groups (MM285)"
# which shows the overlap

#   Probe_ID        design                    
#1 cg00101675_BC21 EPIC;cg00101675           
#2 cg00116289_BC21 EPIC;cg00116289           
#3 cg00211372_TC21 EPIC;cg00211372           
#4 cg00531009_BC21 EPIC;cg00531009           
#5 cg00747726_TC21 EPIC;cg00747726           
#6 cg00896209_TC21 CGI;chr3:89169426-89169681

#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#

library(CETYGO)
library(readr)

source("config.r")

#read in manifest and create probe ID with no suffix
man <- fread(manifest, skip=7, fill=TRUE, data.table=F)
man$probeID <- gsub("_.*$", "", man$IlmnID)

designMan <- as.data.frame(read_tsv("0_metadata/MM285.design.tsv"))
epicMatched <- designMan[grep("EPIC", designMan$design),]

epicIDsinMouse <- #
  ### need to use regexp to extract cg and following 8 digits


#----------------------------------------------------------------------#
# EPIC V1 AND MOUSE ARRAY PROVE OVERLAP
#----------------------------------------------------------------------#


load("pathToMRCNormalisedData/normalised.rdata")

inBoth <- intersect(row.names(celltypeNormbeta), man$probeID) # 141 probes

man[man$probeID %in% inBoth,]

                    