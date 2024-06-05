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

#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#

library(CETYGO)

source("config.r")

#read in manifest and create probe ID with no suffix
man <- fread(manifest, skip=7, fill=TRUE, data.table=F)
man$probeID <- gsub("_.*$", "", man$IlmnID)


#----------------------------------------------------------------------#
# EPIC V1 AND MOUSE ARRAY PROVE OVERLAP
#----------------------------------------------------------------------#


load("pathToMRCNormalisedData/normalised.rdata")

inBoth <- intersect(row.names(celltypeNormbeta), man$probeID) # 141 probes

man[man$probeID %in% inBoth,]

                    