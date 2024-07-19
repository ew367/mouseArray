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

#  col (2) design: contains the annotation of the probes. For the syntenic EPIC probe mapping, search for EPIC prefix in the design column. 

#   Probe_ID        design                    
#1 cg00101675_BC21 EPIC;cg00101675           
#2 cg00116289_BC21 EPIC;cg00116289           
#3 cg00211372_TC21 EPIC;cg00211372           
#4 cg00531009_BC21 EPIC;cg00531009           
#5 cg00747726_TC21 EPIC;cg00747726           
#6 cg00896209_TC21 CGI;chr3:89169426-89169681


# the cg IDs don't appear to have suffixes so it looks to be matched
# to epicv1 data only

#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#

library(readr)
library(dplyr)
library(data.table)

source("config.r")

#read in mouse manifest and create probe ID with no suffix
man <- fread(manifest, skip=7, fill=TRUE, data.table=F)
man$probeID <- gsub("_.*$", "", man$IlmnID)


# read in gene annotation manifest
geneMan <- fread(geneInfo, fill=TRUE, data.table=F)



# read in gene anno file
designMan <- as.data.frame(readr::read_tsv("0_metadata/MM285.design.tsv"))


# read in human manifest
hsManifest <- fread("/lustre/projects/Research_Project-MRC190311/references/EPICArray/MethylationEPIC_v-1-0_B5.csv", skip=7, fill=TRUE, data.table=F)



#----------------------------------------------------------------------#
# CREATE NEW MOUSE MANIFEST
#----------------------------------------------------------------------#

# add col to manifest for gene info
#gene <- c()
#for(i in man$Name){
#  gene <- c(gene, toString(unique(geneMan[geneMan$name == substr(i, 1,10), "Gene"]##)))
#}

#man$gene <- gene

#write.csv(man, "0_metadata/mm10manifestWithGeneInfo.csv")


#----------------------------------------------------------------------#
  # CREATE NEW EPIC MATCHED MANIFEST
#----------------------------------------------------------------------#


# create col with matched epic probe ID

epicMatched <- designMan[grep("EPIC", designMan$design),]
epicMatched$epicProbeID <- substr(epicMatched$design, 6, 15)


# create col with mouse gene

gene <- c()
for(i in epicMatched$Probe_ID){
  gene <- c(gene, toString(unique(geneMan[geneMan$name == substr(i, 1,10), "Gene"])))
}

epicMatched$mouseGene <- gene


# create col with human gene

humanGene <- hsManifest[hsManifest$IlmnID %in% epicMatched$epicProbeID, c("Name", "UCSC_RefGene_Name")]
colnames(humanGene) <- c("epicProbeID", "humanGene")

epicMatched <- left_join(epicMatched, humanGene)


write.csv(epicMatched, "0_metadata/epicMatchedManifest.csv")



                    