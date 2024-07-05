##---------------------------------------------------------------------#
##
## Title: Plot EWAS results
##
## Purpose of script: Plot results from the EWAS analysis
##                    
##                 
##
##                   
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# no significant for pathology neun+

# change outFile for plots to include what coloured by


#----------------------------------------------------------------------#
# LOAD PACKAGES etc
#----------------------------------------------------------------------#

library(dplyr)
library(ggplot2)

source("config.r")



cellType <- "NEUNneg"
#model <- "path"
model <- ""
#var <- "Pathology"
var <- "nullSexM"
colVar <- "Sex"

varP <- paste0(var, "_P")

resOut <- paste0(cellType, "EWAS", model, "out.rdat")
resFile <- file.path("3_analysis/results/", resOut)


normData<-file.path(projDir, "2_normalised/normalisedData.rdat")

#----------------------------------------------------------------------#
# LOAD MANIFESTS
#----------------------------------------------------------------------#

# standard manifest
man <- fread(manifest, skip=7, fill=TRUE, data.table=F)
man <- man[!duplicated(man$Name),] # remove replicate probes

# gene manifest
geneMan <- fread(geneInfo, fill=TRUE, data.table=F)


#----------------------------------------------------------------------#
# LOAD NORMALISED DATA AND SAMPLESHEET
#----------------------------------------------------------------------#

load(normData)

#subset to cell type of interest

## subset to cell type samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell_Type == cellType),]


# subset beta matrix to cell type specific samples
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]



#----------------------------------------------------------------------#
# LOAD RESULTS FILE
#----------------------------------------------------------------------#

# load results

load(resFile)
outtab <- as.data.frame(outtab)
outtab$Name <- substr(rownames(outtab), 1,10)


# calculate bonferroni corrected threshold 

bonfP <- 0.05/nrow(outtab)


#----------------------------------------------------------------------#
# Plot var of interest
#----------------------------------------------------------------------#

# get IDS for 20 most significant dmps for variable of interest
dmps <- row.names(outtab %>% arrange(!!as.symbol(varP)) %>% dplyr::slice(1:20))

outFile <- paste0("3_analysis/plots/", cellType, "Top20", var, "DMPs.pdf")

pdf(outFile)
for(i in dmps){
  
  x <- ifelse(grepl("null", var), substr(var, 5, nchar(var)), var) # if present, strip null from var to match with format from QCmetrics file
  
  x <- ifelse(grepl("M", x), substr(x, 1, nchar(x)-1), x) # if present, strip M from end of var to match with format from QCmetrics file - (problem with Sex)
  
  # get meth and sex data
  dmpBetas <- as.data.frame(celltypeNormbeta[i,])
  dmpBetas$Basename <- rownames(dmpBetas)
  dmpBetas <- left_join(dmpBetas, QCmetrics %>% dplyr::select(Basename, !!as.symbol(x), !!as.symbol(colVar)))
  colnames(dmpBetas)[1] <- "Methylation"

  
  
  # extract info for plot title
  pval <- signif(outtab[i, varP],3)
  gene <- toString(unique(geneMan[geneMan$name == substr(i, 1,10), "Gene"]))
  chr <- unique(geneMan$chrom[geneMan$name == substr(i, 1,10)])
  title <- paste0("Pval = ", pval, ", ", chr, ": " , gene)
  
  
  # if statement for if model is pathology then point not box plot
  
  if(var == "Pathology" | var == "nullAge"){
    
    p <- ggplot(dmpBetas, aes_string(x=x, y="Methylation", colour = colVar))+
      geom_point()+
      #geom_smooth()+
      ggtitle(title)
    
  
  } else {
    
    p <- ggplot(dmpBetas, aes_string(x=colVar, y="Methylation", fill = colVar))+
      geom_boxplot()+
      ggtitle(title)
    
  }
  
  print(p)
  
}

dev.off()

