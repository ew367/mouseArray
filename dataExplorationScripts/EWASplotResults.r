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

#----------------------------------------------------------------------#
# LOAD PACKAGES etc
#----------------------------------------------------------------------#

library(dplyr)
library(ggplot2)

source("config.r")

cellType <- "NEUNneg"
res <- "NEUNnegEWASout.rdat"


resFile <- file.path("3_analysis/results/", res)
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
# LOAD AND ANNOTATE RESULTS FILE
#----------------------------------------------------------------------#

# load results

load(resFile)
outtab <- as.data.frame(outtab)
outtab$Name <- substr(rownames(outtab), 1,10)


# calculate bonferroni corrected threshold 

bonfP <- 0.05/nrow(outtab)

#----------------------------------------------------------------------#
# LOAD NORMALISED DATA AND SAMPLESHEET
#----------------------------------------------------------------------#

load(normData)


#----------------------------------------------------------------------#
# Sex specific DMPs
#----------------------------------------------------------------------#

# take top 20 most sig

sexDmps <- outtab %>% dplyr::select(contains('Sex')) %>% filter(nullSexM_P < 0.05/nrow(outtab))
sexTop20 <- row.names(sexDmps[order(sexDmps$nullSexM_P),][1:20,])


# plot top 20

sexPlotFile <- paste0("3_analysis/plots/", cellType, "Top20SexDMPs.pdf")
pdf(sexPlotFile)

for(i in sexTop20){

# get meth and sex data
sexBetas <- as.data.frame(celltypeNormbeta[i,])
sexBetas$Basename <- rownames(sexBetas)
sexBetas <- left_join(sexBetas, QCmetrics %>% dplyr::select(Basename, Sex))
colnames(sexBetas)[1] <- "Methylation"


# extract info for plot title
pval <- signif(outtab[i, "nullSexM_P"],3)
gene <- toString(unique(geneMan[geneMan$name == substr(i, 1,10), "Gene"]))
title <- paste0("Pval = ", pval, " , gene = ", gene)


# plot
p <- ggplot(sexBetas, aes(x=Sex, y=Methylation, fill = Sex))+
  geom_boxplot()+
  ggtitle(title)

print(p)

}

dev.off()



#----------------------------------------------------------------------#
# Group specific DMPs
#----------------------------------------------------------------------#


groupTop20 <- row.names(outtab[order(outtab$GroupWT_P),][1:20,])


# plot top 20
groupPlotFile <- paste0("3_analysis/plots/", cellType, "Top20groupDMPs.pdf")
pdf(groupPlotFile)


for(i in groupTop20){
  
  # get meth and sex data
  groupBetas <- as.data.frame(celltypeNormbeta[i,])
  groupBetas$Basename <- rownames(groupBetas)
  groupBetas <- left_join(groupBetas, QCmetrics %>% dplyr::select(Basename, Group))
  colnames(groupBetas)[1] <- "Methylation"
  
  
  # extract info for plot title
  pval <- signif(outtab[i, "GroupWT_P"],3)
  gene <- toString(unique(geneMan[geneMan$name == substr(i, 1,10), "Gene"]))
  title <- paste0("Pval = ", pval, " , gene = ", gene)
  
  
  # plot
  p <- ggplot(groupBetas, aes(x=Group, y=Methylation, fill = Group))+
    geom_boxplot()+
    ggtitle(title)
  
  print(p)
  
}

dev.off()

