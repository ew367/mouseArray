##---------------------------------------------------------------------#
##
## Title: App Gene
##
## Purpose of script: look at CpGs localised near App gene
##                    
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# There are 43 probes in the App gene
# 40 in this dataset after QC

# https://www.nature.com/articles/nn.3697#Sec1
# this paper explains where the mutations are


#----------------------------------------------------------------------#
# LOAD PACKAGES etc
#----------------------------------------------------------------------#

source("config.r")

cellType <- "NEUNpos"

window <- 1500 # number of bp either side of APP gene to include
# nb this window doesn't add an extra cpgs

# location of App Gene from USCS - this also matches with what is App in manifest
appStart <- 84954443 
appEnd <- 85173758


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
# GET App CpGs
#----------------------------------------------------------------------#

appCPGs <- man$IlmnID[which(man$CHR == 16 & man$MAPINFO > appStart - window & man$MAPINFO < appEnd + window)]
mutCPGs <- man$IlmnID[which(man$CHR == 16 & man$MAPINFO > (appStart + 206484)-window & man$MAPINFO < appEnd + window)]


#----------------------------------------------------------------------#
# LOAD METHYLATION DATA
#----------------------------------------------------------------------#

load(normData)


#subset to cell type of interest

## subset to cell type samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell_Type == cellType),]

# subset beta matrix to cell type specific samples
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]


#subset to app gene probes only
# note there are 3 probes in app which are not in the QC'd dataset
appBetas <- celltypeNormbeta[rownames(celltypeNormbeta) %in% appCPGs,]
mutBetas <- celltypeNormbeta[rownames(celltypeNormbeta) %in% mutCPGs,]

#----------------------------------------------------------------------#
# PLOT
#----------------------------------------------------------------------#

# extract mapinfo from manifest for plot
mapinfo <- man %>% dplyr::select(IlmnID, MAPINFO)

#check betas and pheno are in same order
identical(colnames(appBetas), QCmetrics$Basename)

# create plot object
plotdf <- cbind(t(appBetas), QCmetrics %>% dplyr::select(Group))
plotdf <- reshape2::melt(plotdf)
colnames(plotdf)[2:3] <- c("IlmnID", "Methylation")
plotdf <- left_join(plotdf, mapinfo)


#plot!

outFile <- paste0("3_analysis/plots/", cellType, "AppGeneCpGs.pdf")

pdf(outFile)
P <- ggplot(plotdf, aes(x=MAPINFO, y=Methylation, color=Group)) +
  geom_point()+
  geom_smooth() +
 ggtitle(paste0(cellType, " - CpGs in App gene"))+
  geom_vline(xintercept = appStart+206484)+
  geom_vline(xintercept = appStart+207984)+
  geom_vline(xintercept = appStart+210922)+
  geom_vline(xintercept = appStart+216617)
  
  
  

print(P)

dev.off()





#### zoomed in on mutation region


# create plot object



plotdf <- cbind(t(mutBetas), QCmetrics %>% dplyr::select(Group))
plotdf <- reshape2::melt(plotdf)
colnames(plotdf)[2:3] <- c("IlmnID", "Methylation")
plotdf <- left_join(plotdf, mapinfo)


#plot!

outFile <- paste0("3_analysis/plots/", cellType, "MutGeneCpGs.pdf")

pdf(outFile)
P <- ggplot(plotdf, aes(x=MAPINFO, y=Methylation, color=Group)) +
  geom_point()+
  geom_smooth() +
  ggtitle(paste0(cellType, " - CpGs in mutation region of App gene"))+
  geom_vline(xintercept = appStart+206484)+
  geom_vline(xintercept = appStart+207984)+
  geom_vline(xintercept = appStart+210922)+
  geom_vline(xintercept = appStart+216617)+
  ylim(0, 0.8)


print(P)

dev.off()




