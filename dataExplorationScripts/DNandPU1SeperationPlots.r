##---------------------------------------------------------------------#
##
## Title: DN / PU1 separation analyses plots
##
## Purpose of script: Creates extra plots in addition to standard QC
##                    output.
##
##
##                    Specifically to investigate the PU1 / DN
##                    lack of separation
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#




#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#

source("config.r")

load(file = file.path(QCDir, "QCmetrics.rdat"))
QCmetrics <- QCmetrics[QCmetrics$IntensityPass,]


load(file = file.path(QCDir, "PCAbetas.rdat"))



#----------------------------------------------------------------------#
# PCA COR PLOTS - TG ONLY 
#----------------------------------------------------------------------#

# subset to cols to cor
colsToKeep <- c(c(c("M.median","U.median","BsCon", "pFiltPass", "Cell_Type", "Batch", "DummyGroup",
                    "Age","Sex", "Pathology"),
                  colnames(QCmetrics)[grep("PC", colnames(QCmetrics))]))

colsToKeep <- colsToKeep[colsToKeep %in% colnames(QCmetrics)]

corDF <- QCmetrics[,colsToKeep]
corDF$Cell_Type <- as.numeric(as.factor(corDF$Cell_Type))
corDF$Batch <- as.numeric(as.factor(corDF$Batch))
corDF$Sex <- as.numeric(as.factor(corDF$Sex))

# subset to TG only
TGcor <- corDF[corDF$DummyGroup == 1,]
corrplot(cor(test, use = "p"))



#----------------------------------------------------------------------#
# PC CLUSTER PLOTS DIFFERENT COLOURS
#----------------------------------------------------------------------#

betas.scores = pca$x
colnames(betas.scores) = paste(colnames(betas.scores), '_betas', sep='')

betas.pca<-pca$sdev^2/sum(pca$sdev^2)
betas.scores<-betas.scores[match(QCmetrics$Basename, QCmetrics$Basename[QCmetrics$IntensityPass]),]
rownames(betas.scores)<-QCmetrics$Basename	


cellCols<-c("navy","turquoise")

plotGroup <- QCmetrics$Batch

par(mfrow = c(1,2))
par(mar = c(4,4,0.5,0.5))
plot(betas.scores[,1], betas.scores[,2],col = cellCols[as.factor(plotGroup)], pch = 16, xlab = "PC1", ylab = "PC2", cex.axis = 1.5, cex.lab = 1.5)

plot(betas.scores[,3], betas.scores[,4],col = cellCols[as.factor(plotGroup)], pch = 16, xlab = "PC3", ylab = "PC4", cex.axis = 1.5, cex.lab = 1.5)

legend("bottomright", pch = 16, col = cellCols, levels(as.factor(plotGroup)))



#----------------------------------------------------------------------#
# PC BOX PLOTS SEPARATED BY GROUP
#----------------------------------------------------------------------#s

ggplot(QCmetrics, aes(x=Full_Cell_Type, y = PC2_betas, fill = Group))+
  geom_boxplot()



#----------------------------------------------------------------------#
# PU1 - DN DIFFERENCE
#----------------------------------------------------------------------#

# subset to celltype
DN <-  QCmetrics[QCmetrics$Full_Cell_Type == "DN",c("Individual_ID", "PC2_betas", "Group")]
PU1 <- QCmetrics[QCmetrics$Full_Cell_Type == "PU1", c("Individual_ID", "PC2_betas", "Group")]

#find samples that have data for both cell types
inBoth <- intersect(DN$Individual_ID, PU1$Individual_ID)
DN <- DN[DN$Individual_ID %in% inBoth,]
PU1 <- PU1[PU1$Individual_ID %in% inBoth,]

#remove sample with ID 30 as it has two PU1 fractions...(?)
DN <- DN[DN$Individual_ID != "30",]
PU1 <- PU1[PU1$Individual_ID != "30",]

# join to single dataframe to plot
colnames(DN)[2] <- "DN_PC2"
colnames(PU1)[2] <- "PU1_PC2"
plotdf <- left_join(DN, PU1)


# calc difference and plot
plotdf$PC2_diff <- plotdf$DN_PC2 - plotdf$PU1_PC2

# boxplot
ggplot(plotdf, aes(x=Group, y = PC2_diff))+
  geom_boxplot()

# histogram
ggplot(plotdf, aes(x=PC2_diff, fill = Group))+
  geom_histogram(position = "dodge")

  
  


