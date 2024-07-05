##---------------------------------------------------------------------#
##
## Title: Look at overlaps in EWAS results from the 2 cell types
##
## Purpose of script: plot Venn diagrams for sex and group DMPS
##                    and correlate effect sizes  
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

library(VennDiagram)
library(dplyr)

var <- "nullSexM"

varP <- paste0(var, "_P")
varCoeff <- paste0(var, "_coeff")


#----------------------------------------------------------------------#
# LOAD RESULTS FILES
#----------------------------------------------------------------------#

# NEUN pos
load("3_analysis/results/NEUNposEWASout.rdat")
npos <- as.data.frame(outtab)
npos$probe <- row.names(npos)


# NEUN neg
load("3_analysis/results/NEUNnegEWASout.rdat")
nneg <- as.data.frame(outtab)
nneg$probe <- row.names(nneg)


#----------------------------------------------------------------------#
# EXTRACT SIG DMPS
#----------------------------------------------------------------------#

nposDmps <- row.names(npos %>% filter(!!as.symbol(varP) < 0.05/nrow(outtab)))
nnegDmps <- row.names(nneg %>% filter(!!as.symbol(varP) < 0.05/nrow(outtab)))


#----------------------------------------------------------------------#
# venn diagram
#----------------------------------------------------------------------#

vennOut <- paste0("3_analysis/plots/", var, "_VennDiagram.jpeg")
  
venn.diagram(
  x = list(nposDmps, nnegDmps),
  category.names = c("NeuN+", "NeuN-"),
  filename = vennOut,
  output=TRUE
)


#----------------------------------------------------------------------#
# correlate effect sizes
#----------------------------------------------------------------------#

# create dataframe with effect sizes for sex DMPs from both cell types

allSig <- unique(c(nnegDmps, nposDmps)) # sig in either cell type


# NeuN neg
plotdf <- nneg[allSig, c("probe", varCoeff)]
colnames(plotdf)[2] <- "NeuNnegCoeff"

# add NeuN pos
plotdf <- left_join(plotdf, npos %>% dplyr::select(probe, varCoeff))
colnames(plotdf)[3] <- "NeuNposCoeff"

# plot
corOut <- paste0("3_analysis/plots/", varCoeff, "Correlation.pdf")

pdf(corOut)
ggplot(plotdf, aes(x=NeuNnegCoeff, y=NeuNposCoeff))+
  geom_point()+
  ggtitle(paste0("Significant Sex DMPs. R=", signif(cor(plotdf$NeuNnegCoeff, plotdf$NeuNposCoeff),2)))+
  geom_abline(colour="red", linetype="dashed")+
  xlim(-0.25, 0.15)+
  ylim(-0.25, 0.15)
dev.off()
