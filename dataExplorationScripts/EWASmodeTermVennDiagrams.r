##---------------------------------------------------------------------#
##
## Title: Look at overlaps in EWAS results from pathology model and age
##        term from the null model
##
## Purpose of script: plot Venn diagram
##                    correlate effect sizes of DMPs in both
##                    
##                 
##
##                   
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# compare pathology DMPs and Age DMPs

#----------------------------------------------------------------------#
# LOAD PACKAGES etc
#----------------------------------------------------------------------#

library(VennDiagram)
library(dplyr)

cellType <- "NEUNneg"

load(paste0("3_analysis/results/", cellType,"EWASpathout.rdat"))
pathRes <- as.data.frame(outtab)

load(paste0("3_analysis/results/", cellType, "EWASout.rdat"))
nullRes <- as.data.frame(outtab)

rm(outtab)

#----------------------------------------------------------------------#
# LOAD RESULTS FILES
#----------------------------------------------------------------------#

# NEUN neg
load("3_analysis/results/NEUNnegEWASout.rdat")
nneg <- as.data.frame(outtab)
nneg$probe <- row.names(nneg)


#----------------------------------------------------------------------#
# EXTRACT SIG DMPS
#----------------------------------------------------------------------#

ageDmps <- row.names(nullRes %>% filter(nullAge_P < 0.05/nrow(nullRes)))
pathDmps <- row.names(pathRes %>% filter(Pathology_P < 0.05/nrow(pathRes))) # n.b. none for neun+
overlap <- intersect(ageDmps, pathDmps) # in both


#----------------------------------------------------------------------#
# venn diagram
#----------------------------------------------------------------------#

vennOut <- "3_analysis/plots/pathologyVsAgeDMPs_VennDiagram.jpeg"
  
venn.diagram(
  x = list(ageDmps, pathDmps),
  category.names = c("Age", "Pathology"),
  filename = vennOut,
  output=TRUE
)


#----------------------------------------------------------------------#
# correlate effect sizes for DMPS
#----------------------------------------------------------------------#

# define plot function
plotFunc <- function(dmpList){
  
  #create plot object
  plotdf <- cbind.data.frame(dmpList,
                             nullRes[dmpList, "nullAge_coeff"],
                             pathRes[dmpList, "Pathology_coeff"])
  colnames(plotdf) <- c("probeID", "Age", "Pathology")
  
  
  # plot params
  corOut <- paste0("3_analysis/plots/", as.list(match.call())$dmpList, "DMP_Correlation.pdf")
  plotTitle <- paste0(" Effect Sizes for ", 
                      as.list(match.call())$dmpList,
                      ": R=", signif(cor(plotdf$Age, plotdf$Pathology),2))
  
  # plot
  pdf(corOut)
  p <- ggplot(plotdf, aes(x=Age, y=Pathology))+
    geom_point()+
    ggtitle(plotTitle)+
    geom_abline(colour="red", linetype="dashed")
  print(p)
  
  dev.off()
  
  return(p)
  
}

# apply function to all DMP lists 
plotFunc(overlap)
plotFunc(pathDmps)  
plotFunc(ageDmps)
  
