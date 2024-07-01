##---------------------------------------------------------------------#
##
## Title: Plot Venn diagrams of EWAS results from the 2 cell types
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

library(VennDiagram)

# NEUN pos
load("3_analysis/results/NEUNposEWASout.rdat")
npos <- as.dataframe(outtab)

nposSexDmps <- npos %>% dplyr::select(contains('Sex')) %>% filter(nullSexM_P < 0.05/nrow(outtab))
nposSexTop20 <- row.names(nposSexDmps[order(nposSexDmps$nullSexM_P),][1:20,])


# NEUN neg
load("3_analysis/results/NEUNnegEWASout.rdat")
nneg <- as.data.frame(outtab)

nnegSexDmps <- nneg %>% dplyr::select(contains('Sex')) %>% filter(nullSexM_P < 0.05/nrow(outtab))
nnegSexTop20 <- row.names(nnegSexDmps[order(nnegSexDmps$nullSexM_P),][1:20,])


#----------------------------------------------------------------------#
# venn diagram
#----------------------------------------------------------------------#

#myCol <- brewer.pal(2, "Pastel2")

venn.diagram(
  x = list(nposSexTop20, nnegSexTop20),
  category.names = c("NeuN+", "NeuN-"),
  filename = '3_analysis/plots/SexDMPs_venn_diagramm.png',
  #fill = myCol,
  output=TRUE
)

