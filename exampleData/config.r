# Example cell sorted data

# project specific parmeters - must include "Sex" and/or "Cell_Type" and/or "Tissue_Type" to run optional checks below
projVar <- c("Sex", "Age", "Plate", "Batch", "Individual_ID", "N_Nuclei", "Pathology", "Cell_Type") 


# optional checks (note cannot have both CT and Tissue checks set to TRUE)
ctCheck = TRUE
sexCheck = TRUE
tissueCheck = FALSE


# calcMouseMethMetrics thresholds
intensThresh <- 2000
pFiltProbeThresh <- 0.05
pFiltSampleThresh <- 5
bsConThresh <- 90

# ctCheck thresholds
neunCT <- "NEUNpos"
predDistinctCT<-c("NEUNpos", "NEUNneg")
studentThres <- 1.5
nSDThres<-3