
# filepaths
manifest <- ("0_metadata/MouseMethylation-12v1-0_A2.csv")

# project specific parmeters - must include "Sex" and "Cell_Type" to run optional checks below
projVar <- c("Genotype", "Chip_Location", "Sex", "Age") 

# empty <- "" # basenames of any sample to exclude here


# optional checks
ctCheck = FALSE
sexCheck = FALSE


# calcMouseMethMetrics thresholds
intensThresh <- 2000
pFiltProbeThresh <- 0.05
pFiltSampleThresh <- 5
bsThresh <- 2
bsConThresh <- 90


