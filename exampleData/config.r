
# filepaths
projDir <- ""
pheno <- file.path(projDir, "0_metadata/sampleSheet.csv")
manifest <- file.path(projDir, "0_metadata/MouseMethylation-12v1-0_A2.csv")
idatPath <- file.path(projDir, "1_raw/")
normDir <- file.path(projDir, "2_normalised/")
QCDir <- file.path(projDir, "2_normalised/QC")


# project specific parmeters - must include "Sex" and "Cell_Type" to run optional checks below
projVar <- c("Genotype", "Sex", "Plate", "Chip_Location") 


# optional checks
ctCheck = FALSE
sexCheck = TRUE


# calcMouseMethMetrics thresholds
intensThresh <- 2000
pFiltProbeThresh <- 0.05
pFiltSampleThresh <- 5
bsThresh <- 2
bsConThresh <- 90


