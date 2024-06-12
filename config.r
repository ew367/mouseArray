
# filepaths
projDir <- "/lustre/projects/Research_Project-191406/cellSortedEWAS"
pheno <- file.path(projDir, "0_metadata/sampleSheet.csv")
manifest <- file.path(projDir, "0_metadata/MouseMethylation-12v1-0_A2.csv")
idatPath <- file.path(projDir, "1_raw/")
normDir <- file.path(projDir, "2_normalised/")
QCDir <- file.path(projDir, "2_normalised/QC")


# dataTypes
array <- "mouse"
cellSorted <- TRUE


MC <- "205333480089_R04C01"   # fully methylated control
empty <- c("205333480089_R05C02", "205333480089_R06C02") # empty wells


# calcMouseMethMetrics thresholds
intensThresh <- 2000
pFiltProbeThresh <- 0.05
pFiltSampleThresh <- 5
bsThresh <- 2
bsConThresh <- 90

neunCT <- "NEUNpos"
predDistinctCT<-c("NEUNpos", "NEUNneg")

studentThres <- 1.5
nSDThres<-3

