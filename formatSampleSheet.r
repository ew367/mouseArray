##---------------------------------------------------------------------#
##
## Title: Merge SampleSheets
##
## Purpose of script: The cell sorted data was originally run in 2 batches:
##                    Apr2022 (24 samples) and
##                    June2022 (180 samples)
##                    This script merges the sampleSheet files and saves
##
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

## The number of neuclei is missing from the n180 file. NAs are used as a
## place holder until it becomes available

## NB need to be in the cellSortedEWAS dir to run

## The sex info was updated manually in EXCEL according to doc from
## Nick (re convo with IC and SP)

## mouse 56 M -> F
## mice 73-75 -> M

#----------------------------------------------------------------------#
# DEFINE PARAMETERS ETC.
#----------------------------------------------------------------------#

library(dplyr)

source("config.r")

n24pheno <- file.path(projDir, "0_metadata/Pheno_N24.csv")
n180pheno <- file.path(projDir, "0_metadata/PhenoALL.csv")

#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#

setwd(projDir)

n24 <- read.csv(n24, stringsAsFactors = F)[1:24,]
n180 <- read.csv(n180, stringsAsFactors = F)

# standardise column names

#colnames(n24)
#[1] "Chip_ID"       "Chip_Position" "Basemane"      "Full_ID"       "Cell_Type"  
#[6] "N.nuclei"      "sample_ID"     "Group"         "Age.months."   "Sex"        
#[11] "Model"         "brain.region" 


#colnames(n180)
#[1] "Plate"           "Plate_Location"  "Chip_Number"     "Chip_ID"        
#[5] "Chip_Position"   "Basename"        "Sample_ID"       "Animal_ID"      
#[9] "Nuclei_Fraction" "Group"           "Age"             "Sex"  


colsNeeded <- c("Basename","Chip_ID","Chip_Position","Plate","Batch",
                "Individual_ID","Sample_ID","Cell_Type","N_Nuclei",
                "Group","Age","Sex")


###############
# n24
###############

# add in new cols
n24$Plate <- 1
n24$Batch <- "n24"
n24$Individual_ID <- n24$sample_ID
n24$Sample_ID <- paste(n24$Individual_ID, n24$Group, n24$Age.months., n24$Cell_Type, sep="_") 

# rename cols
lookup <- c(Basename = "Basemane", N_Nuclei = "N.nuclei", Age = "Age.months.")
n24 <- rename(n24, all_of(lookup))

# select and arrange cols
n24 <- n24 %>% dplyr::select(all_of(colsNeeded))


###############
# n180
###############

# add in new cols
n180$Batch <- "n180"
n180$Individual_ID <- n180$Animal_ID
n180$Cell_Type <- toupper(n180$Nuclei_Fraction)
n180$N_Nuclei <- NA

n180 <- n180 %>% dplyr::select(all_of(colsNeeded))


#----------------------------------------------------------------------#
# COMBINE AND SAVE
#----------------------------------------------------------------------#

sampleSheet <- rbind(n24,n180)

write.csv(sampleSheet, "0_metadata/sampleSheet.csv", row.names = F)
