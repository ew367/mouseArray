# mouseArray
QC and analysis of methylation array data run on the Ilumina Infinium Mouse Methylation BeadChip \
(https://emea.illumina.com/products/by-type/microarray-kits/infinium-mouse-methylation.html)
\
\
This repository contains workflows for the pre-processing of bulk/cell sorted/multiple tissue type data. They are developed for use with a HPC system with the SLURM job scheduler.


# Guidance for using the pipeline scripts

PREQUISITES: 
* A project folder containing the following sub directories:
  * 0_metadata
  * 1_raw
  * 2_normalised
  * logFiles
  
* A config.txt file in the project folder containing file paths specific to the project (see exampleData/config.txt for an example)
* A config.r file containg threshold values etc specific to your project (see exampleData/config.r for an example)

* A file named sampleSheet.csv in the 0_metadata folder which contains the following columns, along with any other (optional) associated phenotype data (see exampleData/sampleSheet.csv):
    "Basename","Chip_ID","Chip_Position", "Individual_ID","Sample_ID".
    
  To run cell type checks a "Cell_Type" column must be included.
  To run sex check a "Sex" column must be included.
  To run tissue check a "Tissue_Type" column must be included.

* A manifest file (e.g. MouseMethylation-12v1-0_A2.csv) in the REFDIR directory (specified in the config.txt file)
    
* raw idat files in the 1_raw folder



# QUALITY CONTROL STAGE

OUTPUT: 
* Quality control objects, metrics and html reports are output to 2_normalised/QCmetrics


The QCjobSubmission.sh script automates the quality control pipeline and can be run on the command line with:

sbatch QCjobSubmission.sh <filepath/to/projectFolder> (slurm job submissions)
or
sh QCjobSubmission.sh <filepath/to/projectFolder> (local machine)


This script will execute:
* Rscript calcMouseMethQCmetrics.r $DATADIR $REFDIR
* Rscript cellTypeChecks.r $DATADIR $REFDIR (if requested in the config.r file)
* Rscript -e "rmarkdown::render('QC.rmd', output_file='QC.html')" --args $1
* Rscript cellTypeChecks.r $DATADIR $REFDIR
* Rscript -e "rmarkdown::render('QC.rmd', output_file='QC.html')" --args $DATADIR $REFDIR
* Rscript normalisation.r $DATADIR $REFDIR

arguments are provided in the config.txt file.





