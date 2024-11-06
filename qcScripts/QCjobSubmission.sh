#!/bin/bash
#SBATCH --export=ALL #export all enviroment variables to the batch job
#SBATCH -p mrcq #submit to the serial queue
#SBATCH --time=24:00:00 ##maximum wall time for the job
#SBATCH -A Research_Project-MRC190311 #research project to submit under
#SBATCH --nodes=1 #specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --error=logFiles/%u/QCDNAdata.err # error file
#SBATCH --output=logFiles/%u/QCDNAdata.log # output file
#SBATCH --job-name=QCDNAdata


#------------------------------------------------------

# 1. command line argument input is your config.txt file

# e.g. to run sbatch QCjobSubmission.sh <filepath/to/projectFolder>

#-----------------------------------------------------


## print start date and time
echo Job started on:
  date -u
JOBNAME="QCDNAdata"

echo Job sumbitted from:
  echo $SLURM_SUBMIT_DIR

#load config.txt file
source $1 || exit 1

# create directry for QC and normalised data
NORMDIR=${DATADIR}/2_normalised

mkdir -p $NORMDIR
mkdir -p ${NORMDIR}/QC

# Move the user to the scripts folder
cd $SCRIPTSDIR


## load modules
module load $RVERS
module load Pandoc


# run the script to load idats to rgset and calculate QC metrics
Rscript calcMouseMethQCmetrics.r $DATADIR $REFDIR

# run cluster cell types script
Rscript cellTypeChecks.r $DATADIR

# create 1st stage QC report
Rscript -e "rmarkdown::render('QC.rmd', output_file='QC.html')" --args $DATADIR



# mv markdown report to correct location
mv QC.html ${NORMDIR}/QC



# create cell types check QC report
#Rscript -e "rmarkdown::render('cellTypeQC.rmd', output_file='cellTypeQC.html')" $DATADIR

# mv markdown report to correct location
mv cellTypeQC.html ${NORMDIR}/QC

# run normalisation script
Rscript normalisation.r $DATADIR


## print finish date and time
echo Job finished on:
  date -u

