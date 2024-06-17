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

# 1. command line argument input is your project folder.

# e.g. to run sbatch QCjobSubmission.sh <filepath/to/projectFolder>

#-----------------------------------------------------


## print start date and time
echo Job started on:
  date -u
JOBNAME="QCDNAdata"

echo Job sumbitted from:
  echo $SLURM_SUBMIT_DIR

# Move the user to the project directory
cd $1


## load modules
module load Pandoc
module load R/4.2.1-foss-2022a

mkdir -p 2_normalised/QC

# run the script to load idats to rgset and calculate QC metrics
Rscript scripts/calcMouseMethQCmetrics.r 

# create 1st stage QC report
Rscript -e "rmarkdown::render('scripts/QC.rmd', output_file='QC.html')" --args $1

# mv markdown report to correct location
mv scripts/QC.html 2_normalised/QC

# run cluster cell types script
Rscript scripts/cellTypeChecks.r

# create cell types check QC report
Rscript -e "rmarkdown::render('scripts/cellTypeQC.rmd', output_file='cellTypeQC.html')" 

# mv markdown report to correct location
mv scripts/cellTypeQC.html 2_normalised/QC

# run normalisation script
Rscript scripts/normalisation.r


## print finish date and time
echo Job finished on:
  date -u

