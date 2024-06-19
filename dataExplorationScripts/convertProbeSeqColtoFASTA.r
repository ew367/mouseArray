##---------------------------------------------------------------------#
##
## Title: Convert probe Sequences to fasta format
##
## Purpose of script: Convert the inputted col to fasta format for use
##                      with blat on the cmd line
##
##
##
#---------------------------------------------------------------------#

#----------------------------------------------------------------------
# Notes
#----------------------------------------------------------------------

# code adapted from a script written by AF

# use Rscript to run while specifying col to use


#----------------------------------------------------------------------
# Define Parameters and set up
#----------------------------------------------------------------------

library(data.table)

args<-commandArgs(trailingOnly = TRUE)
#args <- "SourceSeq"

source("config.r")


# load manifest
man <- fread(manifest, skip=7, fill=TRUE, data.table=F)

sqCol <- which(colnames(man)==args[1])
probeCol <- which(colnames(man)=='IlmnID')

outFile <- paste0("0_metadata/allProbeSeqs_", args[1], "_fasta4.fa")


#----------------------------------------------------------------------
# Define function
#----------------------------------------------------------------------

fasta <- function(x){
	sq <- as.character(x[sqCol])
	probe <- as.character(x[probeCol])
	name <- paste0(">",probe)
	comb <- capture.output(cat(name, sq, sep='\n')) #using capture.output() because cat() prints to screen
	return(comb)
}


#----------------------------------------------------------------------
# Create and save fasta file
#----------------------------------------------------------------------

sq2 <- apply(as.data.frame(man[-which(man[,args[1]] == ""),]), 1, fasta) 

write.table(as.vector(sq2), file=outFile, row.names=F, col.names=F, quote=F)

