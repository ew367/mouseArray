##---------------------------------------------------------------------#
##
## Title: Compare Mouse and Epic array probes
##
## Purpose of script: This script looks at the overlap between probes on
##                    on the mouse array and the Epic Array
##                    
##        
##                   
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# Notes
#----------------------------------------------------------------------#

# https://bioinformatics.stackexchange.com/questions/225/uppercase-vs-lowercase-letters-in-reference-genome

#----------------------------------------------------------------------#
# Import Data
#----------------------------------------------------------------------#

# load manifests

mmManifest <- fread("/lustre/projects/Research_Project-191406/cellSortedEWAS/0_metadata/MouseMethylation-12v1-0_A2.csv", skip=7, fill=TRUE, data.table=F)

mmGene <- fread("/lustre/projects/Research_Project-191406/cellSortedEWAS/0_metadata/MouseMethylation-12v1-0_A1_Annotation_Mus_musculus(1).csv", fill=TRUE, data.table=F)
colnames(mmGene)[1] <- "MouseName"
mmGene <- mmGene[!duplicated(mmGene$MouseName), c("MouseName", "Gene")]


hsManifest <- fread("/lustre/projects/Research_Project-MRC190311/references/EPICArray/MethylationEPIC_v-1-0_B5.csv", skip=7, fill=TRUE, data.table=F)



#----------------------------------------------------------------------#
# Compare mouse/human AlleleA probe sequences
#----------------------------------------------------------------------#

AlleleA_seq <- intersect(mmManifest$AlleleA_ProbeSeq, hsManifest$AlleleA_ProbeSeq)

nrow(mmManifest[mmManifest$AlleleA_ProbeSeq %in% AlleleA_seq,]) # 958 - doesn't give gene name

mmAseq <- mmManifest[mmManifest$AlleleA_ProbeSeq %in% AlleleA_seq[-1], c("Name", "AlleleA_ProbeSeq")]
colnames(mmAseq)[1] <- "MouseName"
mmAseq <- mmAseq[order(mmAseq$AlleleA_ProbeSeq),]

hsAseq <- hsManifest[hsManifest$AlleleA_ProbeSeq %in% AlleleA_seq[-1], c("Name", "AlleleA_ProbeSeq", "UCSC_RefGene_Name")]
colnames(hsAseq)[1] <- "EpicName"
hsAseq <- hsAseq[order(hsAseq$AlleleA_ProbeSeq),]


allAseq <- left_join(mmAseq, hsAseq)
allAseq <- left_join(allAseq, mmGene %>% dplyr::select(MouseName, Gene))  # just take 1st row/genename



#----------------------------------------------------------------------#
# Compare mouse/human AlleleB probe sequences
#----------------------------------------------------------------------#

# find sequences that match in mouse and epic
AlleleB_seq <- intersect(mmManifest$AlleleB_ProbeSeq, hsManifest$AlleleB_ProbeSeq) 

nrow(mmManifest[mmManifest$AlleleB_ProbeSeq %in% AlleleB_seq[-1],]) # remove "" from beginning - 31


# create object for mouse and human data with the name, probe seq, gene info

mmBseq <- mmManifest[mmManifest$AlleleB_ProbeSeq %in% AlleleB_seq[-1], c("Name", "AlleleB_ProbeSeq", "SourceSeq")]
colnames(mmBseq)[1] <- "MouseName"
colnames(mmBseq)[3] <- "MouseSourceSeq"
mmBseq <- mmBseq[order(mmBseq$AlleleB_ProbeSeq),]

hsBseq <- hsManifest[hsManifest$AlleleB_ProbeSeq %in% AlleleB_seq[-1], c("Name", "AlleleB_ProbeSeq", "UCSC_RefGene_Name", "SourceSeq")]
colnames(hsBseq)[1] <- "EpicName"
colnames(hsBseq)[4] <- "SourceSeq"

hsBseq <- hsBseq[order(hsBseq$AlleleB_ProbeSeq),]
  

allBseq <- left_join(mmBseq, hsBseq)
allBseq <- left_join(allBseq, mmGene %>% dplyr::select(MouseName, Gene))  # just take 1st row/genename
  


#----------------------------------------------------------------------#
# Compare mouse/human sourceSeq
#----------------------------------------------------------------------#

sourSeq <- intersect(toupper(mmManifest$SourceSeq), hsManifest$SourceSeq)

# mouse seq is in mix up upper and lower case (dom/rec alleles?)
# still no matches after converting to upper
# https://bioinformatics.stackexchange.com/questions/225/uppercase-vs-lowercase-letters-in-reference-genome 


#----------------------------------------------------------------------#
# Compare mouse AlleleA/AlleleB probe sequences
#----------------------------------------------------------------------#

intersect(mmManifest$AlleleA_ProbeSeq, mmManifest$AlleleB_ProbeSeq) #0 

intersect(hsManifest$AlleleA_ProbeSeq, mmManifest$AlleleB_ProbeSeq) #0

intersect(hsManifest$AlleleA_ProbeSeq, hsManifest$AlleleB_ProbeSeq) #0

intersect(hsManifest$AlleleB_ProbeSeq, mmManifest$AlleleA_ProbeSeq) # 0

