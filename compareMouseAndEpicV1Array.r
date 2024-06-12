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
# Import Data
#----------------------------------------------------------------------#

# load manifests

mmManifest <- fread("/lustre/projects/Research_Project-191406/cellSortedEWAS/0_metadata/MouseMethylation-12v1-0_A2.csv", skip=7, fill=TRUE, data.table=F)

mmGene <- fread("/lustre/projects/Research_Project-191406/cellSortedEWAS/0_metadata/MouseMethylation-12v1-0_A1_Annotation_Mus_musculus(1).csv", fill=TRUE, data.table=F)
colnames(mmGene)[1] <- "MouseName"
mmGene <- mmGene[!duplicated(mmGene$MouseName), c("MouseName", "Gene")]


hsManifest <- fread("/lustre/projects/Research_Project-MRC190311/references/EPICArray/MethylationEPIC_v-1-0_B5.csv", skip=7, fill=TRUE, data.table=F)



#----------------------------------------------------------------------#
# Compare AlleleA probe sequences
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
# Compare AlleleB probe sequences
#----------------------------------------------------------------------#

# find sequences that match in mouse and epic
AlleleB_seq <- intersect(mmManifest$AlleleB_ProbeSeq, hsManifest$AlleleB_ProbeSeq) 

nrow(mmManifest[mmManifest$AlleleB_ProbeSeq %in% AlleleB_seq[-1],]) # remove "" from beginning - 31


# create object for mouse and human data with the name, probe seq, gene info

mmBseq <- mmManifest[mmManifest$AlleleB_ProbeSeq %in% AlleleB_seq[-1], c("Name", "AlleleB_ProbeSeq")]
colnames(mmBseq)[1] <- "MouseName"
mmBseq <- mmBseq[order(mmBseq$AlleleB_ProbeSeq),]

hsBseq <- hsManifest[hsManifest$AlleleB_ProbeSeq %in% AlleleB_seq[-1], c("Name", "AlleleB_ProbeSeq", "UCSC_RefGene_Name")]
colnames(hsBseq)[1] <- "EpicName"
hsBseq <- hsBseq[order(hsBseq$AlleleB_ProbeSeq),]
  

allBseq <- left_join(mmBseq, hsBseq)
allBseq <- left_join(allBseq, mmGene %>% dplyr::select(MouseName, Gene))  # just take 1st row/genename
  


