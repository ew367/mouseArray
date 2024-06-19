##---------------------------------------------------------------------#
##
## Title: analyse blat output
##
## Purpose of script: analyse the blat output from comparing mouse manifest
##                    sourceSeq column to the human genome 
##                    
##                    This scripted is adapted from one written
##                    by AF
##
##                   
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# https://www.ensembl.org/info/website/upload/psl.html


#----------------------------------------------------------------------#
# Set up
#----------------------------------------------------------------------#

library(hiReadsProcessor) # for pslCols()
library(data.table)
library(dplyr)
library(ggplot2)

col = "AlleleA_ProbeSeq"

source("config.r")

'%ni%' <- Negate('%in%')

mmManifest <- fread("/lustre/projects/Research_Project-191406/cellSortedEWAS/0_metadata/MouseMethylation-12v1-0_A2.csv", skip=7, fill=TRUE, data.table=F)


#----------------------------------------------------------------------#
# Define Functions
#----------------------------------------------------------------------#

homology <- function(df){
	df$Perc <- (df$matches/50)*100
	df$Perc90 <- df$Perc>=90
	df$Perc95 <- df$Perc>95
	return(df)
}

annot <- function(df){
	man <- mmManifest[which(mmManifest$IlmnID %in% df$qName),]				# subset EPIC manifest to probes in df
	man <- man[match(df$qName, man$IlmnID),]									# match them
	identical(as.character(df$qName), as.character(man$probeID))					# check identical
	df.anno <- cbind(df, man[,c('CHR','MAPINFO','AlleleA_ProbeSeq', 'SourceSeq','Strand')])	# add chromosome, start, end and Gene to df
	# check if chr name matches between EPIC and df as well as probe position in chr
	#df.anno$posInRegion <- df.anno$tName==df.anno$chr & (df.anno$start >= df.anno$tStart & df.anno$start <= df.anno$tEnd) | (df.anno$end >= df.anno$tStart & df.anno$end <= df.anno$tEnd)	
	return(df.anno)
}


#----------------------------------------------------------------------#
# Load Blat Results
#----------------------------------------------------------------------#

blat <- read.table(paste0("0_metadata/allProbeSeqs_",col,".psl"), sep="\t", skip=5)

colnames(blat) <- names(pslCols()) # rename columns to more standard column names
blat$qName <- as.character(blat$qName)
blat$tName <- gsub("chr","",blat$tName)
blat$Span <- blat$tEnd-blat$tStart


# calculate percentage homology
blat <- homology(blat)

# annotate with EPIC manifest info and check chr and position agreement
blat <- annot(blat)

hist(blat$Perc)

############################################################


blat$Mismatch <- blat$posInRegion==FALSE & blat$Perc90==TRUE #58424. 58424
blat <- blat[order(blat$Perc, decreasing=T),]

#3. Find probes with single or multiple hits ===============================================================================================================
# hit1 <- names(tb)[tb==1] #791578
# cross <- names(tb)[tb>1] #29854

# hit1.blat <- blat[which(blat$qName %in% hit1),] #791578
# cross.blat <- blat[which(blat$qName %in% cross),] #294,320


# hit1.mismatch <- hit1.blat[which(hit1.blat$Mismatch==TRUE),] #22




#3a. Autosomal probes matching to X/Y chromosomes ==========================================================================================================
blat.mis <- blat[which(blat$Mismatch==TRUE),] #58424
write.csv(blat.mis, file="BLAT_cross_perc90.csv")
probes.mis <- unique(blat.mis$qName) #13214
write.csv(probes.mis, file="BLAT_cross_probelist_perc90.csv")
mis.aut <- blat.mis[-which(blat.mis$chrm %in% c('X','Y')),] #57393 - autosomal according to EPIC manifest
mis.aut.cross <- mis.aut[which(mis.aut$tName %in% c('X','Y')),] #3191 - non-autosomal according to BLAT annotation
probes.cross <- unique(mis.aut.cross$qName) #1065 - list of probes with evidence of potential to cross-hybridise to X or Y chr.

# Annotate with number of times probe appears in...
# all high-confidence BLAT hits
blat.conf <- blat[which(blat$Perc90==TRUE),]					# BLAT hits with >90% homology
tb <- table(blat.conf$qName)									# how many times probe appears in blat
tb2 <- tb[which(names(tb) %in% mis.aut.cross$qName)]			# limit table to probes present in mis.aut.cross
tb2 <- tb2[match(mis.aut.cross$qName, names(tb2))]				# match order of table to mis.aut.cross
identical(names(tb2),mis.aut.cross$qName)						# TRUE
mis.aut.cross$nHits <- as.numeric(tb2)							# how many BLAT hits per probe in mis.aut.cross
mis.aut.cross$hitType <- rep('Multiple', nrow(mis.aut.cross))	# multiple = this hit is one of many hits for this probe
mis.aut.cross$hitType[which(mis.aut.cross$nHits==1)] <- 'Single'# single = this is the only hit for the probe

# cross-hybridising autosomal-to-sex probe list
h <- table(mis.aut.cross$qName)									# how many times a probe appears in mis.aut.cross
h <- h[match(mis.aut.cross$qName, names(h))]					# match order of table to mis.aut.cross
identical(names(h),mis.aut.cross$qName)							# TRUE
mis.aut.cross$nXYHits <- as.numeric(h)							# how many hits per probe in mis.aut.cross

write.csv(mis.aut.cross, file="BLAT_crossXY_perc90.csv")		# 3191
write.csv(unique(mis.aut.cross$qName), file="BLAT_crossXY_probelist_perc90.csv") #1065




chrCol <- 'chrm'
df <- mis.aut.cross %>% distinct(qName, .keep_all = TRUE)		# unique probes


df[,chrCol][which(df[,chrCol]=='X')] <- rep('23',length(which(df[,chrCol]=='X')))
df[,chrCol][which(df[,chrCol]=='Y')] <- rep('24',length(which(df[,chrCol]=='Y')))
if(any(df[,chrCol] %ni% c(1:24))){ df <- df[-which(df[,chrCol] %ni% c(1:24)),] }

actualCHR <- table(df[,chrCol])
actualCHR <- actualCHR[order(as.numeric(names(actualCHR)))]
actualCHR <- actualCHR[match(1:24,names(actualCHR))]
names(actualCHR) <- 1:24
actualCHR[which(is.na(actualCHR))] <- 0

pdf("BLAT_crossXY_chrBarplot.pdf", width=10, height=7)
par(mar=c(8, 5, 4.1, 2.1))
barplot(actualCHR, col='lightgreen', xlab='Chromosome', ylab='Number of probes', ylim=c(0,40), main='Probes cross-hybr to X or Y Chr', cex.axis=1.5, cex.names=1, cex.lab=1.5, cex.main=2)
dev.off()


# compare to autosomal sex DMPs
res <- data.frame(readRDS("/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/Bulk/3_analysis/linearRegression/Age/ageReg_fetalBrain_EX3_23pcw_annotAllCols_filtered.rds"))
res <- res[which(res$P.Sex<9e-8),] #14727
res <- res[-which(res$CHR %in% c('X','')),] #1282
res$IlmnID <- as.character(res$IlmnID)
mis.res <- mis.aut.cross
mis.res$Significance <- rep('Not sig',nrow(mis.res))
mis.res$Significance[which(mis.res$qName %in% res$IlmnID)] <- 'Sig'

res$Mismatch <- rep('Not',nrow(res))
res$Mismatch[which(res$IlmnID %in% mis.aut.cross$qName)] <- 'Mismatch' # 48
res.ord <- res[order(res$P.Sex),]

df <- mis.res %>% distinct(qName, .keep_all = TRUE)
chrCol <- 'chrm'
sigCol <- 'Significance'

df[,chrCol][which(df[,chrCol]=='X')] <- rep('23',length(which(df[,chrCol]=='X')))
df[,chrCol][which(df[,chrCol]=='Y')] <- rep('24',length(which(df[,chrCol]=='Y')))
if(any(df[,chrCol] %ni% c(1:24))){ df <- df[-which(df[,chrCol] %ni% c(1:24)),] }

#actualCHR <- table(df[,chrCol],df[,sigCol])
actualCHR <- table(df[,sigCol],df[,chrCol])
actualCHR <- actualCHR[,order(as.numeric(colnames(actualCHR)))]
actualCHR <- actualCHR[,match(1:24,colnames(actualCHR))]
colnames(actualCHR) <- 1:24
actualCHR[which(is.na(actualCHR))] <- 0

pdf("BLAT_crossXY_chrBarplot_colSigFetalAutSexDMP.pdf", width=10, height=7)
par(mar=c(8, 5, 4.1, 2.1))
barplot(actualCHR, col=c('lightgrey','lightblue'), xlab='Chromosome', ylab='Number of probes', ylim=c(0,40), main='Probes cross-hybr to X or Y Chr', cex.axis=1.5, cex.names=1, cex.lab=1.5, cex.main=2)
legend("topright", legend=c('Not an autosomal sex DMP', 'Autosomal sex DMP'), col=c('lightgrey','lightblue'), pch=15)
dev.off()






#####
setwd("/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/Bulk/3_analysis/linearRegression/Sex/crossHybr/")
m <- read.csv("BLAT_crossXY_perc90.csv")
x <- unique(m$qName[which(m$tName=='X')]) #874
y <- unique(m$qName[which(m$tName=='Y')]) #334

res.df <- data.frame(readRDS("/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/Bulk/3_analysis/linearRegression/Age/ageReg_fetalBrain_EX3_23pcw_annotAllCols_filtered.rds"))
res.df$IlmnID <- as.character(res.df$IlmnID)
res.sig <- res.df[which(res.df$P.Sex<9e-8),] #14727 - significant probes
res.sig <- res.sig[-which(is.na(res.sig$IlmnID)),]
res.sig$chrType <- rep('X',nrow(res.sig))
res.sig$chrType[which(res.sig$CHR != 'X')] <- 'Autosomal'
res.sig$mismatchTo <- rep('NA',nrow(res.sig))
res.sig$mismatchTo[which(res.sig$IlmnID %in% x)] <- 'X'
res.sig$mismatchTo[which(res.sig$IlmnID %in% y)] <- 'Y'
res.sig$mismatchTo[which(res.sig$IlmnID %in% x & res.sig$IlmnID %in% y)] <- 'X&Y'

   # NA     X   X&Y     Y
# 14675    36     2    12

res.plot <- res.sig[-which(res.sig$chrType=='X' & res.sig$mismatchTo=='X'),]
res.plot$plotCol <- rep('',nrow(res.plot))
res.plot$plotCol[which(res.plot$chrType=='Autosomal' & res.plot$mismatchTo=='NA')] <- 'Autosomal DMP\nNot mismatch'
res.plot$plotCol[which(res.plot$chrType=='Autosomal' & res.plot$mismatchTo=='X')] <- 'Autosomal DMP\nMismatch to X'
res.plot$plotCol[which(res.plot$chrType=='Autosomal' & res.plot$mismatchTo=='Y')] <- 'Autosomal DMP\nMismatch to Y'
res.plot$plotCol[which(res.plot$chrType=='Autosomal' & res.plot$mismatchTo=='X&Y')] <- 'Autosomal DMP\nMismatch to X & Y'
res.plot$plotCol[which(res.plot$chrType=='X')] <- 'X Chr DMP'

res.plot$plotCol <- factor(res.plot$plotCol, levels=c('Autosomal DMP\nNot mismatch', 'Autosomal DMP\nMismatch to X', 'Autosomal DMP\nMismatch to Y', 'Autosomal DMP\nMismatch to X & Y', 'X Chr DMP'))

pdf("Sex_ES_comparison.pdf", width=12, height=8)
ggplot(res.plot, aes(x=plotCol, y=(Beta.Sex*100))) + 
	geom_violin()+
	xlab('')+
	ylab('Sex effect size (%)')+
	geom_jitter(height = 0, width = 0.1)+
	theme_minimal()+
	theme(axis.text=element_text(size=19), axis.title=element_text(size=22))+
	theme(plot.title = element_text(size=23))
dev.off()


PathToBetas <- "/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/Bulk/2_normalised/"
load(paste0(PathToBetas, "fetalBulk_EX3_23pcw_n91.rdat"))



# investigating the autosomal DMPS mismatching to X
b <- res.plot[which(res.plot$plotCol=='Autosomal DMP\nMismatch to X'),c('IlmnID','Beta.Sex','P.Sex')]
PathToBetas <- "/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/Bulk/2_normalised/"
load(paste0(PathToBetas, "fetalBulk_EX3_23pcw_n91.rdat"))
betas <- betas[b$IlmnID,]
betas <- betas[match(b$IlmnID,rownames(betas)),]
identical(rownames(betas),b$IlmnID)
pdf("checkingAutosomalMisToX.pdf", width=12,height=12)
par(mfrow=c(3,3))
for(i in 1:nrow(betas)){
boxplot(betas[i,] ~ pheno$Sex, main=paste0(rownames(betas)[i], '   ', signif(b$Beta.Sex[i],digits=3),"  p=", signif(b$P.Sex[i], digits=3)))
}
dev.off()

############################################################################################################################################################################################################


# # mis.probes <- c()
# # for(each in unique(cross.blat$qName)){
	# # probe <- cross.blat[which(cross.blat$qName==each),]
	# # if(any(probe$Mismatch==TRUE)){
		# # mis.probes <- c(mis.probes, each)		# any hit for a probe is a mismatch
	# # }
# # }

# probes that have ANY high-confidence mismatch
mis.probes <- readRDS("cross_misProbes.rds")
cross.misprobes <- cross.blat[which(cross.blat$qName %in% mis.probes),] #98125
k <- cross.misprobes[cross.misprobes$Mismatch==TRUE,] #26286 hits. 
length(unique(k$qName)) # 7162 unique probes. 


# cross.mis2 <- cross.misprobes[which(cross.misprobes$Mismatch==TRUE),] #219983. 26286
# cross.mis3 <- cross.mis2[-which(cross.mis2$Perc==100 & cross.mis2$tName==cross.mis2$chrm),] #216106. Span: 48 406412. 162014 hits span<51. (24146)
# cross.mis4 <- cross.mis3[which(cross.mis3$Span<51),] #162014. 23208 probes
#cross.misprobes <- cross.misprobes %>% distinct(qName, .keep_all = TRUE) #27528


to.remove <- c(hit1.mismatch$qName, unique(k$qName))
#write.csv(to.remove, file="crossHybrProbes_allPossProbes.txt") #7184 unique probes
#to.remove <- read.csv("crossHybrProbes_allPossProbes.txt")
#to.remove <- to.remove$x
blat.mismatch <- blat[which(blat$qName %in% to.remove),] #98147 hits



#4. Match to autosomal sex DMPs ============================================================================================================================
res <- data.frame(readRDS("/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/Bulk/3_analysis/linearRegression/Age/ageReg_fetalBrain_EX3_23pcw_annotAllCols_filtered.rds"))
res <- res[which(res$P.Sex<9e-8),] #14727
res <- res[-which(res$CHR %in% c('X','')),] #1282
res$IlmnID <- as.character(res$IlmnID)

res.cross <- res[which(rownames(res) %in% to.remove),] #37 unique sex DMPs are mismatches
cross.res <- blat.mismatch[which(blat.mismatch$qName %in% rownames(res.cross)),] #37 probes equates to 197 mismatch hits
res.cross <- res.cross[match(cross.res$qName, rownames(res.cross)),]								
cross.res <- cbind(cross.res, res.cross[,c('CHR','MAPINFO','Gene','SourceSeq','Strand')])	# add chromosome, start, end and Gene to df

sex.mismatch <- cross.res[which(cross.res$Mismatch==TRUE),] #67 out of 197 are the culprits. All have Perc95==TRUE. 
#saveRDS(sex.mismatch, file="fetalAutSexDMP_mismatches.rds")

par(mar=c(8, 4.1, 4.1, 2.1))
# number of hits per sex DMP mismatch
barplot(table(sex.mismatch$qName), las=2, ylab='Number of hits', main='Distribution of 37 autosomal sex DMP mismatch probes', col='lightblue', cex.names=0.8)
# distribution of percentage homologies for the sex DMP mismatches
barplot(table(sex.mismatch$Perc), las=2, ylab='Number of hits', main='Percentage homologies for \n37 autosomal sex DMP mismatch probes', col='lightblue')
# number of hits split by percentage homology
barplot(table(sex.mismatch$Perc, sex.mismatch$qName), las=2, ylab='Number of hits', main='Number of hits coloured by percentage homology', col=c('green','royalblue','red'), cex.names=0.8)
legend("topright", legend=rownames(table(sex.mismatch$Perc, sex.mismatch$qName)), col=c('green','royalblue','red'), pch=15, title='% homology', bty='n')



res.blat <- res[which(rownames(res) %in% blat$qName),]
blat.res <- blat[which(blat$qName %in% rownames(res.blat)),]
#res.blat2 <- res.blat[match(blat.res$qName, rownames(res.blat)),]
#blat.res <- cbind(blat.res, res.blat2[,c('CHR','MAPINFO','Gene','SourceSeq','Strand')])

epic.blat <- epicManifest[which(epicManifest$probeID %in% blat$qName),]

#df <- res.blat #for 'EPIC_autosomalSexDMPs.pdf'. This is using EPIC b4 
df <-  epic.blat[-which(epic.blat$probeID %in% rownames(res)),] #820152. For 'EPIC_NOTautosomalSexDMPs.pdf'
#df <- sex.mismatch
#chrCol <- 'CHR'
chrCol <- 'chrm'


df[,chrCol][which(df[,chrCol]=='X')] <- rep('23',length(which(df[,chrCol]=='X')))
df[,chrCol][which(df[,chrCol]=='Y')] <- rep('24',length(which(df[,chrCol]=='Y')))
if(any(df[,chrCol] %ni% c(1:24))){ df <- df[-which(df[,chrCol] %ni% c(1:24)),] }

actualCHR <- table(df[,chrCol])
actualCHR <- actualCHR[order(as.numeric(names(actualCHR)))]
actualCHR <- actualCHR[match(1:24,names(actualCHR))]
names(actualCHR) <- 1:24
actualCHR[which(is.na(actualCHR))] <- 0

pdf("EPIC_autosomalSexDMPs.pdf", width=10, height=7)
barplot(actualCHR, col='lightgreen', xlab='EPIC CHR', ylab='Number of hits', ylim=c(0,150))
dev.off()

pdf("EPIC_NOTautosomalSexDMPs.pdf", width=10, height=7)
barplot(actualCHR, col='lightgreen', xlab='EPIC CHR', ylab='Number of hits', ylim=c(0,80000))
dev.off()




#df <- blat.res
df <- blat[-which(blat$qName %in% rownames(res)),] #1083370
chrCol <- 'tName'

#df <- blat[-which(blat$qName %in% rownames(res)),]
#df <- blat.mismatch[-which(blat.mismatch$qName %in% sex.mismatch$qName),]
#chrCol <- 'chrm'

df[,chrCol][which(df[,chrCol]=='X')] <- rep('23',length(which(df[,chrCol]=='X')))
df[,chrCol][which(df[,chrCol]=='Y')] <- rep('24',length(which(df[,chrCol]=='Y')))
if(any(df[,chrCol] %ni% c(1:24))){ df <- df[-which(df[,chrCol] %ni% c(1:24)),] }

predCHR <- table(df[,chrCol])
names(predCHR)[which(names(predCHR) %in% c('X','Y'))] <- c('23','24') #change X and Y to 23 and 24 to allow for ordering numerically
predCHR <- predCHR[order(as.numeric(names(predCHR)))]
predCHR <- predCHR[match(1:24,names(predCHR))]
names(predCHR) <- 1:24
predCHR[which(is.na(predCHR))] <- 0

pdf("BLAT_autosomalSexDMP.pdf", width=10, height=7)
barplot(predCHR, col='lightblue', xlab='BLAT CHR', ylab='Number of hits', ylim=c(0,250))
dev.off()

pdf("BLAT_NOTautosomalSexDMP.pdf", width=10, height=7)
barplot(predCHR, col='lightblue', xlab='BLAT CHR', ylab='Number of hits', ylim=c(0,100000))
dev.off()

pdf("EPICnotSexDMP_crossHybr.pdf", width=8, height=8)
par(mfrow=c(2,1))
barplot(actualCHR, col='lightgreen', xlab='EPIC CHR', ylab='Number of hits', ylim=c(0,max(max(actualCHR),max(predCHR))))
barplot(predCHR, col='lightblue', xlab='BLAT CHR', ylab='Number of hits', ylim=c(0,max(max(actualCHR),max(predCHR))))
dev.off()


# BLAT CHR split by percentage homology
predCHR.perc <- table(df$Perc, df$tName)
colnames(predCHR.perc)[which(colnames(predCHR.perc) %in% c('X','Y'))] <- c('23','24') #change X and Y to 23 and 24 to allow for ordering numerically
predCHR.perc <- predCHR.perc[,order(as.numeric(colnames(predCHR.perc)))]
predCHR.perc <- predCHR.perc[,match(1:24,colnames(predCHR.perc))]
colnames(predCHR.perc) <- 1:24
z <- function(x){ if(all(is.na(x))){ x <- rep(0,length(x)) }else{return(x)} }
g <- apply(predCHR.perc, 2, z)
rownames(g) <- rownames(predCHR.perc)
predCHR.perc <- g

# par(mfrow=c(1,1))
# barplot(predCHR.perc, xlab='BLAT CHR', ylab='Number of hits', col=c('green','royalblue','red'), ylim=c(0,25))
# legend("topleft", legend=rownames(predCHR.perc), col=c('green','royalblue','red'), pch=15, title='% homology', bty='n')


red <- length(which(as.numeric(rownames(predCHR.perc)) <= 80))
orange <- length(which(as.numeric(rownames(predCHR.perc)) > 80 & as.numeric(rownames(predCHR.perc)) <=90))
yellow <- length(which(as.numeric(rownames(predCHR.perc)) > 90 & as.numeric(rownames(predCHR.perc)) <=94))
green <- length(which(as.numeric(rownames(predCHR.perc)) > 94 & as.numeric(rownames(predCHR.perc)) <=98))
darkgreen <- length(which(as.numeric(rownames(predCHR.perc)) > 98))
cols <- c( rep('red',red), rep('orange',orange), rep('yellow',yellow), rep('green',green), rep('darkgreen',darkgreen))
leglabel <- c( '<80%', '80-90%', '90-94%', '94-98%', '100%')

pdf("EPICnotSexDMP_crossHybr_homology.pdf", width=8, height=8)
par(mfrow=c(1,1))
barplot(predCHR.perc, xlab='BLAT CHR', ylab='Number of hits', col=cols, ylim=c(0,max(max(actualCHR),max(predCHR))))
legend("topright", legend=leglabel, col=unique(cols), pch=15, title='Homology', bty='n')
dev.off()


pdf("AutosomalSexDMP_crossHybr_homology.pdf", width=8, height=8)
par(mfrow=c(1,1))
barplot(predCHR.perc, xlab='BLAT CHR', ylab='Number of hits', col=cols, ylim=c(0,250))
legend("topright", legend=leglabel, col=unique(cols), pch=15, title='Homology', bty='n')
dev.off()

pdf("mismatch_AutosomalSexDMP_crossHybr_homology.pdf", width=8, height=8)
par(mfrow=c(1,1))
barplot(predCHR.perc, xlab='BLAT CHR', ylab='Number of hits', col=cols, ylim=c(0,22))
legend("topleft", legend=c('94-98%', '100%'), col=unique(cols), pch=15, title='Homology', bty='n')
dev.off()


res.blat <- res[match(blat$qName, res$IlmnID),] #18944
identical(blat$qName, res2$IlmnID)
