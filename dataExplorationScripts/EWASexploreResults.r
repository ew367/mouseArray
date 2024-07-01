##---------------------------------------------------------------------#
##
## Title: explore EWAS results
##
## Purpose of script: explore output from linear regression within 
##                    each cell type
##                    
##
##                   
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#


#----------------------------------------------------------------------#
# LOAD PACKAGES etc
#----------------------------------------------------------------------#

library(dplyr)

source("config.r")

# load manifest
man <- fread(manifest, skip=7, fill=TRUE, data.table=F)

load("3_analysis/results/EWASout.rdat")

outtab <- as.data.frame(outtab)


#----------------------------------------------------------------------#
# Sex specific DMPs
#----------------------------------------------------------------------#

# sex chrs were removed in normalisation script

#sexDmps <- outtab %>% dplyr::select(contains('Sex')) %>% filter(nullSexM_P < 0.05/nrow(outtab))

#sexDmps$IlmnID <- row.names(sexDmps)

# add chr info
#sexDmps <- left_join(sexDmps, man %>% dplyr::select(IlmnID, CHR))

#table(sexDmps$CHR)



#----------------------------------------------------------------------#
# Age specific DMPs
#----------------------------------------------------------------------#

nullAgeDmps <- outtab %>% dplyr::select(contains('Age')) %>% filter(nullAge_P < 0.05/nrow(outtab))

intAgeDmps <- outtab %>% dplyr::select(contains('Age')) %>% filter(Age_P < 0.05/nrow(outtab))

inboth <- intAgeDmps[intersect(row.names(intAgeDmps), row.names(nullAgeDmps)),] # 3334

sum(inboth$nullAge_coeff > inboth$Age_coeff) # [1] 990 
# int term coef seems to be bigger, so interaction term is doing something, just not passing threshold



# which ones aren't in both?
length(setdiff(row.names(intAgeDmps), row.names(nullAgeDmps))) # 16 - many are v near bonfP threshold

#intAgeDmps[setdiff(row.names(intAgeDmps), row.names(nullAgeDmps)),]




#----------------------------------------------------------------------#
# APP localised DMPs
#----------------------------------------------------------------------#
