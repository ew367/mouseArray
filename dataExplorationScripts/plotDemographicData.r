##---------------------------------------------------------------------#
##
## Title: Plot Demographic data
##
## Purpose of script: plot sex distribution
##                    plot age distribution
##                    
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# this uses QC metrics so is only post QC data


#----------------------------------------------------------------------#
# LOAD DATA
#----------------------------------------------------------------------#

source("config.r")

normData<-file.path(projDir, "2_normalised/normalisedData.rdat")

load(normData)


#----------------------------------------------------------------------#
# PLOTS
#----------------------------------------------------------------------#

ggplot(QCmetrics, aes(x=Sex, fill = Group))+
  geom_histogram(stat="count", position="dodge")


ggplot(QCmetrics, aes(x=Age, fill = Group))+
  geom_histogram(position="dodge")


ggplot(QCmetrics, aes(x=Age, fill = Sex))+
  geom_histogram(position="dodge")


ggplot(QCmetrics, aes(x=Age, fill = Cell_Type))+
  geom_histogram(position="dodge")


