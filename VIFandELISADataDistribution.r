##---------------------------------------------------------------------#
##
## Title: Explore data structure
##
## Purpose of script: explore multicollinearity among predictors
##
##                    plot Elisa distributions                   
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# https://www.r-bloggers.com/2023/12/exploring-variance-inflation-factor-vif-in-r-a-practical-guide/ 
# https://rforpoliticalscience.com/2020/08/03/check-for-multicollinearity-with-the-car-package-in-r/
# above says "As a rule of thumb, a vif score over 5 is a problem. A score over 10 should be remedied (and you should consider dropping the problematic variable from the regression model or creating an index of all the closely related variables)."

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(car)
library(dplyr)
library(ggplot2)

#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#

load("2_normalised/normalisedData.rdat")

elisa <- read.csv("0_metadata/AB42ELISAconcentrations.csv", stringsAsFactors = F)
elisa$Individual_ID <- as.character(elisa$Sample_ID)
elisa$Pathology <- elisa$Mean_N3

QCmetrics <- left_join(QCmetrics, elisa %>% dplyr::select(Individual_ID, Pathology))

QCmetrics$DummyGroup <- ifelse(QCmetrics$Group == "TG", 1, 0)
QCmetrics$InteractionTerm <- QCmetrics$DummyGroup*QCmetrics$Age
#QCmetrics$Group <- as.factor(QCmetrics$Group)




#----------------------------------------------------------------------#
# Plot ELISA distribution
#----------------------------------------------------------------------#


#ggplot(QCmetrics, aes(x = Pathology, fill = Group))+
 # geom_histogram(position = "dodge")+
  #ggtitle("Elisa concentrations (pathology)")

pdf("3_analysis/plots/ElisaConcsByGroupLinePlot.pdf")
ggplot(QCmetrics, aes(x = Age, y = Pathology, colour = Group))+
  geom_line()+
  ggtitle("Elisa concentrations (pathology)")
dev.off()


#----------------------------------------------------------------------#
# DEFINE MODELS
#----------------------------------------------------------------------#

# age as continuous variable (increases possible degrees of freedom compared to factor)
# missing nuclei number at the moment - add once confirmed by wet lab team
# age*group here represents pathology

meanBetas <- colMeans(celltypeNormbeta) 

# group*age and pathology model
fullmodel <- lm(meanBetas ~ QCmetrics$Group + QCmetrics$Age*QCmetrics$Group + QCmetrics$Pathology + QCmetrics$Age + QCmetrics$Sex + QCmetrics$Batch) 

vif(fullmodel, type = "predictor")

# no interaction terms
nullmodel <- lm(meanBetas ~ QCmetrics$Group + QCmetrics$Pathology + QCmetrics$Age + QCmetrics$Sex + QCmetrics$Batch) 

vif(nullmodel)


# interaction model for Group*Age
neunModel<-lm(meanBetas ~ QCmetrics$Group + QCmetrics$Age*QCmetrics$Group + QCmetrics$Age + QCmetrics$Sex + QCmetrics$Batch) 

vif(neunModel, type = "predictor") # need predictor arg to include interaction terms

# null model for Group
nullneunModel <- lm(meanBetas ~ QCmetrics$Group + QCmetrics$Sex + QCmetrics$Batch +  QCmetrics$Age)

vif(nullneunModel)



# with dummy interaction variable

fullDumModel <- lm(meanBetas ~ QCmetrics$Group + QCmetrics$InteractionTerm + QCmetrics$Pathology + QCmetrics$Age + QCmetrics$Sex + QCmetrics$Batch) 

vif(fullDumModel)

noPathDumModel <- lm(meanBetas ~ QCmetrics$Group + QCmetrics$InteractionTerm + QCmetrics$Age + QCmetrics$Sex + QCmetrics$Batch) 

vif(noPathDumModel)

noIntDumModel <- lm(meanBetas ~ QCmetrics$Group + QCmetrics$Pathology + QCmetrics$Age + QCmetrics$Sex + QCmetrics$Batch) 

vif(noIntDumModel)




