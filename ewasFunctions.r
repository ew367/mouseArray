
# functions for EWAS

row <- meanBetas

runNeunPosPathologyEWAS<-function(row,QCmetrics){
  
  pathologyModel <- lm(row ~ QCmetrics$Group + QCmetrics$Pathology + QCmetrics$Age + QCmetrics$Sex + 
                         QCmetrics$Batch)
                
   nullPathology <- lm(row ~ QCmetrics$Group + QCmetrics$Age + QCmetrics$Sex + QCmetrics$Batch)
                
  
  # extract TG/WT main effect and cell proportion effect
  return(c(summary(modelLM)$coefficients["QCmetrics$PhenotypeSchizophrenia",c(1,2,4)],
           summary(modelLM)$coefficients["QCmetrics$Cell.Proportions",c(1,2,4)],
           
           # extract cell specific case control effect
           summary(nullCT)$coefficients["QCmetrics$PhenotypeSchizophrenia",c(1,2,4)]))
}



# group*age: DNAm ~ QCmetrics$Group + QCmetrics$Age*QCmetrics$Group + QCmetrics$Age + QCmetrics$Sex + QCmetrics$Batch
