##---------------------------------------------------------------------#
##
## Title: Normalise the data for downstream analyses
##
## Purpose of script: Filter the data using the QC metrics and cell type
##                    checks
##
##                    Remove flagged probes
##
##                    Normalise the data within celltype
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# holding code from NC


# Separate neuronal and non-neuronal sample betas prior to normalisation

```{r, split cell types, echo = FALSE} 
table(sample_sheet$Nuclei_Fraction)
betas.NeuN <- betas[, colnames(betas) %in% sample_sheet$Basename[sample_sheet$Nuclei_Fraction == "NeuN"]]
betas.DN_PU1 <- betas[, colnames(betas) %in% sample_sheet$Basename[sample_sheet$Nuclei_Fraction != "NeuN"]]


```

# Normalisation with wateRmelon

```{r, normalise, echo=FALSE}

normbetas.NeuN <- betaqn(betas.NeuN)
normbetas.DN_PU1 <- betaqn(betas.DN_PU1)

identical(rownames(normbetas.DN_PU1), rownames(normbetas.NeuN))
normbetas <- cbind(normbetas.DN_PU1, normbetas.NeuN)

save(normbetas, sample_sheet, file = "MouseArray_CellDeconv_FilteredNormalised_Betas_WT.rdat")
# save(normbetas.NeuN, file = "MouseArray_CellDeconv_FilteredNormalised_NeuN-Betas.rdat")
# save(normbetas.DN_PU1, file = "MouseArray_CellDeconv_FilteredNormalised_DN-PU1-Betas.rdat")

# write QC metrics
write.csv(QCmetrics, "MouseArray_CellDeconv_QCmetrics_WT.csv")

# WRITE M-VALUES
save(mVals, file = "MouseArray_CellDeconv_Filtered_mVals.rdat")

## Heatmap using normalised betas
normbetas <- normbetas[, sample_sheet$Basename]
identical(colnames(normbetas), sample_sheet$Basename)
png(file="heatmap_top5000_NORM.png",width=800, height=800)
sigma <- apply(normbetas, 1, sd)
heatmap.2(betas[order(sigma, decreasing = TRUE)[1:5000],], main = "Passed Samples - top 5000", trace = "none", labCol = sample_sheet$Sample_ID, dendrogram = "column", labRow = "", density.info = "none", scale = "none", cexCol = 0.6)
dev.off()

```


```{r density plots normalized, echo = FALSE, message = F}
densityPlot(normbetas, main = "Normalised betas", sampGroups = sample_sheet$Nuclei_Fraction)

# reorder columns
mVals <- mVals[, sample_sheet$Basename]
# density plot of M-values
densityPlot(mVals, main = "M-values", sampGroups = sample_sheet$Nuclei_Fraction)

#multidensity(betas,main="Multidensity") 

```


