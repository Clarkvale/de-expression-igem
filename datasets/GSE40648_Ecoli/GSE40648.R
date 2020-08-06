
#GSE40648 Ecoli MICROARRAY
#R script by BENJAMIN CLARK
#Analysis by Hajar El Mouddene

library(GEOquery)
library(Biobase)
library(limma)
source("microarray_functions.R")
library(arrayQualityMetrics)
library(affy)
library(genefilter)

#This pulls the all the samples from the microarry dataset. This returns a list containing a single expressionSet object. 
list.gse <- getGEO("GSE40648", GSEMatrix =TRUE, AnnotGPL=TRUE)

gse <- list.gse[[1]]

#take a look into the dataset
#head(gse)
View(gse$title)
#dim(gse)

#here we format the feature names. fvarLabels belongs the biobase package and is used to extract features from ExpressionSet Objects
fvarLabels(gse) <- make.names(fvarLabels(gse))



gse <- gse[,-c(4,8)]

exprs(gse) <- normalizeBetweenArrays(exprs(gse))

#DE analysis, removing 4 and 8 
control <- c(1,2,3)
treatment <- c(4,5,6)


ecoli <- de.analysis(gse = gse, microgravity_group = treatment, ground_group = control)

#print out the toptable 
ecoliname <- "datasets/GSE40648_Ecoli/GSE40648_Ecoli.csv"
#write.table(ecoli$TopTable, ecoliname, row.names = FALSE, sep = ",")

filtered.tT <- remove.controls(ecoli$TopTable)
write.table(filtered.tT$TopTable, ecoliname, row.names = FALSE, sep = ",")

#process and extract metadata for all datasets comparisons

metaName <- "datasets/GSE40648_Ecoli/GSE40648_meta"
strain <- "K12 MG1655"
gse_list <- list(ecoli)
labels <- c("ecoli")
extractMetaData(filename = metaName, gse_groups = gse_list, microgravity_type = M.TYPE$RPM, metaLabels = labels, strain = strain)




