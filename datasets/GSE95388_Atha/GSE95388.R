#BENJAMIN CLARK PLANT STUDY SPACE FLOWN

library(GEOquery)
library(Biobase)
library(limma)
source("microarray_functions.R")

#This pulls the all the samples from the microarry dataset. This returns a list containing a single expressionSet object. 
gset <- getGEO("GSE95388", GSEMatrix =TRUE, AnnotGPL=TRUE)

gse <- gset[[1]]
fvarLabels(gse) <- make.names(fvarLabels(gse))



gset <- gset[[1]]



control <- c(6,7,8,9,10)
treatment <- c(1,2,3,4,5)


de <- de.analysis(microgravity_group = treatment,ground_group = control, gse  = gse)
filtered.probes <- remove.controls(de$TopTable)

tT.name <- "datasets/GSE95388_Atha/GSE95388.csv"
write.table(filtered.probes$TopTable, tT.name, row.names = FALSE, sep = ",")

meta.name <- "datasets/GSE95388_Atha/GSE95388_meta"

gse_list <- list(de)
labels <- c("")
extractMetaData(filename = meta.name, gse_groups = gse_list, microgravity_type = M.TYPE$SPACEFLOWN, metaLabels = labels, strain = "COL-0")


