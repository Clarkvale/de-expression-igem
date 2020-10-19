library(GEOquery)
library(limma)
source("microarray_functions.R")
gse <- getGEO("GSE105058", AnnotGPL = T, GSEMatrix = T)[[1]]

limma::plotDensities(gse)
limma::plotMA(gse)

s <- c(1,2,3)
g <- c(4,5,6)

de <- de.analysis(microgravity_group = s, ground_group = g, gse)
de$TopTable <- remove.controls(de$TopTable)
View(de$TopTable)

write.csv(de$TopTable$TopTable, file = "datasets/GSE105058_Atha/GSE105058.csv")

extractMetaData(gse_groups = list(de), metaLabels = "Atha", microgravity_type = M.TYPE$SPACEFLOWN, filename = "datasets/GSE105058_Atha/GSE105058_meta")
