###GSE64468 YEAST MICROARRAY TWO CHANNEL
#BENJAMIN CLARK 

library(GEOquery)
library(Biobase)
library(limma)
library(arrayQualityMetrics)
source("microarray_functions.R")


list.gse <- getGEO("GSE50881", GSEMatrix =TRUE, AnnotGPL=TRUE)

gse <- list.gse[[1]]

#take a look into the dataset
#head(gse)
View(gse$title)
#dim(gse)


fvarLabels(gse) <- make.names(fvarLabels(gse))

design <- cbind(DyeEffect = 1,SpaceVsGround = c(1,-1,1,-1,1,-1,1,-1))

fit <- lmFit(exprs(gse),design)
fit <- eBayes(fit)

topTable.dyes <- topTable(fit, coef = "DyeEffect")
topTable.SvsG <- topTable(fit, adjust.method = "fdr", sort.by = "B", number = length(fit[[1]]))
anno.data <- fData(gse)
