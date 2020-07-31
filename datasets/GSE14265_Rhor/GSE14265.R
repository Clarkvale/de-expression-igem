#DE-ANALYSIS BENJAMIN CLARK RHODOSPIRILLUM RUBRUM S1H

library(GEOquery)
library(Biobase)
library(limma)
library(statmod)
source("microarray_functions.R")

gset <- getGEO("GSE14265", GSEMatrix =TRUE, AnnotGPL=TRUE)

supfiles <- getGEOSuppFiles("GSE14265")


gse <- gset[[1]]
fvarLabels(gse) <- make.names(fvarLabels(gse))

#SpaceFlown
biol.rep <- c(1,3)
reps <- gse[,biol.rep]
space <- simple_2ch(exprs(reps))

corfit <- duplicateCorrelation(reps, ndups = 1, block = biol.rep) 

simple_2ch <- function(ex){
  fit <-lmFit(ex)
  fit <- eBayes(fit, 0.01)
  return(topTable(fit, number = length(fit[[1]]), sort.by = "B", adjust.method = "fdr", confint = FALSE))
}
