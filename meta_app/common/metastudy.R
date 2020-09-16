#pairwise comparisons

source("metastudy_functions.R")
library(limma)
library(GEOquery)
library(affy)
library(arrayQualityMetrics)
library(splines)



heat <- getGEO("GSE132186", GSEMatrix = T, AnnotGPL= T)
h37 <- heat[[1]][,append(1:3,grep(pattern = "WT_37C", heat[[1]]$title))]
targets <- data.frame(Name = h37$title, Target = )



gravity <- getGEO("GSE4136", GSEMatrix =TRUE, AnnotGPL=TRUE)[[1]]
exprs(gravity) <- normalizeCyclicLoess(log2(exprs(gravity)))


#gravity dupcor
design <- cbind(grp1=1, grp2vs1 = c(0,0,0,0,0,0,1,1,1,1,1,1))
block <- c(1,1,1,2,2,2,1,1,1,2,2,2)
dupcor.gravity <- duplicateCorrelation(gravity, design, block = block)

fit <- lmFit(gravity, design, block = block, correlation = dupcor.gravity$consensus.correlation, )
fit <- eBayes(fit)









