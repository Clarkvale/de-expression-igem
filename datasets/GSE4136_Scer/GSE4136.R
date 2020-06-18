
#GSE4136 YEAST MICROARRAY
#BENJAMIN CLARK

library(GEOquery)
library(Biobase)
library(limma)
source("microarray_functions.R")

#This pulls the all the samples from the microarry dataset. This returns a list containing a single expressionSet object. 
list.gse <- getGEO("GSE4136", GSEMatrix =TRUE, AnnotGPL=TRUE)

#if there are multiple datasets just take the GPL file  
if (length(list.gse) > 1) idx <- grep("GPL2529", attr(list.gse, "names")) else idx <- 1

#if not just take the expressionSet object
gse <- list.gse[[idx]]

#take a look into the datase
#head(gse)
#View(gse)
#dim(gse)

#here we format the feature names. fvarLabels belongs the biobase package and is used to extract features from ExpressionSet Objects
fvarLabels(gse) <- make.names(fvarLabels(gse))


#5th gen groups
control <- c(1,2,3)
treatment <- c(7,8,9)
gen5 <- de.analysis(gse = gse, microgravity_group = treatment, ground_group = control)

#print out the toptable 
gen5name <- "datasets/GSE4136_Scer/GSE4136_5thGen.csv"
write.table(gen5$TopTable, gen5name, row.names = FALSE, sep = ",")


metaName <- "datasets/GSE4136_Scer/GSE4136_meta"

#25th gen groups
control <- c(4,5,6)
treatment <- c(10,11,12)

gen25 <- de.analysis(gse = gse, microgravity_group = treatment, ground_group = control)
gen25name <- "datasets/GSE4136_Scer/GSE4136_25thGen.csv"
write.table(gen25$TopTable, gen25name, row.names = FALSE, sep = ",")


#process and extract metadata for all datasets comparisons
gse_list <- list(gen5,gen25)
extractMetaData(filename = metaName, gse_groups = gse_list, microgravity_type = M.TYPE$HARV)




