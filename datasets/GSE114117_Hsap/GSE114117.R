library(GEOquery)
library(edgeR)


gset <- getGEO("GSE114117", GSEMatrix =T, AnnotGPL = T)
#extracting metadata

#organism 
org <- gset$GSE114117_series_matrix.txt.gz$organism_ch1[[1]]

#platform
plat <- gset$GSE114117_series_matrix.txt.gz$platform_id[[1]]

#tissue type
tissue <- gset$GSE114117_series_matrix.txt.gz$characteristics_ch1[[1]]

data.proc <- gset$GSE114117_series_matrix.txt.gz$data_processing[[1]]

gset <- gset[[1]]


#Establishing groups


#2 day osteogenesis 
ground <- c(3)
space <- c(2)

two_day <- append(ground, space)
groups_repr <- c()

groups_repr[which(two_day == ground)] <- "1"
groups_repr[which(two_day== space)] <- "2"
fl <- as.factor(groups_repr)

ex <- exprs(gset)


