
#GSE4136 YEAST MICROARRAY
#BENJAMIN CLARK

library(GEOquery)
library(Biobase)
library(limma)
library(stats)
library(dendextend)
library(gplots)
source("microarray_functions.R")
#This pulls the all the samples from the microarry dataset. This returns a list containing a single expressionSet object. 
list.gse <- getGEO("GSE4136", GSEMatrix =TRUE, AnnotGPL=TRUE)

gse <- list.gse[[1]]

#take a look into the dataset
#head(gse)
View(gse$title)
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



#25th gen groups
control <- c(4,5,6)
treatment <- c(10,11,12)

gen25 <- de.analysis(gse = gse, microgravity_group = treatment, ground_group = control)
gen25name <- "datasets/GSE4136_Scer/GSE4136_25thGen.csv"
write.table(gen25$TopTable, gen25name, row.names = FALSE, sep = ",")


#process and extract metadata for all datasets comparisons

metaName <- "datasets/GSE4136_Scer/GSE4136_meta"
strain <- "BY4743"
gse_list <- list(gen5,gen25)
labels <- c("5thGen", "25thGen")
extractMetaData(filename = metaName, gse_groups = gse_list, microgravity_type = M.TYPE$HARV, metaLabels = labels, strain = strain)





##DE analysis Between generations
micro.5thgen <- c(7,8,9)
micro.25thgen <- c(10,11,12)
gen5.vs.gen25 <- de.analysis(gse = gse, microgravity_group = micro.25thgen, ground_group = micro.5thgen)

failed.probes <- remove.controls(gen5.vs.gen25$TopTable)
##CLUSTERING ANALYSIS


#FILTERING SIGNIFICANT DE GENES
gen5.sig.df.genes <- which(gen5$TopTable$logFC >= 1 & gen5$TopTable$adj.P.Val <= 0.01)
gen5.orfs <- gen5$TopTable$ID[gen5.sig.df.genes]

gen25.sig.df.genes <- which(gen25$TopTable$logFC >= 1 & gen25$TopTable$adj.P.Val <= 0.01)
gen25.orfs <- gen25$TopTable$ID[gen25.sig.df.genes]

grouped.orfs <- union(gen5.orfs,gen25.orfs)

#pulling relevent expression data 
ex <- exprs(gse)
filtered.ex <- ex[rownames(ex) %in% grouped.orfs,]

#take out control columns
filtered.ex <- filtered.ex[,7:12]

#adding our factors
both <- intersect(gen5.orfs,gen25.orfs)
gen5.dif <- setdiff(gen5.orfs, gen25.orfs)
gen25.dif <- setdiff(gen25.orfs,gen5.orfs)


filtered.ex <- log2(filtered.ex)

#Making plots
colnames(filtered.ex) <- c("5gen1", "5gen2", "5gen3", "25gen1", "25gen2", "25gen3")

dev.new(width = 6, height = 256)
plot(heatmap.2(filtered.ex), distfun = function(x) as.dist(1-cor(t(x), method="pearson")))

hc <- as.dendrogram(hclust(as.dist(1-cor(t(filtered.ex), method="pearson")), method="complete"))


#c.index <- which(labels(hc) %in% both)

hc %>% set("labels_col", "white") %>% 
  set("by_labels_branches_col", value = both) %>% 
  plot(main = "Highlighted Union Genes Across Timescales")

#PLOTTING INTERSECTING Genes
filtered.ex.intersect <- ex[rownames(ex) %in% both,][,7:12]
colnames(filtered.ex.intersect) <- c("5gen1", "5gen2", "5gen3", "25gen1", "25gen2", "25gen3")
plot(heatmap(filtered.ex.intersect))


#build annotable
symbols <- gen5$TopTable$Gene.symbol[gen5$TopTable$ID %in% both]
logfc <- gen5$TopTable$logFC[gen5$TopTable$ID %in% both]
p.val <- gen5$TopTable$adj.P.Val[gen5$TopTable$ID %in% both]
