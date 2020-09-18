
#GSE40648 Ecoli MICROARRAY
#R script by BENJAMIN CLARK
#Analysis by Hajar El Mouddene

library(GEOquery)
library(Biobase)
library(limma)
library(org.EcK12.eg.db)
source("microarray_functions.R")
#This pulls the all the samples from the microarry dataset. This returns a list containing a single expressionSet object. 
list.gse <- getGEO("GSE40648", GSEMatrix =TRUE, AnnotGPL=TRUE)

gse <- list.gse[[1]]

#take a look into the dataset
#head(gse)
#View(gse$title)
#dim(gse)

#here we format the feature names. fvarLabels belongs the biobase package and is used to extract features from ExpressionSet Objects
fvarLabels(gse) <- make.names(fvarLabels(gse))


#DE analysis


gse.rm <- gse[,-c(4,8)]
control <- c(1,2,3)
treatment <- c(4,5,6)

ecoli <- de.analysis(gse = gse.rm, microgravity_group = treatment, ground_group = control)

#print out the toptable 
ecoliname <- "datasets/GSE40648_Ecoli/GSE40648_Ecoli.csv"
write.table(ecoli$TopTable, ecoliname, row.names = FALSE, sep = ",")
remove.controls(ecoli$TopTable)
filtered.tT <- remove.controls(ecoli$TopTable)



#adding new GO ids
symbols <- filtered.tT$TopTable$Gene.symbol
symbols.1 <- sapply(symbols, FUN = function(x){return(strsplit(x, "///")[[1]][1])})
names(symbols.1) <- NULL
symbols.1 <- as.vector(symbols.1)

esym <- org.EcK12.egSYMBOL2EG
esym <- as.list(esym)
entrez_ids <- unlist(esym[symbols.1])


gos <- org.EcK12.egGO
gos <- as.list(gos)

mapped_gos <- gos[entrez_ids]

go.db <- as.list(GO.db::GOTERM)
anno.go.genome <- function(go.ids){
  out <- list()
  for(i in 1:length(go.ids)){
    line_ids <- names(go.ids[[i]])
    line.go.items <- list(go.db[line_ids])
    terms <- sapply(names(line.go.items[[1]]),FUN = Term)
    names(terms) <- NULL
    terms <- unlist(terms)
    
    ontologies <- sapply(names(line.go.items[[1]]), FUN = Ontology)
    names(ontologies) <- NULL
    ontologies <- unlist(ontologies)
    
    inner_out <- data.frame(id = names(line.go.items[[1]]), term = terms, ontology = ontologies)
    out[names(go.ids)[i]]<- list(inner_out)
    
    
    
  }
  return(out)
}
go.genome <- anno.go.genome(mapped_gos)
valid_ids <- c()
for(i in 1:length(entrez_ids)){
  if(length(go.genome[[entrez_ids[i]]]$ontology) > 0){  
    GO.Function <- paste(dplyr::filter(go.genome[[entrez_ids[i]]], ontology == "MF")$term , collapse = "///")
    GO.Process <- paste(dplyr::filter(go.genome[[entrez_ids[i]]], ontology == "BP")$term , collapse = "///")
    GO.Component <- paste(dplyr::filter(go.genome[[entrez_ids[i]]], ontology == "CC")$term , collapse = "///")
  }
  else{
    next
  }
  
  ided.rows <- grep(filtered.tT$TopTable$Gene.symbol, pattern = names(entrez_ids[i]))
  correct.probe <-  which.min(filtered.tT$TopTable[ided.rows,]$P.Value)
  
  filtered.tT$TopTable$GO.Function[ided.rows[correct.probe]] <- GO.Function
  
  filtered.tT$TopTable$GO.Component[ided.rows[correct.probe]] <- GO.Component
  
  filtered.tT$TopTable$GO.Process[ided.rows[correct.probe]] <- GO.Process
  
  
  valid_ids <- append(valid_ids, ided.rows[correct.probe])
  
  
}
valid_ids <- valid_ids[-which(duplicated(valid_ids))]
filtered.tT$TopTable <- filtered.tT$TopTable[valid_ids,]



write.table(filtered.tT$TopTable, ecoliname, row.names = FALSE, sep = ",")


#process and extract metadata for all datasets comparisons

metaName <- "datasets/GSE40648_Ecoli/GSE40648_meta"
strain <- "K12 MG1655"
gse_list <- list(ecoli)
labels <- c("ecoli")
extractMetaData(filename = metaName, gse_groups = gse_list, microgravity_type = M.TYPE$RPM, metaLabels = labels, strain = strain)

