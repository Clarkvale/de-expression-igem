###GSE50881 YEAST MICROARRAY TWO CHANNEL
#BENJAMIN CLARK 

library(GEOquery)
library(Biobase)
library(limma)
library(arrayQualityMetrics)
library(FoldGO)
library(GO.db)
source("microarray_functions.R")


list.gse <- getGEO("GSE50881", GSEMatrix =TRUE, AnnotGPL=TRUE)

gse <- list.gse[[1]]

#take a look into the dataset
#head(gse)
#View(gse$title)
#dim(gse)


fvarLabels(gse) <- make.names(fvarLabels(gse))




# design <- cbind(DyeEffect = 1,SpaceVsGround = c(1,-1,1,-1,1,-1,1,-1))
# #block <- fData(gse[,1])$Block
# 
# 
# fit <- lmFit(gse,design)
# fit <- eBayes(fit, 0.01)
# 
# #topTable.dyes <- topTable(fit, coef = "DyeEffect")
# 
# anno.data <- fData(gse)
# 
# f.tT <- remove.candida.controls(topTable.SvsG)
# 
# f.tT <- subset(f.tT, select = c("CGD_Systematic_Name", "Description", "SpaceVsGround", "AveExpr", "F", "P.Value", "adj.P.Val"))


median.blocked.duplicates <- function(GSE){
  blocks <- factor(fData(GSE[,1])$Block)
  spots <- fData(GSE[,1])
  for(i in 1:length(gse)){
    for(block in 1:length(levels(blocks))){
      bi <- which(blocks == as.numeric(levels(blocks)[block]))
      dups <- which(spots$SPOT_ID == spots$SPOT_ID[1])
      num.of.dups <- length(dups)
      dup.length <- dups[2] - dups[1]
      for(j in 1:length(dup.length)){
        med.spot <- which() 
      }
    }
  }
}



get.median.duplicates <- function(gsm, block.num){
  browser()
  block.ids <-  fData(gsm)$SPOT_ID[which(fData(gsm)$Block == block.num)]
  #remove empty spots and autoblanks
  
  block.ids <- block.ids[-which(!grepl("orf19", block.ids))]
  block.index <- which(fData(gsm)$SPOT_ID[which(fData(gsm)$Block == block.num)] %in% block.ids)
  block.ex <- exprs(gsm)[block.index]
  
  dups <- which(block.ids == block.ids[1])
  ndups <- length(dups)
  spacing <- dups[2] - dups[1]
  
  cor <- duplicateCorrelation(block.ex, ndups = ndups, spacing = spacing )
  
  
  unique.spots <- unique(block.ids)
  

  
  counter <- 1
  out <- c()
  added.ar <- c()
  #dupcor <- duplicateCorrelation(gsm[which(block)])
  # while(b.id != unique.spots[length(unique.spots)]){
  #   #find duplicate spots
  #   if(!(b.id %in% added.ar)){
  #     dups <- which(block.ids == b.id)
  #     #pull median values from identical spots
  #     med <- median(na.omit(exprs(gsm[dups])))
  #     #assign to output
  #     out[b.id] <- med
  #     
  #     added.ar <- append(added.ar, b.id)
  #     
  #     
  #   }
  #   #update counter and change id
  #   counter <- counter + 1
  #   b.id <- unique.spots[counter]
  # } 
  return(out)
}

#scrapping the median business and doing duplicatedCorrelation across all slides per block
median.gse <- list()
full.medians <- c()
vblocks <- factor(fData(gse[,1])$Block)
for(gsm in 1:length(colnames(gse))){
  #browser()
  for(b in 1:length(levels(vblocks))){
    
    block.medians <- get.median.duplicates(gsm = gse[,gsm], block.num = as.numeric(b))
    full.medians <- append(full.medians, block.medians) 
  }
  median.gse[[colnames(gse)[gsm]]] <- full.medians 
  full.medians <- c()
}

# median.matrix <- matrix(unlist(median.gse), ncol = 8)
# rownames(median.matrix) <- names(median.gse[[1]])
# colnames(median.matrix) <- names(median.gse)
# median.matrix <- normalizeCyclicLoess(median.matrix)
# 
# fit <- lmFit(median.matrix, design)
# fit <- eBayes(fit, 0.01)
# topTable.SvsG <- topTable(fit, adjust.method = "fdr", sort.by = "B", number = length(fit[[1]]))
# 



##getting GO annotations locally
# setwd("C:/Users/Benja/repos/de-expression-git/de-expression-igem/datasets/GSE50881_Calb")
# go_ids <- read.csv("Candida_albicans_GO.csv")
# ordered.gos <- list()
# for(annos in 1:length(f.tT$CGD_Systematic_Name)){
#   indexes <- which(f.tT$CGD_Systematic_Name[annos] == go_ids$systemic_id)
#   ids <- go_ids$go_id[indexes]  
#   ordered.gos[annos] <- list(ids)
#   
# }
# 
# anno.terms <- list()
# for( i in 1:length(ordered.gos)){
#   GO.PROCESS <- c()
#   GO.FUNCTION <- c()
#   GO.COMPONENT <- c()
#   if( length(ordered.gos[[i]]) == 0){
#   anno.terms[[anno.data$CGD_Systematic_Name[i]]] <- list(GO.FUNCTION = NA, GO.COMPONENT = NA, GO.PROCESS = NA)
#   next
#   }
#   for(j in 1:length(ordered.gos[[i]])){
#     
#     if(Ontology(ordered.gos[[i]][[j]]) == "MF"){
#       GO.FUNCTION <- append(GO.FUNCTION, Term(ordered.gos[[i]][[j]]))
#       
#     }else if(Ontology(ordered.gos[[i]][[j]]) == "CC"){
#       GO.COMPONENT <- append(GO.COMPONENT, Term(ordered.gos[[i]][[j]]))
#     }
#     else if(Ontology(ordered.gos[[i]][[j]]) == "BP"){
#       GO.PROCESS <- append(GO.PROCESS, Term(ordered.gos[[i]][[j]]))
#       
#     }
#   }
#   anno.terms[[anno.data$CGD_Systematic_Name[i]]] <- list(GO.FUNCTION = GO.FUNCTION, GO.COMPONENT = GO.COMPONENT, GO.PROCESS = GO.PROCESS)
# }

