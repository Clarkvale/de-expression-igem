###GSE50881 YEAST MICROARRAY TWO CHANNEL
#BENJAMIN CLARK 

library(GEOquery)
library(Biobase)
library(limma)
library(arrayQualityMetrics)
library(FoldGO)
library(GO.db)
library(dplyr)
source("microarray_functions.R")


list.gse <- getGEO("GSE50881", GSEMatrix =TRUE, AnnotGPL=TRUE)

gse <- list.gse[[1]]

#take a look into the dataset
#head(gse)
#View(gse$title)
#dim(gse)


fvarLabels(gse) <- make.names(fvarLabels(gse))

design <- cbind(DyeEffect = 1,SpaceVsGround = c(1,-1,1,-1,1,-1,1,-1))
#block <- fData(gse[,1])$Block


fit <- lmFit(gse,design)
fit <- eBayes(fit, 0.01)

#topTable.dyes <- topTable(fit, coef = "DyeEffect")

anno.data <- fData(gse)

topTable.SvsG <- topTable(fit, adjust.method = "fdr", sort.by = "B", number = length(fit[[1]]))
f.tT <- remove.candida.controls(topTable.SvsG)

f.tT <- subset(f.tT, select = c("CGD_Systematic_Name", "Description", "SpaceVsGround", "AveExpr", "F", "P.Value", "adj.P.Val", "SPOT_ID", "ID"))



rowid <- c() 
for(i in 1:length(f.tT$ID)){
  rowid <- append(rowid,((f.tT %>% filter(SPOT_ID == f.tT$SPOT_ID[i]) %>% arrange(P.Value))[1,])$ID)
  
}

best.rows.ids <- rowid[-which(duplicated(rowid))]
rm.ftT <- na.omit(f.tT[best.rows.ids,]) %>% arrange(adj.P.Val)

match.id.lowest.pval <- function(id){
  best.p <- as.data.frame(f.tT) %>% dplyr::filter(SPOT_ID == id) %>% dplyr::arrange(P.Value)[1,]
  return(best.p)
}



remove.candida.controls <- function(toptable){
  controls <- which(toptable$CGD_Systematic_Name == "")
  passed.probes <- which(toptable$CGD_Systematic_Name != "")
  return(topTable.SvsG[passed.probes,])
}





##getting GO annotations locally
setwd("C:/Users/Benja/repos/de-expression-git/de-expression-igem/datasets/GSE50881_Calb")
go_ids <- read.csv("Candida_albicans_GO.csv")
ordered.gos <- list()
for(annos in 1:length(rm.ftT$CGD_Systematic_Name)){
  indexes <- which(rm.ftT$CGD_Systematic_Name[annos] == go_ids$systemic_id)
  ids <- go_ids$go_id[indexes]  
  ordered.gos[annos] <- list(ids)
  
}

anno.terms <- list()
full.MF <- c()
full.CC <- c()
full.BP <- c()
for(i in 1:length(ordered.gos)){
  GO.PROCESS <- c()
  GO.FUNCTION <- c()
  GO.COMPONENT <- c()
  if(length(ordered.gos[[i]]) == 0){
   full.MF <- append(full.MF, NA)
   full.CC <-  append(full.CC, NA) 
   full.BP <- append(full.BP, NA)
   next
  }
  for(j in 1:length(ordered.gos[[i]])){
    
    if(Ontology(ordered.gos[[i]][[j]]) == "MF"){
      GO.FUNCTION <- append(GO.FUNCTION, Term(ordered.gos[[i]][[j]]))
      
    }else if(Ontology(ordered.gos[[i]][[j]]) == "CC"){
      GO.COMPONENT <- append(GO.COMPONENT, Term(ordered.gos[[i]][[j]]))
      
    }
    else if(Ontology(ordered.gos[[i]][[j]]) == "BP"){
      GO.PROCESS <- append(GO.PROCESS, Term(ordered.gos[[i]][[j]]))
      
    }
    
  }
  if(length(GO.FUNCTION) != 0){
    full.MF <- append(full.MF, paste(GO.FUNCTION, collapse = "///"))}
  else{
    full.MF <- append(full.MF, NA)
  }
  if(length(GO.PROCESS) != 0){
    full.BP <- append(full.BP, paste(GO.PROCESS, collapse = "///"))}
  else{
    full.BP <- append(full.BP, NA)
  }
  if(length(GO.COMPONENT) != 0){
    full.CC <- append(full.CC, paste(GO.COMPONENT, collapse = "///"))}
  else{
    full.cc <- append(full.cc, NA)
  }
  #anno.terms[[anno.data$CGD_Systematic_Name[i]]] <- list(GO.FUNCTION = GO.FUNCTION, GO.COMPONENT = GO.COMPONENT, GO.PROCESS = GO.PROCESS)
}
names <- sapply(rm.ftT$Description, USE.NAMES = FALSE, FUN = function(x){
  if(startsWith(x , prefix = "|")){
    return(NA)
  }
  else{
    return(stringr::str_extract(strsplit(x, split = "|", fixed = TRUE)[[1]][1], pattern = "\\w+"))
  }
})
names <- unlist(names)

Chromosome_Location <- stringr::str_extract_all(rm.ftT$Description, pattern = "(Contig\\d+:\\D+\\d+..\\d+\\)|Contig\\d+:\\d+..\\d+)")
Chromosome_Location <- sapply(Chromosome_Location, FUN = paste, collapse = "///", USE.NAMES = FALSE)

final.tT <- rm.ftT %>% mutate( GO.Function = full.MF, GO.Component = full.CC, 
                               GO.Process = full.BP, Gene.Name = names, 
                               Chromosome.Location = Chromosome_Location, 
                               ID = NULL, SPOT_ID = NULL, Description = NULL) %>% 
          rename(LogFC = SpaceVsGround, Platform_ORF = CGD_Systematic_Name)


write.csv(final.tT, file = "datasets/GSE50881_Calb/GSE50881.csv")


extractMetaData(gse_groups = list(list(GSE = gse[,c(1,3,5,7)]), list(GSE = gse[,c(2,4,6,8)])), microgravity_type = M.TYPE$SPACEFLOWN, 
                filename = "datasets/GSE50881_Calb/GSE50881_meta", metaLabels = c("dye_swap1", "dye_swap2"), strain = "SC5413")
