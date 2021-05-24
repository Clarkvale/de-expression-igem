# MICROARRAY FUNCTIONS ADAPTED FROM GEO2R SCRIPTS FOR GENERALIZED DE ANALYSIS
#Benjamin Clark

#License

# © Copyright 2020 iGEM Concordia, Maher Hassanain, Benjamin Clark, Hajar Mouddene, Grecia Orlano
# This file is part of AstroBio.
# AstroBio is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or any later version.
# AstroBio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with AstroBio.  If not, see <https://www.gnu.org/licenses/>.

library(enumerations)
library(dplyr)

#An enumeration for Microgravity types
M.TYPE <- create.enum(c("HARV","RPM","HYPERBOLIC","SPACEFLOWN", "RCCS"))

#Confidence intervals to standard error conversion
ci2se <- function(CI.R, CI.L){
  return((CI.R - CI.L)/3.92)
}

#The above function but applied to a topTable 
topTable.SE <- function(topTable){
  return(mapply(CI.R = topTable$CI.R, CI.L = topTable$CI.L, FUN = ci2se))
}

#Row wise OR evaluation for boolean matrices. Returns a vector of True row indices. 
rowOR <- function(matrix){
  return(as.vector(which(apply(matrix, MARGIN = 1, FUN = any))))
}

#Checks if an expression matrix has a logarithmic distribution
logcheck <- function(expression_matrix){
  
  qx <- as.numeric(quantile(expression_matrix, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- !((qx[5] > 100) ||
              (qx[6]-qx[1] > 50 && qx[2] > 0) ||
              (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2))
  
  return(LogC)
  
}

#Extracts metadata from your analysis. Will print a txt file for the whole study
#as well as annotated design matrices in tabulated format. 'contrasts' MUST be in list for this to work 
extractMetaData <- function(gse, design, contrasts,  filename, microgravity_type, 
                            metaLabels, strain = "", cellType = "", description_label = NA){
  if(!is(microgravity_type, "character")){
    e <- simpleError("Not a valid microgravity type")
    stop(e)
  }
  
  
  org <- gse$organism_ch1[[1]]
  if(strain != ""){
    org <- paste(org, strain, sep = " ")
  }
  
  #iterate through the contrasts matrices
  for(c in 1:length(contrasts)){
    #for each contrast matrix pull the columns which we are comparing
    cols <- which(contrasts[[c]] != 0)
    #Get the gsms (rows) pertaining to the columns we selected above, 
    #a gsm is selected on the basis it has at least one non-zero number in the 
    #row within the design matrix 
    gsms <- gse[,rowOR(design[,cols] != 0)]
    
    #pulling general info
    titles <- gsms$title
    if(any(grepl("description", varLabels(gsms), ignore.case = T)) && is.na(description_label)){
      descriptions <- gsms$description
    }
    else{
      descriptions <- description_label
    }
    accessions <-  gsms$geo_accession
    #writing design dataframes
    df <- data.frame(accessions = accessions, treatment = 
                       titles, description = descriptions )
    suppressWarnings(write.csv(df, paste(filename, "_", metaLabels[[c]], ".csv", sep = ""), append = FALSE))
  }
  
  


  #clean the directory
  suppressWarnings(sink())
  unlink(paste(filename, ".txt", sep = ""))
  
  #append and create output
  sink(paste(filename, ".txt", sep = ""), append = TRUE)
  print(paste("ORGANISM:", org))
  print(paste("MICROGRAVITY TYPE:",microgravity_type))
  print(paste("CELL TYPE:", cellType))
  
  
  
  print(experimentData(gse))
  
  sink()
  
}
#pulls columns from affymetrix-styled arrays
pull.output.tT <- function(tT){
  ses <- ci2se(tT$CI.R, tT$CI.L)
  prob <- sapply(tT$B, plogis)
  out_tT <- tT %>% dplyr::select(adj.P.Val, P.Value, t, B, logFC, Gene.symbol, 
                                 Gene.title, Platform_ORF, GO.Function, GO.Function.ID,
                                 GO.Process, GO.Process.ID, GO.Component, GO.Component.ID,
                                 Chromosome.annotation, Gene.ID) %>% 
    
    dplyr::mutate(Standard.Error = ses, Probablity = prob) %>%
    dplyr::rename(Entrez.ID = Gene.ID)
  
  #out_tT <- remove.controls(out_tT)$TopTable
  return(out_tT)
  
}
  
  
#depreciated, only to be used if the experiment contains one contrast (0g vs g) and is a single channel chip
de.analysis <- function(microgravity_group, ground_group, gse){
  
  if (!(is(gse, "ExpressionSet") && is.vector(microgravity_group) && is.vector(ground_group))){
    e <- simpleError("Improper data type(s) in signature")
    stop(e)
  }
  group.set <- append(ground_group, microgravity_group)
  groups_repr <- c()
  
  groups_repr[which(group.set == ground_group)] <- "normal.gravity"
  groups_repr[which(group.set == microgravity_group)] <- "micro.gravity"
  fl <- as.factor(groups_repr)
  
  #filtering out according to groups
  filtered.gse <- gse[, group.set]
  
  #pulling expression data from it 
  ex <- exprs(gse)[,group.set]

  
  #check if the data is log-transformed
  LogC <- logcheck(ex)
  
  
  #if the sample isnt log transformed then we do it ourselves 
  if (!LogC) { 
    ex[which(ex <= 0)] <- NaN
    exprs(filtered.gse) <- log2(ex) }
  
  # set up the data and proceed with analysis
  #sml <- paste("G_", groups_repr, sep="")    # set group names
  filtered.gse$description <- fl
  design <- model.matrix(~ description + 0, filtered.gse)
  colnames(design) <- levels(fl)
  fit <- lmFit(filtered.gse, design)
  cont.matrix <- makeContrasts(micro.gravity-normal.gravity, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2, 0.01)
  
  #pulling the whole fitted microarray dataset and ordering by B value 
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number = Inf, confint = TRUE)
  
  #We can get rid of the "number" argument and replace it with lfc = 2 for a log2 fold change of 2
  #tT <- topTable(fit2, adjust="fdr", sort.by="B", lfc = 2)
  
  #select parameters we want for the output
  #tT <- subset(tT, select=c("adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title",
   #                         "Platform_ORF", "GO.Function", "GO.Process", "GO.Component", "Chromosome.annotation",
    #                        "ID", "Gene.ID", "GO.Function.ID", "GO.Process.ID", "GO.Component.ID"))
  out.list <- list("TopTable" = tT, "GSE" = filtered.gse)
  return(out.list)
  
} 

#Looks up probes in a topTable that do not have an ORF identifier and removes them.
remove.controls <- function(topTable){
  failed.probes <- which(topTable$Platform_ORF == "")
  passed.probes <- which(topTable$Platform_ORF != "")
  f.topTable <- topTable[passed.probes,]
  return(list(TopTable = f.topTable, failed.probes = failed.probes))
}

.compress.ids <- function(entrezid, go.matrix){
  lines <- go.matrix %>% filter(ENTREZID == entrezid)
  cc <- lines %>% filter(ONTOLOGY == "CC")
  bp <- lines %>% filter(ONTOLOGY == "BP")
  mf <- lines %>% filter(ONTOLOGY == "MF")
  
  return(list(cc.go = paste(unique(cc$GO), collapse = "///"),
              cc.term = paste(unique(cc$TERM), collapse = "///"),
              bp.go = paste(unique(bp$GO), collapse = "///"),
              bp.term = paste(unique(bp$TERM), collapse = "///"),
              mf.go = paste(unique(mf$GO), collapse = "///"),
              mf.term = paste(unique(mf$TERM), collapse = "///")))
}

.make.vector <- function(index, list, name){
  return(list[[index]][[name]])
}

get.GOs <- function(TopTable, org.database, Entrez.name){
  require(GO.db)
  
  if(is.null(TopTable[[Entrez.name]])){
    stop("Invalid ID vector. Use a column with Entrez IDs.")
  }
  gos <- AnnotationDbi::select(org.database, keys = as.character(TopTable[[Entrez.name]]), keytype= "ENTREZID", columns= c("GO", "ONTOLOGY"))
  u.ids <- as.numeric(unique(gos$ENTREZID))
  TERM <- apply(FUN = Term, as.array(gos$GO), MARGIN = 1)
  gos <- cbind(gos, TERM)
  
  
  listed.gos <- lapply(u.ids, FUN = .compress.ids, go.matrix = gos)
  vnames <- c("cc.go", "cc.term", "bp.go", "bp.term", "mf.go", "mf.term")
  
  df_out <- data.frame(u.ids)
  
  for(i in 1:length(vnames)){
    df_out <- cbind(df_out,  sapply(1:length(u.ids), FUN = .make.vector, list = listed.gos, name = vnames[i]))
    
  }
  
  colnames(df_out) <- append("Entrez.id", vnames)
  return(df_out)
  
  
}


  

