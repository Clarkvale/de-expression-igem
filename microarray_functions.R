# MICROARRAY FUNCTIONS ADAPTED FROM GEO2R SCRIPTS FOR GENERALIZED DE ANALYSIS
#Benjamin Clark

library(enumerations)
library(arrayQualityMetrics)
M.TYPE <- create.enum(c("HARV","RPM","HYPERBOLIC","SPACEFLOWN", "RCCS"))

logcheck <- function(expression_matrix){
  
  qx <- as.numeric(quantile(expression_matrix, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- !((qx[5] > 100) ||
              (qx[6]-qx[1] > 50 && qx[2] > 0) ||
              (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2))
  
  return(LogC)
  
}

extractMetaData <- function(gse_groups, filename, microgravity_type, metaLabels, strain = ""){
  if(!is(microgravity_type, "character")){
    e <- simpleError("Not a valid microgravity type")
    stop(e)
  }
  
  
  org <- gse_groups[[1]]$GSE$organism_ch1[[1]]
  if(strain != ""){
    org <- paste(org, strain, sep = " ")
  }
  
  
  
  
  #wipe datatable if exists
  #write.table(data.frame(), paste(filename, ".csv", sep = ""), append = FALSE)
  for(i in 1:length(gse_groups)){
   gsms <- gse_groups[[i]]$GSE$geo_accession
   titles <- gse_groups[[i]]$GSE$title
   descriptions <- gse_groups[[i]]$GSE$description
   write.csv(data.frame(accesssions = gsms, treatment = titles, description = descriptions), paste(filename, "_" ,metaLabels[[i]], ".csv", sep = ""), append = FALSE )
  }

  #clean the directory
  sink()
  unlink(paste(filename, ".txt", sep = ""))
  
  #append and create output
  sink(paste(filename, ".txt", sep = ""), append = TRUE)
  print(paste("ORGANISM:", org))
  print(paste("MICROGRAVITY TYPE:",microgravity_type))
  
  
  
  print(experimentData(gse_groups[[1]]$GSE))
  
  sink()
  
}

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
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number = length(fit2[[1]]))
  
  #We can get rid of the "number" argument and replace it with lfc = 2 for a log2 fold change of 2
  #tT <- topTable(fit2, adjust="fdr", sort.by="B", lfc = 2)
  
  #select parameters we want for the output
  tT <- subset(tT, select=c("adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title",
                            "Platform_ORF", "GO.Function", "GO.Process", "GO.Component", "Chromosome.annotation","ID"))
  out.list <- list("TopTable" = tT, "GSE" = filtered.gse)
  return(out.list)
  
} 

remove.controls <- function(topTable){
  failed.probes <- which(topTable$Platform_ORF == "")
  passed.probes <- which(topTable$Platform_ORF != "")
  f.topTable <- topTable[passed.probes,]
  return(list(TopTable = f.topTable, failed.probes = failed.probes))
}


  

