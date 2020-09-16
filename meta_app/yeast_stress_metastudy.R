###Generalized microgravity-stress response Co-Expression Analysis 

library(GEOquery)
library(Biobase)
library(limma)
library(arrayQualityMetrics)
library(FoldGO)
library(MetaVolcanoR)
source("microarray_functions.R")


stress.gse <- getGEO("GSE18", GSEMatrix = TRUE, AnnotGPL = TRUE)
heatshock <- stress.gse$`GSE18-GPL51_series_matrix.txt.gz`
nitro.sorb <- stress.gse$`GSE18-GPL52_series_matrix.txt.gz`
dtt.diamide <- stress.gse$`GSE18-GPL53_series_matrix.txt.gz`
dtt.ade_starvation <- stress.gse$`GSE18-GPL54_series_matrix.txt.gz`
hypoosmotic <- stress.gse$`GSE18-GPL56_series_matrix.txt.gz`
mena.h202 <- stress.gse$`GSE18-GPL64_series_matrix.txt.gz`

#_____________________________________________________________
heat <- heatshock[,c(1,2,3,5,6,7,8)]

hlist <- list()
for(i in 1:length(heat)){
  
}

fvarLabels(heat) <- make.names(fvarLabels(heat))
plot(boxplot((exprs(heat))))



#exprs(early.heat) <- normed.heat
#plot(boxplot((normed.heat)))
#plotMA(((normed.heat)))


heat.tT <- simple_2ch(heat)
early.heat.tT <- remove.controls(early.heat.tT)
write.csv(early.heat.tT$TopTable, file = "data/gse18_heat.csv")
#____________________________________________________________


#___________________________________________________________
early.hypo <- hypoosmotic[,1:3]
fvarLabels(early.hypo) <- make.names(fvarLabels(early.hypo))
#plot(boxplot((exprs(early.hypo))))
#limma::plotMA(early.hypo)
normed.hypo <-  normalizeCyclicLoess(((early.hypo)))

exprs(early.hypo) <- normed.hypo
#plot(boxplot((normed.hypo)))
#plotMA(((normed.hypo)))


early.hypo.tT <- simple_2ch(early.hypo)
early.hypo.tT <- remove.controls(early.hypo.tT)
write.csv(early.hypo.tT$TopTable, file = "data/gse18_hypoosmotic.csv")
#________________________________________________
early.nitro_depletion <- nitro.sorb[,c(20,21,19)]
early.nitro.tT <- do_loess_lmFit(early.nitro_depletion)
write.csv(early.nitro.tT$TopTable, file = "data/gse18_nitrogen_depletion.csv")

#_______________________________
early.dtt <- dtt.diamide[,2:4]
early.dtt.tT <- do_loess_lmFit(early.dtt)
write.csv(early.dtt.tT$TopTable, file = "data/gse18_dtt.csv")

#__________________________________
early.ade <- dtt.ade_starvation[,c(30,31,36)]
early.ade.starve.tT <- do_loess_lmFit(early.ade) 
write.csv(early.ade.starve.tT$TopTable, file = "data/gse18_adenine_starve.csv")
#__________________________________
early.mena <- mena.h202[,c(7,2,3)]
early.mena.tT <- do_loess_lmFit(early.mena)
write.csv(early.mena.tT$TopTable, file = "data/gse18_menadione.csv")

#__________________________________
early.h202 <- mena.h202[,c(14,10,11)]
early.h202.tT <- do_loess_lmFit(early.h202)
#hist(early.h202.tT$TopTable$P.Value)

write.csv(early.h202.tT$TopTable, file = "data/gse18_h202.csv")



#________________________________________
micro.gse <- getGEO("GSE4136", GSEMatrix =TRUE, AnnotGPL=TRUE)[[1]]
fvarLabels(micro.gse) <- make.names(fvarLabels(micro.gse))


#5th gen groups
control <- c(1,2,3)
treatment <- c(7,8,9)
micro.tT <- de.analysis.1ch(gse = micro.gse, microgravity_group = treatment, ground_group = control)

#micro.tT <- de.analysis.1ch(microgravity_group = treatment, ground_group = control, gse = gen5$GSE)
micro.tT <- remove.controls(micro.tT$TopTable)
write.csv(micro.tT$TopTable, file = "common/data/microgravity_5thgen.csv")

#25th Gen groups
control <- c(4,5,6)
treatment <- c(10,11,12)

gen25 <- de.analysis.1ch(gse = micro.gse, microgravity_group = treatment, ground_group = control)
fgen25_rpm <- remove.controls(gen25$TopTable)
write.csv(fgen25_rpm$TopTable, file = "common/data/microgravity_25thgen.csv")




diffList <- list(Microgravity.5thGen = pull.relevant.columns(micro.tT$TopTable),
                 Microgravity.25thGen = pull.relevant.columns(fgen25_rpm$TopTable),
                 Heat.Shock = pull.relevant.columns(early.heat.tT$TopTable),
                 Hypoosmotic.Shock = pull.relevant.columns(early.hypo.tT$TopTable),
                 Nitrogen.Depletion = pull.relevant.columns(early.nitro.tT$TopTable),
                 DTT = pull.relevant.columns(early.dtt.tT$TopTable),
                 Adenine.Starvation = pull.relevant.columns(early.ade.starve.tT$TopTable),
                 Menadione = pull.relevant.columns(early.mena.tT$TopTable),
                 Peroxide = pull.relevant.columns(early.h202.tT$TopTable))



meta_degs_rem <- rem_mv(diffexp=diffList,
                        pcriteria="pvalue",
                        foldchangecol='Log2FC', 
                        genenamecol='Symbol',
                        geneidcol=NULL,
                        collaps=TRUE,
                        llcol='CI.L',
                        rlcol='CI.R',
                        vcol=NULL, 
                        cvar=TRUE,
                        metathr=0.01,
                        jobname="multi-stress-metastudy-volcano-plot",
                        ncores=1,
                        draw = 'HTML')

save(meta_degs_rem, meta_comb, diffList, file = "common/data/appdata.RData")


meta_comb <- MetaVolcanoR::combining_mv(diffexp = diffList, collaps = TRUE)



do_loess_lmFit <- function(gse){
  #plot(hist(exprs(gse)), main = "Sample Distribution Before Normalization")
  fvarLabels(gse) <- make.names(fvarLabels(gse))
  #plot(boxplot((exprs(gse)), main = "Sample Distributions Before LOESS Normalization"))
  #plotMA(gse, main = "Sample MA Plot Before LOESS Normalization")
  normed.gse <-  normalizeCyclicLoess(((gse)))
  
  exprs(gse) <- normed.gse
  #boxplot((normed.gse))
  #plotMA(((normed.gse)))
  
  
  gse.tT <- simple_2ch(gse)
  gse.tT <- remove.controls(gse.tT)
  
}

simple_2ch <- function(ex){
  fit <-lmFit(ex)
  fit <- eBayes(fit, 0.01)
  return(topTable(fit, number = length(fit[[1]]), sort.by = "B", adjust.method = "fdr", confint = TRUE))
}
de.analysis.1ch <- function(microgravity_group, ground_group, gse){
  
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
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number = length(fit2[[1]]), confint = TRUE)
  
  #We can get rid of the "number" argument and replace it with lfc = 2 for a log2 fold change of 2
  #tT <- topTable(fit2, adjust="fdr", sort.by="B", lfc = 2)
  
  #select parameters we want for the output
  #tT <- subset(tT, select=c("adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title",
  #                          "Platform_ORF", "GO.Function", "GO.Process", "GO.Component", "Chromosome.annotation","ID"))
  out.list <- list("TopTable" = tT, "GSE" = filtered.gse)
  return(out.list)
  
}

pull.relevant.columns <- function(tT){
  gene.symbols <- tT$Gene.symbol
  orfs <- tT$Platform_ORF
  nas <- which(gene.symbols == "")
  gene.symbols[nas] <- orfs[nas]
  
  logfc <- as.numeric(tT$logFC)
  CI.L <- as.numeric(tT$CI.L)
  CI.R <- as.numeric(tT$CI.R)
  p.val <- as.numeric(tT$adj.P.Val)
  
  return(data.frame("Symbol" = gene.symbols,"Log2FC" = logfc, "pvalue" = p.val, "CI.L" = CI.L, "CI.R" = CI.R))
  
}

gene.vector <- c(
  "RGI1",
  "TDA10",
  "MPC3",
  "BTN2",
  "OPI10",
  "YMR084W",
  "AFR1",
  "ARG82",
  "RNY1",
  "GYP7",
  "RAV2",
  "TDA1",
  "SOL1",
  "GGA1",
  "FUN19",
  "ATG3",
  "ATG1",
  "AIM17",
  "STP4",
  "IMA3",
  "OSW2",
  "CMK2",
  "VID30",
  "THI4",
  "SDS22",
  "YPI1",
  "RRT8",
  "MDH1",
  "YHR016C",
  "PRM8",
  "IGO2",
  "RPN4",
  "PCL8",
  "TRX3"
)

for( i in 1:length(gene.vector)){
  gene <- gene.vector[i]
  draw_forest(remres = meta_degs_rem, gene = gene, jobname = "", outputfolder = "stress_metaanalysis")
}

meta_degs_comb <- combining_mv(diffexp= diffList,
                               pcriteria='pvalue', 
                               foldchangecol='Log2FC',
                               genenamecol='Symbol',
                               geneidcol=NULL,
                               metafc='Mean',
                               metathr=0.01, 
                               collaps=TRUE,
                               jobname="MetaVolcano",
                               outputfolder="stress_metaanalysis",
                               draw='HTML')



