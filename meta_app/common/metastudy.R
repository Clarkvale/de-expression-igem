#pairwise comparisons

source("metastudy_functions.R")
library(limma)
library(GEOquery)
library(affy)
library(arrayQualityMetrics)
library(splines)
library(MetaVolcanoR)
source("custom_draw.R")

gravity <- getGEO("GSE4136", GSEMatrix =TRUE, AnnotGPL=TRUE)[[1]]
#exprs(gravity) <- normalizeCyclicLoess(log2(exprs(gravity)))



levs <- c(rep("NG.5", 3), rep("NG.25",3),rep("MG.5",3), rep("MG.25",3))
targets <- data.frame(cbind(GSE = gravity$geo_accession, Target = levs))

f <- factor(targets$Target, levels = unique(levs))
design <- model.matrix(~0+f)
colnames(design) <- unique(levs)
fit <- lmFit(gravity, design)
cont.dif <- makeContrasts(
  Dif5 = MG.5 - NG.5,
  
  levels = design
)
contrasts.fit <- contrasts.fit(fit, contrasts = cont.dif)

fit2 <- eBayes(contrasts.fit, 0.01)
g.tT <- topTable(fit2, adjust.method = "fdr", confint = TRUE, number = Inf)


#____________oxidative_stress
peroxide_time <- getGEO("GSE26169", GSEMatrix = T, AnnotGPL = T)
gsms_c <- grep(peroxide_time$GSE26169_series_matrix.txt.gz$title, pattern = "wild type, Control")
gsms_t <- grep(peroxide_time$GSE26169_series_matrix.txt.gz$title, pattern = "wild type, CHP treated")  

gsms <- union(gsms_c,gsms_t)
peroxide_time <- peroxide_time$GSE26169_series_matrix.txt.gz[,gsms]

exprs(peroxide_time) <- normalizeCyclicLoess(exprs(peroxide_time))

lev <- c()

times <- c(0,3,6,12,20)
types <- c("Control", "CHP treated")

for(i in 1:length(times)){
  for(j in 1:length(types)){
    label <- make.names(paste(types[j],times[i], sep = "_"))
    items <- grep(peroxide_time$title, pattern = paste(types[j], ", t=", times[i], " min", sep = ""))
    lev[items] <- label
  }
}
targets <- data.frame(cbind(GSE = peroxide_time$geo_accession, Target = lev))

f <- factor(targets$Target, levels = unique(lev))


design <- model.matrix(~0+f)
#colnames(design) <- unique(lev)
fit <- lmFit(peroxide_time, design)
cont.dif <- makeContrasts(
  Dif6min = (fControl_6-fControl_0) - (fCHP.treated_6-fCHP.treated_0),
  
  levels = design
)
#Dif12min = (fControl_12-fControl_0) - (fCHP.treated_12-fCHP.treated_0),
fit2 <- contrasts.fit(fit, cont.dif)
fit2 <- eBayes(fit2, 0.01)
ox.tT <- topTable(fit2, adjust="fdr", sort.by="B", number = length(fit2[[1]]), confint = TRUE)


grav5.vs.ox6 <- list(Microgravity.5thgen = pull.relevant.columns(g.tT),
                     Oxidative.Stress = pull.relevant.columns(ox.tT))

rem_g5.v.ox6 <- rem_mv(diffexp=grav5.vs.ox6,
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
                        draw = '')

#___________________heat shock__________________________________________


heat <- getGEO("GSE132186", GSEMatrix = T, AnnotGPL = T)
gsms <- grep(heat$GSE132186_series_matrix.txt.gz$title, pattern = "WT_37C")
gsms <- union(c(1,2,3), gsms)
heat_37 <- heat$GSE132186_series_matrix.txt.gz[,gsms]


times <- c(1,3,5,7,10,15,40,80)
temp <- c("25C", "37C")

lev <- c()
for(i in 1:length(times)){
  for(j in 1:length(temp)){
    label <- paste(temp[j],times[i], sep = "_")
    items <- grep(heat_37$title, pattern = paste(temp[j], "_",times[i], sep = ""))
    lev[items] <- label
  }
}

lev[1:3] <- "25C"
targets <- data.frame(cbind(GSE = heat_37$geo_accession, Target = lev))

f <- factor(targets$Target, levels = unique(lev))

design <- model.matrix(~0+f)

fit <- lmFit(heat_37, design)
cont.dif <- makeContrasts(
  Dif7min = f37C_7-f25C,
  
  levels = design
)
fit2 <- contrasts.fit(fit, cont.dif)
fit2 <- eBayes(fit2, 0.01)
h.tT <- topTable(fit2, adjust="fdr", sort.by="B", number = length(fit2[[1]]), confint = TRUE)


grav5.vs.he7 <- list(Microgravity.5thgen = pull.relevant.columns(g.tT),
                     Heat.Stress = pull.relevant.columns(h.tT))

rem_g5.v.he7 <-  rem_mv(diffexp=grav5.vs.he7,
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
                        draw = '')

#__________________________osmotic shock__________________________________
osmotic_shock <- getGEO("GSE13097", GSEMatrix = T, AnnotGPL = T)[[1]]
exprs(osmotic_shock) <- normalizeCyclicLoess(log2(exprs(osmotic_shock)))

times <- c(0,2,4,6,8,10,15)
lev <- c()
for( i in 1:length(times)){
  label <- paste("t", times[i], sep = "")
  print(label)
  lev[grep(osmotic_shock$title, pattern = label)] <- label 
}

targets <- data.frame(cbind(GSE = osmotic_shock$geo_accession, Target = lev))
f <- factor(targets$Target)
design <- model.matrix(~0+f)
fit <- lmFit(osmotic_shock, design)
cont.dif <- makeContrasts(
  Dif6min = ft6-ft0,
  
  levels = design
)
fit2 <- contrasts.fit(fit, cont.dif)
fit2 <- eBayes(fit2, 0.01)
o.tT <- topTable(fit2, adjust="fdr", sort.by="B", number = length(fit2[[1]]), confint = TRUE)
o.tT <- o.tT %>% rename(Gene.symbol = Gene.name , Platform_ORF = ORF )

grav5.vs.osm6 <- list(Microgravity.5thgen = pull.relevant.columns(g.tT),
                      Osmotic.Stress = pull.relevant.columns(o.tT))


rem_g5.v.osm6 <-  rem_mv(diffexp=grav5.vs.osm6,
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
                        draw = '')


diff_all <- list(Microgravity = pull.relevant.columns(g.tT),
                 Osmotic.Stress = pull.relevant.columns(o.tT),
                 Oxidative.Stress = pull.relevant.columns(ox.tT),
                 Heat.Stress = pull.relevant.columns(h.tT))

all_rem <- rem_mv(diffexp=diff_all,
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
                  draw = '')

save(all_rem, rem_g5.v.he7, rem_g5.v.osm6, rem_g5.v.ox6, diff_all, grav5.vs.ox6, grav5.vs.he7, grav5.vs.osm6, file = "data/appdata.RData")



