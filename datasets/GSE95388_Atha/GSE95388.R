library(GEOquery)
library(Biobase)
library(limma)

#This pulls the all the samples from the microarry dataset. This returns a list containing a single expressionSet object. 
gset <- getGEO("GSE95388", GSEMatrix =TRUE, AnnotGPL=TRUE)
fvarLabels(gset) <- make.names(fvarLabels(gset))

#extracting metadata

#organism 
org <- gset$GSE95388_series_matrix.txt.gz$organism_ch1[[1]]

#platform
plat <- gset$GSE95388_series_matrix.txt.gz$platform_id

#tissue type
tissue <- gset$GSE95388_series_matrix.txt.gz$`tissue source:ch1`[[1]]



gset <- gset[[1]]

#View(gset[[1]])

control <- c(6,7,8,9,10)
treatment <- c(1,2,3,4,5)
group.set <- cbind(control, treatment)
groups_repr <- c()

groups_repr[which(group.set == control)] <- "1"
groups_repr[which(group.set == treatment)] <- "2"
fl <- as.factor(groups_repr)

##NEEEDS NORMALIZATION

gset <- gset[, cbind(control,treatment)]

ex <- exprs(gset)[, cbind(control,treatment)]
exprs(gset) <- log2(ex)

labels <- c("control","treated")

# set up the data and proceed with analysis
sml <- paste("G", groups_repr, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G2-G1, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)


tT <- topTable(fit2, adjust="fdr", sort.by="B", number = length(fit2[[1]]))
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
