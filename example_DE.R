library(GEOquery)
library(Biobase)
library(limma)

#This pulls the all the samples from the microarry dataset. This returns a list containing a single expressionSet object. 
gset <- getGEO("GSE4136", GSEMatrix =TRUE, AnnotGPL=TRUE)

#if there are multiple datasets just take the GPL file  
if (length(gset) > 1) idx <- grep("GPL2529", attr(gset, "names")) else idx <- 1

#if not just take the expressionSet object
gset <- gset[[idx]]

#take a look into the datase
head(gset)
View(gset)
dim(gset)

#here we format the feature names. fvarLabels belongs the biobase package and is used to extract features from ExpressionSet Objects
fvarLabels(gset) <- make.names(fvarLabels(gset))

#This is an ugly way of choosing groups but I guess it works. Imagine the 0s are the control replicates and the ones are space-flown.
gsms <- "000XXX111XXX"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# order samples by group
#exprs() is a biobase function for retrieving expression data from eSet objects.
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("g","s")

# log2 transform
ex <- exprs(gset)
#this checks if the data is log-transformed
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)


#pulling the whole microarray dataset and ordering by B value 
tT <- topTable(fit2, adjust="fdr", sort.by="B", number = length(fit2[[1]]))

#We can get rid of the "number" argument and replace it with lfc = 2 for a log2 fold change of 2
#tT <- topTable(fit2, adjust="fdr", sort.by="B", lfc = 2)

tT <- subset(tT, select=c("adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Platform_ORF","Gene.title","GO.Function","GO.Process", "GO.Component","Chromosome.annotation"))

#View the final DE expression table
View(tT)

write.table(tT, file = "GSE4136_YEAST_EXAMPLE_DE.csv", row.names = FALSE, sep = ",")



