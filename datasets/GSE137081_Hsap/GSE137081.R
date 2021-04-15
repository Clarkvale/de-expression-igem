library(edgeR)
library(dplyr)
library(GEOquery)

#Reformatting output from featureCount

counts <- read.delim("raw_counts.tabular", row.names = 1)
genes <- rownames(counts)
counts <- na.omit(sapply(counts, as.numeric))
rownames(counts) <- genes[-1]
counts <- as.data.frame(counts) %>% select(-c(7:12))
d0 <- DGEList(counts)

#bad <- grep(colnames(counts), pattern = "SRR10084987")
#counts <- as.data.frame(counts) %>% select(-11)

#Adding Normalizing factors
d0 <- calcNormFactors(d0)

#filtering low expressed genes
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) 


Treat <- factor(c(rep("MG", 6), rep("NG", 6)))
lines <- factor(c(rep(c(1,2,3),4)))

design <- model.matrix(~lines+Treat)
#grouped <- interaction(flight, lines)


plotMDS(d, col = as.numeric(Treat))

#Voom and modeling
#mm <- model.matrix(~0 + grouped)
y <- voom(d, design, plot = TRUE) #looks okay

fit <- lmFit(y, design)
#contrasts <- makeContrasts(flightmicro.gravity - flightnormal.gravity)

fit2 <- eBayes(fit)
tT <- topTable(fit2, coef = "TreatNG", confint = T, number = Inf)

#setting up metadata extraction
gse <- getGEO("GSE137081")
gsms <- names(gse@gsms)[-c(7:12)]
colnames(counts) <- gsms
miame <- MIAME(name = gse@header$name, lab = gse@header$contact_address, 
               contact = gse@header$contact_name, title = gse@header$title,
               abstract = gse@header$summary, url = gse@header$relation[[1]],
               samples = list(gse@header$sample_id), pubMedIds = gse@header$pubmed_id)

eset <- ExpressionSet(assayData = as.matrix(counts), annotation = gse@header$platform_id, 
                      experimentData = miame)
