# RUV RNA-seq pipeline 1.0
# edgeR pipeline is recommended when there are 3,4 biological replicates per condition
# limma-voom is recommended when there are more than 5,6 biological replicates per condition

library(RUVSeq)
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)
library(RColorBrewer)
library(pamr)

# Retrieve gene quantification matrix
files <- list.files("./genelevel/trim_HSC_only") #retrieve all files in female folder.
# feature Count
x <- readDGE(files, path="./genelevel/trim_HSC_only", columns=c(1,2), sep=" ") #readDGE Reads and merges a set of text files containing gene expression counts.
class(x)
dim(x)


samplenames <- substring(colnames(x), 0, 23)
samplenames
colnames(x) <- samplenames
total_samples = 8

## ----filter lowly expressed genes ##
cpm <- cpm(x)
table(rowSums(x$counts>0)==total_samples)
keep.exprs <- rowSums(cpm>2.0)>=3
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

## ----create sample group for DE analysis---- ##
expo <- as.factor(rep(c("NR","WT"), c(4,4)))
diet <- as.factor(rep(c(""), c(8)))
method <- as.factor(rep(c("HSC"), c(8)))
age <- as.factor(rep(c("O"), 
                     c(8)))
x$samples$expo <- expo
x$samples$diet <- diet
x$samples$method <- method
x$samples$age <- age

group <- factor(paste(method,age,expo, sep=""))
lcpm <- cpm(x, log=FALSE)
genes <- rownames(x)


## ----store_data-------------------------------------------
set <- newSeqExpressionSet(as.matrix(x),
                           phenoData = data.frame(group, row.names=colnames(x)))



## ----rle, fig.cap="No normalization.",fig.subcap=c("RLE plot","PCA plot")----
library(RColorBrewer)
par(mfrow=c(1,2))
colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[group], cex=0.8)
plotPCA(set, col=colors[group], cex=0.8)

## ----uq, fig.cap="Upper-quartile normalization.", fig.subcap=c("RLE plot","PCA plot")----
set <- betweenLaneNormalization(set, which="median")
par(mfrow=c(1,2))
colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[group], cex=0.8)
plotPCA(set, col=colors[group], cex=0.5)

## ----uq, fig.cap="Upper-quartile normalization.", fig.subcap=c("RLE plot","PCA plot")----
design <- model.matrix(~group, data=pData(set))
y <- DGEList(counts=counts(set), group=group)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set))$table

summary(de <- decideTestsDGE(lrt, adjust.method="BH", p.value=0.05, lfc=1))

detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
plotMD(lrt, pch=16, cex=1.2)
abline(h=c(-1.5, 1.5), col="blue")

tab <- topTags(lrt, n=2000, sort.by="p.value", p.value=0.10)
tab_all <- topTags(lrt, n=20000, sort.by="p.value", p.value=1)
summary(de <- decideTestsDGE(lrt, adjust.method="BH", p.value=0.2, lfc=0.583))
tab <- topTags(lrt, n=20000, sort.by="p.value", p.value=1.0)
write.table(tab, file="mygenelist_cohort3_HSC_noRUVg.txt", sep="\t")

empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:10000]))]

## ----emp_ruvg, fig.cap="RUVg normalization based on empirical controls.", fig.subcap=c("RLE plot","PCA plot")----
set2 <- RUVg(set, empirical, k=1)
par(mfrow=c(1,2))
pData(set2)
plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[group])
plotPCA(set2, col=colors[group], cex=0.5)

# example of pair-wise comparison
design <- model.matrix(~group + W_1, data=pData(set2))
y <- DGEList(counts=counts(set2), group=group)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
#  
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)

summary(de <- decideTestsDGE(lrt, adjust.method="BH", p.value=0.05, lfc=0.583))

par(mfrow=c(1,1))
detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
plotMD(lrt, pch=16, cex=0.8)
abline(h=c(-0.381, 0.381), col="blue")

tab <- topTags(lrt, n=2000, sort.by="p.value", p.value=0.10)
tab_all <- topTags(lrt, n=20000, sort.by="p.value", p.value=1)
summary(de <- decideTestsDGE(lrt, adjust.method="BH", p.value=0.2, lfc=0.583))
tab <- topTags(lrt, n=20000, sort.by="p.value", p.value=1.0)
write.table(tab, file="mygenelist_cohort3_HSC.txt", sep="\t")

temp <- set2@assayData[["normalizedCounts"]]
temp2 <- cpm(temp, log=TRUE)
#temp2 <- normalize.quantiles(lcpm_normalized, copy=TRUE)
colnames(temp2) <- colnames(temp)
rownames(temp2) <- rownames(temp)
write.table(temp2, file="mygenelist_pairwise_normalized_logCPM_HSC.txt", sep="\t")
