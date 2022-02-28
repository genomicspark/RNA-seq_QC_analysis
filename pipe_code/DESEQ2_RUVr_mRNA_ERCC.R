# RUVg and RUVr RNA-seq pipeline

library(RUVSeq)
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)
library(RColorBrewer)
library(pamr)

# Retrieve gene quantification matrix
files <- list.files("./fc_isoform") #retrieve all files in female folder.
# feature Count
x <- readDGE(files, path="./fc_isoform", columns=c(1,3), sep="\t") #readDGE Reads and merges a set of text files containing gene expression counts.
class(x)
dim(x)

# Save R object as output (rda)
save(x, file="raw_rnaseq_subset.rda")
# save raw matrix into 

samplenames <- substring(colnames(x), 0, 23)
samplenames
colnames(x) <- samplenames
colnames(x) <- c("FA1","FA2","FA3","FA4","PM1","PM2","PM3","PM4")
total_samples = 8

## ----filter lowly expressed genes ##
cpm <- cpm(x)
table(rowSums(x$counts>0)==total_samples)
keep.exprs <- rowSums(cpm>1.0)>=4
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

## ----create sample group for DE analysis---- ##
expo <- as.factor(rep(c(""), c(8)))
diet <- as.factor(rep(c(""), c(8)))
method <- as.factor(rep(c("",""), c(4,4)))
age <- as.factor(rep(c("FA","PM"), 
                     c(4,4)))
x$samples$expo <- expo
x$samples$diet <- diet
x$samples$method <- method
x$samples$age <- age

group <- factor(paste(age,method, sep=""))
lcpm <- cpm(x, log=FALSE)
genes <- rownames(x)


## ----store_data-------------------------------------------
set <- newSeqExpressionSet(as.matrix(x),
                           phenoData = data.frame(group, row.names=colnames(x)))



## ----rle, fig.cap="No normalization.",fig.subcap=c("RLE plot","PCA plot")----
library(RColorBrewer)
par(mfrow=c(1,2))
colors <- brewer.pal(6, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[group], cex=0.8)
plotPCA(set, col=colors[group], cex=0.8)

## ----uq, fig.cap="Upper-quartile normalization.", fig.subcap=c("RLE plot","PCA plot")----
set <- betweenLaneNormalization(set, which="median")
par(mfrow=c(1,1))
colors <- brewer.pal(6, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[group], cex=0.8)
plotPCA(set, col=colors[group], cex=0.5)

design <- model.matrix(~group, data=pData(set))
y <- DGEList(counts=counts(set), group=group)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set))$table
# RUVr needs res values using deviance
res <- residuals(fit, type="deviance")

summary(de <- decideTestsDGE(lrt, adjust.method="BH", p.value=0.05, lfc=1.0))

detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
#plotMD(lrt, pch=16, cex=1.2)
abline(h=c(-0.5, 0.5), col="blue")

empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]
spikes <- rownames(set)[grep("^ERCC", rownames(set))]
#spikes <- empirical[grep("^ERCC", empirical)]


## ----emp_ruvg, fig.cap="RUVg normalization based on empirical controls.", fig.subcap=c("RLE plot","PCA plot")----
set2 <- RUVg(set, empirical, k=1)
par(mfrow=c(1,1))
pData(set2)
plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[group])
plotPCA(set2, col=colors[group], cex=0.5)

## ----emp_ruvg, fig.cap="RUVg normalization based on empirical controls.", fig.subcap=c("RLE plot","PCA plot")----
set3 <- RUVg(set, spikes, k=1)
par(mfrow=c(1,1))
pData(set3)
plotRLE(set3, outline=FALSE, ylim=c(-4, 4), col=colors[group])
plotPCA(set3, col=colors[group], cex=0.5)

## ----emp_RUVr, fig.cap="RUVg normalization based on empirical controls.", fig.subcap=c("RLE plot","PCA plot")----
set4 <- RUVr(set, genes, k=1, res) #RUVr
par(mfrow=c(1,1))
pData(set4)
colors <- brewer.pal(8, "Set2")
plotRLE(set4, outline=FALSE, ylim=c(-4, 4), col=colors[group])
plotPCA(set4, col=colors[group], cex=0.5)

## DESeq2 implemantation
library(preprocessCore)
cts = counts(set2)
cts_quantile <- normalize.quantiles(cts, copy=TRUE)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = pData(set2),
                              design = ~ group+W_1) 

dds <- DESeq(dds)

# exposure analysis
# group_list1 <- c("14wChowATACFA","14wChowATACFA","14wChowATACFA","14wChowOMNIFA","14wChowOMNIFA","14wChowOMNIPM","22wChowOMNIFA")
# group_list2 <- c("14wChowATACPM","14wHFDATACFA","14wHFDATACPM","14wChowOMNIPM","22wChowOMNIFA","22wChowOMNIRev","22wChowOMNIRev")
#group_list1 <- c("YSH","YSH","YSL","YSL","MSH","MSL","YSH","MSH","OSH")
#group_list2 <- c("MSH","OSH","MSL","OSL","OSH","OSL","YSL","MSL","OSL")
library(preprocessCore)
temp <-   set2@assayData[["normalizedCounts"]]
lcpm_normalized <- cpm(temp, log=TRUE)
temp2 <- normalize.quantiles(lcpm_normalized, copy=TRUE)
colnames(lcpm_normalized) <- colnames(temp)
rownames(lcpm_normalized) <- rownames(temp)
colnames(temp2) <- colnames(temp)
rownames(temp2) <- rownames(temp)
write.table(lcpm_normalized, file="mygenelist_basic_normalized_logCPM.txt", sep="\t")
write.table(temp2, file="mygenelist_quantile_normalized_logCPM.txt", sep="\t")
write.table(temp, file="mygenelist_RUVg_count.txt", sep="\t")

library(plotly)
pca_liver <- prcomp(temp2)
pca_plot_data <- pca_liver$rotation[, c(1,2,3)]
color_code <- c("FA","FA","FA","FA","PM","PM","PM","PM")
name <- rownames(pca_plot_data)
pca_plot_data <- cbind(pca_plot_data, color_code, name)

p <- data.frame(pca_plot_data)
pca_3d_plot <- plot_ly(p, x = ~PC1, y = ~PC2, z = ~PC3, color = ~color_code, colors = c("#FF6633",'#23AB23')) %>%
  add_markers( text = ~name) %>%
  layout(scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3')))

pca_3d_plot

group_list1 <- c("FA")
group_list2 <- c("PM")
column_list1 <- c(1,2,3,4)
column_list2 <- c(5,6,7,8)

cnt <- 0
tot_ind <- c()
for(i in group_list1){
  cnt <- cnt + 1
  filename = paste(group_list1[cnt], group_list2[cnt], sep="-")
  filename
  
  all = paste("DEG2_RUVg",paste(filename,"_all.txt",sep=""),sep="_")
  up = paste("DEG2_RUVg",paste(filename,"_up.txt",sep=""),sep="_")
  down = paste("DEG2_RUVg",paste(filename,"_down.txt",sep=""),sep="_")
  pdf_up = paste("Heatmap2_RUVg",paste(filename,"_up.pdf",sep=""),sep="_")
  pdf_down = paste("Heatmap2_RUVg",paste(filename,"_down.pdf",sep=""),sep="_")
  
  res <- results(dds, contrast=c("group",group_list2[cnt],group_list1[cnt]))
  res$Gene <- rownames(res)
  res$logCPM <- aveLogCPM(counts(set2))
  
  res
  sum(res$padj < 0.05, na.rm=TRUE)
  
  write.table(res, file=all, sep="\t")
  ind_up <- which(res$padj < 0.05 & res$logCPM > 1.5 & res$log2FoldChange > 0.261)
  res_deg <- res[ind_up,]
  write.table(res_deg, file=up, sep="\t")
  tot_ind <- c(tot_ind, ind_up)
  
  ind_dn <- which(res$padj < 0.05 & res$logCPM > 1.5 & res$log2FoldChange < -0.261)
  res_deg <- res[ind_dn,]
  write.table(res_deg, file=down, sep="\t")  
  tot_ind <- c(tot_ind, ind_dn)
  
  library(gplots)
  mypalette <- brewer.pal(11,"RdYlBu")
  morecols <- colorRampPalette(mypalette)
  # Set up colour vector for celltype variable
  
  ind2 <- match(rownames(res)[ind_up], rownames(lcpm))
  sub_set = temp2[ind2, c(column_list1,column_list2)]
  #col.cell <- c("black","red","brown","black","red","brown","black","red","brown")[group]
  #col.cell = col.cell[1:9]
  
  pdf(pdf_up) 
  # Plot the heatmap
  heatmap.2(sub_set,col=rev(morecols(50)),trace="none", 
            main=paste("Differentially Expressed Genes",up,sep="\n"),
            scale="row",
            Rowv=TRUE,Colv=FALSE)
  dev.off()
  ind2 <- match(rownames(res)[ind_dn], rownames(lcpm))
  sub_set = temp[ind2, c(column_list1,column_list2)]
  #col.cell <- c("black","red","brown","black","red","brown","black","red","brown")[group]
  #col.cell = col.cell[1:9]
  
  pdf(pdf_down) 
  # Plot the heatmap
  heatmap.2(sub_set,col=rev(morecols(50)),trace="none", 
            main=paste("Differentially Expressed Genes",down,sep="\n"),
            scale="row",
            Rowv=TRUE,Colv=FALSE)
  dev.off()
}

deg <- read.table("file_list_deg2_up.txt", header=FALSE, sep="\t", row.names=1)
colnames(deg)<-colnames(temp2)
library(gplots)
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
col.cell <- c("black","red","brown","green")[group]
col.cell = col.cell[1:total_samples]
# Plot the heatmap
heatmap.2(as.matrix(deg),col=rev(morecols(50)),trace="none", 
          main="Differentially expressed genes",
          ColSideColors=col.cell,scale="row",
          Rowv=TRUE,Colv=FALSE)

########## data visualization ##############
library(gplots)
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable

ind2 <- match(rownames(res)[tot_ind], rownames(lcpm))
sub_set = lcpm[ind2, 1:total_samples]
col.cell <- c("black","red","brown","green")[group]
col.cell = col.cell[1:total_samples]

# Plot the heatmap
heatmap.2(sub_set,col=rev(morecols(50)),trace="none", 
          main="Differentially expressed genes",
          ColSideColors=col.cell,scale="row",
          Rowv=TRUE,Colv=TRUE)
############################################

