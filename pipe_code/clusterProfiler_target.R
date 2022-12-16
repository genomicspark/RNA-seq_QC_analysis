library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
head(gene.df)

# ont: MF, BP, CC = Go term 
ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)
dotplot(ego)

# Young DEG (up*down)
geneList <- read.delim("blood_DEG_after_consortium_biswal_male.csv", sep=",", header=TRUE)
universe_gene <- unlist(geneList[[2]])
universe_gene.df <- bitr(universe_gene, fromType = "SYMBOL",
                         toType = c("ENSEMBL","ENTREZID"),
                         OrgDb = org.Mm.eg.db)
head(universe_gene.df)
universe_gene <- unlist(universe_gene.df[[3]])

geneList <- read.delim("blood_DEG_after_consortium_biswal_male_up.csv", sep=",", header=TRUE)
gene <- unlist(geneList[[2]])
gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL","ENTREZID"),
                OrgDb = org.Mm.eg.db)
head(gene.df)
gene <- unlist(gene.df[[3]])

ego_up <- enrichGO(gene       = gene,
                    universe      = universe_gene,
                    OrgDb         = org.Mm.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.1,
                    qvalueCutoff  = 0.25,
                    readable      = TRUE)
head(ego_up)
dotplot(ego_up)

geneList <- read.delim("blood_DEG_after_consortium_mutlu_male.csv", sep=",", header=TRUE)
universe_gene <- unlist(geneList[[2]])
universe_gene.df <- bitr(universe_gene, fromType = "SYMBOL",
                toType = c("ENSEMBL","ENTREZID"),
                OrgDb = org.Mm.eg.db)
head(universe_gene.df)
universe_gene <- unlist(universe_gene.df[[3]])

geneList <- read.delim("blood_DEG_after_consortium_mutlu_male_dn.csv", sep=",", header=TRUE)
gene <- unlist(geneList[[2]])
gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL","ENTREZID"),
                OrgDb = org.Mm.eg.db)
head(gene.df)
gene <- unlist(gene.df[[3]])

ego2_up <- enrichGO(gene       = gene,
                universe      = universe_gene,
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.1,
                qvalueCutoff  = 0.25,
                readable      = TRUE)
head(ego2_up)
dotplot(ego2_up)

## common DEG
geneList <- read.delim("blood_DEG_after_consortium_mutlu_male.csv", sep=",", header=TRUE)
universe_gene <- unlist(geneList[[2]])
universe_gene.df <- bitr(universe_gene, fromType = "SYMBOL",
                         toType = c("ENSEMBL","ENTREZID"),
                         OrgDb = org.Mm.eg.db)
head(universe_gene.df)
universe_gene <- unlist(universe_gene.df[[3]])

geneList <- read.delim("jhu-chicaco-common-blood-down-DEGs.csv", sep=",", header=TRUE)
gene <- unlist(geneList[[2]])
gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL","ENTREZID"),
                OrgDb = org.Mm.eg.db)
head(gene.df)
gene <- unlist(gene.df[[3]])

ego2_up <- enrichGO(gene       = gene,
                    universe      = universe_gene,
                    OrgDb         = org.Mm.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
head(ego2_up)
dotplot(ego2_up, showCategory = 5)



