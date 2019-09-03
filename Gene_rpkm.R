setwd("~/Documents/IFN_enhancer/R/")
library(GenomicRanges)
library(GenomicFeatures)
library(gplots)
library(rtracklayer)
source("Script/ref_2_symbol_matrix.R")

##############################################################
##  1. Enhancer gene expression
##############################################################

# Load raw counts
x <- read.table("Data/genes_GM_filtered.counts", header = T)
y <- read.table("Data/genes.bed")
gene <- GRanges(seqnames = y[, 1], strand = y[, 6], IRanges(y[, 2], y[, 3]), id = y[, 4])

gene.length <- width(gene)
names(gene.length) <- mcols(gene)$id

# Normalize by total reads and enhancer length
y <- read.table("Data/total_mapped_reads_GM.txt")
total_ct <- y[, 2]
names(total_ct) <- y[, 1]
x <- x[rownames(x) %in% mcols(gene)$id, ]
gene.rpkm <- t(t(x) / total_ct[colnames(x)]) / gene.length[rownames(x)] * 10^9
save(gene, gene.rpkm, file = "Data/gene_expression_GM_unfiltered.RData")
symbol.rpkm <- ref_2_symbol(gene.rpkm)
write.table(round(symbol.rpkm, 2), file = "Data/Expression_matrix.symbol.csv", quote = F, sep = "\t")

mean <- rowMeans(gene.rpkm)
lab <- which(mean > 1)
gene.rpkm <- gene.rpkm[lab, ]
gene <- gene[mcols(gene)$id %in% rownames(gene.rpkm)]
gene.rpkm <- gene.rpkm[as.character(mcols(gene)$id), ]

save(gene, gene.rpkm, file = "Data/gene_expression_GM.RData")


mean <- rowMeans(gene.rpkm)
group <- cut(mean, quantile(mean, probs = c(0:3) / 3), labels = F)
col <- c("blue", "grey", "red")
lab <- c("Low", "Medium", "High")
par(mfrow = c(3, 1), mai = c(0, 1.5, 0, 1.5))
for (i in 3:1)
{
  boxplot(gene.rpkm[group == i, ], xaxt = "n", border = col[i], ylab = "RPKM")
  legend("topleft", lab[i], text.col = col[i], bty = "n", cex = 1.5)
}
