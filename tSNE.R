# t-SNE of GM eRNA expression profile
setwd("~/Documents/IFN_enhancer/R/")
library(GenomicRanges)
library(GenomicFeatures)
library(gplots)
library(Rtsne)

tSNE_plot <- function(x, seed) {
  set.seed(seed)
  my_palette <- colorRampPalette(c("blue", "red", "orange"), bias = 1)(n = 12)
  # Before t-SNE, normalize features
  x <- t((t(x) - colMeans(x)) / (apply(x, 2, sd)))
  tsne <- Rtsne(x, dims = 2, perplexity = 2, verbose = T, max_iter = 500)
  plot(tsne$Y, col = my_palette, pch = 16, main = "", xaxt = "n", yaxt = "n", xlab = "Dim_1", ylab = "Dim_2")
  lines(tsne$Y, lty = "dashed", col = "grey", lwd = 1.5)
  lab <- substring(rownames(x), first = 3)
  text(tsne$Y, labels = lab, col = my_palette, cex = 1.5)
}

##################################################
# 1. Load data
##################################################

load("Data/gene_expression_GM.RData")
gene.rpkm <- gene.rpkm[, -c(1)]

load("Data/eRNA_expression_GM.RData")
eRNA.rpkm <- eRNA.rpkm[, -c(1)]

##################################################
# 2. t-SNE
##################################################
pdf("Figure/tSNE_mRNA.pdf")
tSNE_plot(t(gene.rpkm[rowMeans(gene.rpkm) > 0.5, ]), 20)
dev.off()
pdf("Figure/tSNE_eRNA.pdf")
tSNE_plot(t(eRNA.rpkm[rowMeans(eRNA.rpkm) > 0.5, ]), 9)
dev.off()
