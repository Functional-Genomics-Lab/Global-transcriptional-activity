# PCA of GM eRNA expression profile
setwd("~/Documents/IFN_enhancer/R/")
library(GenomicRanges)
library(GenomicFeatures)
library(gplots)

PCA_plot <- function(x) {
  # Before PCA, normalize features
  x <- t((t(x) - colMeans(x)) / (apply(x, 2, sd)))

  # PCA
  svd <- svd(x)

  # Decide number of dimensions reduced
  dsquare <- svd$d^2
  par(mfrow = c(1, 1), mai = c(2, 2, 1, 1))
  plot(dsquare / sum(dsquare), xlab = "Principal component", ylab = "Variance explained")
  i <- 1
  while (sum(dsquare[1:i]) / sum(dsquare) < 0.95) {
    i <- i + 1
  }
  vreduce <- svd$v[, 1:i]
  z <- x %*% vreduce

  # Plot principle components
  my_palette <- colorRampPalette(c("blue", "red", "orange"), bias = 1)(n = 12)
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
  plot(z[, 1], z[, 2], col = my_palette, pch = 16, xaxt = "n", yaxt = "n", xlab = "PC1", ylab = "PC2", cex = 1.3)
  lines(z[, 1], z[, 2], lty = "dashed", col = "grey", lwd = 1.5)
  text(z[, 1], z[, 2], substr(colnames(gene.rpkm), 3, 1000L), pos = c(4, 4, 4, 4, 2, 4, 4, 3, 4, 2, 2), col = my_palette, cex = 1.5)
  return(z)
}

##################################################
# 1. Load data
##################################################
load("Data/gene_expression_GM.RData")
gene.rpkm <- gene.rpkm[, -c(1)]

load("Data/eRNA_expression_GM.RData")
eRNA.rpkm <- eRNA.rpkm[, -c(1)]

##################################################
# 2. PCA
##################################################
pdf("Figure/PCA_mRNA.pdf")
mpca <- PCA_plot(t(gene.rpkm[rowMeans(gene.rpkm) > 0.5, ]))
dev.off()
cc <- abs(cor(mpca[, 1], t(gene.rpkm))[1, ])
write.table(names(sort(cc, decreasing = T))[1:500], "Data/PCA_pc1_mRNA.txt", quote = F, row.names = F, col.names = F)

cc <- abs(cor(mpca[, 2], t(gene.rpkm)))[1, ]
write.table(names(sort(cc, decreasing = T))[1:200], "Data/PCA_pc2_mRNA.txt", quote = F, row.names = F, col.names = F)
pdf("Figure/PCA_eRNA.pdf")
epca <- PCA_plot(t(eRNA.rpkm[rowMeans(eRNA.rpkm) > 0.5, ]))
dev.off()
