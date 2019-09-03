# 3C analysis for IFN family

setwd("~/Documents/IFN_enhancer/R")
library(GenomicRanges)
library(GenomicFeatures)
library(gplots)
library(rtracklayer)

load(file = "Data/gene_expression_GM_unfiltered.RData")
load("Data/ref2sym_dedup.RData")
sym_2_ref <- function(x) {
  if (class(x) == "character") {
    return(unique(ref2sym[ref2sym[, 2] %in% x, ]))
  }
}

gene_names <- read.table("Data/Immune_genes_example.txt", colClasses = "character")
gene_refseq <- sym_2_ref(gene_names[, 1])
gene_refseq <- gene_refseq[gene_refseq[, 1] %in% rownames(gene.rpkm), ]
uniq_symbol <- sort(unique(gene_refseq[, 2]))
for (i in 1:length(uniq_symbol)) {
  tp_gene <- gene_refseq[gene_refseq[, 2] == uniq_symbol[i], 1]
  tp <- gene.rpkm[tp_gene, ]
  if (length(tp_gene) > 1) {
    tp <- colMeans(tp)
  }
  if (i == 1) {
    heat.data <- tp
  }
  else {
    heat.data <- rbind(heat.data, tp)
  }
}
rownames(heat.data) <- uniq_symbol
colnames(heat.data) <- c("0h", "30min", "1h", "2h", "4h", "6h", "9h", "12h", "18h", "24h", "48h", "72h")
save(gene_refseq, heat.data, file = "~/Documents/IFN_enhancer/R/Data/Immune_genes_example.RData")

load("~/Documents/IFN_enhancer/R/Data/Immune_genes_example.RData")
heat.data <- heat.data[!(rownames(heat.data) %in%
  c("BCL2", "IFI16", "TBK1", "APOBEC3F", "APOBEC3G")), ] # Clean up data
heat.data <- heat.data[rowMax(heat.data) > 0.1, ]
heat.data <- (heat.data - rowMeans(heat.data)) / apply(heat.data, 1, sd) # Standardize data

par(mar = c(5, 4, 4, 2))
my_palette <- colorRampPalette(c("blue", "white", "red"), bias = 1)(n = 500)
breaks <- c(seq(-3, 3, length.out = 501))
heatmap.2(heat.data,
  col = my_palette, breaks = breaks,
  Colv = F, Rowv = T, dendrogram = "none", labRow = NULL,
  keysize = 1.2, density.info = "density", key.ylab = "", key.xlab = "", key.title = "",
  trace = "none"
)
