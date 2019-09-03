# 3C analysis for IFN family

setwd("~/Documents/IFN_enhancer/R/")
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

time <- c("0h", "30min", "1h", "2h", "4h", "6h", "9h", "12h", "18h", "24h", "48h", "72h")
# time <- c("0h", "6h", "12h", "18h", "24h")
x <- read.table("Data/3C_IFN_Family.txt", header = T)
x <- t(x)
x <- x[1:17, ]
IFN_X <- c(rownames(x), "CDKN2A")
IFN_X <- sym_2_ref(IFN_X)
lab <- which(IFN_X[, 1] %in% rownames(gene.rpkm))
IFN_X <- IFN_X[lab, ]
ifn_exp <- gene.rpkm[IFN_X[, 1], paste("GM", time, sep = "")]
rownames(ifn_exp) <- IFN_X[, 2]
rownames(ifn_exp)[16:18] <- paste("CDKN2A", c(1, 2, 3), sep = "")
colnames(ifn_exp) <- time
save(IFN_X, ifn_exp, file = "~/Documents/IFN_enhancer/R/Data/IFN_family_rpkm.RData")


heat.data <- ifn_exp[rowMax(ifn_exp) > 1, ]
heat.data <- heat.data[sort(rownames(heat.data)), ]
heat.data <- (heat.data - rowMeans(heat.data)) / apply(heat.data, 1, sd)

par(mar = c(5, 4, 4, 2))
my_palette <- colorRampPalette(c("blue", "white", "red"), bias = 1)(n = 500)
breaks <- c(seq(-3, 3, length.out = 501))
heatmap.2(heat.data,
  col = my_palette, breaks = breaks,
  Colv = F, Rowv = T, dendrogram = "none", labRow = NULL,
  keysize = 1.2, density.info = "density", key.ylab = "", key.xlab = "", key.title = "",
  trace = "none"
)
