
setwd("~/Documents/IFN_enhancer/R/")
library(GenomicRanges)
library(GenomicFeatures)
library(gplots)
library(rtracklayer)

# Discordant/concordant enhancer list

# Load RPKM matrix
load("Data/eRNA_expression_GM_unfiltered.RData")
gene.rpkm <- read.table("Data/Expression_matrix.symbol.csv")

rpkm <- rbind(eRNA.rpkm, gene.rpkm)

#################################################################################################################

# Heatmaps
plot_heatmap <- function(l_dc, l_cc) {
  gene_list <- c(l_dc, l_cc)
  heat.data <- rpkm[gene_list, ]
  # heat.data <- log2((heat.data+0.25) /(heat.data[,1]+0.25))
  # heat.data <- (heat.data - rowMeans(heat.data)) / apply(heat.data, 1, sd)
  heat.data <- (heat.data - heat.data[, 1]) / apply(heat.data, 1, sd)
  heat.data <- as.matrix(heat.data)
  for (i in c(1:length(gene_list))) {
    rownames(heat.data)[i] <- substr(gene_list[i], 14, 1000L)
  }
  for (i in c(1:ncol(heat.data))) {
    colnames(heat.data)[i] <- substr(colnames(heat.data)[i], 3, 1000L)
  }
  my_palette <- colorRampPalette(c("blue", "white", "yellow"), bias = 1)(n = 100)
  breaks <- c(seq(-1, 3, length.out = 101))
  heatmap.2(heat.data,
    col = my_palette, breaks = breaks,
    Colv = F, Rowv = F, dendrogram = "none", labRow = "",
    RowSideColors = c(rep("darkgrey", length(l_dc)), rep("lightgrey", length(l_cc))),
    keysize = 1.2, density.info = "density", key.ylab = "", key.xlab = "", key.title = "",
    trace = "none", scale = "none"
  )
}
plot_heatmap(dc[, 2], cc[, 2])

amplitude_index <- function(x, abs = T) {
  fc <- log2((x[, -1] + 0.1) / (x[, 1] + 0.1))
  if (abs) {
    return(apply(abs(fc), 1, max))
  }
  else {
    lab <- apply(abs(fc), 1, function(x) {
      xx <- which(x == max(x))
      return(xx[1])
    })
    tp <- NULL
    for (i in 1:length(lab)) {
      tp <- c(tp, fc[i, lab[i]])
    }
    names(tp) <- rownames(x)
    return(tp)
  }
}

get_idx <- function(pairs, method = 1) {
  # Method 1:
  if (method == 1) {
    res <- c()
    for (i in 1:nrow(pairs)) {
      res <- c(res, cor(t(rpkm[pairs[i, 1], ]), t(rpkm[pairs[i, 2], ]), method = "spearman"))
    }
    return(res)
  }
  # Method 2:
  if (method == 2) {
    res <- c()
    for (i in 1:nrow(pairs)) {
      res <- c(res, max(as.vector(rpkm[pairs[i, 2], 5:9])) / max(as.vector(rpkm[pairs[i, 2], 10:12])))
    }
    return(res)
  }
  # Method 3:
  if (method == 3) {
    res <- c()
    for (i in 1:nrow(pairs)) {
      res <- c(res, apply(rpkm[pairs[i, 2], 5:9], 1, mean) / apply(rpkm[pairs[i, 2], 10:12], 1, mean))
    }
    return(res)
  }
}

idx_dc <- get_idx(dc, 1)
idx_cc <- get_idx(cc, 1)
hist(idx_dc, col = rgb(0.7, 0, 0, 0.5), xlim = range(c(idx_dc, idx_cc)), breaks = 20)
hist(idx_cc, col = rgb(0, 0, 0.7, 0.5), add = T, breaks = 60)

pairs <- read.csv(file = "Data/Induced_EP.csv", colClasses = "character")[, c("symbol", "id_b")]
pairs <- pairs[!(is.na(pairs[, 1])), ]
ai_e <- amplitude_index(rpkm[pairs[, 2], 1:5], F)
ai_g <- amplitude_index(rpkm[pairs[, 1], 1:5], F)
pairs <- pairs[(ai_e > 1) & (ai_g > 1), ]

idx <- get_idx(pairs, 1)
dc <- pairs[idx < quantile(idx, 0.3), ]
cc <- pairs[idx > quantile(idx, 0.7), ]
plot_heatmap(dc[, 1], cc[, 1])
plot_heatmap(dc[, 2], cc[, 2])

idx[pairs[, 1] == "IFNB1"]

# dc <- read.table("~/Documents/eRNA/link/IFN project/discordant_pairs_gene-control-enhancer.pairs.txt", header = F, colClasses = 'character')
# cc <- read.table("~/Documents/eRNA/link/IFN project/concordant_pairs_gene-control-enhancer.pairs.txt", header = F, colClasses = 'character')
pdf("Figure/Discordant_concordant_pairs_gene.pdf")
plot_heatmap(dc[, 1], cc[, 1])
dev.off()
pdf("Figure/Discordant_concordant_pairs_eRNA.pdf")
plot_heatmap(dc[, 2], cc[, 2])
dev.off()
