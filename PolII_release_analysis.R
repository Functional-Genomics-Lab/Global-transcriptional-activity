setwd("~/Documents/IFN_enhancer/R/")
library(GenomicRanges)
library(GenomicFeatures)
library(gplots)
library(rtracklayer)
library("cluster")
source("Script/ref_2_symbol_matrix.R")

######################################################
# 0. functions
######################################################
amplitude_index <- function(x) {
  fc <- log2((x[, -1] + 0.1) / (x[, 1] + 0.1))
  return(apply(fc, 1, max))
}
continuity_index <- function(x) {
  return(diag(cor(t(x[, -ncol(x)]), t(x[, -1]), method = "spearman")))
}
interpolation <- function(xa) {
  xb <- (xa[, 1:(ncol(xa) - 1)] + xa[, 2:ncol(xa)]) / 2
  colnames(xb) <- paste(colnames(xa)[1:(ncol(xa) - 1)], colnames(xa)[2:ncol(xa)], sep = "-")
  x <- cbind(xa, xb)
  lab <- 1
  for (i in 2:ncol(xa)) {
    lab <- c(lab, ncol(xa) + i - 1, i)
  }
  return(x[, lab])
}
# Function to plot color bar
color.bar <- function(lut, min, max = -min, nticks = 11, ticks = seq(min, max, len = nticks), title = "") {
  scale <- (length(lut) - 1) / (max - min)
  dev.new(width = 1.75, height = 5)
  plot(c(0, 10), c(min, max), type = "n", bty = "n", xaxt = "n", xlab = "", yaxt = "n", ylab = "", main = title)
  axis(2, ticks, las = 1)
  for (i in 1:(length(lut) - 1)) {
    y <- (i - 1) / scale + min
    rect(0, y, 10, y + 1 / scale, col = lut[i], border = NA)
  }
}
##############################################################
##  1. Calculate Pol II Release Ratio (PRR)
##############################################################

# Load raw counts
load("Data/gene_expression_GM_unfiltered.RData")
colnames(gene.rpkm) <- substr(colnames(gene.rpkm), 3, 1000L)
gene_list <- rownames(gene.rpkm)[apply(gene.rpkm, 1, max) > 1]

pmt <- read.table("Data/pmt_GM_filtered.counts")
gbd <- read.table("Data/gbd_GM_filtered.counts")

gbd_len <- NULL
tp <- read.table("Data/refseq_groseq.genebody.bed")
for (i in 1:nrow(tp)) {
  gbd_len[as.character(tp[i, 4])] <- (tp[i, 3] - tp[i, 2]) / 1000
}
gbd <- gbd / gbd_len[rownames(gbd)]

gene_list <- gene_list[gene_list %in% rownames(pmt)]
gene_list <- gene_list[apply(pmt[gene_list, ], 1, min) > 10]
pmt <- pmt[gene_list, -c(1, 7)]
gbd <- gbd[gene_list, -c(1, 7)]

prr <- (1 + gbd) / (1 + pmt)
colnames(prr) <- substr(colnames(prr), 3, 1000L)

######################################################
# 2. Analyze PRR dynamics
######################################################

x <- prr[, -c(7, 10)]
# Linear interpolation： increase the number of data point to get robustness
# xa <- prr[,-c(7, 10)]
# x <- interpolation(xa)

ci <- continuity_index(x)
ai <- amplitude_index(x)

rpkm <- gene.rpkm
rpkm <- rpkm[rownames(x), colnames(x)]
rpkm <- interpolation(rpkm)
ai_rpkm <- amplitude_index(rpkm)
ci_rpkm <- continuity_index(rpkm)
ind_rpkm <- ai_rpkm * ci_rpkm

n_bin <- 50
cuts <- log(1:(n_bin + 1), base = 2)
cuts <- cuts / (max(cuts) - min(cuts)) * (floor(max(ind_rpkm)) - min(ind_rpkm)) + min(ind_rpkm)
lab <- as.integer(cut(ind_rpkm, cuts))
lab[ind_rpkm == min(ind_rpkm)] <- 1
lab[isNA(lab)] <- n_bin

col_scale <- lab / n_bin
col_scale <- (col_scale - min(col_scale)) / (max(col_scale) - min(col_scale))
plot(ci, ai,
  col = rgb(col_scale, 0, 1 - col_scale, col_scale),
  main = "Pol II Release Ratio Dynamics", xlab = "Continuity index", ylab = "Amplitude index",
  pch = 20, cex = 1.5
)
show_scale <- c(0:5 / 5)
legend("topleft",
  col = rgb(show_scale, 0, 1 - show_scale, show_scale), pch = 20,
  legend = c("Gene Induction", round(quantile(cuts, show_scale[-1]), 2))
)

a <- names(ci)[ci < (-0.6) & ind_rpkm > 3]
text(ci[a], ai[a] - 0.1, "Annotation issue")

ci_thred <- 0.2
ai_thred <- log2(1.5)

abline(v = ci_thred, lty = "dashed")
abline(h = ai_thred, lty = "dashed")

# Count number of genes in each quadrant
quadrant <- rep(1, length(ci))
names(quadrant) <- names(ci)
quadrant[(ci > ci_thred & ai < ai_thred)] <- 2
quadrant[(ci < ci_thred & ai < ai_thred)] <- 3
quadrant[(ci < ci_thred & ai > ai_thred)] <- 4

prr_stats <- table(quadrant)
prr_stats <- NULL
for (i in 1:4) {
  prr_stats <- c(
    prr_stats,
    sum(quadrant == i),
    sum(quadrant == i & ind_rpkm > 2),
    sum(quadrant == i & ind_rpkm > 3)
  )
}
prr_stats <- matrix(prr_stats, nrow = 4, ncol = 3, byrow = T)

######################################################
# 3. Genes with stimulated elongation
######################################################

genes_elongation <- names(ci)[quadrant == 1 & ind_rpkm > log2(1.5)]
boxplot(prr[genes_elongation, ], ylab = "Pol II Release Ratio", ylim = c(0, 5), col = rainbow(10), outline = F)
boxplot(gene.rpkm[genes_elongation, -c(1, 5, 7)], ylab = "RPKM", ylim = c(0, 10), col = rainbow(10), outline = F)

tp <- read.table("Data/refseq_groseq.bed")
genes_elongation <- tp[tp[, 4] %in% genes_elongation, ]
write.table(genes_elongation, file = "Data/Elong_genes.bed", row.names = F, col.names = F, quote = F, sep = "\t")

# Question: which genes were de-novo activation, without pausing signals in the beginning?

################################################################

# After running Deeptools, we get the metagene profiles for elongated genes

y <- read.table("Heatmap/PRR_early.profile", skip = 2)
n <- nrow(y)
rownames(y) <- substr(as.character(y[, 1]), 3, 1000L)
y <- y[, -c(1, 2)]
tss <- 101
tts <- 600
norm <- apply(y[, 0:49 + tss], 1, max)
y_norm <- y / norm


x_lim <- c(-100, 699)
y_lim <- c(0, 1.25)
plot(-10^5, -1, xlim = x_lim, ylim = y_lim, xlab = "", ylab = "Normalized RPKM", xaxt = "n")
axis(side = 1, at = c(-100, 0, 499, 699), labels = c("-1kb", "TSS", "TES", "2kb"))
use_row <- c(1, 2, 4, 6, 7)
for (i in use_row) {
  lines(c(x_lim[1]:x_lim[2]), y_norm[i, ], col = rainbow(n)[i], lwd = 1.2)
}

legend("topright", col = rainbow(n)[use_row], rownames(y)[use_row], lty = "solid", lwd = 1.5)

















# Heatmap of induced genes identified
# Sort rows by kmeans clustering
# heat.data <- gene.rpkm[names(ci)[quadrant==1 & ind_rpkm>1], ]
heat.data <- gene.rpkm[names(ci)[ci > 0.5 & ind_rpkm > 1], ]
heat.data <- ref_2_symbol(heat.data)
heat.data <- heat.data / apply(heat.data, 1, max)

# Data clustering
dis_matrix <- (1 - cor(t(heat.data), method = "spearman")) / 2
set.seed(0)
centers <- 2
kmedoids <- pam(dis_matrix, k = centers, diss = T) # Method: Partitioning Around Medoids

k_cluster <- kmedoids$clustering
cluster_order <- c(1, 2)

ind <- cluster_order[k_cluster] * 10^4
for (i in 1:centers) {
  lab <- which(k_cluster == i)
  ind[lab] <- ind[lab] - dis_matrix[kmedoids$medoids[i], lab]
}
od <- order(ind, decreasing = F)
heat.data <- heat.data[od, ]
k_cluster <- k_cluster[od]

my_palette <- my_palette <- colorRampPalette(c("blue", "white", "red"), bias = 1)(n = 20)
heatmap.2(heat.data,
  col = my_palette, Colv = F, Rowv = F, dendrogram = "none",
  keysize = 1.5, density.info = "density", labRow = NULL,
  key.ylab = "", key.xlab = "", key.title = "", trace = "none"
)
