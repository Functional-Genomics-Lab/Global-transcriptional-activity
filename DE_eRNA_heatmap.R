setwd("~/Documents/eRNA/link/IFN project/")
library(GenomicRanges)
library(GenomicFeatures)
library(edgeR)
library(DESeq2)
library(cluster)
library(gplots)
library(RColorBrewer)
source("~/Documents/eRNA/Rscript/ref_2_symbol_matrix.R")
source("~/Documents/eRNA/Rscript/sym_2_ref.R")
load("ref2sym_dedup.RData")
find_neighbor <- function(center, target, flanking = 200000) {
  mid <- floor((start(target) + end(target)) / 2)
  start(target) <- mid
  end(target) <- mid + 1
  tp <- center
  center <- floor((start(center) + end(center)) / 2)
  start(tp) <- center - flanking
  end(tp) <- center + flanking
  hits <- findOverlaps(tp, target, ignore.strand = T)
  if (length(hits) > 0) {
    return(as.character(mcols(target[subjectHits(hits)])$id))
  }
  else {
    return(NA)
  }
}
plot_eRNA_gene <- function(erna, mrna) {
  par(mar = c(4, 5, 4, 5))
  a <- eRNA[mcols(eRNA)$id %in% erna]
  b <- gene[mcols(gene)$id %in% sym_2_ref(mrna)]
  tp <- c(start(a), end(a), start(b), end(b))
  region <- paste(as.character(seqnames(a)), min(tp), max(tp), sep = " ")
  x <- eRNA.rpkm[erna, 1:12]
  y <- gene.rpkm[mrna, 1:12]
  xmax <- ceiling(max(x))
  ymax <- ceiling(max(y))
  plot(x / xmax,
    type = "b", col = "red", xaxt = "n", yaxt = "n", main = mrna, xlab = "", ylab = "",
    ylim = c(0, 1), pch = 21
  )
  lines(y / ymax, col = "blue", pch = 22, type = "b")
  axis(side = 1, at = c(1:12), labels = substr(colnames(gene.rpkm)[c(1:12)], 3, 1000L), cex.axis = 0.8)
  axis(side = 2, at = c(0:4) / 4, labels = c(0:4) / 4 * xmax, col.axis = "red", las = 2)
  axis(side = 4, at = c(0:4) / 4, labels = c(0:4) / 4 * ymax, col.axis = "blue", las = 2)
  # legend("topleft", c(erna, "", region), bty="n")
  mtext(text = "eRNA RPKM", side = 2, line = 3, col = "red")
  mtext(text = "Gene RPKM", side = 4, line = 3, col = "blue")
}
enhancer <- function(number, strand = "both") {
  if (strand == "both") {
    lab <- paste("MetaEnhancer", number, c("+", "-"), sep = "_")
  }
  if (strand %in% c("+", "-")) {
    lab <- paste("MetaEnhancer", number, strand, sep = "_")
  }
  return(as.data.frame(eRNA[mcols(eRNA)$id %in% lab]))
}

###############################################################
# 1. Load data
###############################################################

load("gene_expression_GM_unfiltered.RData")
tss <- promoters(gene, upstream = 1, downstream = 1)
load("genes_expression_symbol.RData")
load("eRNA_expression_GM_unfiltered.RData")

load("deg_matrix_eRNA_GM.RData")
# load("deg_matrix_GM_compare_with_30min.RData")
de_genes <- rownames(deg_matrix)[rowSums(deg_matrix[, 1:9]) > 0]
de_mRNA <- grep("MetaEnhancer", de_genes, invert = T, value = T)
de_eRNA <- grep("MetaEnhancer", de_genes, invert = F, value = T)
flanking <- 100000

###############################################################
# 2. Differential expression analysis (You only need to do this once)
###############################################################
# y <- read.table("eRNA_GM_replicate.counts", header=T)
# total <- read.table("total_mapped_reads_GM_replicate.txt")[,2]
# #load("eRNA_GM.RData")
# countData <- y
# rpkmData <- t(t(y/width(eRNA))/total)*10^9
# deseq_padj_cutoff <- 0.01
#
# deg_list <- NULL
# deg_ct <- NULL
# tp <- unlist(strsplit(colnames(countData), "_"))[c(0:11*4+1)]
# groups <- factor(rep(tp, each=2))
# colData <- data.frame(groups=groups, row.names = colnames(countData))
# dds <- DESeqDataSetFromMatrix(countData = countData,
#                               colData = colData,
#                               design = ~ groups)
# dds <- dds[ rowSums(rpkmData > 1) > 2, ] # Pre-filtering
# dds <- DESeq(dds, betaPrior=FALSE )
# norm.counts <- counts(dds, normalized = T)
#
# deg_matrix <- matrix(0, nrow=nrow(countData), ncol=10, dimnames = list(rownames(countData), as.character(tp)[3:12]))
# for(i in 3:12){
#   res <- results(dds, contrast=c("groups",as.character(tp[i]),"GM30min"))
#   lab <- rownames(norm.counts)[which(res$padj < deseq_padj_cutoff)]
#   deg_matrix[lab, i-2] <- 1
# }
# save(norm.counts, countData, deg_matrix, file="deg_matrix_eRNA_GM.RData")
#
# load("deg_matrix_eRNA_GM.RData")
# deg_matrix <- deg_matrix[rownames(deg_matrix) %in% rownames(norm.counts), ]
# deg_matrix_direction <- deg_matrix
# initial <- log2(rowMeans(norm.counts[,1:2])+1)
# for(i in 1:ncol(deg_matrix)){
#   lab <- which(deg_matrix[,i]==1)
#   a <- log2(rowMeans(norm.counts[lab, i*2 + c(1:2)])+1)
#   b <- initial[lab]
#   deg_matrix_direction[lab[a-b<0], i] <- -1
# }
# count_deg <- function(x){return(c(
#   sum(x==1), -sum(x==-1)
# ))}
# tp <- apply(deg_matrix_direction, 2, count_deg)
#
# barplot(tp[1,], ylim=c(-200, 200), col="red", ylab="Number of Differentially Expressed eRNAs",
#         names.arg = substr(colnames(tp), 3, 1000L))
# barplot(tp[2,], add=T, col="blue", names.arg = F)

###############################################################
# 3. Summarize eRNA expression profiles into groups
###############################################################

x <- eRNA.rpkm[, 1:12]
lab <- rowSums(deg_matrix[, 1:9]) > 1
lab <- rownames(deg_matrix)[lab]
lab <- lab[lab %in% rownames(x)]
x <- x[lab, ]
x <- x[, 2:12] - x[, 1]
x <- (x - rowMeans(x)) / apply(x, 1, sd)

# Sort rows by kmeans clustering
zscore <- x
set.seed(0)
centers <- 3
kmeans <- kmeans(zscore, centers, iter.max = 10^5)

k_cluster <- kmeans$cluster
cluster_order <- c(2, 1, 3)
ind <- cluster_order[k_cluster] * 10^4 + rowMeans(zscore)
od <- order(ind, decreasing = F)
zscore <- zscore[od, ]
k_cluster <- k_cluster[od]

# eRNA heatmap
my_palette <- colorRampPalette(c("blue", "black", "yellow"), bias = 1)(n = 500)
heatmap.2(zscore[nrow(zscore):1, ],
  col = my_palette, Colv = F, Rowv = F, dendrogram = "none",
  labRow = NA, keysize = 1.5, density.info = "density",
  key.ylab = "", key.xlab = "", key.title = "", trace = "none"
)
###############################################################
# 4. Summarize gene expression profiles near differentially expressed eRNAs
###############################################################
par(mfcol = c(2, 2), mar = c(4, 4, 4, 5))
for (i in cluster_order) {
  boxplot(zscore[k_cluster == i, ], ylab = "Zscore")
}

flanking <- 500000
for (i in cluster_order) {
  lab <- rownames(zscore)[k_cluster == i]
  tp <- eRNA[mcols(eRNA)$id %in% lab]
  start(tp) <- start(tp) - flanking
  end(tp) <- end(tp) + flanking
  hits <- findOverlaps(tp, tss)
  lab <- as.character(unique(mcols(tss[subjectHits(hits)])$id))
  lab <- ref_2_symbol(lab)
  lab <- lab[lab %in% de_genes]
  if (i == 2) {
    write.table(lab, quote = F, row.names = F, col.names = F, file = "test.txt")
  }
  tp <- gene.rpkm[lab, ]
  tp <- (tp - rowMeans(tp)) / apply(tp, 1, sd)
  print(length(lab))
  # boxplot(tp, ylab="Zscore")
  heatmap.2(tp,
    col = my_palette, Colv = F, Rowv = T, dendrogram = "none",
    keysize = 1.5, density.info = "density",
    key.ylab = "", key.xlab = "", key.title = "", trace = "none"
  )
}

# Some individual cases
par(mfrow = c(1, 2))
my_palette <- colorRampPalette(c("blue", "red", "orange"), bias = 1)(n = 9)
plot(eRNA.rpkm["MetaEnhancer_65_-", 2:10], gene.rpkm["IRF2", 2:10],
  col = my_palette, pch = 16, cex = 1.3,
  xlab = "MetaEnhancer_65_-", ylab = "IRF2"
)
plot(eRNA.rpkm["MetaEnhancer_250_-", 2:10], gene.rpkm["IRF1", 2:10],
  col = my_palette, pch = 16, cex = 1.3,
  xlab = "MetaEnhancer_250_-", ylab = "IRF1"
)

plot(eRNA.rpkm["MetaEnhancer_116_+", 2:10], gene.rpkm["TNFSF10", 2:10],
  col = my_palette, pch = 16, cex = 1.3,
  xlab = "MetaEnhancer_116_+", ylab = "TNFSF10"
)

###############################################################
# 5. Correlation test
###############################################################

###############################################################
# 5.1 A permutation test to see whether nearby pairs tend to have higher correlations
###############################################################

pdf(file = "Correlated_pairs_profile.pdf")
cor_pairs_e <- NULL
cor_pairs_m <- NULL
# par(mfrow=c(3,2))
my_palette <- colorRampPalette(c("blue", "red", "orange"), bias = 1)(n = 9)
for (i in cluster_order[3]) {
  eRNA_list <- rownames(zscore)[k_cluster == i]
  for (j in 1:length(eRNA_list)) {
    if ((mean(eRNA.rpkm[eRNA_list[j], 2:10]) > 10) | (max(eRNA.rpkm[eRNA_list[j], 2:10]) < 1.5)) {
      next
    }
    tp <- eRNA[mcols(eRNA)$id %in% eRNA_list[j]]
    lab <- find_neighbor(center = tp, target = tss, flanking = flanking)
    lab <- ref_2_symbol(lab)
    lab <- lab[lab %in% de_genes]
    lab <- lab[!(substr(lab, 1, 4) %in% "HIST")]
    if (length(lab) > 0) {
      for (k in 1:length(lab)) {
        pvalue <- cor.test(eRNA.rpkm[eRNA_list[j], 2:10], gene.rpkm[lab[k], 2:10])$p.value
        if (pvalue < 0.01) {
          cor_pairs_e <- c(cor_pairs_e, eRNA_list[j])
          cor_pairs_m <- c(cor_pairs_m, lab[k])
          plot_eRNA_gene(eRNA_list[j], lab[k])
          #           plot(eRNA.rpkm[eRNA_list[j], 2:10], gene.rpkm[lab[k], 2:10], col=my_palette, pch=16, cex=1.3,
          #                main=paste(apply(as.data.frame(tp)[1,1:3], 1, as.character), collapse=c(" ")),
          #                xlab=eRNA_list[j], ylab=lab[k])
          legend("bottomright", paste("P <", signif(pvalue, 1)), bty = "n")
        }
      }
    }
  }
}
dev.off()

########################################################################
# Divergence: Output relavent promoter regions for primer design
########################################################################
# 1. eRNAs
tp <- c(cor_pairs_e, cor_pairs_m)
output <- data.frame(chr = tp, start = tp, end = tp, name = tp, stringsAsFactors = F)
extension <- 100000
for (i in 1:length(cor_pairs_e)) {
  tp <- promoters(eRNA[as.character(mcols(eRNA)$id) == cor_pairs_e[i]], upstream = extension, downstream = extension)
  output[i, 1:3] <- apply(as.data.frame(tp)[, 1:3], 1, as.character)
}
# 2. mRNAs
for (i in 1:length(cor_pairs_m)) {
  lab <- sym_2_ref(cor_pairs_m[i])
  lab <- lab[which(lab %in% as.character(mcols(gene)$id))[1]]
  lab <- which(as.character(mcols(gene)$id) %in% lab)[1]
  stopifnot(length(lab) == 1)
  tp <- promoters(gene[lab], upstream = extension, downstream = extension)
  output[length(cor_pairs_e) + i, 1:3] <- apply(as.data.frame(tp)[, 1:3], 1, as.character)
}
write.table(output, file = "correlated_pairs.bed", row.names = F, col.names = F, quote = F)
write.table(cbind(cor_pairs_e, cor_pairs_m), file = "correlated_pairs.txt", row.names = F, col.names = F, quote = F)

par(mfrow = c(3, 2))
flanking <- 200000
my_palette <- colorRampPalette(c("blue", "red", "orange"), bias = 1)(n = 9)
for (i in 1:1) {
  eRNA_list <- rownames(zscore)[k_cluster == i]
  for (j in 1:length(eRNA_list)) {
    tp <- eRNA[mcols(eRNA)$id %in% eRNA_list[j]]
    start(tp) <- start(tp) - flanking
    end(tp) <- end(tp) + flanking
    hits <- findOverlaps(tp, tss, ignore.strand = T)
    lab <- as.character(unique(mcols(tss[subjectHits(hits)])$id))
    lab <- ref_2_symbol(lab)
    lab <- lab[lab %in% de_genes]
    if (length(lab) > 0) {
      for (k in 1:length(lab)) {
        pvalue <- cor.test(eRNA.rpkm[eRNA_list[j], 2:10], gene.rpkm[lab[k], 2:10])$p.value
        if (pvalue < 0.01) {
          a <- eRNA.rpkm[eRNA_list[j], 2:10]
          a <- a / max(a) # Normalization
          b <- gene.rpkm[lab[k], 2:10]
          b <- b / max(b)
          names.arg <- substr(names(a), 3, 1000L)
          barplot(rbind(a, b),
            beside = T, main = paste(apply(as.data.frame(tp)[1, 1:3], 1, as.character), collapse = c(" ")),
            ylim = c(0, 1.4), names.arg = names.arg, col = c("grey10", "grey90")
          )
          legend("topleft", c(sub("MetaEnhancer", "eRNA", eRNA_list[j]), lab[k]),
            pch = 15, col = c("grey10", "grey90"), bty = "n", cex = 1.2
          )
        }
      }
    }
  }
}

###############################################################
# 6. Gene adaptation and their eRNAs
###############################################################

init <- gene.rpkm[de_mRNA, 2]
peak <- rowMax(gene.rpkm[de_mRNA, 3:8])
back <- gene.rpkm[de_mRNA, 10]
adaptation_score <- log2(peak^2 / init / back)
lab <- de_mRNA[order(adaptation_score, decreasing = T)]
ct <- 0
pdf(file = "Asynchromized_pairs_profile.pdf", height = 7, width = 8)
for (i in 1:length(lab)) {
  stopifnot(adaptation_score[lab[i]] > 2)
  ref <- sym_2_ref(lab[i])
  tp <- tss[mcols(tss)$id %in% ref]
  if (length(tp) == 0) {
    next
  }
  tp <- tp[1]
  # neighbor_eRNA <- find_neighbor(tp, eRNA[mcols(eRNA)$id %in% de_eRNA])
  neighbor_eRNA <- find_neighbor(tp, eRNA)
  if (all(is.na(neighbor_eRNA)) == F) {
    if (max(eRNA.rpkm[neighbor_eRNA, ]) < 1) {
      next
    }
    ct <- ct + 1
    for (j in 1:length(neighbor_eRNA)) {
      if (max(eRNA.rpkm[neighbor_eRNA[j], ]) > 1) {
        plot_eRNA_gene(neighbor_eRNA[j], lab[i])
      }
    }
  }
}
dev.off()

cand_mRNA <- c(
  "ACKR1", "TNFSF10", "CD38", "MNDA", "PARP9",
  "PARP14", "SAMD9", "IFITM1", "KIAA0040", "IFI35",
  "IFI16", "IFITM3", "SOCS1", "ZBP1"
)
cand_eRNA <- paste("MetaEnhancer",
  c(
    "1136_+", "2005_+", "811_-", "3795_???", "2972_-",
    "2972_???", "3039_+", "1268_+", "1115_+", "820_???",
    "1136_+", "1268_+", "810_???", "1228_+"
  ),
  sep = "_"
)

pdf(file = "Asynchromized_pairs_profile_candidate.pdf", height = 7, width = 8)
for (i in 1:length(cand_mRNA)) {
  neighbor_eRNA <- cand_eRNA[i]
  plot_eRNA_gene(neighbor_eRNA, cand_mRNA[i])
}
dev.off()

########################################################################
# Divergence: Output relavent promoter regions for primer design
########################################################################
# 1. eRNAs
tp <- c(cand_eRNA, cand_mRNA)
output <- data.frame(chr = tp, start = tp, end = tp, name = tp, stringsAsFactors = F)
extension <- 100000
for (i in 1:length(cand_eRNA)) {
  tp <- promoters(eRNA[as.character(mcols(eRNA)$id) == cand_eRNA[i]], upstream = extension, downstream = extension)
  output[i, 1:3] <- apply(as.data.frame(tp)[, 1:3], 1, as.character)
}
# 2. mRNAs
for (i in 1:length(cand_mRNA)) {
  lab <- sym_2_ref(cand_mRNA[i])
  lab <- lab[which(lab %in% as.character(mcols(gene)$id))[1]]
  lab <- which(as.character(mcols(gene)$id) %in% lab)[1]
  stopifnot(length(lab) == 1)
  tp <- promoters(gene[lab], upstream = extension, downstream = extension)
  output[length(cand_eRNA) + i, 1:3] <- apply(as.data.frame(tp)[, 1:3], 1, as.character)
}
write.table(output, file = "asynchronized_pairs.bed", row.names = F, col.names = F, quote = F)
write.table(cbind(cand_eRNA, cand_mRNA), file = "asynchronized_pairs.txt", row.names = F, col.names = F, quote = F)
