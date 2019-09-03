setwd("~/Documents/IFN_enhancer/R/")
library(RColorBrewer)
library(GenomicRanges)
library(GenomicFeatures)
library(gplots)
library(rtracklayer)
library(cluster)
library(vioplot)
source("Script/ref_2_symbol_matrix.R")
source("Script/sym_2_ref.R")
load("Data/gene_expression_GM_unfiltered.RData")
colnames(gene.rpkm) <- substr(colnames(gene.rpkm), 3, 1000L)

load("~/Documents/IFN_enhancer/R/Data/ref2sym_dedup.RData")
ref_2_symbol <- function(x) {
  stopifnot(all(x %in% ref2sym[, 1]))
  tp <- 1:nrow(ref2sym)
  names(tp) <- ref2sym[, 1]
  return(ref2sym[tp[x], 2])
}

fold_change <- function(x, y, pseudo = 0.1) {
  log_fold <- log2((x + pseudo) / (y + pseudo))
  return(log_fold)
}
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
continuity_index <- function(x) {
  return(diag(cor(t(x[, -ncol(x)]), t(x[, -1]), method = "spearman")))
}
find_pair <- function(a, b, a_mat, b_mat, distance) {
  # Make sure that records of a/b matches with rows of a/b_mat
  id_a <- NULL
  id_b <- NULL
  id_a_all <- rownames(a_mat)
  id_b_all <- rownames(b_mat)
  for (i in 1:length(a)) {
    inner <- flank(a[i], width = distance[1], both = T)
    outer <- flank(a[i], width = distance[2], both = T) # Region flanking TSS
    tp <- setdiff(outer, inner)
    hits <- unique(subjectHits(findOverlaps(tp, b, ignore.strand = T)))
    if (length(hits) == 0) {
      next
    }
    lab <- order(cor(t(b_mat)[, hits], a_mat[i, ]), decreasing = T)[1] # Incase of multiple enhancers, select the one that is the most similar to the gene
    hits <- hits[lab]
    id_a <- c(id_a, id_a_all[i])
    id_b <- c(id_b, id_b_all[hits])
  }
  return(data.frame(id_a, id_b, stringsAsFactors = F))
}
GRange_index <- function(region, id) {
  # call records from GRange objects by the order of ids provided
  tp <- c(1:length(region))
  names(tp) <- as.character(mcols(region)$id)
  return(region[tp[id]])
}
PE_plot <- function(gid, eid, distance, col, black_list = NULL) {
  par(mfrow = c(2, 1), mar = c(3, 4, 1, 1))
  rg <- GRange_index(gene, gid)
  mg <- gene.rpkm[gid, -1]
  re <- GRange_index(eRNA, eid)
  me <- eRNA.rpkm[eid, -1]
  PE <- find_pair(a = rg, b = re, a_mat = mg, b_mat = me, distance = distance)

  # Plot of gene
  lab <- unique(PE[, 1])
  lab <- lab[!(lab %in% black_list)]
  tp <- fold_change(mg[lab, -1], mg[lab, 1])
  plot(-1, -1, xlim = c(0, 11), ylim = c(-3, 4), xaxt = "n", ylab = "Gene Log2 Fold-change", xlab = "")
  axis(side = 1, at = c(1:10), labels = colnames(tp))
  for (j in c(1:10)) {
    vioplot(tp[, j], at = j, add = T, col = col[1])
  }
  text(1, 2.8, paste("n=", length(lab), sep = ""))
  abline(h = 0, lty = "dashed", col = "grey", lwd = 1.5)

  # Plot of eRNA
  lab <- unique(PE[, 2])
  lab <- lab[!(lab %in% black_list)]
  tp <- fold_change(me[lab, -1], me[lab, 1])
  plot(-1, -1, xlim = c(0, 11), ylim = c(-3, 4), xaxt = "n", ylab = "eRNA Log2 Fold-change")
  axis(side = 1, at = c(1:10), labels = colnames(tp))
  for (j in c(1:10)) {
    vioplot(tp[, j], at = j, add = T, col = col[2])
  }
  text(1, 2.8, paste("n=", length(lab), sep = ""))
  abline(h = 0, lty = "dashed", col = "grey", lwd = 1.5)
  return(lab)
}
pair_region <- function(a, b) {
  a <- GRange_index(gene, a)
  b <- GRange_index(eRNA, b)
  tp <- c(start(a), end(a), start(b), end(b))
  tp <- paste(seqnames(a), min(tp), max(tp))
  return(tp)
}
################################################################
# 1. Load data
################################################################

# 1.1 Inducible genes
load(file = "Data/gene_expression_GM_unfiltered.RData")
rpkm <- gene.rpkm[, -1]
rpkm <- rpkm[apply(rpkm, 1, max) > 1, ] # 8484 genes remained
ai_rpkm <- amplitude_index(rpkm[, 1:9], abs = F) # Only consider early time points
ci_rpkm <- continuity_index(rpkm)
index_rpkm <- ai_rpkm * ci_rpkm

deg <- names(index_rpkm)[(ci_rpkm > 0.2 & ai_rpkm > log2(2))]
deg.data <- cbind(rpkm, ai_rpkm, ci_rpkm, index_rpkm)[deg, ]
tp <- c(1:length(gene))
names(tp) <- as.character(gene$id)
deg.region <- gene[tp[deg]] # Select corresponding regions

ind_gene <- deg
rnd_gene <- rownames(gene.rpkm)[apply(gene.rpkm, 1, max) > 1]
rnd_gene <- setdiff(rnd_gene, deg)
rnd_gene <- sample(rnd_gene, 500, replace = F)

write.table(ind_gene, "Data/ISG.list", quote = F, row.names = F, col.names = F)

# 1.2 Expressed enhancers
# Criteria:
# 1. max(RPKM) > 1
# 2. Amplitude index > 1
# 3. Continuity index > 0.2
load(file = "Data/eRNA_expression_GM_unfiltered.RData")
colnames(eRNA.rpkm) <- substr(colnames(eRNA.rpkm), 3, 1000L)
eRNA_rpkm <- eRNA.rpkm[, -1]
eRNA_rpkm <- eRNA_rpkm[apply(eRNA_rpkm, 1, max) > 0.5, ]
ai_eRNA <- amplitude_index(eRNA_rpkm)
ci_eRNA <- continuity_index(eRNA_rpkm)
ind_eRNA <- 1 * ci_eRNA

# Criteria:
enh <- names(ind_eRNA)[ci_eRNA > 0.2]
enh.data <- cbind(eRNA_rpkm, ai_eRNA, ci_eRNA, ind_eRNA)[enh, ]
tp <- c(1:length(eRNA))
names(tp) <- as.character(eRNA$id)
enh.region <- eRNA[tp[enh]]


################################################################
# 3. Optimal EP distance
################################################################

# Make sure that regions match with expression table
ind_r <- GRange_index(gene, ind_gene)
ind_d <- deg.data[ind_gene, c(1:11)]

# Test which is the best distance to assign EP pairs

# Comparison 1
pdf("Figure/Enhancers_induced_near_genes_large_scale.pdf")
par(mfrow = c(2, 1), mar = c(3, 4, 1, 1))
step <- 1 * 10^5
for (i in 1:5) {
  distance <- c(0, step) + (i - 1) * step
  PE <- find_pair(a = ind_r, b = enh.region, a_mat = ind_d, b_mat = enh.data[, 1:11], distance = distance)
  lab <- unique(PE[, 2])
  tp <- log2((eRNA.rpkm[lab, -c(1, 2)] + 0.1) / (eRNA.rpkm[lab, 2] + 0.1))
  plot(-1, -1,
    xlim = c(0, 11), ylim = c(-3, 3), xaxt = "n", ylab = "eRNA Log2 Fold-change",
    main = paste(paste(distance / 1000, collapse = "-"), "kb")
  )
  axis(side = 1, at = c(1:10), labels = colnames(tp))
  for (j in c(1:10)) {
    vioplot(tp[, j], at = j, add = T, col = rainbow(10)[j])
  }
  text(1, 2.8, paste("n=", length(lab), sep = ""))
  abline(h = 0, lty = "dashed", col = "grey", lwd = 1.5)
}
dev.off()

pdf("Figure/Enhancers_induced_near_genes_small_scale.pdf")
par(mfrow = c(2, 1), mar = c(3, 4, 1, 1))
step <- 4 * 10^4
for (i in 1:5) {
  distance <- c(0, step) + (i - 1) * step
  PE <- find_pair(a = ind_r, b = enh.region, a_mat = ind_d, b_mat = enh.data[, 1:11], distance = distance)
  lab <- unique(PE[, 2])
  tp <- log2((eRNA.rpkm[lab, -c(1, 2)] + 0.1) / (eRNA.rpkm[lab, 2] + 0.1))
  plot(-1, -1,
    xlim = c(0, 11), ylim = c(-3, 3), xaxt = "n", ylab = "eRNA Log2 Fold-change",
    main = paste(paste(distance / 1000, collapse = "-"), "kb")
  )
  axis(side = 1, at = c(1:10), labels = colnames(tp))
  for (j in c(1:10)) {
    vioplot(tp[, j], at = j, add = T, col = rainbow(10)[j])
  }
  text(1, 2.8, paste("n=", length(lab), sep = ""))
  abline(h = 0, lty = "dashed", col = "grey", lwd = 1.5)
}
dev.off()



##########################################################################################
# Comparison 2
##########################################################################################

# Violin plot

c_start <- 0
binsize <- 10^5
step <- 10^5
n <- 9

enh_list <- list()
for (i in 1:n) {
  distance <- c(c_start, c_start + binsize)
  PE <- find_pair(
    a = ind_r, b = promoters(enh.region, upstream = 100, downstream = 100),
    a_mat = ind_d, b_mat = enh.data[, 1:11], distance = distance
  )
  lab <- unique(PE[, 2])
  enh_list[[i]] <- lab
  c_start <- c_start + step
}
c_start <- 0
x_lab <- NULL
for (i in 1:n) {
  x_lab <- c(x_lab, paste(c_start, c_start + binsize / 1000, sep = "~"))
  c_start <- c_start + step / 1000
}

pdf(paste("Figure/Distance_effct_P2E.pdf", sep = ""))
par(mfrow = c(2, 1), mar = c(3, 4, 1, 1))
c_start <- 0
for (j in 1:10) {
  timepoint <- colnames(eRNA.rpkm)[-c(1, 2)][j]
  plot(-1, -1, xlim = c(0, n + 1), ylim = c(-4, 4), xaxt = "n", ylab = "eRNA Log2 Fold-change", main = "", xlab = "EP Distance (kb)")
  text(seq(1, n, by = 1), par("usr")[3] - 0.2, labels = x_lab, srt = 45, pos = 1, xpd = TRUE, cex = 0.7)
  for (i in 1:n) {
    lab <- unlist(enh_list[[i]])
    tp <- log2((eRNA.rpkm[lab, timepoint] + 0.1) / (eRNA.rpkm[lab, 2] + 0.1))
    vioplot(tp, at = i, add = T, col = rainbow(n)[i])
    text(i, -3.5, paste("n=", length(lab), sep = ""), cex = 0.7)
  }
  abline(h = 0, lty = "dashed", col = "grey", lwd = 1.5)
  text(0.5, 4, timepoint, cex = 1)
}
dev.off()

##########################################################################################
# Curve

c_start <- 0
binsize <- 10^5
step <- 10^5
n <- 9

enh_list <- list()
for (i in 1:n) {
  distance <- c(c_start, c_start + binsize)
  PE <- find_pair(
    a = ind_r, b = promoters(enh.region, upstream = 100, downstream = 100),
    a_mat = ind_d, b_mat = enh.data[, 1:11], distance = distance
  )
  lab <- unique(PE[, 2])
  enh_list[[i]] <- lab
  c_start <- c_start + step
}

c_start <- 0
x_lab <- NULL
for (i in 1:n) {
  x_lab <- c(x_lab, paste(c_start, c_start + binsize / 1000, sep = "~"))
  c_start <- c_start + step / 1000
}

fc <- NULL
nn <- NULL
c_start <- 0
for (i in 1:n) {
  lab <- unlist(enh_list[[i]])
  tp <- log2((eRNA.rpkm[lab, -1] + 0.1) / (eRNA.rpkm[lab, 1] + 0.1))
  tp <- colMeans(tp)
  if (is.null(fc)) {
    fc <- tp
  }
  else {
    fc <- rbind(fc, tp)
  }
  nn <- c(nn, length(lab))
}
rownames(fc) <- x_lab

pdf("Figure/Enhanchers_near_genes_curve.pdf")
par(mfrow = c(1, 1), mar = c(3, 4, 1, 1))
nx <- 8
lw <- 1.5
plot(-1, -1, xlim = c(0.5, nx + 1.5), ylim = c(-0.25, 0.8), xaxt = "n", ylab = "eRNA Log2 Fold-change", main = "", xlab = "EP Distance (kb)")
text(seq(1, nx, by = 1), par("usr")[3] - 0.01, labels = colnames(eRNA.rpkm)[2:(nx + 1)], pos = 1, xpd = TRUE, cex = 1)
for (k in 1:2) {
  tp <- fc[k, 1:nx]
  lines(tp, col = rainbow(n)[k], lw = lw)
  points(tp, col = rainbow(n)[k])
  text(nx + 0.8, tp[nx], x_lab[k], cex = 0.7)
}
plot(-1, -1, xlim = c(0.5, nx + 1.5), ylim = c(-0.25, 0.8), xaxt = "n", ylab = "eRNA Log2 Fold-change", main = "", xlab = "EP Distance (kb)")
text(seq(1, nx, by = 1), par("usr")[3] - 0.01, labels = colnames(eRNA.rpkm)[2:(nx + 1)], pos = 1, xpd = TRUE, cex = 1)
for (k in 3:4) {
  tp <- fc[k, 1:nx]
  lines(tp, col = rainbow(n)[k], lw = lw)
  points(tp, col = rainbow(n)[k])
  text(nx + 0.8, tp[nx], x_lab[k], cex = 0.7)
}
plot(-1, -1, xlim = c(0.5, nx + 1.5), ylim = c(-0.4, 0.8), xaxt = "n", ylab = "eRNA Log2 Fold-change", main = "", xlab = "EP Distance (kb)")
text(seq(1, nx, by = 1), par("usr")[3] - 0.01, labels = colnames(eRNA.rpkm)[2:(nx + 1)], pos = 1, xpd = TRUE, cex = 1)
for (k in 5:6) {
  tp <- fc[k, 1:nx]
  lines(tp, col = rainbow(n)[k], lw = lw)
  points(tp, col = rainbow(n)[k])
  text(nx + 0.8, tp[nx], x_lab[k], cex = 0.7)
}
plot(-1, -1, xlim = c(0.5, nx + 1.5), ylim = c(-0.4, 0.8), xaxt = "n", ylab = "eRNA Log2 Fold-change", main = "", xlab = "EP Distance (kb)")
text(seq(1, nx, by = 1), par("usr")[3] - 0.01, labels = colnames(eRNA.rpkm)[2:(nx + 1)], pos = 1, xpd = TRUE, cex = 1)
for (k in 7:8) {
  tp <- fc[k, 1:nx]
  lines(tp, col = rainbow(n)[k], lw = lw)
  points(tp, col = rainbow(n)[k])
  text(nx + 0.8, tp[nx], x_lab[k], cex = 0.7)
}
dev.off()
##########################################################################################

################################################################
# 4. Associate enhancers to genes
################################################################

# Now we decided to use 100kb as the cutoff
# Show that induced genes had induced neighbor genes
paired_color <- brewer.pal(12, "Paired")
pdf("Figure/Enhanchers_near_genes_all_types.pdf")
ind_enh <- PE_plot(ind_gene, rownames(enh.data), c(0, 2 * 10^5), col = paired_color[c(6, 5)])
rnd_enh <- PE_plot(rnd_gene, rownames(enh.data), c(0, 2 * 10^5), col = paired_color[c(10, 9)], black_list = ind_enh)
dev.off()

pdf("Figure/Enhanchers_near_genes_all_types_100kb.pdf")
ind_enh <- PE_plot(ind_gene, rownames(enh.data), c(0, 10^5), col = paired_color[c(6, 5)])
rnd_enh <- PE_plot(rnd_gene, rownames(enh.data), c(0, 10^5), col = paired_color[c(10, 9)], black_list = ind_enh)
dev.off()














###################################################################################
# Optional
###################################################################################
# par(mfrow=c(2,2))
# deg.up1 <- unique(rownames(heat.deg)[k_cluster==3])
# deg.up2 <- unique(rownames(heat.deg)[k_cluster==2]); deg.up2 <- setdiff(deg.up2, deg.up1)
# deg.dn1 <- unique(rownames(heat.deg)[k_cluster==1]); deg.dn1 <- setdiff(deg.dn1, deg.up1)
# deg.rnd <- setdiff(rownames(deg.data), c(deg.up1, deg.up2, deg.dn1))
# deg.rnd <- sample(deg.rnd, 200, replace = F)
# deg.list <- list(deg.up1, deg.up2, deg.dn1, deg.rnd)
# col <- c("red", "gold", "blue", "grey")
# for(i in 1:length(deg.list)){
#   tp <- unlist(deg.list[i])
#   tp <- log2((deg.data[tp, -1]+0.1) / (deg.data[tp, 1]+0.1))
#   plot(-1,-1, xlim=c(0, 11), ylim=quantile(tp, probs=c(0.05,0.99)),
#        xlab="Time Point", ylab="Log2 Fold Change", xaxt="n")
#   for(j in c(1:10)){vioplot(tp[,j], at=j, add=T, col=col[i])}
#   abline(h=0, lty="dashed")
# }
#
# par(mfrow=c(2,2))
# enh.up1 <- unique(rownames(heat.enh)[k_cluster==3])
# enh.up2 <- unique(rownames(heat.enh)[k_cluster==2]); enh.up2 <- setdiff(enh.up2, enh.up1)
# enh.dn1 <- unique(rownames(heat.enh)[k_cluster==1]); enh.dn1 <- setdiff(enh.dn1, enh.up1)
# enh.rnd <- setdiff(rownames(enh.data), c(enh.up1, enh.up2, enh.dn1))
# enh.rnd <- sample(enh.rnd, 200, replace = F)
# enh.list <- list(enh.up1, enh.up2, enh.dn1, enh.rnd)
# col <- c("red", "gold", "blue", "grey")
# for(i in 1:length(enh.list)){
#   tp <- unlist(enh.list[i])
#   tp <- log2((enh.data[tp, -1]+0.1) / (enh.data[tp, 1]+0.1))
#   plot(-1,-1, xlim=c(0, 11), ylim=quantile(tp, probs=c(0.01,0.99)),
#        xlab="Time Point", ylab="Log2 Fold Change", xaxt="n")
#   for(j in c(1:10)){vioplot(tp[,j], at=j, add=T, col=col[i])}
#   abline(h=0, lty="dashed")
# }
#
# par(mfrow=c(3,2))
# for(i in 1:10){
#   vioplot(heat.enh[k_cluster==3, i], heat.enh[k_cluster==2, i], heat.enh[k_cluster==1, i])
#   title(colnames(heat.enh)[i])
# }
#
# for(i in 1:10){
#   vioplot(heat.deg[k_cluster==3, i], heat.deg[k_cluster==2, i], heat.deg[k_cluster==1, i])
#   title(colnames(heat.enh)[i])
# }
