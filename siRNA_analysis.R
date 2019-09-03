# Analyze siRNA data

setwd("~/Documents/IFN_enhancer/R/")
source("Script/ref_2_symbol_matrix.R")
library(GenomicRanges)
library(GenomicFeatures)
library(gplots)
library(rtracklayer)

GRange_index <- function(region, id) {
  # call records from GRange objects by the order of ids provided
  tp <- c(1:length(region))
  names(tp) <- as.character(mcols(region)$id)
  return(region[tp[id]])
}

x <- read.table("Data/siRNA_efficiency.txt", header = T)
efold <- x[, 1]
mfold <- x[, 2]
stats <- c(
  sum(efold < 1),
  sum(efold < 0.7),
  sum(efold < 0.5),
  sum((efold < 0.7) & (mfold < 0.7)),
  sum((efold < 0.5) & (mfold < 0.5)) # Number of successful KD experiments
)
tp <- round(sum((efold < 0.7) & (mfold < 1)) / sum(efold < 0.7) * 100, 1)
print(paste(tp, "% mRNA repressed when eRNA KD by >30%", sep = ""))
tp <- round(sum((efold < 0.5) & (mfold < 1)) / sum(efold < 0.5) * 100, 1)
print(paste(tp, "% mRNA repressed when eRNA KD by >50%", sep = ""))


pairs_id <- rownames(x) # Each EP pair may correspond to several siRNA designs
lab <- grep("[A,B,C,D,E]", pairs_id)
pairs_id[lab] <- substr(pairs_id[lab], 1, nchar(pairs_id[lab]) - 1)
uid <- unique(pairs_id)

ids <- read.table("Data/siRNA_EP_ids", colClasses = "character")
lab <- grep("-[1,2,3,4,5]", ids[, 1])
ids[lab, 1] <- substr(ids[lab, 1], 1, nchar(ids[lab, 1]) - 2)
rownames(ids) <- uid
ids_each_row <- ids[pairs_id, ]
names(efold) <- ids_each_row[, 2]
names(mfold) <- ids_each_row[, 1]

load("Data/gene_expression_GM_unfiltered.RData")
gene.rpkm <- ref_2_symbol(gene.rpkm)
gdata <- gene.rpkm[ids_each_row[, 1], ]
load("Data/eRNA_expression_GM_unfiltered.RData")
edata <- eRNA.rpkm[ids_each_row[, 2], ]

pairs_stats <- NULL
for (i in uid) {
  ii <- which(pairs_id == i)
  ee <- efold[ii]
  mm <- mfold[ii]
  pairs_stats <- c(
    pairs_stats,
    sum(ee < 1),
    sum(ee < 0.7),
    sum(ee < 0.5),
    sum((ee < 0.7) & (mm < 0.7)),
    sum((ee < 0.5) & (mm < 0.5))
  )
}
pairs_stats <- matrix(as.integer(pairs_stats > 0), ncol = 5, byrow = T)
rownames(pairs_stats) <- uid

pdf(file = "Figure/siRNA_summary.pdf")
par(mfrow = c(1, 1))
barplot(rbind(stats, colSums(pairs_stats)),
  beside = T, ylim = c(0, 50),
  col = c("navy", "darkorange"),
  ylab = "Number of EP pairs",
  names.arg = c(
    "eFold<1",
    "eFold<0.7",
    "eFold<0.5",
    "both<0.7",
    "both<0.5"
  )
)
legend("topright", c("# of siRNA", "# of EP pairs"), col = c("navy", "darkorange"), pch = 15)
dev.off()


#############################################################################################
pdf("Figure/mfold_boxplot.pdf")
par(mar = c(4, 5, 5, 4))
boxplot(list(log2(mfold[efold < 0.7] + 0.001), log2(mfold[efold > 0.7] + 0.001)),
  notch = F, at = c(1, 2), xlim = c(0.5, 2.5), ylim = c(-3, 3), boxwex = 0.4,
  names = c("eRNA<70%", "eRNA>70%"), col = c("#2c7fb8", "#d95f0e"),
  ylab = "Log2 mRNA FC", outline = F
)
abline(h = 0, lty = "dashed", col = "grey")
dev.off()

# heat.data <- log2(cbind(efold, mfold)+0.001)
# heat.data <- heat.data[order(efold),]
# my_palette <- colorRampPalette(c("blue", "white", "yellow"), bias=1)(n = 100)
# breaks <- c(seq(-3, 3, length.out=101))
# pdf("Figure/mfold_efold_all.pdf")
# heatmap.2(heat.data, col=my_palette, breaks=breaks,
#           Colv=F, Rowv=F, dendrogram="none", labRow=F,
#           keysize=1.2, density.info="density", key.ylab="", key.xlab="", key.title="",
#           trace="none")
# dev.off()
#
# xx <- cbind(x, ids_each_row[,2:1])
# t1 <- xx[xx[,1]<0.7, ]
# t1 <- t1[order(t1[,2]), ]
# t2 <- xx[xx[,1]>0.7, ]
# t2 <- t2[order(t2[,2]), ]
# xx <- rbind(t1, t2)
# heat.data <- as.matrix(log2(xx[,c(1,2)]+0.001))
# breaks <- c(seq(-3, 3, length.out=101))
# pdf("Figure/mfold_efold_all_sort_mfold.pdf")
# heatmap.2(heat.data, col=my_palette, breaks=breaks,
#           Colv=F, Rowv=F, dendrogram="none", labRow=xx[,3],
#           keysize=1.2, density.info="density", key.ylab="", key.xlab="", key.title="",
#           trace="none")
# dev.off()

# We selected inducible enhancers to match with inducible genes (foreground, fg), which are <200kb from each others
# The rest serve as negative controls (background, bg).
xx <- cbind(x, ids_each_row[, 2:1])
selected_pairs <- read.table("Data/PE.ulist", colClasses = "character")
fg <- which((xx[, 4] %in% selected_pairs[, 1]) & (xx[, 3] %in% selected_pairs[, 2]))
bg <- which(!(c(1:nrow(xx)) %in% fg))

# ##########################################################
# # When eRNAs were effectively repressed
# ##########################################################
# u_rows <- function(x){
#   rn <- unique(rownames(x))
#   for(i in 1:length(rn)){
#     tp <- x[rownames(x)==rn[i],]
#     if(!is.null(nrow(tp))){tp <- colMeans(tp)}
#     if(i==1){out <- tp}
#     else{out <- rbind(out, tp)}
#   }
#   rownames(out) <- rn
#   return(out)
# }
# # 1. Foreground
# tm <- mfold[fg]
# te <- efold[fg]
# tm <- tm[te<0.7]
# te <- te[te<0.7]
# heat.data <- log2(cbind(tm, te)+2^(-4))
# rownames(heat.data) <- paste(names(tm), names(te), sep="*")
# heat.data <- u_rows(heat.data)
# lab <- order(heat.data[,1])
# heat.data <- heat.data[lab,]
# fc_fg <- heat.data[,1]
# heat.eRNA <- NULL
# heat.gene <- NULL
# for(i in 1:nrow(heat.data)){
#   tp <- strsplit(rownames(heat.data)[i], "[*]")
#   heat.gene <- c(heat.gene, tp[[1]][1])
#   heat.eRNA <- c(heat.eRNA, tp[[1]][2])
# }
# tp_list <- cbind(heat.gene, heat.eRNA)
# write.table(tp_list, file = "Data/mfold_eRNA_repressed_fg.list",
#             quote=F, row.names = F, col.names = F, sep="\t")
# my_palette <- colorRampPalette(c("blue", "white", "yellow"), bias=1)(n = 100)
# breaks <- c(seq(-3, 3, length.out=101))
# pdf("Figure/mfold_eRNA_repressed_fg.pdf")
# heatmap.2(heat.data, col=my_palette, breaks=breaks,
#           Colv=F, Rowv=F, dendrogram="none",
#           keysize=1.2, density.info="density", key.ylab="", key.xlab="", key.title="",
#           trace="none")
# dev.off()
#
# # Show siRNA and GROseq results case by case
# pdf("Figure/siRNA_and_GROseq_fg.pdf")
# par(mfrow=c(2,2), mai=c(0.4,.8,0.5,0.1))
# color <- c("navy", "darkorange")
# pch <- c(15, 16)
# timepoints <- substr(colnames(gene.rpkm), 3, 1000L)
# x_at <- c(1,3,5,7,9,11)
# for(i in 1:nrow(heat.data)){
#   barplot(c(2^heat.data[i, 2], 2^heat.data[i, 1]), space=0.8, xlim=c(0.2,4),
#           ylim=c(0, max(c(1, 2^heat.data[i,1], 2^heat.data[i,2]))),
#           ylab="siRNA FC", col=color, names.arg = tp_list[i,c(2,1)])
#   ta <- edata[tp_list[i,2],]
#   tb <- gdata[tp_list[i,1],]
#   plot(ta, ylim=c(min(ta, tb), max(ta, tb)), pch=pch[1], xaxt="n",
#        xlab="", ylab="RPKM", main="", col=color[1])
#   axis(1, at=x_at,labels=timepoints[x_at], col.axis="black", las=1)
#   lines(ta, col=color[1])
#   points(tb, pch=pch[2], col=color[2])
#   lines(tb, col=color[2])
#   legend("topleft", c(tp_list[i,2], tp_list[i,1]), bty="n", col=color, lty="solid")
# }
# dev.off()

##########################################################
# When eRNAs were effectively repressed
##########################################################
u_rows <- function(x) {
  rn <- unique(rownames(x))
  for (i in 1:length(rn)) {
    tp <- x[rownames(x) == rn[i], ]
    if (!is.null(nrow(tp))) {
      tp <- c(
        colMeans(tp),
        sd(tp[, 1]), # Error bar data for mRNA
        sd(tp[, 2]) # Error bar data for eRNA
      )
    }
    else {
      tp <- c(tp, -1, -1)
    } # If not replicates, set error to -1
    if (i == 1) {
      out <- tp
    }
    else {
      out <- rbind(out, tp)
    }
  }
  rownames(out) <- rn
  return(out)
}
# 1. Foreground
tm <- mfold[fg]
te <- efold[fg]
tm <- tm[te < 0.7]
te <- te[te < 0.7]
fold.data <- cbind(tm, te)
rownames(fold.data) <- paste(names(tm), names(te), sep = "*")
heat.data <- log2(fold.data + 2^(-4))
heat.data <- u_rows(heat.data)
heat.data <- heat.data[, 1:2]
lab <- order(heat.data[, 1])
heat.data <- heat.data[lab, ]
fc_fg <- heat.data[, 1]
heat.eRNA <- NULL
heat.gene <- NULL
for (i in 1:nrow(heat.data)) {
  tp <- strsplit(rownames(heat.data)[i], "[*]")
  heat.gene <- c(heat.gene, tp[[1]][1])
  heat.eRNA <- c(heat.eRNA, tp[[1]][2])
}
tp_list <- cbind(heat.gene, heat.eRNA)
write.table(tp_list,
  file = "Data/mfold_eRNA_repressed_fg.list",
  quote = F, row.names = F, col.names = F, sep = "\t"
)
my_palette <- colorRampPalette(c("blue", "white", "yellow"), bias = 1)(n = 100)
breaks <- c(seq(-3, 3, length.out = 101))
pdf("Figure/mfold_eRNA_repressed_fg.pdf")
heatmap.2(heat.data[, 1:2],
  col = my_palette, breaks = breaks,
  Colv = F, Rowv = F, dendrogram = "none",
  keysize = 1.2, density.info = "density", key.ylab = "", key.xlab = "", key.title = "",
  trace = "none"
)
dev.off()

# Show siRNA and GROseq results case by case
bar.data <- u_rows(fold.data)
bar.data <- bar.data[rownames(heat.data), ]
pdf("Figure/siRNA_and_GROseq_fg.pdf")
par(mfrow = c(2, 2), mai = c(0.4, 0.8, 0.5, 0.1))
color <- c("darkorange", "navy")
pch <- c(15, 16)
timepoints <- substr(colnames(gene.rpkm), 3, 1000L)
x_at <- c(1, 3, 5, 7, 9, 11)
for (i in 1:nrow(bar.data)) {
  tp <- barplot(
    # c(bar.data[i, 2], bar.data[i, 1]),
    c(2^(heat.data[i, 1]), 2^(heat.data[i, 2])),
    space = 0.8, xlim = c(0.2, 4),
    # ylim=c(0, max(c(1, bar.data[i,1], bar.data[i,2]))),
    ylim = c(0, max(c(1, 2^(heat.data[i, 1]), 2^(heat.data[i, 2])))),
    ylab = "siRNA FC", col = color, names.arg = tp_list[i, c(1, 2)]
  )
  if (bar.data[i, 3] != -1) { # If we have replicates, plot error bar
    # arrows(tp[,1], bar.data[i,c(2,1)]-bar.data[i,c(4,3)]/2,
    #        tp[,1], bar.data[i,c(2,1)]+bar.data[i,c(4,3)]/2,
    #        length=0.15, angle=90, code=3)
    arrows(tp[, 1], 2^heat.data[i, c(1, 2)] - bar.data[i, c(3, 4)] / 2,
      tp[, 1], 2^heat.data[i, c(1, 2)] + bar.data[i, c(3, 4)] / 2,
      length = 0.15, angle = 90, code = 3
    )
  }
  ta <- gdata[tp_list[i, 1], ]
  tb <- edata[tp_list[i, 2], ]
  plot(ta,
    ylim = c(min(ta, tb), max(ta, tb)), pch = pch[1], xaxt = "n",
    xlab = "", ylab = "RPKM", main = "", col = color[1]
  )
  axis(1, at = x_at, labels = timepoints[x_at], col.axis = "black", las = 1)
  lines(ta, col = color[1])
  points(tb, pch = pch[2], col = color[2])
  lines(tb, col = color[2])
  legend("topleft", c(tp_list[i, 1], tp_list[i, 2]), bty = "n", col = color, lty = "solid")
}
dev.off()

# 2. Background
tm <- mfold[bg]
te <- efold[bg]
tm <- tm[te < 0.7]
te <- te[te < 0.7]
heat.data <- log2(cbind(tm, te) + 2^(-4))
rownames(heat.data) <- paste(names(tm), names(te), sep = "*")
heat.data <- u_rows(heat.data)
lab <- order(heat.data[, 1])
heat.data <- heat.data[lab, ]
fc_bg <- heat.data[, 1]
heat.eRNA <- NULL
heat.gene <- NULL
for (i in 1:nrow(heat.data)) {
  tp <- strsplit(rownames(heat.data)[i], "[*]")
  heat.gene <- c(heat.gene, tp[[1]][1])
  heat.eRNA <- c(heat.eRNA, tp[[1]][2])
}
tp_list <- cbind(heat.gene, heat.eRNA)
write.table(tp_list,
  file = "Data/mfold_eRNA_repressed_bg.list",
  quote = F, row.names = F, col.names = F, sep = "\t"
)
my_palette <- colorRampPalette(c("blue", "white", "yellow"), bias = 1)(n = 100)
breaks <- c(seq(-3, 3, length.out = 101))
pdf("Figure/mfold_eRNA_repressed_bg.pdf")
heatmap.2(heat.data,
  col = my_palette, breaks = breaks,
  Colv = F, Rowv = F, dendrogram = "none",
  keysize = 1.2, density.info = "density", key.ylab = "", key.xlab = "", key.title = "",
  trace = "none"
)
dev.off()

# Show siRNA and GROseq results case by case
pdf("Figure/siRNA_and_GROseq_bg.pdf")
par(mfrow = c(2, 2), mai = c(0.4, .8, 0.5, 0.1))
color <- c("navy", "darkorange")
pch <- c(15, 16)
timepoints <- substr(colnames(gene.rpkm), 3, 1000L)
x_at <- c(1, 3, 5, 7, 9, 11)
for (i in 1:nrow(heat.data)) {
  barplot(c(2^heat.data[i, 2], 2^heat.data[i, 1]),
    space = 0.8, xlim = c(0.2, 4),
    ylim = c(0, max(c(1, 2^heat.data[i, 1], 2^heat.data[i, 2]))),
    ylab = "siRNA FC", col = color, names.arg = tp_list[i, c(2, 1)]
  )
  ta <- edata[tp_list[i, 2], ]
  tb <- gdata[tp_list[i, 1], ]
  plot(ta,
    ylim = c(min(ta, tb), max(ta, tb)), pch = pch[1], xaxt = "n",
    xlab = "", ylab = "RPKM", main = "", col = color[1]
  )
  axis(1, at = x_at, labels = timepoints[x_at], col.axis = "black", las = 1)
  lines(ta, col = color[1])
  points(tb, pch = pch[2], col = color[2])
  lines(tb, col = color[2])
  legend("topleft", c(tp_list[i, 2], tp_list[i, 1]), bty = "n", col = color, lty = "solid")
}
dev.off()

##########################################################
# Statistical test
##########################################################

pdf("Figure/siRNA_ks-test.pdf")
boxplot(list(fc_bg, fc_fg),
  notch = F, at = c(1, 2), xlim = c(0.5, 2.5), ylim = c(-3, 3), boxwex = 0.4,
  names = c("Control", "Selected"), col = c("#2c7fb8", "#d95f0e"),
  ylab = "Log2 mRNA FC", outline = F
)
abline(h = 0, lty = "dashed", col = "grey")
pv <- ks.test(fc_fg, fc_bg, alternative = "greater")$p.value
text(x = 2, y = 2.5, paste("P-value:", pv, sep = "\t"))
dev.off()

# heat.data <- log2(mfold[efold<0.7])
# heat.data <- heat.data[order(heat.data)]
# my_palette <- colorRampPalette(c("blue", "white", "yellow"), bias=1)(n = 100)
# breaks <- c(seq(-3, 3, length.out=101))
# pdf("Figure/mfold_eRNA_repressed.pdf")
# heatmap.2(cbind(heat.data, heat.data), col=my_palette, breaks=breaks,
#           Colv=F, Rowv=F, dendrogram="none",
#           keysize=1.2, density.info="density", key.ylab="", key.xlab="", key.title="",
#           trace="none")
# dev.off()
#
# pdf("Figure/mfold_eRNA_activated.pdf")
# heat.data <- log2(mfold[efold>1])
# heat.data <- heat.data[order(heat.data)]
# my_palette <- colorRampPalette(c("blue", "white", "yellow"), bias=1)(n = 100)
# breaks <- c(seq(-3, 3, length.out=101))
# heatmap.2(cbind(heat.data, heat.data), col=my_palette, breaks=breaks,
#           Colv=F, Rowv=F, dendrogram="none",
#           keysize=1.2, density.info="density", key.ylab="", key.xlab="", key.title="",
#           trace="none")
# dev.off()


#############################################################################################
# Distance affect siRNA results?
#############################################################################################

# smb <- NULL
# rfs <- NULL
# for(i in 1:length(gene)){
#   tp <- ref_2_symbol(as.character(mcols(gene[i])$id))
#   if(length(tp)==1){
#     smb <- c(smb, tp)
#     rfs <- c(rfs, as.character(mcols(gene[i])$id))
#   }
# }
# names(smb) <- rfs
#
# ds <- NULL
# lab <- which(efold<0.7)
# for(i in 1:length(lab)){
#   te <- GRange_index(eRNA, names(efold)[lab[i]])
#   tm <- rfs[smb == names(mfold)[lab[i]]]
#   tm <- GRange_index(gene, tm[1])
#   tp <- abs(start(promoters(te, 1, 0)) - start(promoters(tm, 1,0)))
#   ds <- c(ds, tp)
#   print(ds)
# }
# lab <- lab[order(ds)]
# barplot(mfold[lab])
# cor(mfold[lab], ds, method="spearman")
