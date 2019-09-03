# Gro-seq data from Kim Lab
# Analysis of enhancer RNA expression dynamics of GM12878 cells
# Build related data files for potential usage

setwd("~/Documents/IFN_enhancer/")

library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)

###################################################
### Load Gene Annotation
###################################################
# hg18.Refseq <- makeTxDbFromUCSC(genome = "hg18", tablename = "refGene")
# gene <- transcripts(hg18.Refseq); gene <- unique(gene); gene_length <- width(gene)
# save(gene, gene_length, file = "C:/Users/PengXie/Documents/eRNA/link/UCSC/hg18_refseq.RData")
# table <- data.frame(chr=seqnames(gene),
#                     start=start(gene),
#                     end=end(gene),
#                     id=mcols(gene)[,"tx_name"],
#                     depth=rep(0, length(gene)),
#                     strand=as.vector(strand(gene))
# )
# table <- table[grep("NM", table[,"id"]),]
# table <- table[-grep("_", table[,"chr"]),]
# write.table(table, file="genes.bed", col.names=F, row.names=F, quote=F, sep="\t")
load(file = "~/Documents/eRNA/link/UCSC/hg18_refseq.RData")

###################################################
### Load H3K4me1 peaks for GM
###################################################
x <- read.table("~/Documents/eRNA/link/UCSC/wgEncodeBroadChipSeqPeaksGm12878H3k4me1.broadPeak")
hist(x[, 8], breaks = 100)
k4me1 <- GRanges(x[, 1], IRanges(x[, 2], x[, 3]), strand = "*")
hist(log10(width(k4me1)))
k4me1 <- k4me1[width(k4me1) < 10^4]
save(k4me1, file = "~/Documents/eRNA/link/UCSC/H3K4me1_GM.RData")

###################################################
### Load H3K27ac peaks for GM
###################################################
x <- read.table("~/Documents/eRNA/link/UCSC/wgEncodeBroadChipSeqPeaksGm12878H3k27ac.broadPeak")
hist(x[, 8], breaks = 100)
k27ac <- GRanges(x[, 1], IRanges(x[, 2], x[, 3]), strand = "*")
hist(log10(width(k27ac)))
k27ac <- k27ac[width(k27ac) < 10^4]
save(k27ac, file = "~/Documents/eRNA/link/UCSC/H3K27ac_GM.RData")

###################################################
### Load NFKB peaks for GM
###################################################
x <- read.table("~/Documents/eRNA/link/UCSC/GM12878_NFKB_narrowPeak.bed")
hist(x[, 8], breaks = 100)
nfkb <- GRanges(x[, 1], IRanges(x[, 2], x[, 3]), strand = "*", level = x[, 8])
hist(log10(width(nfkb)))
nfkb <- nfkb[width(nfkb) < 10^4]
save(nfkb, file = "~/Documents/eRNA/link/UCSC/NFKB_GM.RData")

###################################################
### Load DNase-seq peaks for GM
###################################################
x <- read.table("~/Documents/eRNA/link/UCSC/GM12878-DS9432.hotspot.obs.twopass.merge150.wgt10.zgt2.bed")
hist(x[, 8], breaks = 100)
dnase <- GRanges(x[, 1], IRanges(x[, 2], x[, 3]), strand = "*", level = x[, 8])
hist(log10(width(dnase)))
dnase <- dnase[width(dnase) < 10^4]
save(dnase, file = "~/Documents/eRNA/link/UCSC/DNase_GM.RData")

###################################################
### Load transcribed regions by GRO-seq HOMER
###################################################
lab <- paste("GM", c("0h", "30min", "1h", "2h", "4h", "6h", "9h", "12h", "18h", "24h", "48h", "72h"), sep = "")
save_homer <- function(lab) {
  infile <- paste("R/Data/homer_region/", lab, ".txt", sep = "")
  outfile <- paste("R/Data/HOMER_", lab, ".RData", sep = "")
  y <- read.table(file = infile)
  txn <- GRanges(y[, 2], IRanges(y[, 3], y[, 4]), strand = y[, 5], read.depth = y[, 6])
  txn <- txn[order(y[, 6], decreasing = T)]
  txn.plus <- txn[strand(txn) == "+"]
  txn.minus <- txn[strand(txn) == "-"]
  txn.tss <- list(
    plus = promoters(txn.plus, upstream = 0, downstream = 1),
    minus = promoters(txn.minus, upstream = 0, downstream = 1)
  )
  save(txn, txn.plus, txn.minus, txn.tss, file = outfile)
}
i <- 0
for (i in c(1:length(lab))) {
  save_homer(lab = lab[i])
}
save_homer(lab = "meta_GM")
