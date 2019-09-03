# Find the IRF binding sites in the enhancer sequences

setwd("~/Documents/IFN_enhancer/R/")
library(GenomicRanges)
library(GenomicFeatures)
library(gplots)
library(rtracklayer)

eRNA <- read.table("Data/eRNA_GM.tss.bed")
rownames(eRNA) <- as.character(eRNA[, 4])

colclass <- rep("character", 6)
colclass[2] <- "integer"
irf <- read.table("Data/Inducible_eRNA.IRF.positions", colClasses = colclass)
elist <- unique(irf[, 1])
irf_pos <- NULL
for (i in elist) {
  tp <- irf[irf[, 1] == i, ]
  lab <- which(tp[, 6] == max(tp[, 6]))[1]
  irf_pos <- c(irf_pos, tp[lab, 2])
}

irf_bed <- eRNA[elist, ]
for (i in 1:nrow(irf_bed)) {
  irf_bed[i, 2] <- irf_bed[i, 2] - irf_pos[i]
}
irf_bed[, 3] <- irf_bed[, 2] + 1
write.table(irf_bed, file = "Data/Inducible_eRNA.IRF.bed", col.names = F, row.names = F, quote = F, sep = "\t")
