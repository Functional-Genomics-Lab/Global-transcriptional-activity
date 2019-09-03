# Prepare data for differential expression analysis
setwd("~/Documents/IFN_enhancer/R/")
library(GenomicRanges)
library(GenomicFeatures)
source("Script/ref_2_symbol_matrix.R")

###############################################################
# 1. Load data
###############################################################
load("Data/eRNA_GM.RData") # Loaded "eRNA"
load("Data/gene_expression_GM_unfiltered.RData") # Loaded "gene" and "gene.rpkm".

# Process raw-count table of genes
total <- read.table("Data/total_mapped_reads_GM_replicate.txt")[, 2]
y <- read.table("Data/genes_GM_replicate_filtered.counts", header = T) # This "filtered" means records with the same Refseq id removed
gene <- gene[mcols(gene)$id %in% rownames(y)]
y <- y[as.character(mcols(gene)$id), ]
# We prefer to analyze at the level of genes instead of transcripts
countData <- ref_2_symbol(y)
countData <- t(apply(countData, 1, function(x) {
  return(as.integer(x))
}))
colnames(countData) <- colnames(y)
rpkmData <- t(t(y / width(gene)) / total) * 10^9
rpkmData <- ref_2_symbol(rpkmData)

# Process raw-count table of eRNAs
y <- read.table("Data/eRNA_GM_replicate.counts", header = T)
countData <- rbind(countData, y)
rpkmData <- rbind(rpkmData, t(t(y / width(eRNA)) / total) * 10^9)

save(countData, rpkmData, file = "Data/DEseq_input.RData")

# gene.rpkm matrix with rows representing a gene symbol
# load("Data/gene_expression_GM_unfiltered.RData")
# gene.rpkm.symbol <- ref_2_symbol(gene.rpkm)
# save(gene.rpkm.symbol, file="Data/gene.rpkm.symbol.RData")
