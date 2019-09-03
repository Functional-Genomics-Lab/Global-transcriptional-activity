setwd("~/Documents/IFN_enhancer/")
library(GenomicRanges)
library(GenomicFeatures)
library(gplots)
library(rtracklayer)

###################################################
### Load Data
###################################################
load(file = "~/Documents/eRNA/link/UCSC/hg18_refseq.RData")

# extend all the genes by 10kb downstream
x <- read.table("R/Data/chromInfo.txt", header = F)
chr.ends <- x[, 2] - 1
names(chr.ends) <- x[, 1]
gene.ext <- gene
extend <- 10^4
upstream <- 10^3

lab <- which(strand(gene.ext) == "+")
end(gene.ext[lab]) <- apply(rbind(end(gene.ext[lab]) + extend, chr.ends[as.character(seqnames(gene.ext[lab]))]), 2, min)
start(gene.ext[lab]) <- apply(rbind(start(gene.ext[lab]) - upstream, rep(1, length(lab))), 2, max)

lab <- which(strand(gene.ext) == "-")
start(gene.ext[lab]) <- apply(rbind(start(gene.ext[lab]) - extend, rep(1, length(lab))), 2, max)
end(gene.ext[lab]) <- apply(rbind(end(gene.ext[lab]) + upstream, chr.ends[as.character(seqnames(gene.ext[lab]))]), 2, min)

load(file = "~/Documents/eRNA/link/UCSC/H3K27ac_GM.RData")
load(file = "~/Documents/eRNA/link/UCSC/H3K4me1_GM.RData")

###################################################
### Number of new eRNAs
###################################################
lab <- paste("GM", c("0h", "30min", "1h", "2h", "4h", "6h", "9h", "12h", "18h", "24h", "48h", "72h"), sep = "")

load("R/Data/HOMER_GM0h.RData")
gm_0h <- txn

newer_eRNA <- function(homer_file) { # Enhancers even without markers before treatment
  load(file = homer_file)
  # Intergenic only
  hits <- findOverlaps(txn, gene.ext, ignore.strand = T)
  ta <- txn[-unique(queryHits(hits))]
  # non-overlap with enhancer markers
  markers <- c(k4me1, k27ac)
  ta <- ta[ta %outside% markers, ]
  # Remove pre-existing eRNA
  hits <- findOverlaps(ta, gm_0h, ignore.strand = T)
  tb <- ta[-unique(queryHits(hits))]
  print(c("Enhancers without neither marker nor eRNA:", length(tb)))
  return(tb)
}

new_eRNA <- function(homer_file) { # Enhancers with markers but no eRNA before treatment
  load(file = homer_file)
  # Intergenic only
  hits <- findOverlaps(txn, gene.ext, ignore.strand = T)
  ta <- txn[-unique(queryHits(hits))]
  # Remove pre-existing enhancer
  markers <- c(k4me1, k27ac)
  # tb <- ta[ta %outside% markers, ]
  ta <- ta[ta %over% markers, ] # All expressed enhancers with markers
  # Remove pre-existing eRNA
  hits <- findOverlaps(ta, gm_0h, ignore.strand = T)
  tb <- ta[-unique(queryHits(hits))] # Newly activated enhancers
  print(c("Enhancers with marker:", length(ta)))
  print(c("Enhancers with marker but no eRNA:", length(tb)))
  return(tb)
}
ea <- new_eRNA("R/Data/HOMER_meta_GM.RData")
eb <- newer_eRNA("R/Data/HOMER_meta_GM.RData")
