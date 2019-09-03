setwd("~/Documents/IFN_enhancer/R/Data/4C_Nov_2017/")
library(GenomicRanges)

#########################################################################
# Extract genomic fragments that may represent true interactions
#########################################################################
get_useful_fragments <- function(files, chrom = "all", l_motif = c(6, 4)) {
  if (chrom == "all") {
    chrom <- paste("chr", c(c(1:22), "X"), sep = "")
  }

  x <- read.table(files[1])
  x[, 4] <- 1
  y <- read.table(files[2])
  y[, 4] <- 2

  for (i in c(1:length(chrom))) {
    z <- rbind(x[as.character(x[, 1]) == chrom[i], ], y[as.character(y[, 1]) == chrom[i], ])
    z <- sort(GRanges(seqnames = z[, 1], IRanges(z[, 2], z[, 3]), scores = z[, 4]))

    lab <- mcols(z)[, "scores"]
    if (length(lab) == 0) {
      next
    }
    lab <- lab[1:(length(lab) - 1)] + lab[2:length(lab)]
    lab <- which(lab == 3)
    print(c(chrom[i], length(lab)))

    if (i == 1) {
      fragments_4c <- GRanges(rep(chrom[i], length(lab)), IRanges(start(z[lab]), start(z[lab + 1]) + l_motif[mcols(z)[lab + 1, "scores"]]))
    }
    else {
      fragments_4c <- c(fragments_4c, GRanges(rep(chrom[i], length(lab)), IRanges(start(z[lab]), start(z[lab + 1]) + l_motif[mcols(z)[lab + 1, "scores"]])))
    }
  }
  return(fragments_4c)
}

files <- c("../hg19/GAATTC.bed", "../hg19/GATC.bed")
interaction_fragments <- get_useful_fragments(files, "all")
# We need to exclude primer fragments for later analysis
primer_region <- GRanges(seqnames = "chr9", IRanges(21095600, 21095790))
hits <- findOverlaps(primer_region, interaction_fragments)
primer_fragments <- interaction_fragments[subjectHits(hits)]
interaction_fragments <- interaction_fragments[-subjectHits(hits)]

#########################################################################
# Load reads (pair-end)
#########################################################################
# reads <- read.table('merged.condordant.pairs')  # Slow step!
# reads <- GRanges(seqnames=reads[,2], IRanges(reads[,3], reads[,4]+75))

chrom_all <- paste("chr", c(c(1:22), "X"), sep = "")
rpm <- rep(0, length(interaction_fragments))
stats <- data.frame(chr = chrom_all, total_concordant = rep(0, length(chrom_all)), interaction_reads = rep(0, length(chrom_all)))
for (i in 1:length(chrom_all)) {
  qr <- reads[as.character(seqnames(reads)) == chrom_all[i]] # query
  hits <- findOverlaps(qr, primer_fragments)
  if (length(hits) > 0) {
    qr <- qr[-queryHits(hits)]
  }
  hits <- findOverlaps(qr, interaction_fragments, type = "within")

  stats[i, "total_concordant"] <- length(qr)
  stats[i, "interaction_reads"] <- length(hits)
  rpm <- rpm + countRnodeHits(hits)
}

rpm <- rpm / sum(rpm) * 10^6
mcols(interaction_fragments)["rpm"] <- rpm
# Find read pairs within fragments

barplot(stats[, 3] / stats[, 2], names.arg = chrom_all)

# Enrichment of intrachromosomal interaction
chrom_size <- read.table("../chromInfo.hg19.txt", row.names = 1) / 10^6
inter_reads <- sum(stats[stats[, 1] != "chr9", 2]) / sum(chrom_size[chrom_all[-9], ])
intra_reads <- sum(stats[stats[, 1] == "chr9", 2]) / sum(chrom_size[chrom_all[9], ])

inter_reads <- sum(stats[stats[, 1] != "chr9", 3]) / sum(chrom_size[chrom_all[-9], ])
intra_reads <- sum(stats[stats[, 1] == "chr9", 3]) / sum(chrom_size[chrom_all[9], ])
