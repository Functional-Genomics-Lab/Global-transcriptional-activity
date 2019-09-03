# Last update: Feb 6 2018

# source("https://bioconductor.org/biocLite.R")
# biocLite("BSgenome.Hsapiens.UCSC.hg19")
# biocLite("BSgenome.Hsapiens.UCSC.hg19.masked")
# biocLite("BSgenome.Mmusculus.UCSC.mm9")
# biocLite("BSgenome.Mmusculus.UCSC.mm9.masked")
# biocLite("r3Cseq")

setwd("~/Documents/IFN_enhancer/R/Data/4C_01292018/")
library(r3Cseq)
library(BSgenome.Hsapiens.UCSC.hg19)
library(RColorBrewer)
#######################################################
# Method_1
#######################################################
track_plot <- function(my_obj, tag, scope_size = 200, xaxt = FALSE, ylim = c(0, 4500)) {
  viewpoint <- getViewpoint(my_obj)
  anchor_point <- as.integer((end(viewpoint) + start(viewpoint)) / 2)
  scope <- GRanges(viewpoint)
  start(scope) <- anchor_point - scope_size * 10^3
  end(scope) <- anchor_point + scope_size * 10^3
  frags <- GRanges(expRPM(my_obj))
  hits <- findOverlaps(frags, scope)
  frags <- frags[queryHits(hits)]
  x <- as.integer((end(frags) + start(frags)) / 2) - anchor_point
  x <- x / 1000
  y <- mcols(frags)["RPMs"][, 1]
  if (xaxt) {
    plot(x, y, type = "l", ylim = ylim, ylab = "RPM", xlab = "")
  }
  else {
    plot(x, y, type = "l", xaxt = "n", ylim = ylim, ylab = "RPM", xlab = "")
  }
  legend("topright", tag, bty = "n")
}

method_1 <- function(experiment, viewpoint_chromosome, viewpoint_primer_forward, viewpoint_primer_reverse, restrictionEnzyme, distance = 2e5, ylim = c(0, 4000)) {
  outfile <- paste(experiment, as.character(distance), sep = "")
  outfile <- paste(outfile, ".pdf", sep = "")
  pdf(outfile)
  tags <- paste(experiment, c("0h", "6h", "12h", "18h", "24h"), sep = "")
  # tags <- paste(experiment, c("0h", "6h"), sep="")
  par(mfrow = c(length(tags), 1), mai = c(0, 1, 0, 1), omi = c(1, 0, 1, 0))
  i <- 0
  obj.list <- list()
  for (tag in tags) {
    expFile <- paste(tag, ".trim.bam", sep = "")
    expFile <- paste("bam/", expFile, sep = "")
    my_obj <- new("r3Cseq",
      organismName = "hg19", alignedReadsBamExpFile = expFile,
      isControlInvolved = F, viewpoint_chromosome = viewpoint_chromosome,
      viewpoint_primer_forward = viewpoint_primer_forward,
      viewpoint_primer_reverse = viewpoint_primer_reverse,
      expLabel = tag, restrictionEnzyme = restrictionEnzyme
    )
    getRawReads(my_obj)
    getReadCountPerRestrictionFragment(my_obj)
    calculateRPM(my_obj)
    getInteractions(my_obj)
    expResults <- expInteractionRegions(my_obj)
    # export3Cseq2bedGraph(my_obj)
    i <- i + 1
    obj.list[[i]] <- my_obj
    if (i == length(tags)) {
      track_plot(my_obj, tag = tag, scope_size = distance, xaxt = TRUE, ylim = ylim)
    }
    else {
      track_plot(my_obj, tag = tag, scope_size = distance, xaxt = FALSE, ylim = ylim)
    }
  }
  dev.off()
  return(obj.list)
}

local_interaction <- function(my_obj, scope_size = 1000) {
  viewpoint <- getViewpoint(my_obj)
  anchor_point <- as.integer((end(viewpoint) + start(viewpoint)) / 2)
  scope <- GRanges(viewpoint)
  start(scope) <- anchor_point - scope_size * 10^3
  end(scope) <- anchor_point + scope_size * 10^3
  frags <- GRanges(expRPM(my_obj))
  hits <- findOverlaps(frags, scope)
  frags <- frags[queryHits(hits)]
  x <- as.integer((end(frags) + start(frags)) / 2) - anchor_point
  x <- x / 1000
  y <- mcols(frags)["RPMs"][, 1]
  return(list(x, y))
}
#############################################################################################
timepoints <- c("0h", "6h", "12h", "18h", "24h")

IFNb1_ED <- method_1(
  experiment <- "IFNb1-ED_",
  viewpoint_chromosome = "chr9",
  viewpoint_primer_forward = "CACTCAGTCTATTAGCATTAGATC",
  viewpoint_primer_reverse = "AATATGTTTAGGTTTTATGGCTGGCTTTGGGGAAG",
  restrictionEnzyme = "EcoRI",
  distance = 500
)
save(IFNb1_ED, file = "IFNB1_ED.Robj")
track_plot(IFNb1_ED, "IFNb1-ED_", scope_size = 200, xaxt = FALSE, ylim = c(0, 4500))

IFNb1_HD <- method_1(
  experiment <- "IFNb1-HD_",
  viewpoint_chromosome = "chr9",
  viewpoint_primer_forward = "AACGCAGATAATAAGATC",
  viewpoint_primer_reverse = "AAGCTTCAGAGGTCAAAGGGCTAAT",
  restrictionEnzyme = "HindIII",
  distance = 500
)
save(IFNb1_HD, file = "IFNB1_HD.Robj")

TNFSF10 <- method_1(
  experiment <- "TNFSF10-1_",
  viewpoint_chromosome = "chr3",
  viewpoint_primer_forward = "CTCCATGGGTCAGGTTGTTTTCTTGATGAA",
  viewpoint_primer_reverse = "CCAATGCTTCTCCCATGGATC",
  restrictionEnzyme = "EcoRI",
  distance = 1000,
  ylim = c(0, 5000)
)
save(TNFSF10, file = "TNFSF10.Robj")

####################################################################################################
load("IFNB1_ED.Robj")
load("TNFSF10.Robj")
bar.col <- brewer.pal(3, "Blues")[3:1]

pdf("Total_local_ineractions.pdf")
local.sum <- c()
dist.list <- c(10000, 1000, 200)
dist.names <- c("10M", "1M", "200kb")
for (my_obj in IFNb1_ED) {
  tp <- local_interaction(my_obj, scope_size = max(dist.list))
  for (d in dist.list) {
    local.sum <- append(local.sum, sum(tp[[2]][abs(tp[[1]]) < d]))
  }
}
local.sum <- matrix(local.sum, nrow = 3, byrow = F, dimnames = list(dist.names, timepoints))
barplot(local.sum,
  names.arg = timepoints, beside = T, ylab = "Interactions (RPM)", ylim = c(0, 55000),
  col = bar.col, main = "IFNb1"
)
legend("topright", dist.names, col = bar.col, pch = 15)


local.sum <- c()
for (my_obj in TNFSF10) {
  tp <- local_interaction(my_obj, scope_size = max(dist.list))
  for (d in dist.list) {
    local.sum <- append(local.sum, sum(tp[[2]][abs(tp[[1]]) < d]))
  }
}
local.sum <- matrix(local.sum, nrow = 3, byrow = F, dimnames = list(dist.names, timepoints))
barplot(local.sum,
  names.arg = timepoints, beside = T, ylab = "Interactions (RPM)", ylim = c(0, 55000),
  col = bar.col, main = "TNFSF10"
)
legend("topright", dist.names, col = bar.col, pch = 15)

dev.off()


# method_1(
#   experiment <- "IFNb1-HD_",
#   viewpoint_chromosome='chr9',
#   viewpoint_primer_forward='AACGCAGATAATAAGATC',
#   viewpoint_primer_reverse='AAGCTTCAGAGGTCAAAGGGCTAAT',
#   restrictionEnzyme='HindIII',
#   distance=50
# )
# method_1(
#   experiment <- "PARP94_",
#   viewpoint_chromosome='chr3',
#   viewpoint_primer_forward='TTTACCCTATTGTGGCATATATGAATTC',
#   viewpoint_primer_reverse='GATCTGCACCTGAGCCAAAGAAAC',
#   restrictionEnzyme='EcoRI',
#   distance=50
# )
