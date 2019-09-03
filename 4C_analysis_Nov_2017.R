setwd("~/Documents/IFN_enhancer/R/Data/4C_Nov_2017/")
library(r3Cseq)
library(BSgenome.Hsapiens.UCSC.hg19)

#######################################################
# Method_1
#######################################################

tags <- paste("L2_", c("0h", "6h", "12h", "18h", "24h"), sep = "")

for (tag in tags) {
  expFile <- paste(tag, ".sorted.bam", sep = "")

  my_obj <- new("r3Cseq",
    organismName = "hg19", alignedReadsBamExpFile = expFile,
    isControlInvolved = F, viewpoint_chromosome = "chr9",
    viewpoint_primer_forward = "gccccattcccagagtaaggaatgttgcacaaacaaacccagaagaatctag",
    viewpoint_primer_reverse = "ggaggaggttgaagtaatatgtttaggttttatggctggctttggggaag",
    expLabel = tag, restrictionEnzyme = "EcoRI"
  )
  getRawReads(my_obj)
  getReadCountPerRestrictionFragment(my_obj)
  calculateRPM(my_obj)
  getInteractions(my_obj)

  expResults <- expInteractionRegions(my_obj)
  export3Cseq2bedGraph(my_obj)
  pdf(paste(tag, ".pdf", sep = ""))
  par(mfrow = c(1, 1))
  plotOverviewInteractions(my_obj)
  plotInteractionsNearViewpoint(my_obj, distance = 2e5)
  dev.off()
  tp <- as.data.frame(expResults)
  tp <- tp[tp["space"] == "chr9", ]
  bed_data <- data.frame(chr = tp[, 1], start = tp[, 2], end = tp[, 3], name = row.names(tp), score = tp[, 6] / max(tp[, 6]) * 1000)
  write.table(bed_data, file = paste(tag, ".chr9_interactions.bed", sep = ""), quote = F, row.names = F, col.names = F)
}


#######################################################
# Method_1: time points as replicates
#######################################################
# 'r3CSeqInBatch' kept asking control bam files, which doesn't make sense for our purpose. So I decided to go with a merged bam file

my_obj <- new("r3Cseq",
  organismName = "hg19", alignedReadsBamExpFile = "merged.bam",
  isControlInvolved = F, viewpoint_chromosome = "chr9",
  viewpoint_primer_forward = "gccccattcccagagtaaggaatgttgcacaaacaaacccagaagaatctag",
  viewpoint_primer_reverse = "ggaggaggttgaagtaatatgtttaggttttatggctggctttggggaag",
  expLabel = tag, restrictionEnzyme = "EcoRI"
)

getRawReads(my_obj)
getReadCountPerRestrictionFragment(my_obj)
calculateRPM(my_obj)
getInteractions(my_obj)

expResults <- expInteractionRegions(my_obj)
export3Cseq2bedGraph(my_obj)
pdf(paste("merged.pdf", sep = ""))
par(mfrow = c(1, 1))
plotOverviewInteractions(my_obj)
plotInteractionsNearViewpoint(my_obj, distance = 2e5)
dev.off()
tp <- as.data.frame(expResults)
tp <- tp[tp["space"] == "chr9", ]
bed_data <- data.frame(chr = tp[, 1], start = tp[, 2], end = tp[, 3], name = row.names(tp), score = tp[, 6] / max(tp[, 6]) * 1000)
write.table(bed_data, file = "merged.chr9_interactions.bed", quote = F, row.names = F, col.names = F)
#######################################################
# Method_2
#######################################################

expFile <- "4C12h.bam"
contrFile <- "4C0h.bam"

my_obj <- new("r3Cseq",
  organismName = "hg19", alignedReadsBamExpFile = expFile,
  alignedReadsBamContrFile = contrFile, isControlInvolved = TRUE, viewpoint_chromosome = "chr9",
  viewpoint_primer_forward = "tgttgcacaaacaaacccaga",
  viewpoint_primer_reverse = "TATGGCTGGCTTTGGGGAAG",
  expLabel = "GM12h", contrLabel = "GM0h", restrictionEnzyme = "EcoRI"
)


getRawReads(my_obj)
getReadCountPerRestrictionFragment(my_obj)
calculateRPM(my_obj)
getInteractions(my_obj)

expResults <- expInteractionRegions(my_obj)
contrResults <- contrInteractionRegions(my_obj)
export3CseqRawReads2bedGraph(my_obj)
export3Cseq2bedGraph(my_obj)
plotOverviewInteractions(my_obj)
plotInteractionsNearViewpoint(my_obj, distance = 2e5)
