# source("https://bioconductor.org/biocLite.R")
# biocLite("BSgenome.Hsapiens.UCSC.hg19")
# biocLite("BSgenome.Hsapiens.UCSC.hg19.masked")
# biocLite("BSgenome.Mmusculus.UCSC.mm9")
# biocLite("BSgenome.Mmusculus.UCSC.mm9.masked")
# biocLite("r3Cseq")

setwd("~/Documents/eRNA/link/IFN project/4C/")
library(r3Cseq)
library(BSgenome.Mmusculus.UCSC.mm9)
library(BSgenome.Hsapiens.UCSC.hg19)

expFile <- "Myb_prom_FL_rep1.bam"
contrFile <- "Myb_prom_FB_rep1.bam"

my_obj <- new("r3Cseq",
  organismName = "mm9", alignedReadsBamExpFile = expFile,
  alignedReadsBamContrFile = contrFile, isControlInvolved = TRUE, viewpoint_chromosome = "chr10",
  viewpoint_primer_forward = "TCTTTGTTTGATGGCATCTGTT",
  viewpoint_primer_reverse = "AAAGGGGAGGAGAAGGAGGT",
  expLabel = "Myb_prom_FL_rep1", contrLabel = "MYb_prom_FB_rep1", restrictionEnzyme = "HindIII"
)


getRawReads(my_obj)
getReadCountPerRestrictionFragment(my_obj)
calculateRPM(my_obj)
getInteractions(my_obj)

expResults <- expInteractionRegions(my_obj)
contrResults <- contrInteractionRegions(my_obj)
export3CseqRawReads2bedGraph(my_obj) # Raw reads
export3Cseq2bedGraph(my_obj) # RPM
plotOverviewInteractions(my_obj)
plotInteractionsNearViewpoint(my_obj)
