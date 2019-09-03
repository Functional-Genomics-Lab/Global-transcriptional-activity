liftOver_function <- function(regions, from = "hg18", to = "hg19", ch = NULL) {
  library(rtracklayer)
  if (is.null(ch)) {
    if (to == "hg19") {
      ch <- import.chain("~/Documents/eRNA/link/UCSC/hg18ToHg19.over.chain")
    }
    if (to == "hg18") {
      ch <- import.chain("~/Documents/eRNA/link/UCSC/hg19ToHg18.over.chain")
    }
  }
  seqlevelsStyle(regions) <- "UCSC" # necessary
  regions <- liftOver(regions, ch)
  regions <- unlist(regions)
  if (length(regions) == 0) {
    return(NULL)
  }
  names(regions) <- c(1:length(regions))
  genome(regions) <- to
  return(regions)
}
