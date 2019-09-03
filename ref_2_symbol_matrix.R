# Converte rownames of a matrix from Refseq to gene symbol
load("~/Documents/IFN_enhancer/R/Data/ref2sym_dedup.RData")
ref_2_symbol <- function(x) {
  if (class(x) == "matrix" | class(x) == "data.frame") {
    refseq <- rownames(x)
    stopifnot(any(refseq %in% ref2sym[, 1]))
    symbol <- unique(ref2sym[which(ref2sym[, 1] %in% refseq), 2])
    y <- matrix(0, ncol = ncol(x), nrow = length(symbol))
    for (i in 1:length(symbol)) {
      lab <- ref2sym[which(ref2sym[, 2] == symbol[i]), 1]
      lab <- lab[lab %in% refseq]
      if (length(lab) > 1) {
        y[i, ] <- colMeans(x[lab, ])
      }
      if (length(lab) == 1) {
        y[i, ] <- as.numeric(x[lab, ])
      }
    }
    rownames(y) <- symbol
    colnames(y) <- colnames(x)
    return(y)
  }
  if (class(x) == "character") {
    return(unique(ref2sym[ref2sym[, 1] %in% x, 2]))
  }
}
