load("Data/ref2sym_dedup.RData")
sym_2_ref <- function(x) {
  if (class(x) == "character") {
    return(unique(ref2sym[ref2sym[, 2] %in% x, 1]))
  }
}
