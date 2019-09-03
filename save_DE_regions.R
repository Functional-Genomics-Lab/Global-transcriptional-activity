save_bed <- function(region, ids, tag, liftover = TRUE) {
  source("Script/liftover_function.R")
  x <- region[mcols(region)$id %in% ids]
  x <- promoters(x, upstream = 2000, downstream = 2000)
  if (liftover) {
    x <- liftOver_function(x)
  }
  output <- data.frame(x)
  output <- cbind(output[, 1:3], output[, 6], rep(0, nrow(output)), output[, 5])
  file_name <- paste("Data/", tag, sep = "")
  write.table(output, file = file_name, row.names = F, col.names = F, quote = F, sep = "\t")
}

de_2_bed <- function(x, cols, class, tag) {
  source("Script/sym_2_ref.R")
  if (length(cols) == 1) {
    cols <- rep(cols, 2)
  }
  x <- x[, cols]
  ids <- rownames(x)
  any_pos <- apply(x == 1, 1, any)
  any_neg <- apply(x == -1, 1, any)
  # Case 1
  if (class == "stable") {
    tp <- apply(x == 0, 1, all)
  }
  # Case 2
  if (class == "up") {
    tp <- cbind(any_pos, !any_neg)
    tp <- apply(tp, 1, all)
  }
  # Case 3
  if (class == "down") {
    tp <- cbind(!any_pos, any_neg)
    tp <- apply(tp, 1, all)
  }
  ids <- ids[tp]
  print(c("Number of candidates:", length(ids)))
  save_bed(gene, sym_2_ref(ids), paste(tag, ".gene_tss.bed", sep = ""))
  save_bed(eRNA, ids, paste(tag, ".eRNA_tss.bed", sep = ""))
  return(ids)
}

test <- de_2_bed(deg_matrix_direction, 1:11, class = "stable", tag = "Stable")
test <- de_2_bed(deg_matrix_direction, 4:7, class = "up", tag = "Early_up")
test <- de_2_bed(deg_matrix_direction, 4:7, class = "down", tag = "Early_down")
test <- de_2_bed(deg_matrix_direction, 10:11, class = "up", tag = "Late_up")
test <- de_2_bed(deg_matrix_direction, 10:11, class = "up", tag = "Late_down")
