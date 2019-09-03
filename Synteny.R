library(GenomicRanges)

################################################################################################
ref_2_symbol_local <- function(x) {
  output <- NULL
  for (i in 1:length(x)) {
    if (x[i] %in% ref2sym[, 1]) {
      output <- c(output, ref2sym[which(ref2sym[, 1] == x[i])[1], 2])
    }
    else {
      output <- c(output, NA)
    }
  }
  return(output)
}

synteny_background <- function() {
  # Define background pairs
  load("Data/gene_expression_GM_unfiltered.RData")
  pmt <- promoters(gene, 250, 250)

  # extend all the genes by 10kb downstream
  x <- read.table("Data/chromInfo.txt", header = F)
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

  # random sequences in the promoter flanking region
  n_random <- 10^4
  ep <- read.csv("Data/Induced_EP_for_synteny.csv", header = T, colClasses = "character")
  tp <- as.integer(sub("kb", "", ep$EP_d)) * 1000
  ext <- sample(tp, n_random, replace = T)
  hist(as.integer(sub("kb", "", ep$EP_d)), breaks = 20, col = rgb(0.5, 0.0, 0.0, 0.5), freq = F)
  hist(ext / 1000, breaks = 20, col = rgb(0.5, 0.5, 0.5, 0.5), freq = F, add = T)

  par(mfcol = c(3, 2))
  pmt_sel <- NULL
  rnd_sel <- NULL
  ext_sel <- NULL
  all_lab <- c(1:length(ext))
  lab <- all_lab # Randomly selsect promoters and distal regions, use "lab" to keep labels of distances that have not been used

  while (length(lab) > 10) {
    re_lab <- sample(c(1:length(pmt)), replace = F, length(lab)) # Randomly selsect promoters
    t_pmt <- pmt[re_lab]
    str_lab <- which(strand(t_pmt) == "-")
    t_shift <- ext
    t_shift[str_lab] <- ext[str_lab] * (-1)
    t_rnd <- shift(pmt[re_lab], t_shift)
    # Remove cases where random regions were out of chromosome
    chr_max <- chr.ends[as.character(seqnames(t_rnd))]
    chr_lab <- which((start(t_rnd) < 0) | (end(t_rnd) > chr_max))
    # Remove cases where random regions were not intergenic
    hits <- findOverlaps(t_rnd, gene.ext, ignore.strand = T)
    gene_lab <- unique(queryHits(hits))
    lab <- unique(c(chr_lab, gene_lab))
    t_pmt <- t_pmt[-lab]
    t_rnd <- t_rnd[-lab]
    t_ext <- ext[-lab]
    if (is.null(pmt_sel)) {
      pmt_sel <- t_pmt
      rnd_sel <- t_rnd
      ext_sel <- t_ext
    } else {
      pmt_sel <- c(pmt_sel, t_pmt)
      rnd_sel <- c(rnd_sel, t_rnd)
      ext_sel <- c(ext_sel, t_ext)
    }
    if (length(lab) == 0) {
      break
    }
    ext <- ext[lab]
    hist(as.integer(sub("kb", "", ep$EP_d)), breaks = 40, col = rgb(0.5, 0.0, 0.0, 0.5), freq = F, main = length(lab))
    hist(ext_sel / 1000, breaks = 40, col = rgb(0.5, 0.5, 0.5, 0.5), freq = F, add = T)
  }

  # Output files
  out_pmt <- as.data.frame(pmt_sel)[, c(1, 2, 3, 6, 4, 5)]
  out_rnd <- as.data.frame(rnd_sel)[, c(1, 2, 3, 6, 4, 5)]
  out_rnd$id <- paste("rnd", c(1:nrow(out_rnd)), sep = "_")

  # Output a table for the background pairs
  out_table <- data.frame(id_a = out_pmt[, 4], id_b = out_rnd[, 4], symbol = ref_2_symbol_local(out_pmt[, 4])) # "ref_2_symbol_local" is time-consuming
  # out_table <- data.frame(id_a=out_pmt[,4], id_b=out_rnd[,4], symbol=NA)
  start <- apply(cbind(out_pmt[, 2], out_rnd[, 2]), 1, min) # Record the spanning region
  end <- apply(cbind(out_pmt[, 3], out_rnd[, 3]), 1, max)
  out_table$coord <- paste(out_rnd[, 1], start, end)
  out_table$EP_d <- out_rnd[, 2] - out_pmt[, 2] # Record the distance
  str_lab <- which(out_pmt[, 6] == "-")
  out_table$EP_d[str_lab] <- out_table$EP_d[str_lab] * (-1)
  out_table$EP_d <- paste(as.integer(out_table$EP_d / 1000), "kb", sep = "")
  write.csv(out_table, "Data/Synteny_bg_pairs.csv", row.names = F)

  out_pmt <- as.data.frame(pmt_sel)[, c(1, 2, 3, 6, 4, 5)]
  out_rnd <- as.data.frame(rnd_sel)[, c(1, 2, 3, 6, 4, 5)]
  out_rnd$id <- paste("rnd", c(1:nrow(out_rnd)), sep = "_")
  write.table(out_pmt, "Data/Synteny_bg_promoter.bed", col.names = F, row.names = F, sep = "\t", quote = F)
  write.table(out_rnd, "Data/Synteny_bg_random200kb.bed", col.names = F, row.names = F, sep = "\t", quote = F)
}

ortholog_number <- function(ep, enh, pmt) {
  out <- c(
    nrow(ep),
    sum(ep[, 1] %in% pmt[, 4]),
    sum(ep[, 2] %in% enh[, 4]),
    sum((ep[, 1] %in% pmt[, 4]) & (ep[, 2] %in% enh[, 4]))
  )
  return(out)
}

synteny_analysis <- function(ep, enh, pmt) {
  ep <- ep[((ep[, 1] %in% pmt[, 4]) & (ep[, 2] %in% enh[, 4])), ]
  str_sign <- c(1, -1)
  names(str_sign) <- c("+", "-")
  synteny <- c(ep[1, 1], ep[1, 2], ep[1, 3], ep[1, 4], ep[1, 5], 1e+20)
  for (i in 1:nrow(ep)) {
    te <- enh[as.character(enh[, 4]) == ep[i, 2], ]
    tp <- pmt[which(as.character(pmt[, 4]) == ep[i, 1]), ]
    if (max(c(nrow(te), nrow(tp))) > 10) {
      next
    }
    for (j in 1:nrow(te)) {
      for (k in 1:nrow(tp)) {
        temp <- c(ep[i, 1], ep[i, 2], ep[i, 3], ep[i, 4], ep[i, 5], 1e+20)
        if (te[j, 1] == tp[k, 1]) { # If EP in the same chromosome, calculate distance
          temp[6] <- (te[j, 2] + te[j, 3] - tp[k, 2] - tp[k, 3]) / 2000 * str_sign[tp[k, 6]]
        }
        synteny <- rbind(synteny, temp)
      }
    }
  }
  if (nrow(synteny) <= 2) {
    return(NULL)
  }
  synteny <- synteny[-1, ]
  rownames(synteny) <- NULL
  colnames(synteny) <- c("Refseq", "Enhancer", "Symbol", "Human_coord", "Human_distance", "Mouse_distance")
  synteny <- as.data.frame(synteny)
  synteny$Human_distance <- as.integer(sub("kb", "", synteny$Human_distance))
  synteny$Mouse_distance <- as.numeric(as.character(synteny$Mouse_distance))

  # Select for the best match when liftover results in multiple regions
  synteny_sel <- synteny[1, ]
  pairs <- paste(synteny$Refseq, synteny$Enhancer)
  u_pairs <- unique(pairs)
  for (i in 1:length(u_pairs)) {
    tp <- synteny[pairs == u_pairs[i], ]
    lab <- abs(tp[, 6] - tp[, 5])
    min_dist <- tp[which(lab == min(lab))[1], 6]
    tp <- tp[1, ]
    tp[6] <- min_dist
    synteny_sel <- rbind(synteny_sel, tp)
  }
  if (nrow(synteny_sel) <= 2) {
    return(NULL)
  }
  synteny_sel <- synteny_sel[-1, ]
  rownames(synteny_sel) <- NULL
  colnames(synteny_sel) <- c("Refseq", "Enhancer", "Symbol", "Human_coord", "Human_distance", "Mouse_distance")
  synteny_sel <- as.data.frame(synteny_sel)

  # How many of the pairs were in the same chromosome?
  print(paste("Same chromosome ratio:", sum(abs(synteny_sel$Mouse_distance) < 1e+20) / nrow(synteny_sel)))
  # lab <- abs(synteny_sel$Mouse_distance)<500
  # plot(synteny_sel$Human_distance[lab], synteny_sel$Mouse_distance[lab], xlab="Human EP Distance (kb)", ylab="Mouse EP Distance (kb)")
  return(synteny_sel)
}

ep_dis_compare <- function(fg, bg, name) {
  fg_dis <- fg[, 5] - fg[, 6]
  fg_dis <- fg_dis[abs(fg_dis) < 1e+7]
  bg_dis <- bg[, 5] - bg[, 6]
  bg_dis <- bg_dis[abs(bg_dis) < 1e+7]
  plot(ecdf(abs(bg_dis)),
    do.points = F, verticals = T, xlim = c(0, 100), col = "blue",
    xlab = "Cross-species difference of EP distance (kb)", ylab = "Cumulative probability", main = ""
  )
  legend("bottomright", c("Human", name), col = c("red", "blue"), lwd = 1)
  if (length(fg_dis) > 10) {
    plot(ecdf(abs(fg_dis)), do.points = F, verticals = T, add = T, col = "red")
    pv <- ks.test(abs(bg_dis), abs(fg_dis), alternative = "less")$p.value
  }
  else {
    pv <- NA
  }
  text(80, 0.8, pv)
  return(pv)
}

ep_dis_wrapper <- function(tag, name) {
  out <- NULL

  ep <- read.csv("Data/Induced_EP_for_synteny.csv", header = T, colClasses = "character")
  enh <- read.table(paste("Data/Induced_EP.enhancer.", tag, ".bed", sep = ""), colClasses = c("character", "integer", "integer", "character", "integer", "character"))
  pmt <- read.table(paste("Data/Induced_EP.promoter.", tag, ".bed", sep = ""), colClasses = c("character", "integer", "integer", "character", "integer", "character"))
  fg <- synteny_analysis(ep, enh, pmt)
  out <- c(out, ortholog_number(ep, enh, pmt))
  out <- c(out, sum(abs(fg$Mouse_distance) < 500)) # Define EP within 500kb as conserved

  ep <- read.csv("Data/Synteny_bg_pairs.csv", header = T, colClasses = "character")
  enh <- read.table(paste("Data/Synteny_bg_random200kb.", tag, ".bed", sep = ""), colClasses = c("character", "integer", "integer", "character", "integer", "character"))
  pmt <- read.table(paste("Data/Synteny_bg_promoter.", tag, ".bed", sep = ""), colClasses = c("character", "integer", "integer", "character", "integer", "character"))
  bg <- synteny_analysis(ep, enh, pmt)
  out <- c(out, ortholog_number(ep, enh, pmt))
  out <- c(out, sum(abs(bg$Mouse_distance) < 1e+20))

  out <- c(out, ep_dis_compare(fg, bg, name))
  return(out)
}

conserved_ep <- function(tag) {
  ep <- read.csv("Data/Induced_EP_for_synteny.csv", header = T, colClasses = "character")
  output <- rep(0, nrow(ep))
  enh <- read.table(paste("Data/Induced_EP.enhancer.", tag, ".bed", sep = ""), colClasses = c("character", "integer", "integer", "character", "integer", "character"))
  pmt <- read.table(paste("Data/Induced_EP.promoter.", tag, ".bed", sep = ""), colClasses = c("character", "integer", "integer", "character", "integer", "character"))
  fg <- synteny_analysis(ep, enh, pmt)
  for (i in c(1:nrow(fg))) {
    tp <- fg[i, ]
    if (abs(tp[6]) < 500) {
      output[which((ep[, 1] == as.character(tp$Refseq)) & (ep[, 2] == as.character(tp$Enhancer)))] <- 1
    }
  }
  return(output)
}

################################################################################################

synteny_background()
# Show that fg and bg have the same distribution
pdf("Figure/Synteny_EP_Distance_histogram.pdf")
fg <- read.csv("Data/Induced_EP_for_synteny.csv", colClasses = "character")
bg <- read.csv("Data/Synteny_bg_pairs.csv", colClasses = "character")
par(mfrow = c(1, 1))
hist(as.integer(sub("kb", "", fg$EP_d)), breaks = 40, col = rgb(0.5, 0.0, 0.0, 0.5), freq = F, xlab = "EP distance (kb)", main = "")
hist(as.integer(sub("kb", "", bg$EP_d)), breaks = 40, col = rgb(0.5, 0.5, 0.5, 0.5), freq = F, add = T)
legend("topright", c("Foreground", "Background"), col = c(rgb(0.5, 0.0, 0.0, 0.5), rgb(0.5, 0.5, 0.5, 0.5)), pch = 15)
dev.off()

# Run liftOver to find ortholog regions

################################################################################################

genome_list <- c("PanTro4", "CalJac1", "mm9", "Rn4", "CavPor3", "OryCun1", "BosTau2", "CanFam2", "LoxAfr1", "DasNov1", "AnoCar1")
species_list <- c("Chimp", "Marmoset", "Mouse", "Rat", "Guinea Pig", "Rabbit", "Cow", "Dog", "Elephant", "Armadillo", "Lizard")

pdf("Figure/Synteny_distance_test.pdf")
par(mfrow = c(3, 2))
out <- NULL
for (i in 1:length(genome_list)) {
  print(species_list[i])
  out <- c(out, ep_dis_wrapper(genome_list[i], species_list[i]))
}
dev.off()


out_table <- matrix(out, byrow = T, ncol = 11)
colnames(out_table) <- c(
  "Pairs_fg", "Pairs_w_pmt_fg", "Pairs_w_enh_fg", "Pairs_w_pair_fg", "Same_chr_fg",
  "Pairs_bg", "Pairs_w_pmt_bg", "Pairs_w_enh_bg", "Pairs_w_pair_bg", "Same_chr_bg", "P-value"
)
rownames(out_table) <- species_list

pdf("Figure/Synteny_number_of_ortholog.pdf")
names <- species_list[-c(2, 10)]
# Fig 1:
# Number of all EP pairs with orthologs / Number of orthologous promoters
tp <- cbind(out_table[, 4] / out_table[, 2], out_table[, 9] / out_table[, 7])[-c(2, 10), ]
od <- order(tp[, 1], decreasing = T)
plot(tp[od, 1], pch = 8, col = "red", cex = 2, xlim = c(0, length(od) + 1), ylim = c(0, 1.1), xlab = "Species", ylab = "Percentage of Orthologous EP pairs")
points(tp[od, 2], pch = 16, col = "blue", cex = 1.5)
text(tp[od, 1] + 0.05, names[od], cex = 0.75)
legend("topright", c("Selected Pairs", "Background"), col = c("red", "blue"), pch = c(8, 16))

# Fig 2:
# Number of all EP pairs whose ortholog distance < 500kb / Number of orthologous promoters
tp <- cbind(out_table[, 5] / out_table[, 2], out_table[, 10] / out_table[, 7])[-c(2, 10), ]
od <- order(tp[, 1], decreasing = T)
plot(tp[od, 1], pch = 8, col = "red", cex = 2, xlim = c(0, length(od) + 1), ylim = c(0, 1.1), xlab = "Species", ylab = "Percentage of Orthologous EP pairs (<500kb)")
points(tp[od, 2], pch = 16, col = "blue", cex = 1.5)
text(tp[od, 1] + 0.05, names[od], cex = 0.75)
legend("topright", c("Selected Pairs", "Background"), col = c("red", "blue"), pch = c(8, 16))

dev.off()

################################################################################################

# Find highly conserved EP pairs

out <- NULL
for (i in 1:length(genome_list)) {
  out <- c(
    out,
    conserved_ep(genome_list[i])
  )
}

out_table <- matrix(out, byrow = F, ncol = length(genome_list))
colnames(out_table) <- species_list
out_table <- out_table[, order(colSums(out_table), decreasing = T)]
lab <- order(rowSums(out_table), decreasing = T)
out_table <- out_table[lab, ]

ep <- read.csv("Data/Induced_EP_for_synteny.csv", header = T, colClasses = "character")
ep <- ep[lab, ]

image(t(out_table), col = c("purple", "black"))
