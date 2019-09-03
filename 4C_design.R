setwd("~/Documents/eRNA/link/IFN project/")
source("~/Documents/eRNA/Rscript/sym_2_ref.R")
load("gene_expression_GM_unfiltered.RData")
load("eRNA_GM.RData")

primary_fragment <- function(x, center = 100000, motif_len) {
  # Restriction sites were mapped to eRNA/mRNA promoter regions
  # Promoter regions = TSS +/- 'center'
  name_list <- unique(as.character(x[, 1]))
  region_list <- NULL
  i <- 1
  for (name in name_list) {
    tp <- x[x[, 1] == name, 2]
    left <- max(tp[tp < center]) - center # Extract the fragment that covers the central portion
    right <- min(tp[tp > center]) - center
    if (substr(name, 1, 12) == "MetaEnhancer") {
      region <- eRNA[mcols(eRNA)$id == name]
    }
    else {
      lab <- sym_2_ref(name)
      lab <- lab[which(lab %in% as.character(mcols(gene)$id))[1]]
      lab <- which(as.character(mcols(gene)$id) %in% lab)[1]
      region <- gene[lab]
    }
    region <- promoters(region, upstream = 1, downstream = 1)
    if (as.character(strand(region[1])) == "+") {
      ref_pos <- start(region[1])
    }
    else {
      ref_pos <- start(region[1]) - 1
    }
    start(region) <- ref_pos + left
    end(region) <- ref_pos + right + motif_len
    region <- as.data.frame(region)
    if (is.null(region_list)) {
      region_list <- region
      i <- i + 1
    }
    else {
      region_list[i, ] <- region
      i <- i + 1
    }
  }
  return(region_list[, c(1, 2, 3, 7)])
}

secondary_fragment <- function(x, fragments_1, motif_len) {
  # Loop through all fragments
  name_list <- unique(as.character(x[, 1]))
  region_list <- NULL
  i <- 1
  better_size <- 2000
  for (name in name_list) {
    # Take the left/right-most ends, compare which one is closer to the TSS
    tp <- x[x[, 1] == name, 2]
    left <- min(tp)
    right <- max(tp)
    # Convert to genomic coordinates
    f1 <- fragments_1[fragments_1[, 4] == name, ]
    left <- left + f1[2] + motif_len - 1
    right <- right + f1[2] - 1

    # Get the genomic coordinates of TSS
    if (substr(name, 1, 12) == "MetaEnhancer") {
      region <- eRNA[mcols(eRNA)$id == name]
    }
    else {
      lab <- sym_2_ref(name)
      lab <- lab[which(lab %in% as.character(mcols(gene)$id))[1]]
      lab <- which(as.character(mcols(gene)$id) %in% lab)[1]
      region <- gene[lab]
    }
    tss <- start(promoters(region, upstream = 1, downstream = 1))

    # Case 1: when left and right ends are equally close to tss, choose a better fragment size
    region <- as.data.frame(region)
    if (abs(left - tss) == abs(right - tss)) {
      l_left <- left - f1[2]
      l_right <- f1[3] - right
      if (abs(l_left - better_size) < abs(l_right - better_size)) {
        start(region) <- f1[2]
        end(region) <- left
      }
      else {
        region[2] <- right
        region[3] <- f1[3]
      }
    }
    else { # Case 2: which one is closer to the TSS
      if (abs(left - tss) < abs(right - tss)) {
        region[2] <- f1[2]
        region[3] <- left
      }
      else {
        region[2] <- right
        region[3] <- f1[3]
      }
    }

    if (is.null(region_list)) {
      region_list <- region
    }
    else {
      region_list[i, ] <- region
    }
    i <- i + 1
  }
  return(region_list[, c(1, 2, 3, 7)])
}

secondary_fragment_both <- function(x, fragments_1, motif_len) {
  # Loop through all fragments
  name_list <- unique(as.character(x[, 1]))
  region_list <- NULL
  i <- 1
  better_size <- 2000
  for (name in name_list) {
    # 1. Get record for the whole region
    if (substr(name, 1, 12) == "MetaEnhancer") {
      region <- eRNA[mcols(eRNA)$id == name]
    }
    else {
      lab <- sym_2_ref(name)
      lab <- lab[which(lab %in% as.character(mcols(gene)$id))[1]]
      lab <- which(as.character(mcols(gene)$id) %in% lab)[1]
      region <- gene[lab]
    }
    region <- as.data.frame(region)
    # 2. Take the left and right-most ends
    tp <- x[x[, 1] == name, 2]
    left <- min(tp)
    right <- max(tp)
    # Convert to genomic coordinates
    f1 <- fragments_1[fragments_1[, 4] == name, ]
    # Left end fragment
    left <- left + f1[2] + motif_len - 1
    l_region <- region
    l_region[2] <- f1[2]
    l_region[3] <- left
    l_region[7] <- paste(region[7], "left", sep = "_")
    if (is.null(region_list)) {
      region_list <- l_region
    }
    else {
      region_list[i, ] <- l_region
    }
    i <- i + 1

    # Right end fragment
    right <- right + f1[2] - 1
    r_region <- region
    r_region[2] <- right
    r_region[3] <- f1[3]
    r_region[7] <- paste(region[7], "right", sep = "_")
    if (is.null(region_list)) {
      region_list <- r_region
    }
    else {
      region_list[i, ] <- r_region
    }
    i <- i + 1
  }
  return(region_list[, c(1, 2, 3, 7)])
}

################################################################################################
# Induced EP
################################################################################################
# 1. EcoRI-DpnII
x <- read.table("Induced_EP.enhancer_EcoRI_position.txt", header = F)
f_1 <- primary_fragment(x, motif_len = 6)
write.table(f_1, file = "Induced_EP.enhancer_EcoRI_pos.bed", row.names = F, col.names = F, quote = F)

x <- read.table("Induced_EP.enhancer_EcoRI_DpnII_position.txt", header = F)
# f_2 <- secondary_fragment(x, f_1, motif_len = 4)
f_2 <- secondary_fragment_both(x, f_1, motif_len = 4)
write.table(f_2, file = "Induced_EP.enhancer_EcoRI_DpnII_pos.bed", row.names = F, col.names = F, quote = F)



# 2. HindIII-DpnII
x <- read.table("Induced_EP.enhancer_HindIII_position.txt", header = F)
f_1 <- primary_fragment(x, motif_len = 6)
write.table(f_1, file = "Induced_EP.enhancer_HindIII_pos.bed", row.names = F, col.names = F, quote = F)

x <- read.table("Induced_EP.enhancer_HindIII_DpnII_position.txt", header = F)
# f_2 <- secondary_fragment(x, f_1, motif_len = 4)
f_2 <- secondary_fragment_both(x, f_1, motif_len = 4)
write.table(f_2, file = "Induced_EP.enhancer_HindIII_DpnII_pos.bed", row.names = F, col.names = F, quote = F)

################################################################################################
# Control EP
################################################################################################
# 1. EcoRI-DpnII
x <- read.table("Control_EP.enhancer_EcoRI_position.txt", header = F)
f_1 <- primary_fragment(x, motif_len = 6)
write.table(f_1, file = "Control_EP.enhancer_EcoRI_pos.bed", row.names = F, col.names = F, quote = F)

x <- read.table("Control_EP.enhancer_EcoRI_DpnII_position.txt", header = F)
# f_2 <- secondary_fragment(x, f_1, motif_len = 4)
f_2 <- secondary_fragment_both(x, f_1, motif_len = 4)
write.table(f_2, file = "Control_EP.enhancer_EcoRI_DpnII_pos.bed", row.names = F, col.names = F, quote = F)



# 2. HindIII-DpnII
x <- read.table("Control_EP.enhancer_HindIII_position.txt", header = F)
f_1 <- primary_fragment(x, motif_len = 6)
write.table(f_1, file = "Control_EP.enhancer_HindIII_pos.bed", row.names = F, col.names = F, quote = F)

x <- read.table("Control_EP.enhancer_HindIII_DpnII_position.txt", header = F)
# f_2 <- secondary_fragment(x, f_1, motif_len = 4)
f_2 <- secondary_fragment_both(x, f_1, motif_len = 4)
write.table(f_2, file = "Control_EP.enhancer_HindIII_DpnII_pos.bed", row.names = F, col.names = F, quote = F)
