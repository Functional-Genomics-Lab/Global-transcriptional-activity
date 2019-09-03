# Co-occurence of motifs in EP pairs
library(gplots)
setwd("~/Documents/IFN_enhancer/R/")
ee <- read.table("Data/Inducible_eRNA.motif.positions", colClasses = "character", header = T, sep = "\t")
pp <- read.table("Data/Inducible_genes.motif.positions", colClasses = "character", header = T, sep = "\t")

####################################################################

get_expressed_tf <- function() {
  gene.exp <- read.csv("Data/Expression_matrix.symbol.csv", sep = "\t")
  tf_list <- unique(c(ee$Motif.Name, pp$Motif.Name))
  active_tf <- NULL
  for (i in 1:length(tf_list)) {
    tp <- strsplit(tf_list[i], "_HUMAN")[[1]][1]
    if (tp %in% active_tf) {
      next
    } # Avoid duplicates
    if (tp %in% rownames(gene.exp)) {
      if (max(gene.exp[tp, ]) > 1) {
        active_tf <- c(active_tf, tf_list[i])
      }
    }
  }
  return(active_tf)
}
occurence_table <- function(xx, motif_list) {
  region_list <- unique(xx$PositionID)
  zeros <- rep(0, length(region_list))
  names(zeros) <- region_list
  output <- NULL
  for (i in 1:length(motif_list)) {
    temp <- zeros
    tp <- xx[xx$Motif.Name == motif_list[i], 1]
    tp <- unique(tp)
    temp[tp] <- 1
    output <- c(output, temp)
  }
  output <- matrix(output, nrow = length(motif_list), byrow = T, dimnames = list(motif_list, region_list))
  return(output)
}
motif_dist <- function(xx) {
  motif_list <- rownames(xx)
  output <- matrix(rep(0, length(motif_list)^2), nrow = length(motif_list), dimnames = list(motif_list, motif_list))
  for (i in 1:(length(motif_list) - 1)) {
    for (j in (i + 1):length(motif_list)) {
      output[i, j] <- sum((xx[i, ] - xx[j, ])^2) / ncol(xx)
      output[j, i] <- output[i, j]
    }
  }
  return(output)
}
motif_cooccur <- function(xx) {
  motif_list <- rownames(xx)
  output <- matrix(rep(1, length(motif_list)^2), nrow = length(motif_list), dimnames = list(motif_list, motif_list))
  for (i in 1:(length(motif_list) - 1)) {
    for (j in (i + 1):length(motif_list)) {
      output[i, j] <- sum((xx[i, ] + xx[j, ]) == 2) / ncol(xx)
      output[j, i] <- output[i, j]
    }
  }
  return(output)
}
####################################################################

motif_list <- get_expressed_tf() # Only use motifs of expressed TFs
ee <- ee[ee$Motif.Name %in% motif_list, ]
pp <- pp[pp$Motif.Name %in% motif_list, ]

##################################################
# Occurnece table of enhancers
##################################################

motif_sel <- NULL # Only use motifs with enough occurence
short_id <- NULL
ee_motif <- occurence_table(ee, motif_list)
for (i in 1:length(motif_list)) {
  if ((sum(ee_motif[i, ]) / ncol(ee_motif)) > 0.2) {
    motif_sel <- c(motif_sel, motif_list[i])
    short_id <- c(short_id, strsplit(motif_list[i], "_HUMAN")[[1]][1])
  }
}

ee_motif <- occurence_table(ee, motif_sel)
pp_motif <- occurence_table(pp, motif_sel)

n <- length(short_id)
pv_ee <- matrix(1, nrow = n, ncol = n, dimnames = list(short_id, short_id))
pv_pp <- matrix(1, nrow = n, ncol = n, dimnames = list(short_id, short_id))
for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    pv_ee[i, j] <- fisher.test(ee_motif[i, ], ee_motif[j, ])$p.value
    pv_ee[j, i] <- pv_ee[i, j]
    pv_pp[i, j] <- fisher.test(pp_motif[i, ], pp_motif[j, ])$p.value
    pv_pp[j, i] <- pv_pp[i, j]
  }
}

heat.data <- (-1) * log10(pv_ee)
heatmap.2(heat.data, symm = TRUE, trace = "none")

heat.data <- (-1) * log10(pv_pp)
heatmap.2(heat.data, symm = TRUE, trace = "none")
