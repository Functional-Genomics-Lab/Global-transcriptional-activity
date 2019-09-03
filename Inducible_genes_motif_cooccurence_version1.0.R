# Co-occurence of motifs in EP pairs

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
short_id <- function(xx) {
  short <- xx
  for (i in 1:length(xx)) {
    short[i] <- strsplit(xx[i], "_HUMAN")[[1]][1]
  }
  return(short)
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
ep_cooccur <- function(ee = ee, pp = pp, ep = ep, tf1, tf2) {
  ne <- NULL
  np <- NULL
  for (i in 1:nrow(ep)) {
    ne <- c(ne, sum((ee$PositionID == ep$id_b[i]) &
      (ee$Motif.Name == tf1)) > 0)
    np <- c(np, sum((pp$PositionID == ep$id_a[i]) &
      (pp$Motif.Name == tf2)) > 0)
  }
  ne <- as.integer(ne)
  np <- as.integer(np)
  return(cbind(ne, np))
}

####################################################################

motif_list <- get_expressed_tf()
ee <- ee[ee$Motif.Name %in% motif_list, ]
pp <- pp[pp$Motif.Name %in% motif_list, ]
ep <- read.csv("Data/Induced_EP.csv", colClasses = "character")

##################################################
# Co-occurnece of motifs in EP pairs
##################################################
motif_sel <- NULL # Only use motifs with enough occurence
short_id <- NULL
for (i in 1:length(motif_list)) {
  tp <- ep_cooccur(ee, pp, ep, motif_list[i], motif_list[i])
  if (colSums(tp)[1] > 20 & colSums(tp)[2] > 20) {
    motif_sel <- c(motif_sel, motif_list[i])
    short_id <- c(short_id, strsplit(motif_list[i], "_HUMAN")[[1]][1])
  }
}
ee_sel <- ee[ee$Motif.Name %in% motif_sel, ]
pp_sel <- pp[pp$Motif.Name %in% motif_sel, ]

n <- length(motif_sel)
cooccur_table <- matrix(ncol = n, nrow = n, dimnames = list(short_id, short_id))
pv_pos <- matrix(ncol = n, nrow = n, dimnames = list(short_id, short_id))
pv_neg <- matrix(ncol = n, nrow = n, dimnames = list(short_id, short_id))
for (i in 1:n) {
  for (j in 1:n) {
    tp <- ep_cooccur(ee_sel, pp_sel, ep, motif_sel[i], motif_sel[j])
    cooccur_table[i, j] <- sum(rowSums(tp) == 2) / nrow(tp)
    ta <- factor(tp[, 1], levels = c(0, 1))
    tb <- factor(tp[, 2], levels = c(0, 1))
    pv_pos[i, j] <- (-1) * log10(fisher.test(ta, tb, alternative = "greater")$p.value) # TF1 and TF2 are possitively associated
    pv_neg[i, j] <- (-1) * log10(fisher.test(ta, tb, alternative = "less")$p.value) # TF1 and TF2 are mutually exclusive
  }
}

heatmap.2(cooccur_table, symm = F, symbreaks = F, trace = "none")
heatmap.2(pv_pos, symm = F, symbreaks = F, trace = "none")
# heatmap.2(pv_neg, symm=F, symbreaks=F, trace="none")

pv_org <- 10^(-pv_pos)
pv <- sort(as.vector(pv_org))
pv_adj <- p.adjust(pv, method = "BH")
cutoff <- pv[sum(pv_adj < 0.05)] # Adjust P values by BH procedure and find significance cutoff
out_table <- NULL
for (i in 1:n) {
  for (j in 1:n) {
    if (pv_org[i, j] < cutoff) {
      out_table <- c(out_table, short_id[i], short_id[j]) # Fields: Enhancer_TF, Promoter_TF
      tp <- ep_cooccur(ee_sel, pp_sel, ep, motif_sel[i], motif_sel[j])
      out_table <- c(out_table, colSums(tp), sum(rowSums(tp) == 2)) # Fields: Enhancer_occurrence, Promoter_occurrence, Co_occurrence
      out_table <- c(
        out_table, 10^(-pv_pos[i, j]), # Fields: P_val
        pv_adj[which(pv == pv_org[i, j])[1]]
      ) # Fields: P_val_adj
    }
  }
}
out_table <- matrix(out_table, ncol = 7, byrow = T)
colnames(out_table) <- c("Enhancer_TF", "Promoter_TF", "Enhancer_occurrence", "Promoter_occurrence", "Co_occurrence", "P-val", "P-val.adj")
out_table <- as.data.frame(out_table)
out_table <- out_table[order(out_table$"P-val.adj"), ]
rownames(out_table) <- c(1:nrow(out_table))

# Control experiment: TF pairs that are exclusive
pv_org <- 10^(-pv_neg)
pv <- sort(as.vector(pv_org))
pv_adj <- p.adjust(pv, method = "BH")
cutoff <- pv[sum(pv_adj < 0.05)] # Adjust P values by BH procedure and find significance cutoff
out_table <- NULL
for (i in 1:n) {
  for (j in 1:n) {
    if (pv_org[i, j] < cutoff) {
      out_table <- c(out_table, short_id[i], short_id[j]) # Fields: Enhancer_TF, Promoter_TF
      tp <- ep_cooccur(ee_sel, pp_sel, ep, motif_sel[i], motif_sel[j])
      out_table <- c(out_table, colSums(tp), sum(rowSums(tp) == 2)) # Fields: Enhancer_occurrence, Promoter_occurrence, Co_occurrence
      out_table <- c(
        out_table, 10^(-pv_pos[i, j]), # Fields: P_val
        pv_adj[which(pv == pv_org[i, j])[1]]
      ) # Fields: P_val_adj
    }
  }
}
##################################################
# Occurnece table of enhancers
##################################################
ee_motif <- occurence_table(ee, motif_list)
# Motif similarity measured by Euclidian distance
ee_motif_dist <- motif_dist(ee_motif)
print(summary(ee_motif_dist[upper.tri(ee_motif_dist)]))
for (i in 1:(length(motif_list) - 1)) {
  for (j in (i + 1):length(motif_list)) {
    if (ee_motif_dist[i, j] < 0.1) {
      print(motif_list[c(i, j)])
    }
  }
}

# Motif co-occurence ratio
ee_motif_cooccur <- motif_cooccur(ee_motif)
print(summary(ee_motif_cooccur[upper.tri(ee_motif_cooccur)]))

for (i in 1:(length(motif_list) - 1)) {
  for (j in (i + 1):length(motif_list)) {
    if (ee_motif_cooccur[i, j] > 0.5) {
      print(c(motif_list[c(i, j)], ee_motif_dist[i, j], ee_motif_cooccur[i, j]))
    }
  }
}

##################################################
# Occurnece table of promoters
##################################################
pp_motif <- occurence_table(pp, motif_list)
# Motif similarity measured by Euclidian distance
pp_motif_dist <- motif_dist(pp_motif)
print(summary(pp_motif_dist[upper.tri(pp_motif_dist)]))
for (i in 1:(length(motif_list) - 1)) {
  for (j in (i + 1):length(motif_list)) {
    if (pp_motif_dist[i, j] < 0.1) {
      print(c(motif_list[c(i, j)], pp_motif_dist[i, j], pp_motif_cooccur[i, j]))
    }
  }
}

# Motif co-occurence ratio
pp_motif_cooccur <- motif_cooccur(pp_motif)
print(summary(pp_motif_cooccur[upper.tri(pp_motif_cooccur)]))

for (i in 1:(length(motif_list) - 1)) {
  for (j in (i + 1):length(motif_list)) {
    if (pp_motif_cooccur[i, j] > 0.5) {
      print(c(motif_list[c(i, j)], pp_motif_dist[i, j], pp_motif_cooccur[i, j]))
    }
  }
}
