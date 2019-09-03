# Find motifs that are highly enriched in interesting genomic regions, e.g. inducible enhancer

setwd("~/Documents/IFN_enhancer/R/")

####################################################################
enrichment <- function(x, y, region) {
  core <- which(x > region[1] & x < region[2])
  flank <- which(x < region[1] | x > region[2])
  pseudo <- 10^(-20)
  return(mean(y[core] + pseudo) / mean(y[flank] + pseudo))
}
continuity_index <- function(x) {
  return(diag(cor(t(x[, -ncol(x)]), t(x[, -1]), method = "spearman")))
}
#####################################################################

gene.exp <- read.csv("Data/Expression_matrix.symbol.csv", sep = "\t")

motif_analysis <- function(input, title) {
  n_factor <- (ncol(input) - 5) / 3
  x <- input[, 1]
  y <- input[c(1:n_factor) * 3 - 1]
  short_id <- colnames(y)
  for (i in 1:ncol(y)) {
    short_id[i] <- strsplit(short_id[i], "_HUMAN")[[1]][1]
  }
  print(paste("Total TF analyzed:", n_factor))

  # Match motif profile matrix with TF expression
  tf.list <- NULL
  tf.exp <- NULL # Shape of matrix: n_TF x n_timepoints
  tf.profile <- NULL # Shape of matrix: n_TF x n_bins
  for (i in 1:ncol(y)) {
    if (!(short_id[i] %in% rownames(gene.exp))) {
      next
    }
    lab <- which(rownames(gene.exp) == short_id[i])[1]
    if (max(gene.exp[lab, ]) < 1) {
      next
    } # Ignore TFs that were not expressed
    tf.list <- c(tf.list, short_id[i])
    if (length(tf.list) == 1) {
      tf.profile <- y[, i]
      tf.exp <- gene.exp[lab, ]
    }
    else {
      tf.profile <- rbind(tf.profile, y[, i])
      tf.exp <- rbind(tf.exp, gene.exp[lab, ])
    }
  }
  print(paste("TF expressed:", length(tf.list)))

  # 1) Fetures of motifs
  # Relative enrichment
  fold <- NULL
  region <- c(-100, 100)
  for (i in 1:nrow(tf.profile)) {
    fold <- c(fold, enrichment(x, tf.profile[i, ], region))
  }
  # Absolute enrichemnt
  peak <- apply(tf.profile, 1, max)

  # 2) Features of expression
  # Maximum fold increase
  pseudo <- 1e-1
  increase <- (apply(tf.exp[, -1], 1, max) + pseudo) / (tf.exp[, 1] + pseudo)
  # Continuity
  ci <- continuity_index(tf.exp)
  ai <- amplitude_index(tf.exp, abs = F)
  index_lab <- which((ci > 0.2) & (ai > log2(1.5)))
  print(paste("TF induced:", length(index_lab)))

  # Timing
  timing <- rep(0, length(tf.list))
  for (i in 1:length(tf.list)) {
    tp <- tf.exp[i, ]
    cutoff <- as.numeric(max(tp) / 2 + tp[1] / 2)
    timing[i] <- which(tp >= cutoff)[1]
  }

  col_scale <- timing / 12
  col_face <- rep(rgb(0.5, 0.5, 0.5, 0.0), length(tf.list))
  col_edge <- rep(rgb(0.5, 0.5, 0.5, 0.0), length(tf.list))
  col_face[index_lab] <- rgb(0, 1 - col_scale[index_lab], col_scale[index_lab], 0.5)
  col_edge[index_lab] <- rgb(0, 1 - col_scale[index_lab], col_scale[index_lab], 1.0)

  plot(peak, fold,
    bg = col_face, col = col_edge, cex = increase, pch = 21, xlim = c(0, 0.01), ylim = c(0, 4),
    xlab = "Absolute motif enrichment", ylab = "Relative motif enrichment", main = title
  )

  lab <- which(peak > 0.004 & fold > 2 & ci > 0.2 & ai > log2(2))
  if (length(lab) > 0) {
    text(x = cbind(peak[lab], fold[lab] - 0.2), tf.list[lab])
  }

  return(list(tf.list, tf.exp, tf.profile, peak, fold, ci, ai, increase))
}

ely <- read.table("Data/Motif_HOCO_Early_enhancer.profile", header = T, sep = "\t")
lat <- read.table("Data/Motif_HOCO_Late_enhancer.profile", header = T, sep = "\t")

pdf("Figure/Motif_enrichment_and_TF_expression.pdf")
res_ely <- motif_analysis(ely, "Early enhancers")
res_lat <- motif_analysis(lat, "Late enhancers")
# Plot figure legend
plot(-1, -1, xlim = c(0, 0.01), ylim = c(0, 5), xlab = "Absolute motif enrichment", ylab = "Relative motif enrichment")
text(0.005, 5, "TF induction time")
text(0.001, 4.7, "Early")
text(0.009, 4.7, "Late")
points(
  x = seq(0.002, 0.008, length.out = 12), y = rep(4.7, 12), pch = 21,
  bg = rgb(0, 1 - c(1:12) / 12, c(1:12) / 12, 0.5),
  col = rgb(0, 1 - c(1:12) / 12, c(1:12) / 12, 1), cex = 2
)
text(0.005, 4.4, "TF maximum increase")
text(0.001, 4.1, "1-Fold")
text(0.009, 4.1, "4-Fold")
points(
  x = seq(0.002, 0.008, length.out = 7), y = rep(4.1, 7),
  pch = 21, bg = "grey", col = "darkgrey", cex = seq(1, 4, length.out = 7)
)
dev.off()

x <- ely[, 1]
y <- ely[, c(0:639) * 3 + 2]
fold <- NULL
region <- c(-100, 100)
for (i in 1:ncol(y)) {
  fold <- c(fold, enrichment(x, y[, i], region))
}
short_id <- colnames(y)
for (i in 1:ncol(y)) {
  short_id[i] <- strsplit(short_id[i], "_HUMAN")[[1]][1]
}

pdf("Figure/Motif_inducible_enhancer.pdf")
par(mfrow = c(2, 2))
for (i in which(fold > 2)) {
  if (max(y[, i]) > 10^(-4)) {
    plot(x, y[, i], type = "l", main = short_id[i], xlab = "", ylab = "Sites per base per peak")
    # legend("topright", as.character(round(fold[i],2)))
  }
}
dev.off()


par(mfrow = c(3, 2))
for (i in which(short_id == "IRF1_HUMAN")) {
  if (max(y[, i]) > 10^(-3)) {
    plot(x, y[, i], type = "l", main = substr(colnames(y)[i], 1, 20), xlab = "", ylab = "Sites per base per peak")
    # legend("topright", as.character(round(fold[i],2)))
  }
}
