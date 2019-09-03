library(RColorBrewer)
setwd("~/Documents/IFN_enhancer/R/")
source("Script/ref_2_symbol_matrix.R")
marker <- list(color = brewer.pal(6, "Paired"))
marker <- marker$color

time <- c("0h", "6h", "12h", "18h", "24h")
##################################################################################################
# 1. Meta-gene profiles
##################################################################################################
plot_profile <- function(k, factor, type, ylim, legend = F) {
  x <- read.table(paste(c("~/Documents/IFN_enhancer/R/Heatmap/", factor, "_near_", type, "_tss.profile"), collapse = ""), skip = 2)
  xx <- x[, 3:402]

  main_txt <- paste(c(factor, "Near", type, "TSS"), collapse = " ")
  ylab_txt <- paste(c("RPKM (", as.character(x[k, 2]), " - Stable)"), collapse = " ")
  plot(c(-199:200) * 10 + 5, xx[k, ] - xx[2, ],
    main = main_txt, xlab = "Distance to TSS", ylab = ylab_txt,
    ylim = ylim,
    col = marker[1], type = "l", lwd = 2
  )
  for (i in 1:4) {
    lines(c(-199:200) * 10 + 5, xx[i * 3 + k, ] - xx[i * 3 + 2, ], col = marker[i + 1], lwd = 2)
  }
  if (legend) {
    legend("topright", legend = x[c(0:4) * 3 + 1, 1], lwd = 2, col = marker[1:5])
  }
}

##################################################################################################
par(mfrow = c(2, 2), mar = c(4, 4, 4, 5))

plot_profile(k = 1, factor = "Med1f", type = "eRNA", ylim = c(-0.5, 4.5))
plot_profile(k = 1, factor = "SMC1F", type = "eRNA", ylim = c(-0.5, 4.5), legend = T)
plot_profile(k = 1, factor = "Med1f", type = "gene", ylim = c(-0.5, 2.5))
plot_profile(k = 1, factor = "SMC1F", type = "gene", ylim = c(-0.5, 2.5))

##################################################################################################
par(mfrow = c(2, 2), mar = c(4, 4, 4, 5))

plot_profile(k = 3, factor = "Med1f", type = "eRNA", ylim = c(-0.5, 4.5))
plot_profile(k = 3, factor = "SMC1F", type = "eRNA", ylim = c(-0.5, 4.5), legend = T)
plot_profile(k = 3, factor = "Med1f", type = "gene", ylim = c(-0.5, 2.5))
plot_profile(k = 3, factor = "SMC1F", type = "gene", ylim = c(-0.5, 2.5))


##################################################################################################
# 2. Correlation of ChIP-seq and GRO-seq
##################################################################################################
pmt <- read.table("Data/Early_up.gene_tss.bed")
pmt <- GRanges(seqnames = pmt[, 1], IRanges(pmt[, 2], pmt[, 3]), strand = pmt[, 6], id = pmt[, 4])
group_num <- c(780, 4980, 806) # Three groups of TSS: early_up, stable, early_down
n_col <- 400 # number of columns for each sample
n_sam <- 5 # number of samples
load("Data/gene_expression_GM_unfiltered.RData")

cor_analysis <- function(x, y, tag) {
  ##### Get the Average ChIP-seq signal and fold changes #####
  for (i in 1:n_sam) {
    start <- (i - 1) * n_col
    tp <- apply(x[, (1:n_col) + start], 1, mean)
    if (i == 1) {
      avg_rpkm <- tp
    }
    else {
      avg_rpkm <- cbind(avg_rpkm, tp)
    }
  }

  ##### Now focus on the "Early_up" group #####
  xx <- avg_rpkm[1:group_num[1], ]
  region <- y[1:group_num[1], ]
  fc <- xx[, 2:n_sam] / xx[, 1]
  lab <- order(rowMax(fc), decreasing = T) # Reorder by fold change
  tp <- cbind(rowMax(fc)[lab], region[lab, 1:3], xx[lab, ])
  xx <- xx[lab, ]
  region <- region[lab, ]

  ##### Find correlation of ChIPseq and expression #####

  tp <- GRanges(seqnames = region[, 1], IRanges(region[, 2], region[, 3]))
  hits <- findOverlaps(tp, pmt, type = "equal", ignore.strand = T, select = "first")
  id <- as.character(mcols(pmt)$id[hits])

  gro_rpkm <- gene.rpkm[id, paste("GM", time, sep = "")]
  xx_norm <- t(t(xx) - apply(avg_rpkm, 2, mean))

  cc_18h <- NULL
  for (i in 1:100) {
    cc_18h <- c(cc_18h, cor(gro_rpkm[i, -c(3, 4)], xx_norm[i, -c(3, 4)], method = "pearson"))
  }
  print(median(cc_18h))

  # Optional: plot some cases
  par(mfrow = c(4, 2), mar = c(4, 4, 4, 5))
  for (i in 1:8) {
    plot(gro_rpkm[i, -c(3, 4)], xx_norm[i, -c(3, 4)],
      col = marker, pch = 16, cex = 2,
      xlab = "GROseq RPKM", ylab = paste(tag, " Norm RPKM"),
      main = ref_2_symbol(id[i])
    )
  }
  return(cc_18h)
}

x <- read.table("Heatmap/Med1f_near_gene_tss.matrix", skip = 2)
y <- read.table("Heatmap/Med1f_near_gene_tss.region", skip = 1)
cc_med <- cor_analysis(x, y, tag = "Med1f")

x <- read.table("Heatmap/SMC1F_near_gene_tss.matrix", skip = 2)
y <- read.table("Heatmap/SMC1F_near_gene_tss.region", skip = 1)
cc_coh <- cor_analysis(x, y, tag = "SMC1a")
