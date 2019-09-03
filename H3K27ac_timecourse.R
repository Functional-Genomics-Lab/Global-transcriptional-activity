setwd("~/Documents/IFN_enhancer/R/Heatmap/")
library(RColorBrewer)

tag <- c("H3K27ac", "H3K27ac-2")
for (k in c(1:2)) {
  x <- read.table(paste(tag[k], "_eRNA_tss.profile", sep = ""), skip = 2)
  xx <- x[c(1:5) * 2 - 1, 3:ncol(x)]
  norm <- rowMeans(x[c(1:5) * 2, 3:ncol(x)])
  xx <- xx / norm

  pdf(paste(tag[k], "_eRNA_tss_normalized.profile.pdf", sep = ""))
  col <- brewer.pal(6, "Spectral")[-3]
  lwd <- 2
  plot(1, -1000,
    xlim = c(-200, 200) / 100, ylim = c(0.2, 2.4),
    xlab = "Distance from TSS (kb)", ylab = "Normalized RPKM"
  )
  for (i in 1:5) {
    lines(c(-199.5:199.5) / 100, xx[i, ], col = col[i], lwd = lwd)
  }
  legend("topright", legend = as.character(x[c(1:5) * 2, 1]), col = col, lty = "solid", lwd = lwd, bty = "n")
  dev.off()
}
