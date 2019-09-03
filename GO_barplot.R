setwd("~/Documents/IFN_enhancer/R/Data/GO_analysis/HOMER/Inducible_eRNA/")

pdf("~/Documents/IFN_enhancer/R/Figure/Inducible_eRNA_HOMER_GO_BP_barplot.pdf")
x <- read.table("Inducible_eRNA_HOMER_GO_BP_for_plot.txt", sep = "\t", header = T)
bp <- barplot(x[, 9], horiz = T, xlim = c(0, 20), col = "darkorange", space = 0.4)
for (i in 1:10) {
  text(5, bp[i], x[i, 1])
}
dev.off()
