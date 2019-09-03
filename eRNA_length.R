setwd("~/Documents/IFN_enhancer/R/")
library(GenomicRanges)

load("Data/eRNA_expression_GM_unfiltered.RData")

ll <- width(eRNA)
print(median(ll))
pdf("Figure/eRNA_length.pdf")
hist(log10(ll), breaks = 40, xlim = c(2, 5), ylim = c(0, 2000), col = "grey", xlab = "Log10(eRNA length)", main = "")
dev.off()
