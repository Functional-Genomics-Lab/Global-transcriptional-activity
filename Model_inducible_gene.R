setwd("~/Documents/IFN_enhancer/R/")

x <- c(0:13)
y <- c(0, 0.5, 1, 2, 7, 9, 10, 9.5, 7, 4, 3, 2.5, 2.3, 2.2)

pdf("Figure/Model_inducible_gene.pdf")
plot(x[-1], y[-1], type = "l", lwd = 2, col = "grey", xaxt = "n", yaxt = "n", xlab = "Time", ylab = "Expression level", title = "Inducible gene")
points(x[-1], y[-1], pch = 16, col = rgb(0, 0, 1, ), cex = 2)

lines(x[-1], y[-length(x)], lty = "dashed", lwd = 2, col = "grey")
points(x[-1], y[-length(x)], pch = 1, col = rgb(0, 0, 1, ), cex = 2)
legend("topright", c("Real", "One-step delay"), pch = c(16, 1), col = rgb(0, 0, 1), bty = "n")
dev.off()
