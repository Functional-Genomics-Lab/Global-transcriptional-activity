

setwd("~/Documents/IFN_enhancer/R/")
custom_bar <- function(mu, sd, names, ymax, col, ylab, title) {
  tp <- barplot(as.matrix(mu),
    beside = T, space = 0.7, ylim = c(0, ymax),
    col = col,
    ylab = ylab,
    names.arg = names,
    cex.names = 0.7,
    main = title
  )
  for (j in 1:length(tp)) {
    lines(x = rep(tp[j], 2), y = c(mu[j] - sd[j] / 2, mu[j] + sd[j] / 2))
    lines(x = c(tp[j] - bw, tp[j] + bw), y = rep(mu[j] - sd[j] / 2, 2))
    lines(x = c(tp[j] - bw, tp[j] + bw), y = rep(mu[j] + sd[j] / 2, 2))
  }
  abline(h = 1, col = "grey", lty = "dashed", lwd = 1.5)
}

bw <- 0.1 # 1/2 of error bar width

# 1. eRNA KD and gene expression
x <- read.table("Data/Figure7_1.txt", header = F)

pdf(file = "Figure/Figure7_1.pdf")
ymax <- c(4, 2.5, 2.5, 1)
title <- c("#5 KD", "#30 KD", "#38 KD", "Gene KD")
names <- c("eRNA#5", "eRNA#30", "eRNA#38", "TNFSF10")
col <- c("navy", "navy", "navy", "darkorange")
ylab <- "Expression FoldChange"

par(mfrow = c(2, 2))
for (i in 1:4) {
  mu <- x[, i * 2 - 1]
  sd <- x[, i * 2]
  custom_bar(mu = mu, sd = sd, names = names, ymax = ymax[i], col = col, ylab = ylab, title = title[i])
}
dev.off()

# 2. eRNA KD and EP interaction
x <- read.table("Data/Figure7_2.txt", header = F)

pdf(file = "Figure/Figure7_2.pdf")
ymax <- rep(1.25, 3)
title <- c("#5 KD", "#30 KD", "#38 KD")
names <- c("Neg_A", "Interaction_A", "Neg_B", "Interaction_B")
col <- "grey"
ylab <- "Interaction FoldChange"

par(mfrow = c(2, 2))
for (i in 1:3) {
  mu <- x[, i * 2 - 1]
  sd <- x[, i * 2]
  custom_bar(mu = mu, sd = sd, names = names, ymax = ymax[i], col = col, ylab = ylab, title = title[i])
}
dev.off()
# 3. eRNA KD and EE interaction

x <- read.table("Data/Figure7_3.txt", header = F)

pdf(file = "Figure/Figure7_3.pdf")
ymax <- rep(1.25, 1)
title <- c("#5 KD")
names <- c("Neg_C", "Interaction_C")
col <- "grey"
ylab <- "Interaction FoldChange"

par(mfrow = c(2, 2))
for (i in 1:1) {
  mu <- x[, i * 2 - 1]
  sd <- x[, i * 2]
  custom_bar(mu = mu, sd = sd, names = names, ymax = ymax[i], col = col, ylab = ylab, title = title[i])
}
dev.off()
