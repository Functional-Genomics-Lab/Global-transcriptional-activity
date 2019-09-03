setwd("~/Documents/IFN_enhancer/R/")

bin_smooth <- function(x, bin, step = 1) {
  out <- NULL
  for (i in c(1:(length(x) + 1 - bin))) {
    out <- c(out, mean(x[i:(i + bin - 1)]))
  }
  return(out)
}

#########################################################################################################################

plt <- function(tag, xlim, ylim, bin) {
  fg <- read.csv(paste("Data/phastCons44way", tag, "all_Inducible_eRNA.profile", sep = "_"), sep = "\t", header = T)
  fg_1 <- read.csv(paste("Data/phastCons44way", tag, "Inducible_eRNA.profile", sep = "_"), sep = "\t", header = T)
  # fg_1 <- read.csv(paste("Data/phastCons44way", tag, "Early_enhancer_1.profile", sep="_"), sep="\t", header=T)

  bg <- read.csv(paste("Data/phastCons44way", tag, "Synteney_bg_random200kb_1.profile", sep = "_"), sep = "\t", header = T)

  plot(bin_smooth(fg[, 1], bin), bin_smooth(fg[, 2], bin),
    type = "l", col = "blue", lwd = 2, xlim = xlim, ylim = ylim,
    main = tag, xlab = "Distance to center (bp)", ylab = "PhastCons Score"
  )
  lines(bin_smooth(bg[, 1], bin), bin_smooth(bg[, 2], bin), col = "grey", lwd = 2)
  legend("topright", c("Inducible Enh", "Background"), col = c("blue", "grey"), lwd = 2)

  # lines(bin_smooth(fg_1[,1], bin), bin_smooth(fg_1[,2], bin), col="red", lwd=2)
  # legend("topright", c("EP Enh", "Inducible Enh", "Background"), col=c("red", "blue", "grey"), lwd=2)
}

pdf("Figure/PhastCons_of_Enhancer.pdf", width = 21)
par(mfrow = c(1, 3))
plt("primate", c(-5000, 5000), c(0.1, 0.18), bin = 20)
plt("mammal", c(-5000, 5000), c(0.06, 0.125), bin = 20)
plt("vertebrate", c(-5000, 5000), c(0.06, 0.125), bin = 20)
dev.off()
