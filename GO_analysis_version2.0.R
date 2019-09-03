setwd("~/Documents/IFN_enhancer/R/")

############################################################################################
# 1. DE genes
############################################################################################
load(file = "Data/deg_matrix_direction.RData")
lab <- which(substring(rownames(deg_matrix_direction), 1, 12) == "MetaEnhancer")
mat <- deg_matrix_direction[-lab, ]
# Convert gene symbol to Entrez ID
library(org.Hs.eg.db)
eg <- as.list(org.Hs.egALIAS2EG)
eg <- eg[!is.na(eg)]

rm_list <- NULL
eg_list <- NULL
for (i in 1:nrow(mat)) {
  tp <- which(names(eg) == rownames(mat)[i])
  if (length(tp) == 1) {
    eg_list <- c(eg_list, eg[[tp]][1])
  }
  else {
    rm_list <- c(rm_list, i)
  }
}
mat <- mat[-rm_list, ]

time_points <- colnames(mat)
for (i in 1:ncol(mat)) {
  tp <- mat[, i]
  up <- eg_list[tp == 1]
  dn <- eg_list[tp == -1]
  outfile <- paste("Data/GO_analysis/", "Up_", time_points[i], ".txt", sep = "")
  write.table(up, outfile, quote = F, row.names = F, col.names = F)
  outfile <- paste("Data/GO_analysis/", "Dn_", time_points[i], ".txt", sep = "")
  write.table(dn, outfile, quote = F, row.names = F, col.names = F)
}
write.table(eg_list, "Data/GO_analysis/GM_background.txt", row.names = F, col.names = F, quote = F)

############################################################################################
# 2. DE enhancers
############################################################################################

load(file = "Data/deg_matrix_direction.RData")
lab <- which(substring(rownames(deg_matrix_direction), 1, 12) == "MetaEnhancer")
mat <- deg_matrix_direction[lab, ]

time_points <- colnames(mat)
for (i in 1:ncol(mat)) {
  tp <- mat[, i]
  up <- rownames(mat)[tp == 1]
  dn <- rownames(mat)[tp == -1]
  outfile <- paste("Data/GO_analysis/", "Up_enh", time_points[i], ".txt", sep = "")
  write.table(up, outfile, quote = F, row.names = F, col.names = F)
  outfile <- paste("Data/GO_analysis/", "Dn_enh", time_points[i], ".txt", sep = "")
  write.table(dn, outfile, quote = F, row.names = F, col.names = F)
}
