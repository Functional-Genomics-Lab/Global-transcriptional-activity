# Not working

setwd("~/Documents/IFN_enhancer/R/")

load(file = "Data/deg_matrix_direction.RData")
lab <- which(substring(rownames(deg_matrix_direction), 1, 12) == "MetaEnhancer")
mat <- deg_matrix_direction[-lab, ]

# Convert gene symbol to Entrez ID
library(org.Hs.eg.db)
eg <- as.list(org.Hs.egALIAS2EG)
eg <- eg[!is.na(eg)]
eg_all <- NULL
for (i in 1:length(eg)) {
  for (j in eg[[i]]) {
    eg_all <- c(eg_all, j)
  }
}
eg_all <- unique(eg_all)

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

## Create table
tp <- mat[, 5]
tp[tp != 1] <- 0
tp <- as.factor(tp)
names(tp) <- eg_list
up <- tp[tp == 1]
names(up) <- eg_list[tp == 1]
g <- rep(0, length(eg_all))
bg <- factor(bg, levels = c(0, 1))
names(bg) <- eg_all
bg[names(up)] <- 1

sampleGOdata <- new("topGOdata",
  description = "Simple session", ontology = "BP",
  allGenes = bg, geneSel = up,
  nodeSize = 10, # Min. no. of genes annotated to a GO
  annot = annFUN.org, mapping = "org.Hs.eg.db"
)

## Run test
resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")

## Summarize data in table
allRes <- GenTable(sampleGOdata,
  classicFisher = resultFisher,
  classicKS = resultKS, elimKS = resultKS.elim,
  orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10
)
allRes
