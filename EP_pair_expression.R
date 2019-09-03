setwd("~/Documents/IFN_enhancer/R/")
library(GenomicRanges)
library(GenomicFeatures)
library(gplots)
library(rtracklayer)
library(cluster)

amplitude_index <- function(x) {
  fc <- log2((x[, -1] + 0.1) / (x[, 1] + 0.1))
  fc <- abs(fc)
  return(apply(fc, 1, max))
}
continuity_index <- function(x) {
  return(diag(cor(t(x[, -ncol(x)]), t(x[, -1]), method = "spearman")))
}

# EP paired heatmap
################################################################
# 1. Load data
################################################################

# 1.1 DE genes
# Criteria:
# 1. max(RPKM) > 1
# 2. Amplitude index > 1
# 3. Continuity index > 0.2
load(file = "Data/gene_expression_GM_unfiltered.RData")
rpkm <- gene.rpkm
rpkm <- rpkm[apply(rpkm, 1, max) > 1, ]
colnames(rpkm) <- substr(colnames(rpkm), 3, 1000L)
ai_rpkm <- amplitude_index(rpkm)
ci_rpkm <- continuity_index(rpkm)
ind_rpkm <- ai_rpkm * ci_rpkm

# Criteria: RPKM index > 0.5
deg <- names(ind_rpkm)[(ci_rpkm > 0.5 & ai_rpkm > 1)]
deg.data <- cbind(rpkm, ai_rpkm, ci_rpkm, ind_rpkm)[deg, ]

tp <- c(1:length(gene))
names(tp) <- as.character(gene$id)
deg.region <- gene[tp[deg]] # Select corresponding regions


# 1.2 Expressed enhancers
load(file = "Data/eRNA_expression_GM_unfiltered.RData")
rpkm <- eRNA.rpkm
rpkm <- rpkm[apply(rpkm, 1, max) > 0.5, ]
colnames(rpkm) <- substr(colnames(rpkm), 3, 1000L)
ai_rpkm <- amplitude_index(rpkm)
ci_rpkm <- continuity_index(rpkm)
ind_rpkm <- 1 * ci_rpkm

# Criteria: RPKM index > 0.1
enh <- names(ind_rpkm)[ind_rpkm > 0.5]
enh.data <- cbind(rpkm, ai_rpkm, ci_rpkm, ind_rpkm)[enh, ]

tp <- c(1:length(eRNA))
names(tp) <- as.character(eRNA$id)
enh.region <- eRNA[tp[enh]]

################################################################
# 2. Associate genes with enhancers
################################################################

deg.lab <- NULL
enh.lab <- NULL

r <- 100000 # Find enhancers within this distance
for (i in 1:length(deg)) {
  tp <- flank(deg.region[i], width = r, both = T)
  hits <- subjectHits(findOverlaps(tp, enh.region, ignore.strand = T))
  if (length(hits) == 0) {
    next
  }
  hits <- hits[order(ind_rpkm[hits], decreasing = T)[1]]
  deg.lab <- c(deg.lab, i)
  enh.lab <- c(enh.lab, hits)
}

deg.pr <- deg.region[deg.lab] # Paired regions
enh.pr <- enh.region[enh.lab]
deg.pd <- deg.data[deg.lab, ] # Paired data
enh.pd <- enh.data[enh.lab, ]


deg.pd[order(deg.pd[, "ai_rpkm"], decreasing = T)[1:30], ]
enh.pd[order(deg.pd[, "ai_rpkm"], decreasing = T)[1:30], ]


################################################################
# 3. Heatmap
################################################################

# Heatmap of induced genes identified
# Sort rows by kmeans clustering
heat.deg <- log2((deg.pd[, -c(1, 13:15)] + 0.1) / (deg.pd[, 1] + 0.1))
heat.enh <- log10((enh.pd[, -c(1, 13:15)] + 0.1) / (enh.pd[, 1] + 0.1))

# Data clustering
dis_matrix <- (1 - cor(t(heat.deg), method = "spearman")) / 2
set.seed(0)
centers <- 3
kmedoids <- pam(dis_matrix, k = centers, diss = T) # Method: Partitioning Around Medoids

k_cluster <- kmedoids$clustering
cluster_order <- c(1, 2, 3)

ind <- cluster_order[k_cluster] * 10^4
for (i in 1:centers) {
  lab <- which(k_cluster == i)
  # ind[lab] <- ind[lab] - dis_matrix[kmedoids$medoids[i], lab]
  ind[lab] <- ind[lab] - apply(heat.deg[lab, ], 1, max)
}
od <- order(ind, decreasing = F)
heat.deg <- heat.deg[od, ]
heat.enh <- heat.enh[od, ]
k_cluster <- k_cluster[od]

my_palette <- my_palette <- colorRampPalette(c("blue", "white", "red"), bias = 1)(n = 20)
heatmap.2(heat.deg,
  col = my_palette, Colv = F, Rowv = F, dendrogram = "none",
  keysize = 1.5, density.info = "density", labRow = F,
  key.ylab = "", key.xlab = "", key.title = "", trace = "none"
)

heatmap.2(heat.enh,
  col = my_palette, Colv = F, Rowv = F, dendrogram = "none",
  keysize = 1.5, density.info = "density", labRow = F,
  key.ylab = "", key.xlab = "", key.title = "", trace = "none"
)
