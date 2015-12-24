stop()
q('no')

library(DESeq2)
library(NMF)
library(dplyr)

epitSE2 <- readRDS('epitSE2.RDS')
rw <- assays(epitSE2)$fpkm
genes <- rownames(rw)

w <- list()
for (i in 1:30) {
  w[[i]] <- read.table(paste0('A/', i, '.csv'), sep=',', header=F) %>%
    data.matrix
}


# feature extraction ------------------------------------------------------

leadingGenesSemiNMF <- genes[tail(order(featureScore(w[[1]])), 500)]
saveRDS(leadingGenesSemiNMF, '../s2-leading_genes/leadingGenesSemiNMF.RDS')

#   -----------------------------------------------------------------------




h <- list()
for (i in 1:30) {
  h[[i]] <- read.table(paste0('Y/', i, '.csv'), sep=',', header=F) %>%
    data.matrix
}

n <- ncol(h[[1]])



# consensus clustering ----------------------------------------------------

pairwise_matrix <- matrix(0, nrow=n, ncol=n)
for (h_idx in 1:length(h)) {
  clusters <- apply(h[[h_idx]], 2, which.max)
  for (i in 1:n) {
    for (j in 1:n) {
      pairwise_matrix[i, j] <- pairwise_matrix[i, j] + (clusters[i] == clusters[j])
    }
  }
}

cutree(hclust(dist(pairwise_matrix)), 2)

#   -----------------------------------------------------------------------


