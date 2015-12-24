stop()
q('no')

library(DESeq2)
library(NMF)
library(dplyr)

glioSE5 <- readRDS('glioSE5.RDS')
rw <- assays(glioSE5)$fpkm
genes <- rownames(rw)


h <- list()
for (i in 1:30) {
  h[[i]] <- read.table(paste0('Y_glio/', i, '.csv'), sep=',', header=F) %>%
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

#==== RDS ====#
#saveRDS(cutree(hclust(dist(pairwise_matrix)), 5), 'glioClustersSemiNMF.RDS')
#-------------#
glioClustersSemiNMF <- readRDS('glioClustersSemiNMF.RDS')
#=============#

#   -----------------------------------------------------------------------


