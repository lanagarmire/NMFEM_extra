stop()
q('no')

library(DESeq2)

glioSE5 <- readRDS('glioSE5.RDS')
rw <- assays(glioSE5)$fpkm
rwl <- zx::log_trans(rw)
phe <- colData(glioSE5)$phe

library(mclust)

mclust_results <- Mclust(t(rwl), G=5)

#==== RDS ====#
#saveRDS(mclust_results, 'mclust_results.RDS')
#-------------#
#mclust_results <- readRDS('mclust_results.RDS')
#=============#

mclust_results
summary(mclust_results)
str(mclust_results)
cls <- apply(mclust_results$z, 1, which.max)
rand_measure_zx::rand_measure(cls, phe)
