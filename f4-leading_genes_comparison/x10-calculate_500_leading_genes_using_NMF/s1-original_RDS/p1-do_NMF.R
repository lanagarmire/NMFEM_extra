stop()
q('no')

library(DESeq2)
library(NMF)
library(zx)
library(ggplot2)


epitSE2 <- readRDS('epitSE2.RDS')
rw <- assays(epitSE2)$fpkm
phe <- colData(epitSE2)$phe
rwl <- log_trans(rw)

pr <- prcomp(t(rwl))
ggdat <- data.frame(x=pr$x[,1], y=pr$x[,2], phe=phe)
ggplot() +
  geom_point(aes(x=x, y=y, color=phe), ggdat)


#==== permutation test ====#
leadingGenesNMFp_list <- list()
for (i in 1:1000) {
  message(i)
  rwlp <- apply(rwl, 2, function(x)sample(x))
  rownames(rwlp) <- rownames(rwl)
  renp <- nmf(rwlp, rank=2, nrun=30, method="brunet", seed=12345, .options="p30v3")
  leadingGenesNMFp <- names(tail(sort(featureScore(renp)), 500))
  leadingGenesNMFp_list[[i]] <- leadingGenesNMFp
}
saveRDS(leadingGenesNMFp_list, 'leadingGenesNMFp_list.RDS')

p_values <- sapply(rownames(rwl), function(g)mean(sapply(leadingGenesNMFp_list, function(l)g %in% l)))
#==========================#

#==== jacknife test ====#
leadingGenesNMFj_list <- list()
for (i in 1:71) {
  message(i)
  rwlj <- rwl[, -i]
  rownames(rwlj) <- rownames(rwl)
  rwlj <- rwlj[apply(rwlj, 1, function(x)any(x > 0)), ]
  renj <- nmf(rwlj, rank=2, nrun=30, method="brunet", seed=12345, .options="p30v3")
  leadingGenesNMFj <- names(tail(sort(featureScore(renj)), 500))
  leadingGenesNMFj_list[[i]] <- leadingGenesNMFj
}
#==== RDS ====#
#saveRDS(leadingGenesNMFj_list, 'leadingGenesNMFj_list.RDS')
#-------------#
leadingGenesNMFj_list <- readRDS('leadingGenesNMFj_list.RDS')
#=============#

length(Reduce(intersect, leadingGenesNMFj_list))
#=======================#

ren <- nmf(rwl, rank=2, nrun=30, method="brunet", seed=12345, .options="p30v3")

leadingGenesNMF <- names(tail(sort(featureScore(ren)), 500))

#==== how many of these genes appear in the jackknife? ====#
p_values_j <- sapply(leadingGenesNMF, function(x)mean(sapply(leadingGenesNMFj_list, function(y)x %in% y)))
p_values_background <- sapply(rownames(rwl), function(x)mean(sapply(leadingGenesNMFj_list, function(y)x %in% y)))

ggdat <- data.frame(x = c(p_values_j, p_values_background), color = c(rep("Top 500", length(p_values_j)), rep("Background", length(p_values_background))))

ggplot() +
  geom_density(aes(x=x, fill=color), ggdat) +
  labs(x="frequency of the gene appear in the Jackknife", y="density") +
  theme_gray(base_size = 24)
ggsave('../s3-plots/jackknife.pdf')
#==========================================================#

p_values[leadingGenesNMF]

#==== RDS ====#
saveRDS(leadingGenesNMF, '../s2-leading_genes_RDS/leadingGenesNMF.RDS')
#-------------#
leadingGenesNMF <- readRDS('../s2-leading_genes_RDS/leadingGenesNMF.RDS')
#=============#



