stop()
q('no')

library(NMF)
library(DESeq2)
library(ggplot2)
library(zx)
library(hash)
library(scales)
library(readr)
library(stringr)
library(reshape2)
library(tidyr)
library(dplyr)
library(Rtsne)
library(NMF)

#==== choose ====#
data_set_name <- 'glioSE5'
num_clusters <- 5
#----------------#
data_set_name <- 'epitSE'
num_clusters <- 4
#----------------#
data_set_name <- 'hscmSE2'
num_clusters <- 2
#================#

se <- readRDS(paste0('data_sets/', data_set_name, '.RDS'))
sample_ids <- rownames(colData(se))
phes <- colData(se)$phe %>% unique
rw <- assays(se)$fpkm
rw <- rw[apply(rw, 1, function(x)any(x>0)),]
rwl <- zx::log_trans(rw)
phe <- colData(se)$phe

#==== quality control ====#
#sample_meta <- colData(se) %>% 
#  as.data.frame %>%
#  add_rownames('sample') %>%
#  mutate(sample=factor(sample)) %>%
#  mutate(phenotype=factor(phenotype))
#
#tb <- rwl %>%
#  melt(c('gene', 'sample'), value.name='log_fpkm') %>%
#  tbl_df
#
#ggplot() + 
#  geom_line(aes(x=log_fpkm, group=sample, color=phenotype), tb %>% filter(log_fpkm>0) %>% inner_join(sample_meta), stat='density') + 
#  scale_color_brewer(type='qual', palette=6)
#=========================#


# PCA ---------------------------------------------------------------------
#comp <- prcomp(t(rwl))
#   -----------------------------------------------------------------------



# t-SNE -------------------------------------------------------------------
ret <- Rtsne(t(rwl))

rwt <- t(ret$Y)
colnames(rwt) <- colnames(rwl)

#==== RDS ====#
#saveRDS(rwt, paste0('RDS/rwt_', data_set_name, '.RDS'))
#-------------#
rwt <- readRDS(paste0('RDS/rwt_', data_set_name, '.RDS'))
#=============#


#==== plot ====#
ggdat <- data.frame(x=ret$Y[,1], y=ret$Y[,2])

ggplot() +
  geom_point(aes(x=x, y=y, color=phe), ggdat)

ggplot() +
  geom_point(aes(x=x, y=y, color=predict(ren)), ggdat)
#==============#

#   -----------------------------------------------------------------------

rand_measure_results <- list()
rand_measure_results_t <- list()

# NMF run on tsne result --------------------------------------------------
ren <- nmf(rwl, rank=num_clusters, nrun=30, method="brunet", .options="p32v3", seed=12345)
#==== RDS ====#
#saveRDS(ren, paste0('RDS/ren_', data_set_name, '.RDS'))
#-------------#
ren <- readRDS(paste0('RDS/ren_', data_set_name, '.RDS'))
#=============#

rand_measure_results[['NMF']] <- rand_measure(predict(ren), phe)

ren_t <- nmf(rwt - min(rwt), rank=num_clusters, nrun=30, method="brunet", .options="p32v3", seed=12345)
#==== RDS ====#
#saveRDS(ren_t, paste0('RDS/ren_t_', data_set_name, '.RDS'))
#-------------#
ren_t <- readRDS(paste0('RDS/ren_t_', data_set_name, '.RDS'))
#=============#

rand_measure_results_t[['NMF']] <- rand_measure(predict(ren_t), phe)
#--------------------------------------------------------------------------



# tSNE + kmeans ------------------------------------------------------------------

#==== choose ====#
rwtt <- rw
rwtt <- rwt
#================#


# 30 times average
kms <- list()
for (i in 1:30) {
  message(i)
  kms[[i]] <- kmeans(t(rwtt), num_clusters)
}

km_clusters <- lapply(kms, function(x)x$cluster)

clustering_consensus <- function(clusterings, num_types = length(unique(clusterings[[1]]))) {
  n <- length(clusterings[[1]])
  pairwise_matrix <- matrix(0, nrow=n, ncol=n)
  for (c_idx in 1:length(clusterings)) {
    clusters <- clusterings[[c_idx]]
    for (i in 1:n) {
      for (j in 1:n) {
        pairwise_matrix[i, j] <- pairwise_matrix[i, j] + (clusters[i] == clusters[j])
      }
    }
  }
  pairwise_matrix / length(clusterings)
}

km_cons_mat <- clustering_consensus(km_clusters)

km_cons <- cutree(hclust(dist(km_cons_mat)), num_clusters)

clusterings[['tSNE + Kmeans']] <- km_cons %>% factor %>% as.numeric

#==== choose ====#
rand_measure_results[['Kmeans']] <- rand_measure(km_cons %>% factor %>% as.numeric, phe)
rand_measure_results_t[['Kmeans']] <- rand_measure(km_cons %>% factor %>% as.numeric, phe)
#================#

#   -----------------------------------------------------------------------


# Hierarchical clustering -------------------------------------------------
reh <- hclust(dist(t(rwl)))
rand_measure_results[['Hclust']] <- rand_measure(cutree(reh, num_clusters), phe)

reh_t <- hclust(dist(t(rwt)))
rand_measure_results_t[['Hclust']] <- rand_measure(cutree(reh_t, num_clusters), phe)

clusterings[['tSNE + Hclust']] <- as.numeric(factor(cutree(reh_t, num_clusters)))
#   -----------------------------------------------------------------------


# get clusterings ---------------------------------------------------------

library(combinat)

phes <- levels(factor(phe))
perms <- phes %>% length %>% permn

best_clusterings <- lapply(clusterings, function(c) {
  calculate_accuracy <- . %>% sum(factor(c, labels=phes[.])==phe)
  accuracies <- perms %>% sapply(calculate_accuracy)
  best_perm <- perms[[which.max(accuracies)]]
  factor(c, labels=phes[best_perm])
})

best_clusterings <- c(list(labels = factor(phe)), best_clusterings)
#   -----------------------------------------------------------------------


# scatter plot ------------------------------------------------------------

ggdat_list <- list()
for (i in 1:length(best_clusterings)) {
  method <- names(best_clusterings)[i]
  ggdat_list[[method]] <- data.frame(method=method, 
                                x=rwt[1,], y=rwt[2,], 
                                cluster=best_clusterings[[i]], 
                                correct=(best_clusterings[[i]]==phe))
}
ggdat <- bind_rows(ggdat_list) %>% mutate(method=factor(method, names(ggdat_list)))

ggplot() + 
  geom_point(aes(x=x, y=y, color=cluster, shape=correct), ggdat, size=3) +
  facet_wrap(~ method) + 
  scale_color_brewer(type='qual', palette=6) + 
  scale_shape_manual(values=c(4, 16)) + 
  labs(x='x', y='y') +
  theme_gray(base_size=18)
ggsave(paste0('../s3-plots/scatter_plot_', data_set_name, '.pdf'), width = 11, height = 4)
 
#   -----------------------------------------------------------------------


methods <- names(rand_measure_results)

ggdat <- rbind(data.frame(after_tsne=F, rand_measure=unlist(rand_measure_results)) %>% add_rownames('method'),
               data.frame(after_tsne=T, rand_measure=unlist(rand_measure_results_t)) %>% add_rownames('method')) %>%
  mutate(method=factor(method, methods))

#==== RDS ====#
saveRDS(ggdat, paste0('RDS/ggdat_', data_set_name, '.RDS'))
#-------------#
ggdat <- readRDS(paste0('RDS/ggdat_', data_set_name, '.RDS'))
#=============#

ggplot() + 
  geom_bar(aes(x=method, y=rand_measure, fill=after_tsne), ggdat, stat='identity', position='dodge') + 
  scale_fill_discrete(name='after t-SNE') +
  scale_x_discrete(name='method', labels=c('NMF', 'K-means', 'Hclust')) +
  scale_y_continuous(name='Rand measure', labels=percent_format()) +
  theme_grey(base_size=18)
ggsave(paste0('../s2-bar_plot/bar_plot_', data_set_name, '.pdf'), width=6.5, height=4)


# 