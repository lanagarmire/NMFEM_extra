stop()
q('no')

library(DESeq2)
library(ggplot2)
library(zx)
library(hash)
library(combinat)
library(readr)
library(stringr)
library(reshape2)
library(tidyr)
library(mclust)
library(dplyr)
library(NMF)


#==== choose ====#
data_set_name <- 'glioSE5'
data_set_name <- 'boneSE'
data_set_name <- 'hscmSE'
#================#

se <- readRDS(paste0(data_set_name, '.RDS'))

phe <- factor(colData(se)$phe)
phes <- levels(phe)
num_types <- length(phes)

rw <- assays(se)$rpkm
rw <- rw[apply(rw, 1, function(x)any(x>0)),]

rwl <- log(rw+1)


clusterings <- list()

# NMF ---------------------------------------------------------------------

ren <- nmf(rwl, rank=num_types, nrun=30, method="brunet", .options="p32v3", seed=12345)

#==== saveRDS ====#
saveRDS(ren, paste0('RDS/ren_', data_set_name, '.RDS'))
#-----------------#
ren <- readRDS(paste0('RDS/ren_', data_set_name, '.RDS'))
#=================#

nmf_cons_mat <- consensus(ren)

clusterings[['NMF']] <- predict(ren) %>% factor %>% as.numeric
#   -----------------------------------------------------------------------

# SemiNMF -----------------------------------------------------------------
#glioClustersSemiNMF <- readRDS('glioClustersSemiNMF.RDS')
#clusterings[['SemiNMF']] <- glioClustersSemiNMF %>% factor %>% as.numeric
#   -----------------------------------------------------------------------

# kmeans ------------------------------------------------------------------

# 30 times average
kms <- list()
for (i in 1:30) {
  message(i)
  kms[[i]] <- kmeans(t(rwl), num_types)
}

#==== saveRDS ====#
saveRDS(kms, paste0('RDS/kms_', data_set_name, '.RDS'))
#-----------------#
kms <- readRDS(paste0('RDS/kms_', data_set_name, '.RDS'))
#=================#

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

km_cons <- cutree(hclust(dist(km_cons_mat)), num_types)


#==== saveRDS ====#
saveRDS(km_cons, paste0('RDS/km_cons_', data_set_name, '.RDS'))
#-----------------#
km_cons <- readRDS(paste0('RDS/km_cons_', data_set_name, '.RDS'))
#=================#

clusterings[['Kmeans']] <- km_cons %>% factor %>% as.numeric

#   -----------------------------------------------------------------------


#==== plot consensus map ====#
dendro_order <- function(mat) {order.dendrogram(as.dendrogram(hclust(dist(mat))))}

ggdat_list <- list()

ggdat_list[['NMF']] <- nmf_cons_mat[dendro_order(nmf_cons_mat), dendro_order(nmf_cons_mat)] %>%
  melt(c('s1', 's2'), value.name='consensus') %>%
  mutate(method='NMF')
  
ggdat_list[['Kmeans']] <- km_cons_mat[dendro_order(km_cons_mat), dendro_order(km_cons_mat)] %>%
  melt(c('s1', 's2'), value.name='consensus') %>%
  mutate(method='Kmeans')

ggdat <- bind_rows(ggdat_list) %>% mutate(method=factor(method, names(ggdat_list)))

theme_empty <- theme_bw()
theme_empty$line <- element_blank()
theme_empty$rect <- element_blank()
theme_empty$axis.text <- element_blank()
theme_empty$plot.title <- element_blank()
theme_empty$axis.title <- element_blank()
theme_empty$plot.margin <- structure(c(0, 0, -1, -1), unit = "lines", valid.unit = 3L, class = "unit")

ggplot() +
  geom_raster(aes(x=s1, y=s2, fill=consensus), ggdat) +
  facet_wrap(~ method) +
  theme_empty
ggsave('../s2-plots/consensus.pdf', width=7, height=3)

ggdat_color_bar_list <- list()
ggdat_color_bar_list[['NMF']] <- data.frame(method='NMF', phe=phe[dendro_order(nmf_cons_mat)]) %>% add_rownames('x') %>% mutate(x=as.double(x))
ggdat_color_bar_list[['Kmeans']] <- data.frame(method='Kmeans', phe=phe[dendro_order(km_cons_mat)]) %>% add_rownames('x') %>% mutate(x=as.double(x))
ggdat_color_bar <- bind_rows(ggdat_color_bar_list) %>% mutate(method=factor(method, names(ggdat_list)))

ggplot() +
  geom_raster(aes(x=x-1, y=1, fill=phe), ggdat_color_bar) +
  facet_wrap(~ method) + 
  theme_empty
ggsave('../s2-plots/consensus_color_bar.pdf', width=7, height=3)

#============================#


# Mclust ------------------------------------------------------------------

# this takes a long time and the generated object is HUGE
mclust_results <- Mclust(t(rwl), G=num_types)

# so I save this instead
mclust_zs <- mclust_results$z

#==== RDS ====#
saveRDS(mclust_zs, paste0('RDS/mclust_zs_', data_set_name, '.RDS'))
#-------------#
mclust_zs <- readRDS(paste0('RDS/mclust_zs_', data_set_name, '.RDS'))
#=============#

cls <- apply(mclust_zs, 1, which.max)


clusterings[['Mclust']] <- cls %>% factor %>% as.numeric
#   -----------------------------------------------------------------------


# Hclust ------------------------------------------------------------------
hc <- cutree(hclust(dist(t(rwl))), num_types)
clusterings[['Hclust']] <- hc %>% factor %>% as.numeric
#   -----------------------------------------------------------------------


# get clusterings ---------------------------------------------------------


phes <- levels(phe)
perms <- phes %>% length %>% permn

best_clusterings <- lapply(clusterings, function(c) {
  calculate_accuracy <- . %>% sum(factor(c, labels=phes[.])==phe)
  accuracies <- perms %>% sapply(calculate_accuracy)
  best_perm <- perms[[which.max(accuracies)]]
  factor(c, labels=phes[best_perm])
})

best_clusterings <- c(list(Labels=phe), best_clusterings)

#==== RDS ====#
saveRDS(best_clusterings, paste0('RDS/best_clusterings_', data_set_name, '.RDS'))
#-------------#
best_clusterings <- readRDS(paste0('RDS/best_clusterings_', data_set_name, '.RDS'))
#=============#

#   -----------------------------------------------------------------------


# get rand measures -------------------------------------------------------

rand_measures <- sapply(best_clusterings, function(c){
  rand_measure(c, phe)
})

#==== RDS ====#
saveRDS(rand_measures, paste0('RDS/rand_measures_', data_set_name, '.RDS'))
#-------------#
rand_measures <- readRDS(paste0('RDS/rand_measures_', data_set_name, '.RDS'))
#=============#

#   -----------------------------------------------------------------------


# get chisq ---------------------------------------------------------------

chisq_measures <- sapply(best_clusterings, function(c) {
  chisq.test(table(phe, c))$p.value
})

#   -----------------------------------------------------------------------



# barplot -----------------------------------------------------------------

library(scales)

ggdat <- data.frame(method=factor(names(rand_measures), names(sort(-rand_measures))), rand_measure=rand_measures, chisq_measure=chisq_measures)
ggdat <- ggdat[-1, ]

ggplot() +
  geom_bar(aes(x=method, y=rand_measure), ggdat, stat='identity') +
  scale_y_continuous(name='rand measure', label=percent_format()) +
  scale_x_discrete(name='method') +
  theme_gray(base_size=18)
ggsave('../s2-plots/bar_plot.pdf', width=5.5, height=4)

#   -----------------------------------------------------------------------


methods <- names(best_clusterings)

# scatter plot ------------------------------------------------------------

pr <- prcomp(t(rwl))

ggdat_list <- list()
for (i in 1:length(best_clusterings)) {
  method <- names(best_clusterings)[i]
  ggdat_list[[method]] <- data.frame(method=method, 
                                pc1=pr$x[,1], pc2=pr$x[,2], 
                                cluster=best_clusterings[[i]], 
                                correct=(best_clusterings[[i]]==phe))
}
ggdat <- bind_rows(ggdat_list) %>% mutate(method=factor(method, names(ggdat_list)))

ggdat %>% sample_n(20)

ggplot() + 
  geom_point(aes(x=pc1, y=pc2, color=cluster, shape=correct), ggdat, size=3) +
  facet_wrap(~ method) + 
  scale_color_brewer(type='qual', palette=6) + 
  scale_shape_manual(values=c(4, 16)) + 
  labs(x='PC1', y='PC2') +
  theme_gray(base_size=18)
ggsave('../s2-plots/scatter_plot.pdf', height=8, width=12)
  
#   -----------------------------------------------------------------------



# confusion matrix --------------------------------------------------------

ggdat_list <- list()
for (i in 1:length(best_clusterings)) {
  name <- names(best_clusterings)[i]
  ggdat_list[[name]] <- table(phe, best_clusterings[[name]]) %>%
    melt(c('label', 'method_results'), value.name='num_occurences') %>%
    mutate(label=factor(label, phes)) %>%
    mutate(method_results=factor(method_results, phes)) %>%
    mutate(method=names(best_clusterings)[i])
}

ggdat <- bind_rows(ggdat_list) %>%
  mutate(method=factor(method, methods))

ggplot() +
  geom_tile(aes(x=label, y=method_results, fill=num_occurences), ggdat) +
  geom_text(aes(x=label, y=method_results, label=num_occurences), ggdat, color='white') +
  facet_wrap(~ method) +
  labs(x='original label', y='results from methods') +
  scale_fill_continuous(name=wrap_format(10)('number of occurrences')) +
  theme_gray(base_size=18) + 
  theme(axis.text.x=element_text(hjust=1, vjust=1, angle=45))
ggsave('../s2-plots/confusion_matrices.pdf', height=7, width=10)

#   -----------------------------------------------------------------------








