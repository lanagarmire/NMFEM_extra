stop()
q('no')

library(NMF)
library(DESeq2)
library(ggplot2)
library(hash)
library(scales)
library(readr)
library(stringr)
library(reshape2)
library(tidyr)
library(dplyr)
library(zx)

rw <- read.table('RPKM_table_HSCMPP.csv', sep=',', header=T, row.names=1)
phe <- scan('phe.csv', what="", sep=',')

se <- SummarizedExperiment(SimpleList(fpkm=data.matrix(rw)), colData=DataFrame(row.names=colnames(rw), phe=phe))

se2 <- se[,grepl('1 cell', colData(se)$phe)]
colData(se2)$phe <- sub('(.*) 1 cell', '\\1', colData(se2)$phe)
se2 <- se2[apply(assay(se2), 1, function(x)any(x>0)),]

saveRDS(se2, '../s2-data_RDS/hscmSE2.RDS')

rwl <- log_trans(assay(se2))


# PCA ---------------------------------------------------------------------

rep <- prcomp(t(rwl))

phe <- colData(se2)$phe

ggdat <- data.frame(x=rep$x[,1], y=rep$x[,2])
ggplot() + 
  geom_point(aes(x=x, y=y, color=phe), ggdat)

#   -----------------------------------------------------------------------



# t-SNE -------------------------------------------------------------------

library(Rtsne)

ret <- Rtsne(t(rwl))

phe <- colData(se2)$phe

ggdat <- data.frame(x=ret$Y[,1], y=ret$Y[,2])
ggplot() + 
  geom_point(aes(x=x, y=y, color=phe), ggdat)

rwt <- t(ret$Y)
rwt <- rwt - min(rwt)
colnames(rwt) <- colnames(rwl)

#   -----------------------------------------------------------------------




# NMF ---------------------------------------------------------------------
ren <- nmf(rwl, rank=2, nrun=30, method="brunet", seed=12345, .options="p32v3")
#==== save RDS ====#
#saveRDS(ren, 'ren.RDS')
#------------------#
ren <- readRDS('ren.RDS')
#==================#
rand_measure_results_nmf <- data.frame(method='nmf', after_tsne=F, rand_measure=rand_measure(predict(ren), phe))

ren_t <- nmf(rwt, rank=2, nrun=30, method="brunet", seed=12345, .options="p32v3")
rand_measure_results_nmf_t <- data.frame(method='nmf', after_tsne=T, rand_measure=rand_measure(predict(ren_t), phe))
#   -----------------------------------------------------------------------

# kmeans ------------------------------------------------------------------
rek <- kmeans(t(rwl), 2)
rand_measure_results_kmeans <- data.frame(method='kmeans', after_tsne=F, rand_measure=rand_measure(unname(rek$cluster), phe))

rek_t <- kmeans(t(rwt), 2)
rand_measure_results_kmeans_t <- data.frame(method='kmeans', after_tsne=T, rand_measure=rand_measure(unname(rek_t$cluster), phe))
#--------------------------------------------------------------------------


# Hierarchical clustering -------------------------------------------------
reh <- hclust(dist(t(rwl)))
rand_measure_results_hclust <- data.frame(method='hclust', after_tsne=F, rand_measure=rand_measure(cutree(reh, 2), phe))

reh_t <- hclust(dist(t(rwt)))
rand_measure_results_hclust_t <- data.frame(method='hclust', after_tsne=T, rand_measure=rand_measure(cutree(reh_t, 2), phe))
#   -----------------------------------------------------------------------

ggdat <- bind_rows(list(rand_measure_results_nmf,
                        rand_measure_results_nmf_t,
                        rand_measure_results_kmeans,
                        rand_measure_results_kmeans_t,
                        rand_measure_results_hclust,
                        rand_measure_results_hclust_t)) %>% 
  mutate(method=factor(method, c('nmf', 'kmeans', 'hclust')),
         after_tsne=factor(after_tsne, c(F, T)))

ggplot() + 
  geom_bar(aes(x=method, y=rand_measure, fill=after_tsne), ggdat, stat='identity', position='dodge') + 
  scale_fill_discrete(name='after t-SNE') +
  scale_x_discrete(name='method', labels=c('NMF', 'K-means', 'Hierarchical Clustering')) +
  scale_y_continuous(name='Rand measure', labels=percent_format()) +
  theme_grey(base_size=18)



