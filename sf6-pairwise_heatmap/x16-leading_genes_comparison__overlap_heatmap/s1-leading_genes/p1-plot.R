stop()
q('no')

library(org.Mm.eg.db)
library(foreach)
library(doMC)
library(topGO)
library(DESeq2)
library(limma)
library(igraph)
library(dplyr)

se <- readRDS('epitSE2.RDS')
rw <- assays(se)$fpkm
rwl <- zx::log_trans(rw)
phe <- colData(se)$phe

leading_genes <- list()
leading_genes[['NMF']] <- readRDS('leading_genes/leadingGenesNMF.RDS')
leading_genes[['DESeq2']] <- readRDS('leading_genes/leadingGenesDESeq2.RDS')
leading_genes[['EdgeR']] <- readRDS('leading_genes/leadingGenesEdgeR.RDS')
leading_genes[['Monocle']] <- readRDS('leading_genes/leadingGenesMonocle.RDS')
leading_genes[['B2013']] <- readRDS('leading_genes/leadingGenesB2013.RDS')
leading_genes[['SCDE']] <- readRDS('leading_genes/leadingGenesSCDE.RDS')
leading_genes[['SemiNMF']] <- readRDS('leading_genes/leadingGenesSemiNMF.RDS')
leading_genes[['MAST']] <- readRDS('leading_genes/leadingGenesMAST.RDS')

# intersection heatmap ----------------------------------------------------

library(reshape2)
library(ggplot2)
library(scales)

methods <- names(leading_genes)

self_cross <- function(v, f) {
  n <- length(v)
  mat <- matrix(nrow=n, ncol=n)
  for (i in 1:n) {
    for (j in 1:n) {
      mat[i, j] <- f(v[[i]], v[[j]])
    }
  }
  dimnames(mat) <- list(c1=names(v), c2=names(v))
  mat
}

overlap_table <- self_cross(leading_genes, function(x, y)length(intersect(x, y)))
methods_dendro <- as.dendrogram(hclust(as.dist(1/overlap_table)))
method_cluster_order <- methods[order.dendrogram(methods_dendro)]

ggdat <- overlap_table %>% melt %>%
  mutate(c1=factor(c1, levels=method_cluster_order)) %>%
  mutate(c2=factor(c2, levels=method_cluster_order))

ggplot() +
  geom_tile(aes(x=c1, y=c2, fill=value), ggdat) +
  geom_text(aes(x=c1, y=c2, label=value), ggdat, color='#55B1F7') + 
  labs(x='method', y='method') +
  scale_fill_continuous(name=wrap_format(10)('number of overlap genes'))
ggsave('../s2-plots/overlap_plot.pdf')

pdf('../s2-plots/dendrogram.pdf', height=3)
plot(methods_dendro)
dev.off()

#   -----------------------------------------------------------------------