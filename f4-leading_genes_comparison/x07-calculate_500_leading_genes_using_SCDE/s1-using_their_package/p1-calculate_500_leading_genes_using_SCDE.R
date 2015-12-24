stop()
q('no')

library(scde)
library(DESeq2)

N_CORES <- 10

epitSE2 <- readRDS('epitSE2.RDS')

rw <- assay(epitSE2)
phe <- colData(epitSE2)$phe

fpkm_matrix <- assays(epitSE2)$fpkm

o.ifm <- scde.error.models(counts=rw, groups=phe, n.cores=N_CORES, threshold.segmentation=T,
                           save.crossfit.plots=F, save.model.plots=F, verbose=1)

o.prior <- scde.expression.prior(models=o.ifm, counts=rw, length.out=400, show.plot=F)

ediff <- scde.expression.difference(o.ifm, rw, o.prior, groups=factor(phe), n.randomizations=100, n.cores=N_CORES, verbose=1)

saveRDS(rownames(head(ediff[order(-abs(ediff$Z)),], 500)), '../s2-leading_genes/leadingGenesSCDE.RDS')
