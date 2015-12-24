stop()
q('no')

library(DESeq2)
library(reshape2)
library(dplyr)
library(zx)

epitSE2 <- readRDS('epitSE2_B2013.RDS')
rw <- assays(epitSE2)$fpkm
rwl <- log_trans(rw)
phe <- colData(epitSE2)$phe
phe_tb <- data.frame(sample = colnames(rwl), phe = phe)

tbl <-rwl %>%
  melt(varnames = c('gene', 'sample')) %>%
  inner_join(phe_tb, 'sample')

tbl %>% head

# MAST --------------------------------------------------------------------

library(MAST)

sca <- SingleCellAssay(tbl, 'sample', 'gene', 'value', 'epit', phenovars = 'phe')
scaf <- MAST::filter(sca)
re_zlm <- zlm.SingleCellAssay(~ phe, scaf, method = 'glm', ebayes = T)

relr <- lrTest(re_zlm, 'phe')

p_values <- relr[,'hurdle', 'Pr(>Chisq)']

leadingGenesMAST <- names(p_values)[p_values < 0.05]
length(leadingGenesMAST)

saveRDS(leadingGenesMAST, '../s2-leading_genes/leadingGenesMAST_B2013.RDS')

#   -----------------------------------------------------------------------
