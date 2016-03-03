stop()
q('no')

library(dplyr)
library(DESeq2)
library(edgeR)
library(RDAVIDWebService)
library(org.Mm.eg.db)
library(monocle)


epitSE2 <- readRDS('in-epitSE2_B2013.RDS')
rwf <- assays(epitSE2)$fpkm
phe <- colData(epitSE2)$phe

##### monocle #####
pd <- new("AnnotatedDataFrame", data=as.data.frame(colData(epitSE2)))
DataFrame(pd@data)
fd <- new("AnnotatedDataFrame", data=data.frame(row.names=rownames(rwf)))
fd@data
mnd <- newCellDataSet(rwf, phenoData=pd, featureData=fd)
mnd <- detectGenes(mnd, min_expr = 0.1)
expressed_genes <- rownames(subset(fData(mnd), num_cells_expressed >= 30))
lmnd <- log(exprs(mnd))

redd <- differentialGeneTest(mnd, fullModelFormulaStr="expression~phenotype")

#==== RDS ====#
saveRDS(redd, 'redd.RDS')
#-------------#
saveRDS(redd, 'redd.RDS')
#=============#


leadingGenesMonocle <- rownames(redd)[redd$qval < 0.05]
saveRDS(leadingGenesMonocle, "leadingGenesMonocle_B2013.RDS")
# length = 129


