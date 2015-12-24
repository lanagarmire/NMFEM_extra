stop()
q('no')

library(rbenchmark)
library(zx)
library(DESeq2)
library(ggplot2)
library(readr)
library(hash)
library(purrr)
library(stringr)
library(reshape2)
library(tidyr)
library(dplyr)

#==== read all samples ====#
srr_to_srx <- read_csv('Runinfo.csv') %>%
  transmute(srr=factor(Run), srx=factor(Experiment)) %>%
  sample_frac()

rw <- read.table('count_table.csv', sep=',', header=T, row.names=1)
rw <- rw %>% data.matrix

tb_all <- rw %>%
  melt(c('gene', 'srr'), value.name='counts') %>%
  tbl_df %>%
  mutate(srx=factor(srr, levels=srr_to_srx$srr, labels=srr_to_srx$srx))

tb_all_meta <- read_csv('Summary.csv') %>%
  transmute(srx=factor(`Experiment Accession`), phe=sub('^.*mRNA-seq_(.*?)[_;].*$','\\1', `Experiment Title`))
#==========================#


#==== pick only single-cell samples ====#
tb_meta <- tb_all_meta %>% 
  filter(phe %>% str_detect('^MGH..$'))

tb <- tb_all %>%
  semi_join(tb_meta)
#=======================================#


#==== prepare data for DESeq2 ====#
# remove lowly expressed samples so that DESeq2 can run
num_zeros_in_each_sample <- tb %>%
  group_by(srx) %>%
  summarize(num_zeros=sum(counts==0)) %>%
  arrange(num_zeros)

tb_filtered <- tb %>%
  semi_join(num_zeros_in_each_sample %>% head(446))  # this number is chosen by trial-and-error starting from num of all samples,
                                                     # keep reducing it until it passes the test
tb_filtered %>%
  group_by(gene) %>%
  summarize(all_non_zero=all(counts > 0)) %>%
  filter(all_non_zero)

tb <- tb_filtered %>%
  group_by(gene) %>%
  filter(any(counts>0)) %>%
  ungroup

tb_meta <- tb_meta %>%
  semi_join(tb_filtered)

tb_meta %>% summarize(n_distinct(srx))
#=================================#

#==== get rid of lowly expressed samples ====#
rw_temp <- tb %>% acast(gene ~ srx, value.var='counts')
col_sums <- colSums(rw_temp)

ggdat <- data.frame(col_sums = col_sums)

cutoff <- quantile(col_sums, 0.05)

ggplot() +
  geom_density(aes(x = col_sums), ggdat) +
  scale_x_log10() +
  geom_vline(aes(xintercept = cutoff), linetype = "longdash") +
  labs(x = "Total Counts", y = "Distribution")

rw_temp_f <- rw_temp[, col_sums >= cutoff]
tb <- rw_temp_f %>%
  melt(c('gene', 'srx'), value.name='value')

tb %>% print(n = 10)
#============================================#

#==== normalize counts using DESeq2 ====#
rw_temp <- tb %>% acast(gene ~ srx, value.var='counts')
col_data_temp <- data.frame(row.names=tb_meta$srx, phe=factor(tb_meta$phe))
se <- SummarizedExperiment(SimpleList(counts=rw_temp), 
                           colData=DataFrame(col_data_temp))
dds <- DESeqDataSet(se, ~ phe)
dds <- estimateSizeFactors(dds)
geneSizes <- read.table('~/rna/___resources/human_gene_sizes.txt', sep='\t', header=T, row.names=1)
mcols(dds)$basepairs <- geneSizes[rownames(dds),'gene_size']


tb_counts <- rw_temp %>%
  melt(c('gene', 'srx'), value.name='value')

tb_normalized <- counts(dds, normalized=T) %>%
  melt(c('gene', 'srx'), value.name='value')

tb_fpkm <- fpkm(dds, robust=T) %>%
  melt(c('gene', 'srx'), value.name='value')

tb <- bind_tbls(list(counts=tb_counts, normalized=tb_normalized, fpkm=tb_fpkm), 'value_type')

tb <- tb %>%
  spread(value_type, value)
#============================================#

#--------------------------------------------------------------------------

rw <- read.table('count_table.csv', sep=',', header=T, row.names=1)

srrTable <- read.table('Runinfo.csv', sep=',', header=T, row.names=1)
srxTable <- read.table('Summary.csv', sep=',', header=T, row.names=1)

colnames(rw) <- srrTable[colnames(rw),]$Experiment
# remove the samples with the most zeros until OK for DESeq2 (there exists a nonzero row).
any(apply(rw[,tail(order(apply(rw, 2, function(x)sum(x>0))),751)], 1, function(x)all(x>0)))
rw <- rw[,tail(order(apply(rw, 2, function(x)sum(x>0))),751)]
# remove zero rows
rw <- rw[apply(rw, 1, function(x)any(x>0)),]

glioSE <- SummarizedExperiment(SimpleList(counts=as.matrix(rw)), colData=DataFrame(srxTable[colnames(rw),]))

##### normalize counts #####
dds <- DESeqDataSet(glioSE, ~ Experiment.Title)
dds <- estimateSizeFactors(dds)
geneSizes <- read.table('~/ta.sk/___resources/human_gene_sizes.txt', sep='\t', header=T, row.names=1)
mcols(dds)$basepairs <- geneSizes[rownames(dds),'gene_size']
mediFpkm <- fpkm(dds, robust=T)
mediNormalized <- counts(dds, normalized=T)
assays(glioSE)$normalized <- mediNormalized
assays(glioSE)$fpkm <- mediFpkm
############################

mediNormalizedLogged <- log(mediNormalized+1)

phe <- colData(glioSE)$Experiment.Title
phe <- sub('^.*mRNA-seq_(.*?)[_;].*$','\\1',phe)
colData(glioSE) <- DataFrame(phenotype=phe, row.names=colnames(rw))


##### save RDS #####
saveRDS(glioSE, '../s2.SE_RDS_file/glioSE.RDS')
####################

##### draw CV plot #####
me <- apply(mediNormalized, 1, mean)
co <- apply(mediNormalized, 1, var)

ggplot() +
  geom_point(aes(x=me, y=co/me^2), alpha=0.2) +
  scale_x_log10() +
  scale_y_log10()

barplot(colSums(mediNormalized))
########################
##### PCA #####
pr <- prcomp(t(mediNormalizedLogged))

ggplot() +
  geom_point(aes(x=pr$x[,1],y=pr$x[,2],color=phe)) +
  scale_color_manual(values=1:8)
###############


##### select 2 phe to compare #####
glioSE3 <- glioSE[,(phe=="MGH26")|(phe=="MGH30")|(phe=="MGH31")]

#==== export comparison table ====#
saveRDS(glioSE3, '../s2.SE_RDS_file/glioSE3.RDS')
#=================================#

glioSE5 <- glioSE[,grep('^MGH..$', colData(glioSE)$phenotype)]

#==== export comparison table ====#
saveRDS(glioSE5, '../s2.SE_RDS_file/glioSE5.RDS')
#=================================#

#==== normalize counts ====#
dds2 <- DESeqDataSet(glioSE2, ~ Phenotype)
dds2 <- estimateSizeFactors(dds2)
#mediNormalized2 <- counts(dds2, normalized=T)
mediNormalized2 <- apply(assay(glioSE2), 1, function(x)x/sqrt(sum(x^2)))
#==========================#

#mediNormalizedLogged2 <- log(mediNormalized2+1)
mediLogged2 <- log(assay(glioSE2)+1)

mediLoggedNormalized2 <- apply(mediLogged2, 2, function(x)x/sqrt(sum(x^2)))

phe <- colData(glioSE2)$phenotype

#pr <- prcomp(t(mediNormalizedLogged2))
pr <- prcomp(t(mediLoggedNormalized2))

ggplot() +
  geom_point(aes(x=pr$x[,1],y=pr$x[,2],color=phe)) +
  scale_color_manual(values=1:8)

###################################
