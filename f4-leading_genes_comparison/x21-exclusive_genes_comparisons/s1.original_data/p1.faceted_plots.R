stop()
q('no')

library(DESeq2)
library(ggplot2)
library(plyr)
library(dplyr)

leading_genes_list <- list()
leading_genes_list[['NMF']] <- readRDS('leading_genes/leadingGenesNMF.RDS')
#leading_genes_list[['SemiNMF']] <- readRDS('leading_genes/leadingGenesSemiNMF.RDS')
#leading_genes_list[['B2013']] <- readRDS('leading_genes/leadingGenesB2013.RDS')
leading_genes_list[['SCDE']] <- readRDS('leading_genes/leadingGenesSCDE.RDS')
leading_genes_list[['MAST']] <- readRDS('leading_genes/leadingGenesMAST_B2013.RDS')
leading_genes_list[['Monocle']] <- readRDS('leading_genes/leadingGenesMonocle_B2013.RDS')
leading_genes_list[['DESeq2']] <- readRDS('leading_genes/leadingGenesDESeq2.RDS')
leading_genes_list[['EdgeR']] <- readRDS('leading_genes/leadingGenesEdgeR.RDS')


# exclusiveNMF <- setdiff(leadingGenesNMF, c(leadingGenesDESeq2, leadingGenesEdgeR, leadingGenesMonocle))
# exclusiveDESeq2 <- setdiff(leadingGenesDESeq2, c(leadingGenesNMF, leadingGenesEdgeR, leadingGenesMonocle))
# exclusiveEdgeR <- setdiff(leadingGenesEdgeR, c(leadingGenesDESeq2, leadingGenesNMF, leadingGenesMonocle))
# exclusiveMonocle <- setdiff(leadingGenesMonocle, c(leadingGenesDESeq2, leadingGenesEdgeR, leadingGenesNMF))

se <- readRDS('epitSE2.RDS')
phe <- colData(se)$phe
rwn <- assays(se)$normalized

##### ma-plot #####
a <- apply(rwn, 1, function(x)mean(log(x+1)))
m <- apply(rwn, 1, function(x){mean(log(x[phe=='E16.5']+1)) - mean(log(x[phe=='E14.5']+1))})

ma_ggdat <- data.frame(gene=rownames(se),
                       a=a,
                       m=m)

gen_data_frame <- function(method_name) ma_ggdat %>% mutate(is_exclusive=gene %in% leading_genes_list[[method_name]], method=method_name)

ggdat <- bind_rows(lapply(names(leading_genes_list), function(x)gen_data_frame(x))) %>%
  mutate(method=factor(method, levels=names(leading_genes_list))) %>%
  mutate(is_exclusive=factor(is_exclusive, labels=c('No', 'Yes')))

ggplot(ggdat) +
  geom_point(aes(x=a, y=m, color=is_exclusive, alpha=is_exclusive), ggdat %>% filter(is_exclusive=='No'), size=1.5) +
  geom_point(aes(x=a, y=m, color=is_exclusive, alpha=is_exclusive), ggdat %>% filter(is_exclusive=='Yes'), size=1.5) +
  geom_segment(aes(x=-0.05, y=0, xend=0.5-0.05, yend=(sum(phe=='E14.5')/sum(phe=='E16.5')+1)*0.5), linetype=3, size=0.5) +
  geom_segment(aes(x=-0.05, y=0, xend=0.95, yend=-sum(phe=='E16.5')/sum(phe=='E14.5')-1), linetype=3, size=0.5) +
  facet_wrap(~ method, ncol=2) +
  scale_color_manual(name='Exclusively\nFound', labels=c('No', 'Yes'), values=c('gray', 'blue')) +
  scale_alpha_manual(name='Exclusively\nFound', labels=c('No', 'Yes'), values=c(0.5, 0.5)) +
  labs(x='A', y='M') +
  ggtitle('MA-plot')

ggsave('../s2.plots/maplot.pdf', width=10, height=10)
ggsave('../s2.plots/maplot.png', width=10, height=10)

###################

##### tally num zeros #####
non_zero_dat <- data.frame(gene=rownames(se),
                       num_non_zeros_14=apply(rwn, 1, function(x)sum(x[phe=='E14.5']!=0)),
                       num_non_zeros_16=apply(rwn, 1, function(x)sum(x[phe=='E16.5']!=0)))

gen_data_frame <- function(method_name) non_zero_dat %>% filter(gene %in% leading_genes_list[[method_name]]) %>% mutate(method=method_name)

ggdat <- bind_rows(lapply(names(leading_genes_list), function(x)gen_data_frame(x))) %>%
  group_by(method, num_non_zeros_14, num_non_zeros_16) %>%
  summarise(count=length(gene)) %>%
  mutate(zero_group=(num_non_zeros_14==0|num_non_zeros_16==0))
ggdat$method <- factor(ggdat$method, levels=names(leading_genes_list))

ggplot(ggdat) +
  geom_point(aes(x=num_non_zeros_14, y=num_non_zeros_16, size=count, color=zero_group)) +
  facet_wrap(~ method, ncol=2) +
  labs(x='num nonzeros in E14.5', y='num nonzeros in E16.5') +
  scale_size_continuous(name='count') +
  scale_color_manual(name='one of the \ngroups has no \npositive reads', labels=c('No', 'Yes'), values=c('black', 'red')) +
  ggtitle('Comparison between num nonzeros in each group')

ggsave('../s2.plots/maplot.pdf', width=10, height=7)
########################

##### mean expression level density plot #####
ggdat <- names(leading_genes_list) %>% lapply(function(x)data.frame(method=x, a=a[leading_genes_list[[x]]])) %>% bind_rows
ggdat$method <- factor(ggdat$method, levels=names(leading_genes_list))

ggplot(ggdat) +
  geom_line(aes(x=a, color=method), stat='density') +
  scale_color_brewer(type='qual', palette=6) +
  labs(x='average log expression', y='density')

ggsave('../s2.plots/mean_expression_level_density_plot.pdf', width=8, height=5.5)

NMF_m <- ggdat %>%
  filter(method=='NMF')
NMF_m_density <- density(NMF_m$m)
exp(NMF_m_density$x[which.max(NMF_m_density$y)])
##############################################


##### compare # cases lowly expressed group has all zeros #####
# helper <- function(method, gene_list) {
#   data.frame(method=method, 
#              num_non_zeros_14=apply(rwNormalized[gene_list,phe=="E14.5"], 1, function(x)sum(x!=0)),
#              num_non_zeros_16=apply(rwNormalized[gene_list,phe=="E16.5"], 1, function(x)sum(x!=0))) %>%
#     mutate(min_non_zero=pmin(num_non_zeros_14, num_non_zeros_16))
# }
# 
# helper('Monocle', exclusiveMonocle)
# 
# ggdat <- rbind(helper('NMF', exclusiveNMF),
#                helper('DESeq2', exclusiveDESeq2),
#                helper('EdgeR', exclusiveEdgeR),
#                helper('Monocle', exclusiveMonocle)) %>%
#   group_by(method) %>%
#   summarize(all_zero_weak_side_cases=sum(min_non_zero==0))
# ggdat$method <- factor(ggdat$method, levels=c('NMF', 'DESeq2', 'EdgeR', 'Monocle'))
# 
# ggplot(ggdat) +
#   geom_bar(aes(x=method, y=all_zero_weak_side_cases), stat='identity') +
#   geom_text(aes(x=method, y=all_zero_weak_side_cases, label=all_zero_weak_side_cases), vjust=-0.3) +
# labs(x='method', y='#cases that lowly expressed group has all zeros')
###############################################################



##### cv plot #####
me <- apply(rwNormalized, 1, mean)
co <- apply(rwNormalized, 1, var)

ggdat <- data.frame(m=rep(me, 4), cv2=rep(co/me^2, 4), 
                    is.exclusive=factor(c(rownames(rwNormalized) %in% exclusiveNMF,
                                          rownames(rwNormalized) %in% exclusiveDESeq2,
                                          rownames(rwNormalized) %in% exclusiveEdgeR,
                                          rownames(rwNormalized) %in% exclusiveMonocle),
                                        levels=c(T, F)),
                    method=factor(rep(c('NMF', 'DESeq2', 'EdgeR', 'Monocle'), each=length(m)),
                                  levels=c('NMF', 'DESeq2', 'EdgeR', 'Monocle')))

ggplot(ggdat) +
  geom_point(aes(x=m, y=cv2, color=is.exclusive, alpha=is.exclusive), size=1.5, subset=.(is.exclusive==F)) +
  geom_point(aes(x=m, y=cv2, color=is.exclusive, alpha=is.exclusive), size=1.5, subset=.(is.exclusive==T)) +
  geom_segment(aes(x=5e-3, y=length(phe)*1.1, xend=5, yend=length(phe)*1.1), linetype=3, size=0.5) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ method) +
  scale_color_manual(name='Exclusively\nFound', labels=c('No', 'Yes'), values=c('gray', 'blue')) +
  scale_alpha_manual(name='Exclusively\nFound', labels=c('No', 'Yes'), values=c(0.5, 0.5)) +
  labs(x='Mean', y='CV2') +
  ggtitle('CV-plot')
ggsave('../s2.plots/cvplot.pdf', width=10, height=7)
###################
