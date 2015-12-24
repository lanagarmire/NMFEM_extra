#stop()
#q('no')

library(org.Mm.eg.db)
library(foreach)
library(doMC)
library(topGO)
library(DESeq2)
library(limma)
library(igraph)
library(dplyr)

source('NMFEM.R')

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

for (method_index in 1:length(leading_genes)) {
  final_tb <- run_workflow(rwl, phe, leading_genes[[method_index]], 
                           '~/res/x01-mouse_ppi/s3-processed_tables/mppi.txt', 'org.Mm.eg.db',
                           verbose_level_ = 2, n_threads_ = 50)
  
  write.table(final_tb, paste0('../s2-tables/', names(leading_genes)[method_index], '_modules.txt'), sep='\t', row.names=F, col.names=T, quote=F)
}



# gen graphs --------------------------------------------------------------
#
#expr_matrix_ <- rwl
#grouping_vec_ <- phe
#selected_genes_ <- leading_genes[['NMF']]
#ppi_edges_file_ <- '~/res/x01-mouse_ppi/s3-processed_tables/mppi.txt'
#ppi_edges_file_header_ <- T
#verbose_level_ <- 2
#n_threads_ <- 50
#min_size_ <- 10
#max_size_ <- 100
#n_top_modules_ <- 5
#
#ppi_edges <- read.table(ppi_edges_file_, sep = '\t', header = ppi_edges_file_header_)
#
#ppi <- construct_weighted_graph_from_edge_list(ppi_edges, expr_matrix_, grouping_vec_)
#
#genes <- rownames(expr_matrix_)
#
#
#genes <- intersect(genes, V(ppi)$name)  # 14865 -> 8247
#ppi <- induced_subgraph(ppi, genes)
#ppi <- induced.subgraph(ppi, V(ppi)$name[clusters(ppi)$membership==1])  # 8247 -> 8174
#genes <- V(ppi)$name
#expr_matrix_ <- expr_matrix_[genes, ]
#
#selected_genes_ <- intersect(selected_genes_, genes)
#
#spinglass_results <- generate_spinglass_communities(ppi, selected_genes_, 
#                                                    min_size_ = min_size_, 
#                                                    max_size_ = max_size_, 
#                                                    n_threads_ = n_threads_, 
#                                                    verbose_level_ = verbose_level_)
#
#seed_p_values <- get_monte_carlo_simulation_p_values(spinglass_results$graphs, n_threads_ = n_threads_, verbose_level_ = verbose_level_)
#
#top_modules_tb <- spinglass_results$info %>%
#  inner_join(seed_p_values, 'seed') %>%
#  arrange(p_values) %>%
#  head(n_top_modules_)
#
#top_modules <- top_modules_tb$seed
#
#
#library(scales)
#library(ggplot2)
#
#rescale_t_stat <- function(t, half_range = 2) {
#  t[t < -half_range] <- -half_range
#  t[t > half_range] <- half_range
#  t[t < -half_range] <- -half_range
#  t[between(t, -half_range/2, half_range/2)] <- 0
#  t <- t/half_range
#  return(t/2+0.5)
#}
#
#t_stat_to_color <- . %>% rescale_t_stat %>% colour_ramp(c('green', 'gray',  'red'))(.)
#
#
#lmFit_results <- lmFit(expr_matrix_, model.matrix(~ 0 + grouping_vec_))
#contrast_fit_results <- contrasts.fit(lmFit_results, c(-1, 1))
#contrast_fit_results <- eBayes(contrast_fit_results)
#node_weights <- abs(contrast_fit_results$t[,1])
#
#
#graph_to_plot <- graph.union(spinglass_results$graphs[top_modules])
#V(graph_to_plot)$label.color <- ifelse(V(graph_to_plot)$name %in% top_modules, 'red', 'dodgerblue4')
#V(graph_to_plot)$label.cex <- ifelse(V(graph_to_plot)$name %in% top_modules, 1.3, 1)
#V(graph_to_plot)$color <- contrast_fit_results$t[,1][V(graph_to_plot)$name] %>% t_stat_to_color
#
#pdf(paste0('../s3-graphs/', top_modules[i],'.pdf'), width=20, height=20)
#plot(graph_to_plot,
#     vertex.size = 5,
#     vertex.frame.color = NA,
#     vertex.label.family = 'sans',
#     edge.width = 4,
#     edge.color = 'gray',
#     edge.arrow.mode = '-')
#
#library(shape)
#colorlegend(t_stat_to_color(seq(-2, 2, 0.5)), seq(-2, 2, 0.5), ratio.colbar = 0.3, 
#            xlim = c(-1.55, -1.4), ylim = c(0.5, 1), align = "r", cex = 1)
#text(-1.5, 0.4, c("t(mRNA)"), cex = 1)
#dev.off()
#
#   -----------------------------------------------------------------------
