#' Contruct weighted graph using node weights
#'
#' @param edge_list_ A \code{data.frame}. First two columns will be parsed as a edge list.
construct_weighted_graph_from_edge_list <- 
  function(edge_list_, expr_matrix_, grouping_vec_) {
    library(dplyr)
    library(limma)
    
    lmFit_results <- lmFit(expr_matrix_, model.matrix(~ 0 + grouping_vec_))
    contrast_fit_results <- contrasts.fit(lmFit_results, c(-1, 1))
    contrast_fit_results <- eBayes(contrast_fit_results)
    node_weights <- abs(contrast_fit_results$t[,1])
    
    
    v1_tb <- data.frame(v1 = names(node_weights),
                        v1_weight = node_weights)
    
    v2_tb <- data.frame(v2 = names(node_weights),
                        v2_weight = node_weights)
    
    edge_list <- edge_list_[, 1:2]
    colnames(edge_list) <- c('v1', 'v2')
    
    # filter edge_list so that genes not defined in the node_weights are discarded
    # also calculate the edge weights
    edge_list <- edge_list %>% 
      inner_join(v1_tb, 'v1') %>% 
      inner_join(v2_tb, 'v2') %>% 
      mutate(weight = (v1_weight + v2_weight)/2) %>% 
      select(v1, v2, weight)
    
    
    v_tb <- data.frame(v = names(node_weights),
                       weight = node_weights,
                       t_stat = contrast_fit_results$t[,1])
    
    weighted_graph <- graph_from_data_frame(edge_list, vertices = v_tb)
    
    return(weighted_graph)
  }




#' Generate communities using Spin Glass algorithm
generate_spinglass_communities <-
  function(weighted_graph_, seed_genes_, 
           n_threads_ = 1, gamma_ = 0.5, min_size_ = 1, max_size_ = 100, verbose_level_ = 1) {
    library(igraph)
    library(foreach)
    library(doMC)
    registerDoMC(n_threads_)
    
    if (verbose_level_ >= 1) message('* Running Spinglass algorithm ...')
    if (verbose_level_ >= 2) message('num of threads for parallel computing = ', n_threads_)
    if (verbose_level_ >= 2) message('num of seed genes = ', length(seed_genes_))
    
    modules <- foreach(i = 1:length(seed_genes_)) %dopar% {
      module <- cluster_spinglass(weighted_graph_, 
                                  weights = NULL,
                                  spins = 25, 
                                  start.temp = 1,
                                  stop.temp = 0.1,
                                  cool.fact = 0.99,
                                  update.rule = c("config"),
                                  gamma = gamma_,
                                  vertex = seed_genes_[i])
      module$seed <- seed_genes_[i]
      module
    }
    
    if (verbose_level_ >= 1) message('  - Generating module reports ...')
    
    modules <- modules[sapply(modules, function(x) between(length(x$community), min_size_, max_size_))]
    if (verbose_level_ >= 2) message('num of modules = ', length(modules))
    
    results <- list()
    
    genelists <- list()
    graphs <- list()
    info <- list()
    for (i in 1:length(modules)) {
      seed <- modules[[i]]$seed
      graph <- modules[[i]]$community %>% 
        induced_subgraph(weighted_graph_, .)
      vertex_indices <- modules[[i]]$community
      
      graphs[[seed]] <- graph
      
      info$seed[i] <- seed
      info$size[i] <- vertex_indices %>% length
      info$connectivity[i] <- graph %>% 
        degree_distribution() %>% 
        sum(. * 1:length(.))
      
      genelists[[seed]] <- V(weighted_graph_)$name[vertex_indices]
    }
    info <- data.frame(info)
    
    results$info <- info
    results$genelists <- genelists
    results$graphs <- graphs
    
    if (verbose_level_ >= 2) message('num of modules = ', nrow(results$info))
    if (verbose_level_ >= 2) message('module seeds = ', paste(results$info$seed, collapse = ', '))
    
    
    if (verbose_level_ >= 1) message('  - Done.')
    
    if (verbose_level_ >= 1) message('* Done.')
    
    return(results)
  }



#' Do Monte Carlo simulation
#' 
#' @param gene_weights_ A numeric vector. Should be the same one the user supplied to \code{construct_weighted_graph}.
get_monte_carlo_simulation_p_values <-
  function(module_graphs_, 
           n_monte_carlo_ = 1000, n_threads_ = 1, verbose_level_ = 1) {
    
    library(doMC)
    library(foreach)
    registerDoMC(n_threads_)
    
    if (verbose_level_ >= 1) message('* Doing Monte Carlo simulations to calculate p-values ... ')
    
    n_modules <- length(module_graphs_)
    
    if (verbose_level_ >= 2) message('num of modules = ', n_modules)
    if (verbose_level_ >= 2) message('num of monte carlo runs = ', n_monte_carlo_)
    if (verbose_level_ >= 2) message('num of threads for parallel computing = ', n_threads_)
    
    mean_weights <- vector(length = n_modules)
    for (i in 1:n_modules) {
      mean_weights[i] <- mean(E(module_graphs_[[i]])$weight)
    }
    
    monte_carlo_result_matrix <- matrix(nrow = n_modules, ncol = n_monte_carlo_)
    for (modules_index in 1:n_modules) {
      edgelist <- as_edgelist(module_graphs_[[modules_index]], name = F)
      vertex_weights <- V(module_graphs_[[modules_index]])$weight
      mean_list <- foreach(run_index = 1:n_monte_carlo_) %dopar% {
        vertex_weights <- sample(vertex_weights)
        mean(vertex_weights[edgelist])
      }
      monte_carlo_result_matrix[modules_index, ] <- unlist(mean_list)
    }
    
    p_values <- apply(monte_carlo_result_matrix > mean_weights, 1, mean)
    
    results <- list()
    results$seed <- names(module_graphs_)
    results$p_values <- p_values
    results <- data.frame(results)
    
    if (verbose_level_ >= 1) message('* Done. ')
    
    return(results)
  }


topGO_analysis <-
  function(top_module_genelists_, all_genes_, eg_db_, verbose_level_ = 1) {
    library(topGO)
    library(eg_db_, character.only = T)
    
    if (verbose_level_ >= 1) message('* Doing GO enrichment analysis ... ')
    if (verbose_level_ >= 2) message('num of modules = ', length(top_module_genelists_))
    if (verbose_level_ >= 2) message('module seeds = ', paste(names(top_module_genelists_), collapse = ', '))
    
    go_analysis_df <- list()
    for (i in 1:length(top_module_genelists_)) {
      seed <- names(top_module_genelists_)[i]
      
      if (verbose_level_ >= 1) message('  - Doing module with seed gene: ', seed, ' (', i, '/', length(top_module_genelists_), ') ...')
      
      go_analysis_df$seed[i] <- seed
      genelist <- top_module_genelists_[[i]]
      
      temp_genelist <- top_module_genelists_[[i]]
      interested_genes <- as.integer(all_genes_ %in% temp_genelist)
      names(interested_genes) <- all_genes_
      
      sink("/dev/null")
      GOdata <- invisible(
        new("topGOdata",
            ontology = "BP",
            allGenes = interested_genes,
            geneSel = function(x){x==1},
            nodeSize = 5,
            annot = annFUN.org,
            mapping = eg_db_, 
            ID = 'symbol'))
      
      
      resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
      gen_table <- GenTable(GOdata, classicFisher = resultFisher, numChar=999)
      
      sink()
      
      go_analysis_df$first_term[i] <- gen_table$Term[1]
      go_analysis_df$first_fisher[i] <- gen_table$classicFisher[1]
      go_analysis_df$second_term[i] <- gen_table$Term[2]
      go_analysis_df$second_fisher[i] <- gen_table$classicFisher[2]
      
      if (verbose_level_ >= 1) message('  - Done.')
    }
    go_analysis_df <- data.frame(go_analysis_df)
    
    if (verbose_level_ >= 1) message('* Done.')
    
    return(go_analysis_df)
  }





#' Run the NMFEM workflow.
#' 
#' This function depends on all other functions in this package.
#' 
#' @param output_file_ A string. Path and filename of the output table.
run_workflow <-
  function(expr_matrix_, grouping_vec_, selected_genes_, ppi_edges_file_, eg_db_,
           ppi_edges_file_header_ = T, verbose_level_ = 1, n_top_modules_ = 5, 
           n_threads_ = 1, min_size_ = 10, max_size_ = 100) {
    message('---------------------------- NMFEM ----------------------------')
    
    ppi_edges <- read.table(ppi_edges_file_, sep = '\t', header = ppi_edges_file_header_)
    
    ppi <- construct_weighted_graph_from_edge_list(ppi_edges, expr_matrix_, grouping_vec_)
    
    genes <- rownames(expr_matrix_)
    
    
    genes <- intersect(genes, V(ppi)$name)  # 14865 -> 8247
    ppi <- induced_subgraph(ppi, genes)
    ppi <- induced.subgraph(ppi, V(ppi)$name[clusters(ppi)$membership==1])  # 8247 -> 8174
    genes <- V(ppi)$name
    expr_matrix_ <- expr_matrix_[genes, ]
    
    selected_genes_ <- intersect(selected_genes_, genes)
    
    spinglass_results <- generate_spinglass_communities(ppi, selected_genes_, 
                                                        min_size_ = min_size_, 
                                                        max_size_ = max_size_, 
                                                        n_threads_ = n_threads_, 
                                                        verbose_level_ = verbose_level_)
    
    seed_p_values <- get_monte_carlo_simulation_p_values(spinglass_results$graphs, n_threads_ = n_threads_, verbose_level_ = verbose_level_)
    
    top_modules_tb <- spinglass_results$info %>%
      inner_join(seed_p_values, 'seed') %>%
      arrange(p_values) %>%
      head(n_top_modules_)
    
    top_modules <- top_modules_tb$seed
    
    go_analysis_df <- topGO_analysis(spinglass_results$genelists[top_modules], genes, eg_db_, verbose_level_ = verbose_level_)
    
    final_tb <- top_modules_tb %>%
      inner_join(go_analysis_df, 'seed') %>%
      arrange(p_values)
    
    return(final_tb)
  }
