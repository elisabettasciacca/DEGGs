#' Generate subnetworks
#'
#' Generate subgroup specific gene-gene interaction networks with interaction
#' pvalues
#'
#' @param normalised_counts a data frame containig the normalised sequencing
#' count data. Sample IDs must be in columns and gene/miRNA/TFs in rows.
#' Objects of class `matrix` are not allowed.
#' @param metadata a data frame of sample information matching the sample IDs in
#' `normalised_counts`
#' @param subgroup_variable column name in `metadata` which contains the
#' subgroup definition for each sample in `normalised_counts`
#' @param regression_method whether to use robust linear modeling for the
#' interaction p-values. Options are 'rlm' (default) or 'lm'
#' @param subgroups optional character vector indicating which subgroups
#' are used for comparison. If not specified, all subgroups in
#' `subgroup_variable` will be considered
#' @param entrezIDs Logical whether gene ids in `normalised_counts` are in
#' entrez id (TRUE) or gene symbols (FALSE). Default FALSE
#' @param convert_to_gene_symbols Logical to be used when using entrez ids.
#' If TRUE (default), the output will show gene symbols
#' @return a DEGGs object with subgroup specific networks incorporating pvalues
#' for each interaction
#' @export
generate_subnetworks <- function(normalised_counts, metadata, subgroup_variable,
                                 regression_method = 'rlm', subgroups = NULL,
                                 entrezIDs = FALSE,
                                 convert_to_gene_symbols = TRUE,
                                 cores = parallel::detectCores()/2){

  if(is.data.frame(normalised_counts) == FALSE){
    message(paste0("normalised_counts is not a dataframe"))
    return
  }

  if(entrezIDs == TRUE && convert_to_gene_symbols == TRUE){
    num_entrez_rows <- nrow(normalised_counts)
    normalised_counts$genesymbol <- suppressMessages(
      AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, rownames(normalised_counts),
                           'SYMBOL', 'ENTREZID'))
    normalised_counts <- subset(normalised_counts,
                                !is.na(normalised_counts$genesymbol))
    rownames(normalised_counts) <- normalised_counts$genesymbol
    normalised_counts$genesymbol <- NULL

    message(paste0(num_entrez_rows - nrow(normalised_counts),
                   " genes had no matching gene symbol."))

    # generate graph of the main network (gene symbols)
    edges <- metapathway_gene_symbols
    main_graph <- igraph::graph.data.frame(edges, directed = FALSE)
  } else if(entrezIDs == FALSE){
    # generate graph of the main network (gene symbols)
    edges <- metapathway_gene_symbols
    main_graph <- igraph::graph.data.frame(edges, directed = FALSE)
  } else {
    # generate graph of the main network (entrezIDs)
    edges <- metapathway_entrez_IDs
    main_graph <- igraph::graph.data.frame(edges, directed = FALSE)
  }

  metadata <- tidy_metadata(subgroups = subgroups, metadata = metadata,
                            subgroup_variable = subgroup_variable)

  # align metadata and count data
  normalised_counts <- normalised_counts[, rownames(metadata)]
  metadata <- metadata[colnames(normalised_counts), , drop = FALSE]

  if(is.null(subgroups)) {
    subgroups <- levels(metadata[, subgroup_variable])
  }

  permutations <- as.data.frame(t(gtools::permutations(n = length(subgroups), r = 2,
                                               v = subgroups,
                                               repeats.allowed = FALSE)))

  # create lists of samples for each subgroup
  subgroups_df_list <- lapply(subgroups, function(one_subgroup){
    metadata_subset_subgroup <- subset(metadata,
                                       metadata[, subgroup_variable] == one_subgroup)
    subgroup_df <- normalised_counts[,rownames(metadata_subset_subgroup)]
    subgroup_df <- apply(subgroup_df, 1, mean, na.rm = TRUE) #to check
    rownames(subgroup_df) <- NULL
    return(subgroup_df)
  })
  names(subgroups_df_list) <- subgroups

  # calculate pvalues list (with parallelisation)
  percentile_vector <- seq(0.7, 0.98, by = 0.05)
  sig_edges_count <- 0

  if (Sys.info()["sysname"] == "Windows") {
    cl <- parallel::makeCluster(cores)
    parallel::clusterExport(cl, c("percentile_vector", "subgroups_df_list",
                        "permutations", "regression_method", "edges",
                        "subgroup_variable", "subgroups", "calc_pvalues_percentile",
                        "normalised_counts", "calc_pvalues_network", "metadata",
                        "sig_edges_count"), envir = environment())
    parallel::clusterEvalQ(cl, {
      library("gtools")
      library("igraph")
      library("dplyr")
      library("edgeR")
      library("networkD3")
    })
    pvalues_list <- pbapply::pblapply(cl = cl, percentile_vector, function(percentile){
      calc_pvalues_percentile(normalised_counts = normalised_counts,
                              metadata = metadata,
                              percentile = percentile,
                              subgroups_df_list = subgroups_df_list,
                              permutations = permutations,
                              regression_method = regression_method,
                              edges = edges,
                              subgroup_variable = subgroup_variable,
                              subgroups_length = length(subgroups))
    })
    parallel::stopCluster(cl)
  } else{
    # parallelisation for mac/linux
    pvalues_list <- pbmcapply::pbmclapply(mc.cores = cores, percentile_vector,
                               function(percentile)(
                                 calc_pvalues_percentile(normalised_counts = normalised_counts,
                                                         metadata = metadata,
                                                         percentile = percentile,
                                                         subgroups_df_list = subgroups_df_list,
                                                         permutations = permutations,
                                                         regression_method = regression_method,
                                                         edges = edges,
                                                         subgroup_variable = subgroup_variable,
                                                         subgroups_length = length(subgroups))
                               ))
  }
  names(pvalues_list) <- percentile_vector

  # extract the network filtered on the best threshold percentile
  # (highest number of statistically significant interactions)
  sig_pvalues <- unlist(lapply(pvalues_list,
                               function(networks) (networks$sig_pvalues_count)))

  if(is.null(sig_pvalues)){
    stop("Any significant difference between categories.")
  }

  best_percentile <- sig_pvalues[which(sig_pvalues == max(sig_pvalues))]
  if(length(best_percentile) > 1){
    best_percentile <- best_percentile[1]
  }
  message(paste0("Percolation analysis: genes whose expression is below the ",
                 as.numeric(names(best_percentile)) * 100,
                 "th percentile are removed from networks."))
  final_networks <- pvalues_list[[as.character(names(best_percentile))]]

  return(list("subnetworks" = final_networks, "metadata" = metadata,
              "normalised_counts" = normalised_counts))


}


#' Tidy metadata by removal of unwanted subgroup samples (if specified),
#' subgroups with less than five samples, and NAs
#'
#' @param subgroups optional character vector indicating which subgroups
#' are used for comparison. If not specified, all subgroups in
#' `subgroup_variable` will be considered
#' @param metadata a data frame of sample information matching the sample IDs in
#' `normalised_counts`
#' @param subgroup_variable column name in `metadata` which contains the
#' subgroup definition for each sample in `normalised_counts`
#' @return tidied and aligned metadata
tidy_metadata <- function(subgroups, metadata, subgroup_variable){

  # select specified subgroups (optional)
  if(!is.null(subgroups)) {
    metadata <- subset(metadata, metadata[, subgroup_variable] %in% subgroups)
    if(is.factor(metadata[, subgroup_variable])) (
      metadata[, subgroup_variable] <- droplevels(metadata[, subgroup_variable])
    )
  }

  # if subgroup_variable is not a factor, conver to factor
  if(!is.factor(metadata[, subgroup_variable])){
    metadata[, subgroup_variable] <- as.factor(metadata[, subgroup_variable])
    message(paste0(subgroup_variable, " was converted to factor."))
  }

  # remove subgroups of less than five observations
  # (regression would not be reliable enough)
  tbl <- table(metadata[, subgroup_variable])
  if (length(names(tbl)[tbl < 5]) > 0) {
    message(paste0("The ", names(tbl)[tbl < 5],
                   " subgroup did not contain enough samples
                   (less than five observations). ",
                   "This subgroup will be removed.\n"))
    metadata <- subset(metadata, metadata[, subgroup_variable]
                       %in% names(tbl)[tbl > 5])
  }

  # remove NAs
  NAs_number <- length(which(is.na(metadata[, subgroup_variable])))
  if(NAs_number > 0) {
    metadata <- metadata[which(!is.na(metadata[, subgroup_variable])),]
    metadata[, subgroup_variable] <- droplevels(metadata[, subgroup_variable])
    message(paste0(subgroup_variable, " column contained NAs. ",
                   NAs_number, " samples were removed."))
  }

  return(metadata)
}


#' Compute interaction pvalues for a percentile
#'
#' @param percentile A float number that represents the percentile
#' @param subgroups_df_list list of subgroup dataframes
#' @param permutations dataframe containing the subgroups permutations in rows
#' @param regression_method whether to use robust linear modeling for the
#' interaction p-values. Options are 'rlm' (default) or 'lm'
#' @param edges A dataframe of the meta pathway edges
#' @param subgroup_variable A string that represent the column name of the subgroup
#' inside metadata the dataframe
#' @param subgroups_length An integer number that represent the number
#' of subgroups inside the metadata dataframe
#' @importFrom magrittr %>%
#' @return The list of float numbers of the significant pvalues
#' for a specific percentile
calc_pvalues_percentile <- function(normalised_counts, metadata, percentile,
                                    subgroups_df_list, permutations,
                                    regression_method = "rlm", edges,
                                    subgroup_variable,
                                    subgroups_length){

  # 1st filtering step (remove low expressed genes,
  # i.e. genes under the percentile threshold)
  cut_off <- stats::quantile(as.matrix(normalised_counts), prob = percentile)
  subgroups_df_list <- lapply(subgroups_df_list, function(subgroup_df)(
    subgroup_df[subgroup_df > cut_off]
  ))

  # subgroup network list creation
  subgroups_network_list <- lapply(subgroups_df_list, function(subgroup_df){
    subgroupEdges <- edges %>%
      dplyr::filter(from %in% names(subgroup_df)) %>%
      dplyr::filter(to   %in% names(subgroup_df))
    return(igraph::graph.data.frame(subgroupEdges, directed=FALSE))
  })

  # intersections
  network_intersection_list <- lapply(permutations, function(permutation){
    return(igraph::intersection(subgroups_network_list[[permutation[1]]],
                        subgroups_network_list[[permutation[2]]]))
  })

  # union
  union_graph <- base::do.call(igraph::union, network_intersection_list)

  # remove unions (i.e. remove overlapping edges across networks)
  subgroups_network_list <- lapply(subgroups_network_list, function(subgroup_network){
    return(igraph::difference(subgroup_network, union_graph))
  })

  # count tot edges left
  tot_edges_list <- unlist(lapply(subgroups_network_list, function(subgroup_network){
    igraph::gsize(subgroup_network)
  }))
  tot_edges <- sum(tot_edges_list)

  # calculate interaction p values
  pvalues_list <- lapply(subgroups_network_list, function(subgroup_network){

    return(calc_pvalues_network(subgroup_network = subgroup_network,
                                normalised_counts = normalised_counts,
                                metadata = metadata,
                                subgroup_variable = subgroup_variable,
                                regression_method = regression_method,
                                subgroups_length = subgroups_length))

  })

  # count significant p values
  if(tot_edges > sig_edges_count){
    # num tot edges greater than previous sig edges count
    p_values.sig.count <- unlist(lapply(pvalues_list, function(subgroup_network){
      length(which(subgroup_network$p.value < 0.05))
    }))
    num.sig_pvalues <- sum(p_values.sig.count)
    if(num.sig_pvalues > sig_edges_count){
      sig_edges_count <<- num.sig_pvalues
    }
    pvalues_list <- append(pvalues_list, num.sig_pvalues)
    names(pvalues_list)[length(pvalues_list)] <- "sig_pvalues_count"

    return(pvalues_list)

  } else {
    # previous sig edges count greater than total num of edges
    return(NULL)

  }

}


#' Calculate the pvalues for specific subgroup network samples
#'
#' @param subgroup_network A igraph network related to a specific subgroup samples
#' @param normalised_counts A dataframe that represent the gene expressions
#' of samples
#' @param metadata A dataframe that represent the subgroup of samples
#' @param subgroup_variable A string that represent the column name of the subgroup
#' inside metadata the dataframe
#' @param regression_method whether to use robust linear modeling for the
#' interaction p-values. Options are 'rlm' (default) or 'lm'
#' @param subgroups_length An integer number that represent the number
#' of subgroups inside the metadata dataframe
#' @return The list of pvalues
calc_pvalues_network <- function(subgroup_network, normalised_counts, metadata,
                                 subgroup_variable, regression_method = 'rlm',
                                 subgroups_length){

  genes <- as_edgelist(subgroup_network, names = TRUE)
  genes <- base::unique(genes)

  # prepare data
  df_list <- mapply(function(gene_B, gene_A){
    return(data.frame(t(normalised_counts[gene_A, ]),
                      t(normalised_counts[gene_B, ]),
                      metadata[, subgroup_variable],
                      check.names = FALSE))
  }, gene_B = genes[, 2], gene_A = genes[, 1], SIMPLIFY = F)

  # calculate regressions with interaction term
  p_values <- lapply(df_list, function(df){
    if(subgroups_length == 2){
      if(regression_method == "lm"){
        # gene_B ~ gene_A * subgroup
        lmfit <- stats::lm(df[,2] ~ df[,1] * df[,3])
        p_interaction <- stats::coef(summary(lmfit))[4,4]
      }
      if(regression_method == "rlm"){
        # gene_B ~ gene_A * subgroup
        robustfit <- MASS::rlm(df[,2] ~ df[,1] * df[,3])
        p_interaction <- sfsmisc::f.robftest(robustfit, var=3)$p.value
      }
      output <- data.frame(from = colnames(df)[1], to = colnames(df)[2],
                           p.value = p_interaction)
    }
    if(subgroups_length >= 3){
      # one-way ANOVA
      # gene_B ~ gene_A * subgroup
      res_aov <- stats::aov(df[,2] ~ df[,1] * df[,3], data = df)
      p_interaction <- summary(res_aov)[[1]][["Pr(>F)"]][3]
      output <- data.frame(from = colnames(df)[1], to = colnames(df)[2],
                           p.value = p_interaction)
    }
    return(output)
  })
  p_values <- do.call("rbind", p_values)
  p_values$from <- as.character(p_values$from)
  p_values$to   <- as.character(p_values$to)

  return (p_values)
}
