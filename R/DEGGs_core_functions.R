#' An S4 class to define the Deggs output
#'
#' @slot subnetworks a list of specific networks in the form of tables (data.frame),
#' and the count of significant interactions
#' @slot normalised_counts a data frame containig the normalised tidied
#' count data.
#' @slot metadata a tidied data frame of sample information matching sample IDs
#' in `normalised_counts`.
#' @slot subgroup_variable column name in `metadata` which contains the
#' subgroup definition for each sample in `normalised_counts`.
#' @slot regression_method the chosen regression method used to calculate
#' interaction p-values. It can be either 'rlm' or 'lm'.
#' @slot subgroups character vector indicating which subgroups
#' are used for comparison.
setClass("Deggs", slots = list(
  subnetworks       = "list",
  normalised_counts = "data.frame",
  metadata          = "data.frame",
  subgroup_variable = "character",
  regression_method = "character",
  subgroups         = "character",
  use_qvalues       = "logical"
))



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
#' @param use_qvalues Whether to use Storey's q-values. If FALSE, unadjusted
#' p-values will be used.
#' @param cores Number of cores to use.
#' @return a Deggs object with subgroup specific networks incorporating pvalues
#' for each interaction
#' @export
generate_subnetworks <- function(normalised_counts, metadata, subgroup_variable,
                                 regression_method = 'rlm', subgroups = NULL,
                                 network = NULL,
                                 entrezIDs = FALSE,
                                 convert_to_gene_symbols = TRUE,
                                 use_qvalues = TRUE,
                                 cores = parallel::detectCores()/2){

  sig_var <- ifelse(use_qvalues, "q.value", "p.value")

  if(is.data.frame(normalised_counts) == FALSE){
    message(paste0("normalised_counts is not a dataframe"))
    return
  }

  if(!subgroup_variable %in% colnames(metadata)){
    stop("subgroup_variable must be %in% colnames metadata")
  }

  if(entrezIDs == TRUE && convert_to_gene_symbols == TRUE){
    num_entrez_rows <- nrow(normalised_counts)
    normalised_counts$genesymbol <- suppressMessages(
      AnnotationDbi::mapIds(x = org.Hs.eg.db::org.Hs.eg.db,
                            keys = rownames(normalised_counts),
                            keytype = 'ENTREZID',
                            column = 'SYMBOL'))
    normalised_counts <- subset(normalised_counts,
                                !is.na(normalised_counts$genesymbol))
    rownames(normalised_counts) <- normalised_counts$genesymbol
    normalised_counts$genesymbol <- NULL

    message(paste0(num_entrez_rows - nrow(normalised_counts),
                   " genes had no matching gene symbol."))

    # main network (gene symbols)
    if(is.null(network))(
      edges <- metapathway_gene_symbols %>%
        dplyr::filter(from %in% rownames(normalised_counts)) %>%
        dplyr::filter(to   %in% rownames(normalised_counts))
    ) else (
      edges <- network %>%
        dplyr::filter(from %in% rownames(normalised_counts)) %>%
        dplyr::filter(to   %in% rownames(normalised_counts))
    )
  } else if (entrezIDs == FALSE){
    # main network (gene symbols)
    if(is.null(network))(
    edges <- metapathway_gene_symbols %>%
      dplyr::filter(from %in% rownames(normalised_counts)) %>%
      dplyr::filter(to   %in% rownames(normalised_counts))
    ) else (
      edges <- network %>%
        dplyr::filter(from %in% rownames(normalised_counts)) %>%
        dplyr::filter(to   %in% rownames(normalised_counts))
    )

  } else {
    # main network (entrezIDs)
    if(is.null(network))(
    edges <- metapathway_entrez_IDs %>%
      dplyr::filter(from %in% rownames(normalised_counts)) %>%
      dplyr::filter(to   %in% rownames(normalised_counts))
    ) else (
      edges <- network %>%
        dplyr::filter(from %in% rownames(normalised_counts)) %>%
        dplyr::filter(to   %in% rownames(normalised_counts))
    )

  }
  edges$edge_ID <- paste(edges$from, edges$to)
  edges <- edges[!duplicated(edges$edge_ID),]
  edges$edge_ID <- NULL

  metadata <- tidy_metadata(subgroups = subgroups, metadata = metadata,
                            subgroup_variable = subgroup_variable)

  # align metadata and count data
  nodes <- c(edges[,1], edges[,2]) %>%
    unique()
  normalised_counts <- normalised_counts[rownames(normalised_counts) %in% nodes,
                                         rownames(metadata)]
  metadata <- metadata[colnames(normalised_counts), , drop = FALSE]

  if(is.null(subgroups)) {
    subgroups <- levels(metadata[, subgroup_variable])
  }

  combinations <- combn(subgroups, m = 2) %>%
    as.data.frame()

  # create subgroups (duplicating count data for each subgroup)
  subgroups_df_list <- lapply(subgroups, function(one_subgroup){
    metadata_subset_subgroup <- subset(metadata,
                                       metadata[, subgroup_variable] == one_subgroup)

    subgroup_df <- normalised_counts[,rownames(metadata_subset_subgroup)]
    subgroup_df <- apply(subgroup_df, 1, mean, na.rm = TRUE)
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
                        "combinations", "regression_method", "edges",
                        "subgroup_variable", "subgroups", "calc_pvalues_percentile",
                        "normalised_counts", "calc_pvalues_network", "metadata",
                        "sig_edges_count", "sig_var"), envir = environment())

    parallel::clusterEvalQ(cl, {
      library("igraph")
      library("dplyr")
      library("edgeR")
      library("networkD3")
    })
    pvalues_list <- pbapply::pblapply(cl = cl, percentile_vector, function(percentile){
      calc_pvalues_percentile(normalised_counts = normalised_counts,
                              sig_var = sig_var,
                              metadata = metadata,
                              percentile = percentile,
                              subgroups_df_list = subgroups_df_list,
                              combinations = combinations,
                              regression_method = regression_method,
                              edges = edges,
                              subgroup_variable = subgroup_variable,
                              subgroups_length = length(subgroups),
                              sig_edges_count = sig_edges_count)
    })
    parallel::stopCluster(cl)
  } else{
    # Parallelisation for mac/linux
    pvalues_list <- pbmcapply::pbmclapply(mc.cores = cores, percentile_vector,

                               function(percentile)(
                                 calc_pvalues_percentile(normalised_counts = normalised_counts,
                                                         sig_var = sig_var,
                                                         metadata = metadata,
                                                         percentile = percentile,
                                                         subgroups_df_list = subgroups_df_list,
                                                         combinations = combinations,
                                                         regression_method = regression_method,
                                                         edges = edges,
                                                         subgroup_variable = subgroup_variable,
                                                         subgroups_length = length(subgroups),
                                                         sig_edges_count = sig_edges_count)
                               ))
  }
  names(pvalues_list) <- percentile_vector

  # extract the network filtered on the best threshold percentile
  # (highest number of statistically significant interactions)
  sig_pvalues <- unlist(lapply(pvalues_list,
                               function(networks) (networks$sig_pvalues_count)))

  if(is.null(sig_pvalues)){
    stop("Any significant difference across subgroups.")
  }

  best_percentile <- sig_pvalues[which(sig_pvalues == max(sig_pvalues))]
  if(length(best_percentile) > 1){
    best_percentile <- best_percentile[1]
  }
  message(paste0("Percolation analysis: genes whose expression is below the ",
                 as.numeric(names(best_percentile)) * 100,
                 "th percentile are removed from networks."))
  final_networks <- pvalues_list[[as.character(names(best_percentile))]]

  degg <- new("Deggs",
              subnetworks = final_networks,
              normalised_counts = normalised_counts,
              metadata = metadata,
              subgroup_variable = subgroup_variable,
              regression_method = regression_method,
              subgroups = subgroups,
              use_qvalues = use_qvalues)

  return(degg)
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
#' @param combinations dataframe containing the subgroups combinations in rows
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
calc_pvalues_percentile <- function(normalised_counts,
                                    sig_var,
                                    metadata,
                                    percentile,
                                    subgroups_df_list,
                                    combinations,
                                    regression_method = "rlm",
                                    edges,
                                    subgroup_variable,
                                    subgroups_length,
                                    sig_edges_count){

  # 1st filtering step (remove low expressed genes,
  # i.e. genes under the percentile threshold)
  cut_off <- stats::quantile(as.matrix(normalised_counts), prob = percentile)
  user_message <- paste0("No gene above the threshold (", percentile*100, "th percentile).")

  subgroups_df_list <- lapply(subgroups_df_list, function(subgroup_df){
    subgroup_df <- subgroup_df[subgroup_df > cut_off]
    if(length(subgroup_df) == 0) (
      subgroup_df <- user_message
    )
    return(subgroup_df)
  })

  # format to string vectors to allow easy detection of overlapping edges
  networks_to_string <- lapply(subgroups_df_list, function(subgroup_df){
    if(class(subgroup_df) != "character") {
      subgroupEdges <- edges %>%
        dplyr::filter(from %in% names(subgroup_df)) %>%
        dplyr::filter(to   %in% names(subgroup_df))
      network_to_string <- do.call(paste, subgroupEdges)
    } else {
      network_to_string <- user_message
    }
    return(network_to_string)
  })

  # find overlapping edges across subgroups
  common_links <- lapply(combinations, function(combination){
    return(intersect(networks_to_string[[combination[1]]],
                     networks_to_string[[combination[2]]]))
  }) %>%
    unlist() %>%
    unique()

  if(user_message %in% common_links)(
    common_links <- common_links[!(common_links %in% user_message)]
  )

  # remove overlapping edges and format as data frame
  subgroups_network_list <- lapply(networks_to_string, function(subnetwork){
    subnetwork <- subnetwork[!subnetwork %in% common_links]
    if(!(user_message %in% subnetwork))(
      subnetwork <- data.frame('from' = unlist(lapply(strsplit(subnetwork, " "), `[[`, 1)),
                               'to'   = unlist(lapply(strsplit(subnetwork, " "), `[[`, 2)) )
    )
    return(subnetwork)
  })

  # count tot edges left
  tot_edges <- unlist(lapply(subgroups_network_list, function(subgroup_network){
    if(class(subgroup_network) != "character")(
      nrow(subgroup_network)
    )
  })) %>%
    sum()

  # calculate interaction p values
  pvalues_list <- lapply(subgroups_network_list, function(subgroup_network){

    return(calc_pvalues_network(subgroup_network = subgroup_network,
                                normalised_counts = normalised_counts,
                                metadata = metadata,
                                subgroup_variable = subgroup_variable,
                                regression_method = regression_method,
                                subgroups_length = subgroups_length,
                                sig_var = sig_var))
  })

  # count significant p-values
  if(tot_edges > sig_edges_count){
    # num tot edges greater than previous sig edges count
    p_values_sig_count <- unlist(lapply(pvalues_list, function(subgroup_network){
      if(!is.null(dim(subgroup_network)))(
        return(length(which(subgroup_network[, sig_var] < 0.05))) # these are either pvalues or qvalues
      )
      else (
        return(0)
      )
    }))
    num_sig_pvalues <- sum(p_values_sig_count) # these are either pvalues or qvalues
    if(num_sig_pvalues > sig_edges_count){
      sig_edges_count <<- num_sig_pvalues
    }
    pvalues_list <- append(pvalues_list, num_sig_pvalues)
    names(pvalues_list)[length(pvalues_list)] <- "sig_pvalues_count"

    return(pvalues_list)

  } else {
    # previous sig edges count greater than num tot edges
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
#' of subgroups
#' @return The list of pvalues
calc_pvalues_network <- function(subgroup_network, normalised_counts, sig_var,
                                 metadata, subgroup_variable,
                                 regression_method = 'rlm', subgroups_length){

  if(class(subgroup_network) != "character"){
    if(nrow(subgroup_network) > 0){
      # prepare data
      df_list <- mapply(function(gene_B, gene_A){
        return(data.frame(t(normalised_counts[gene_A, ]),
                          t(normalised_counts[gene_B, ]),
                          metadata[, subgroup_variable],
                          check.names = FALSE))
      }, gene_A = subgroup_network$from, gene_B = subgroup_network$to, SIMPLIFY = F)

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
    } else {
      p_values <- "No specific links for this subgroup."
      }
    } else {
      p_values <- "No specific links for this subgroup."
    }

  if(class(p_values) == "list"){
    # make a data frame with all values
    p_values <- as.data.frame(do.call("rbind", p_values))
    p_values$from <- as.character(p_values$from)
    p_values$to <- as.character(p_values$to)

    if(sig_var == "q.value"){
      # adding Storey's q-values
      q.values <- try(qvalue::qvalue(p_values[, "p.value"])$qvalues)
      if (class(q.values) == "try-error") {
        if(nrow(p_values > 1))(
          q.values <-  qvalue::qvalue(p = p_values[, "p.value"], pi0 = 1)$qvalues
        ) else (
          q.values <- NA
        )
      }
      p_values$q.value <- q.values
    }
  }

  if(class(p_values) != "character"){
    rownames(p_values) <- paste(p_values$from, p_values$to, sep = "-")
  }

  return(p_values)
}
