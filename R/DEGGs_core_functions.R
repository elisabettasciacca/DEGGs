#' An S4 class to define the deggs output
#'
#' @slot subnetworks a list of data frames containing the subgroup networks and
#' a number indicating the total count of links that reached significance in all
#' subgroups (considering p or adj p < 0.05)
#' @slot normalised_counts a data frame containig the normalised count data
#' given in input.
#' @slot metadata a data frame of sample data given in input and matching
#' the sample IDs in `normalised_counts`.
#' @slot subgroup_variable column name in `metadata` that contains the
#' subgroup definition for each sample in `normalised_counts`.
#' @slot regression_method the chosen regression method used to calculate
#' interaction p values. It can be either 'rlm' or 'lm'.
#' @slot subgroups character vector indicating which subgroups
#' are used for comparison.
#' @slot use_qvalues logical. Indicates whether p values are adjusted via
#' Storey's q value method
methods::setClass("deggs", slots = list(
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
#' p values
#'
#' @param normalised_counts a data frame containing normalised counts from an
#' high throughput sequencing experiment.
#' Sample IDs must be in columns and gene/miRNA/TFs in rows.
#' Objects of class `matrix` are not allowed.
#' @param metadata a data frame of sample data with rownames matching the
#' sample IDs in `normalised_counts` colnames
#' @param subgroup_variable column name in `metadata` that contains the
#' subgroup identifier for each sample in `normalised_counts`
#' @param regression_method whether to use robust linear modelling to calculate
#' link p values. Options are 'rlm' (default) or 'lm'.
#' @param subgroups optional character vector indicating which subgroups
#' must be used for comparison. If not specified, all subgroups in
#' `subgroup_variable` will be used.
#' @param network network of biological interactions provided by the user. The 
#' network must be provided in the form of a table of class data.frame with two 
#' columns named "from" and "to". 
#' If NULL (default) a network of 10,537 molecular interactions obtained from
#' KEGG, mirTARbase, miRecords and transmiR will be used.
#' This has been obtained via the `exportgraph` function of the MITHrIL tool
#' (Alaimo et al., 2016). 
#' @param entrezIDs logical (default FALSE) used to define whether gene ids in
#'  `normalised_counts` are entrez ids (TRUE) or gene symbols (FALSE).
#' @param convert_to_gene_symbols logical to be used when using entrez ids.
#' If TRUE (default), the output will show gene symbols.
#' @param use_qvalues whether to use Storey's q values for multiple test
#' adjustment. If FALSE (default), unadjusted p values will be used and shown
#' in the output.
#' @param cores number of cores to use.
#' @importFrom rlang .data
#' @return a `deggs` object containing subgroup specific networks incorporating
#' p values or adjusted p values for each link. 
#' @seealso [`deggs-class`].
#' @export
generate_subnetworks <- function(normalised_counts,
                                 metadata,
                                 subgroup_variable,
                                 regression_method = 'rlm',
                                 subgroups = NULL,
                                 network = NULL,
                                 entrezIDs = FALSE,
                                 convert_to_gene_symbols = TRUE,
                                 use_qvalues = FALSE,
                                 cores = parallel::detectCores()/2){

  sig_var <- ifelse(use_qvalues, "q.value", "p.value")

  if(is.data.frame(normalised_counts) == FALSE){
    stop(paste0("normalised_counts is not a dataframe"))
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
    
    if(num_entrez_rows - nrow(normalised_counts) > 0) (
      message(paste0(num_entrez_rows - nrow(normalised_counts),
                     " genes had no matching gene symbol."))
    )
    
    # main network (gene symbols)
    if(is.null(network))(
      edges <- metapathway_gene_symbols %>%
        dplyr::filter(.data$from %in% rownames(normalised_counts)) %>%
        dplyr::filter(.data$to   %in% rownames(normalised_counts))
    ) else (
      edges <- network %>%
        dplyr::filter(.data$from %in% rownames(normalised_counts)) %>%
        dplyr::filter(.data$to   %in% rownames(normalised_counts))
    )
  } else if (entrezIDs == FALSE){
    # main network (gene symbols)
    if(is.null(network))(
    edges <- metapathway_gene_symbols %>%
      dplyr::filter(.data$from %in% rownames(normalised_counts)) %>%
      dplyr::filter(.data$to   %in% rownames(normalised_counts))
    ) else (
      edges <- network %>%
        dplyr::filter(.data$from %in% rownames(normalised_counts)) %>%
        dplyr::filter(.data$to   %in% rownames(normalised_counts))
    )

  } else {
    # main network (entrezIDs)
    if(is.null(network))(
    edges <- metapathway_entrez_IDs %>%
      dplyr::filter(.data$from %in% rownames(normalised_counts)) %>%
      dplyr::filter(.data$to   %in% rownames(normalised_counts))
    ) else (
      edges <- network %>%
        dplyr::filter(.data$from %in% rownames(normalised_counts)) %>%
        dplyr::filter(.data$to   %in% rownames(normalised_counts))
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

  combinations <- utils::combn(subgroups, m = 2) %>%
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

  # calculate p values list (with parallelisation)
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
      library("dplyr")
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

  degg <- methods::new("deggs",
                   subnetworks = final_networks,
                   normalised_counts = normalised_counts,
                   metadata = metadata,
                   subgroup_variable = subgroup_variable,
                   regression_method = regression_method,
                   subgroups = subgroups,
                   use_qvalues = use_qvalues)

  return(degg)
}


#' Tidying up of the metadata table. Samples belonging to unwanted subgroups 
#' (if specified) will be removed as well as subgroups with less than five samples,
#' and NAs. 
#'
#' @param subgroups optional character vector indicating which subgroups
#' are used for comparison. If not specified, all subgroups in
#' `subgroup_variable` will be used.
#' @param metadata a data frame of sample data with rownames matching the
#' sample IDs in `normalised_counts` colnames.
#' @param subgroup_variable column name in `metadata` that contains the
#' subgroup identifier for each sample in `normalised_counts`.
#' @return tidied and aligned metadata.
tidy_metadata <- function(subgroups,
                          metadata,
                          subgroup_variable){

  # select specified subgroups (optional)
  if(!is.null(subgroups)) {
    metadata <- subset(metadata, metadata[, subgroup_variable] %in% subgroups)
    if(is.factor(metadata[, subgroup_variable])) (
      metadata[, subgroup_variable] <- droplevels(metadata[, subgroup_variable])
    )
  }

  # if subgroup_variable is not a factor, convert to factor
  if(!is.factor(metadata[, subgroup_variable])){
    metadata[, subgroup_variable] <- as.factor(metadata[, subgroup_variable])
    message(paste0(subgroup_variable, " was converted to factor."))
  }

  # remove subgroups with less than five observations
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


#' Compute interaction p values for a single percentile value
#'
#' @inheritParams generate_subnetworks
#' @param subgroups_length integer number indicating the number of subgroups
#' @param subgroups_df_list list of subgroup data frames
#' @param sig_var Inherited from `generate_subnetworks`. It can be 
#' `q.value` or `p.value` depending on how `use_qvalues` was set in the 
#' `generate_subnetworks` function (default `FALSE`).
#' @param percentile a float number indicating the percentile to use. 
#' @param combinations data frame containing the subgroups combinations in rows
#' @param regression_method whether to use robust linear modelling to obtain 
#' p value of the interactions. Options are 'rlm' (default) or 'lm'
#' @param edges network of biological interactions in the form of a table of 
#' class data.frame with two columns: "from" and "to". 
#' @param sig_edges_count number of significant edges (p < 0.05)
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @return The list of float numbers of the significant pvalues
#' for a specific percentile
calc_pvalues_percentile <- function(normalised_counts,
                                    metadata,
                                    subgroup_variable,
                                    subgroups_length,
                                    subgroups_df_list,
                                    sig_var,
                                    percentile,
                                    combinations,
                                    regression_method = "rlm",
                                    edges,
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
    if(!is.character(subgroup_df)) {
      subgroupEdges <- edges %>%
        dplyr::filter(.data$from %in% names(subgroup_df)) %>%
        dplyr::filter(.data$to   %in% names(subgroup_df))
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
    if(!is.character(subgroup_network))(
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

  # count significant p values
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
#' @inheritParams calc_pvalues_percentile
#' @param subgroup_network network table for a specific subgroup
#' @importFrom methods is
#' @return a list of p values
calc_pvalues_network <- function(normalised_counts,
                                 metadata,
                                 sig_var,
                                 subgroup_variable,
                                 subgroups_length,
                                 regression_method = 'rlm',
                                 subgroup_network){

  if(!is.character(subgroup_network)){
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

  if(is.list(p_values)){
    # make a data frame with all values
    p_values <- as.data.frame(do.call("rbind", p_values))
    p_values$from <- as.character(p_values$from)
    p_values$to <- as.character(p_values$to)

    if(sig_var == "q.value"){
      # adding Storey's q values
      q.values <- try(qvalue::qvalue(p_values[, "p.value"])$qvalues)
      if (is(q.values, "try-error")) {
        if(nrow(p_values > 1))(
          q.values <-  qvalue::qvalue(p = p_values[, "p.value"], pi0 = 1)$qvalues
        ) else (
          q.values <- NA
        )
      }
      p_values$q.value <- q.values
    }
  }

  if(!is.character(p_values)){
    rownames(p_values) <- paste(p_values$from, p_values$to, sep = "-")
  }

  return(p_values)
}


#' Output a table of all the significant gene-gene interactions across subgroups
#'
#' @param deggs_object an object of class `deggs` generated from
#' `generate_subnetworks`
#' @return a `data.frame` listing all the significant gene-gene interactions 
#' found across subgroups.
#' This can be used as features selection method when building machine learning
#' models for the prediction of the subgroups.
#' @export
#' @seealso [`deggs-class`].
extract_sig_deggs <- function(deggs_object){
  use_qvalues <- deggs_object@use_qvalues
  sig_var <- ifelse(use_qvalues, "q.value", "p.value")

  counts   <- deggs_object@normalised_counts
  metadata <- deggs_object@metadata
  subgroup_variable <- deggs_object@subgroup_variable
  model <- deggs_object@regression_method
  # extracting subnetworks (we just need to exclude sig_pvalues_count from the list)
  condition <- lapply(deggs_object@subnetworks, is.list)
  subnetworks_list <- deggs_object@subnetworks[unlist(condition)]

  # exrtact all significant gene pairs (from any subgroup),
  # these will be tested for prediction
  sig.edges <- lapply(subnetworks_list, function(subnetwork){
    subnetwork <- subnetwork[which(subnetwork[, sig_var] < 0.05), ]
  })
  sig.edges <- do.call(rbind, sig.edges)
}
