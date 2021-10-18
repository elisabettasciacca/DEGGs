#' View subnetworks in basic igraph style
#'
#' View igraph style networks for a specific subgroup
#'
#' @param DEGGs_object an object of class DEGGs generated from
#' `generate_subnetworks`
#' @param subgroup the subgroup of interest. The associated network will be
#' visualised
#' @return an igraph plot of the selected subnetwork
#' @export
print_simple_network <- function(DEGGs_object, subgroup){
  network <- subnetworks_object[["subnetworks"]][[subgroup]]
  df <- as.data.frame(t(data.frame(matrix(unlist(network),
                                          nrow = length(network),
                                          byrow = TRUE))))
  colnames(df)[1] <- "from"
  colnames(df)[2] <- "to"
  colnames(df)[3] <- "weight"
  # ES: TODO: add vertex names
  df <- unique(df)
  df$color <- ifelse(df$weight >= 0.05, "gray", "blue")
  g <- igraph::graph.data.frame(d = df, directed = FALSE)
  plot(g, vertex.size = 3, vertex.label = NA, edge.width = igraph::E(g)$width)
}


#' View subnetworks in network3D style
#'
#' View network3D style networks for a specific subgroup
#'
#' @param DEGGs_object an object of class DEGGs generated from
#' `generate_subnetworks`
#' @param subgroup the subgroup of interest. The associated network will be
#' visualised
#' @return a network3D plot of the selected subnetwork
#' @export
print_force_network <- function(DEGGs_object, subgroup){
  network <- subnetworks_object[["subnetworks"]][[subgroup]]
  nodes <- data.frame("name" = unique(c(network$from, network$to)))
  nodes$ID <- 0:(nrow(nodes)-1) # assign sequential IDs to nodes
  network <- base::merge(nodes, network, by.x = "name", by.y = "from")
  network <- dplyr::rename(network, source = ID)
  network <- base::merge(nodes, network, by.x = "name", by.y = "to")
  network <- dplyr::rename(network, target = ID)
  network <- dplyr::rename(network, value = p.value)
  network$name   <- NULL
  network$name.y <- NULL
  nodes$group <- "node"

  networkD3::forceNetwork(Links = network, Nodes = nodes,
                          Source = 'source',
                          Target = 'target',
                          NodeID = 'name',
                          Group = "group",
                          linkColour = ifelse(network$value >= 0.05, "gray", "red"),
                          linkWidth = networkD3::JS("function(d) { return d; }"),
                          opacity = 1,
                          opacityNoHover = 0.8,
                          zoom = TRUE)
}


#' Print the regression graph of a specific genes
#'
#' @param DEGGs_object an object of class DEGGs generated from
#' `generate_subnetworks`
#' @param subgroup_variable column name in `metadata` which contains the
#' subgroup definition for each sample in `normalised_counts`
#' @param regression_method whether to use robust linear modeling for the
#' interaction p-values. Options are 'rlm' (default) or 'lm'
#' @param gene_A A string of gene selected
#' @param gene_B A string of gene selected
#' @param subgroups optional character vector indicating which subgroups
#' are used for comparison. If not specified, all subgroups in
#' `subgroup_variable` will be considered
#' @param legend_offset optional numeric offset to horizontally move the legend
#' @return base graphics plot showing gene-gene regressions for each subgroup
#' while the pvalue of the interaction term of
#' *gene A ~ gene B \* subgroup* is reported on top
#' @export
print_regressions <- function (DEGGs_object, gene_A, gene_B,
                               subgroup_variable,
                               regression_method = 'rlm',
                               subgroups = NULL,
                               legend_offset = -0.33){

  metadata <- subnetworks_object$metadata
  normalised_counts <- subnetworks_object$normalised_counts

  if(is.null(subgroups)) {
    subgroups <- levels(metadata[, subgroup_variable])
  }

  subgroups_length <- length(subgroups)

  # prepare data frame
  df <- data.frame(t(normalised_counts[gene_A, ]),
                   t(normalised_counts[gene_B, ]),
                   metadata[, subgroup_variable],
                   check.names = FALSE)
  colnames(df)[3] <- subgroup_variable

  # compute gene-gene regression
  if(subgroups_length == 2){
    if(regression_method == "lm"){
      # gene A ~ gene B * subgroup
      lmfit <- stats::lm(df[,2] ~ df[,1] * df[,3])
      p_interaction <- stats::coef(summary(lmfit))[4,4]
      fit <- lapply(subgroups, function(i) {
        x <- df[df[, subgroup_variable] == i, 1]
        y <- df[df[, subgroup_variable] == i, 2]
        if(length(x) > 0 || length(y) > 0) return(stats::lm(y ~ x))
      })
    }
    if(regression_method == "rlm"){
      # gene A ~ gene B * subgroup
      robustfit <- MASS::rlm(df[,2] ~ df[,1] * df[,3])
      p_interaction <- sfsmisc::f.robftest(robustfit, var = 3)$p.value
      fit <- lapply(subgroups, function(i) {
        x <- df[df[, subgroup_variable] == i, 1]
        y <- df[df[, subgroup_variable] == i, 2]
        if(length(x) > 0 || length(y) > 0) return(MASS::rlm(y ~ x))
      })
    }
  }
  if(subgroups_length >= 3){
    # one-way ANOVA
    # gene A ~ gene B * subgroup
    res_aov <- stats::aov(df[,2] ~ df[,1] * df[,3], data = df)
    p_interaction <- summary(res_aov)[[1]][["Pr(>F)"]][3]
    fit <- lapply(subgroups, function(i) {
      x <- df[df[, subgroup_variable] == i, 1]
      y <- df[df[, subgroup_variable] == i, 2]
      if(regression_method == "lm"){
        if(length(x) > 0 || length(y) > 0) return(stats::lm(y ~ x));
      }
      if(regression_method == "rlm"){
        if(length(x) > 0 || length(y) > 0) return(MASS::rlm(y ~ x))
      }

    })
  }

  # Plot
  col <- viridis::viridis(n = subgroups_length)
  x_adj <- (max(df[,1], na.rm = T) - min(df[,1], na.rm = T)) * 0.05
  new_x <- seq(min(df[,1], na.rm = T) - x_adj,
               max(df[,1], na.rm=T) + x_adj,
               length.out=100)

  # prediction of the fitted model
  preds <- lapply(fit,
                  function(i) stats::predict(i,
                                             newdata = data.frame(x = new_x),
                                             interval = 'confidence'))

  par(mar = c(5.2, 6, 3.3, 10.7), xpd = FALSE)
  plot(df[,1], df[,2], type = 'n', bty = 'l', las = 1, cex.axis = 1.1,
       font.main = 1, cex.lab = 1.3, xlab = colnames(df)[1],
       ylab = colnames(df)[2])

  cols <- col[df[, subgroup_variable]]
  pch <- c(16 : (16 + subgroups_length -1))[df[, subgroup_variable]]
  for(i in 1:subgroups_length) {

    # plot confidence intervals
    polygon(c(rev(new_x), new_x), c(rev(preds[[i]][ ,3]), preds[[i]][ ,2]),
            col = adjustcolor(col[i], alpha.f = 0.15), border = NA)

    # plot regression lines
    abline(fit[[i]], col = col[i], lwd = 1.5)
    row <- points(df[df[, subgroup_variable] == subgroups[i], 1],
                  df[df[, subgroup_variable] == subgroups[i], 2], cex = 1.5,
                  pch = pch, col = adjustcolor(cols, alpha.f = 0.7))
  }

  mtext(bquote(paste("P"["interaction"]*"=",
                     .(format(p_interaction, digits = 2)))),
        cex = 1.3, side = 3, adj = 0.04)

  par(xpd=TRUE)
  legend("topright", bty = "n", inset = c(legend_offset, 0), legend = subgroups,
         col = col, lty = 1, cex = 1.2)
}




#' View interactive shiny app style subnetworks
#'
#' Explore subnetworks and interactively select regression plots
#'
#' @param DEGGs_object an object of class DEGGs generated from
#' `generate_subnetworks`
#' @param subgroup the subgroup of interest: the associated network will be
#' visualised
#' @importFrom magrittr %>%
#' @return a network plot with selectable nodes and links
#' @export
View_interactive_subnetwork <- function(DEGGs_object, subgroup,
                                        subgroup_variable){

  edges <- subnetworks_object[["subnetworks"]][[subgroup]]
  edges$id <- paste(edges$from, edges$to, sep = "-")
  rownames(edges) <- edges$id
  colnames(edges)[3] <- "width"
  nodes <- data.frame("id" = unique(c(edges$from, edges$to)),
                      "label" = unique(c(edges$from, edges$to)))

  # sever
  server <- shiny::shinyServer(function(input, output) {

    output$network <- visNetwork::renderVisNetwork({
      visNetwork::visNetwork(nodes, edges) %>%
        visNetwork::visEdges(arrows ="to") %>%
        visNetwork::visIgraphLayout() %>%
        visNetwork::visLayout(randomSeed = 12) %>% # to have always the same network
        visNetwork::visOptions(highlightNearest=TRUE,
                               nodesIdSelection = TRUE)  %>%
        visNetwork::visEvents(select = "function(data) {
                Shiny.onInputChange('current_nodes_selection', data.nodes);
                Shiny.onInputChange('current_edges_selection', data.edges);
                ;}")

    })

    #render data table restricted to selected nodes
    output$tbl <- DT::renderDT(
      edges %>%
        dplyr::filter(id %in% input$current_edges_selection),
      options = list(lengthChange = FALSE)
    )

    output$regressionPlot <-  shiny::renderPlot({
      req(input$current_edges_selection !="")
      print_regressions(gene_A = edges[input$current_edges_selection, "from"],
                        gene_B = edges[input$current_edges_selection, "to"],
                        DEGGs_object = subnetworks_object,
                        subgroup_variable = subgroup_variable)
    })

  })


  # user interface
  ui <- shiny::fluidPage(
    # dim 100vh , display = flex, type = column, grow = auto (in CSS)
    # in ogni riga puoi impostare un flexdraw (un peso) (see flexbox)
    # generate two rows,
    shiny::fluidRow(
      shiny::column(width = 12, visNetwork::visNetworkOutput("network",
                                                             height = "50vh"), # column widths in a fluidRow should sum to 12
                    shiny::column(width = 6, DT::DTOutput('tbl')),
                    shiny::column(width = 6, shiny::plotOutput('regressionPlot'))
      )
    )
  )

  # ui <- shinyUI(
  #     fluidPage(visNetworkOutput("network", height = "1000px"))
  # )
  #
  shiny::shinyApp(ui = ui, server = server)

}

