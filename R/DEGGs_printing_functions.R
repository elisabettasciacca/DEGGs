#' View subnetworks in network3D style
#'
#' View network3D style networks for a specific subgroup
#'
#' @param deggs_object an object of class DEGGs generated from
#' `generate_subnetworks`
#' @param subgroup the subgroup of interest. The associated network will be
#' visualised
#' @return a network3D plot of the selected subnetwork
#' @export
print_force_network <- function(deggs_object, subgroup){
  network <- deggs_object@subnetworks[[subgroup]]
  nodes <- data.frame("name" = unique(c(network$from, network$to)))
  nodes$ID <- 0:(nrow(nodes)-1) # assign sequential IDs to nodes
  network <- base::merge(nodes, network, by.x = "name", by.y = "from")
  network <- dplyr::rename(network, source = ID)
  network <- base::merge(nodes, network, by.x = "name", by.y = "to")
  network <- dplyr::rename(network, target = ID)
  #network <- dplyr::rename(network, value = p.value)
  network$value <- 1
  network$name   <- NULL
  network$name.y <- NULL
  nodes$group <- 1
  
  
  networkD3::forceNetwork(Links = network, Nodes = nodes, 
                          Source = 'source', 
                          Target = 'target', 
                          NodeID = 'name', 
                          Group = "group", 
                          Value = "value",
                          arrows = TRUE,
                          linkColour = ifelse(network$p.value >= 0.05, "gray", "red"),
                          #linkWidth = networkD3::JS("function(d) { return d; }"),
                          opacity = 1,
                          opacityNoHover = 0.8, 
                          zoom = TRUE)
}


#' Print the regression graph of a specific genes
#'
#' @param deggs_object an object of class DEGGs generated from
#' `generate_subnetworks`
#' @param gene_A A string of gene selected
#' @param gene_B A string of gene selected
#' @param legend_offset optional numeric offset to horizontally move the legend
#' @return base graphics plot showing gene-gene regressions for each subgroup
#' while the pvalue of the interaction term of
#' *gene A ~ gene B \* subgroup* is reported on top
#' @export
print_regressions <- function (deggs_object, 
                               gene_A, 
                               gene_B, 
                               legend_offset = -0.4){
  
  metadata <- deggs_object@metadata
  normalised_counts <- deggs_object@normalised_counts
  subgroups = deggs_object@subgroups
  subgroup_variable <- deggs_object@subgroup_variable
  regression_method <- deggs_object@regression_method
  
  if(class(deggs_object) != "Deggs"){
    stop("The deggs_object must be of deggs class")
  }
  
  
  if(is.null(subgroups)) {
    subgroups <- levels(metadata[, subgroup_variable])
  }
  
  subgroup_length <- length(subgroups)
  
  # check if the genes selected exist
  is_gene_A_present <- FALSE
  is_gene_B_present <- FALSE
  
  lapply(subgroups, function(subgroup){
    if(gene_A %in% deggs_object@subnetworks[
      which(names(deggs_object@subnetworks) == subgroup)]
      [[1]]$from){
      
      is_gene_A_present <<- TRUE
      
    }
    if(gene_B %in% deggs_object@subnetworks[
      which(names(deggs_object@subnetworks) == subgroup)]
      [[1]]$to){
      
      is_gene_B_present <<- TRUE
      
    }
  })
  
  
  if(is_gene_A_present == FALSE){
    stop("The gene A selected is not present in the network")
  }
  else if(is_gene_B_present == FALSE){
    stop("The gene B selected is not present in the network")
  }
  
  # prepare data frame
  df <- data.frame(t(normalised_counts[gene_A, ]), 
                   t(normalised_counts[gene_B, ]), 
                   metadata[, subgroup_variable], 
                   check.names = FALSE)
  colnames(df)[3] <- subgroup_variable 
  
  # compute gene-gene regression
  if(subgroup_length == 2){
    if(regression_method == "lm"){
      lmfit <- lm(df[,2] ~ df[,1] * df[,3])
      # i.e.: gene_A ~ gene_B * subgroup
      p_interaction <- coef(summary(lmfit))[4,4]
      fit <- lapply(subgroups, function(i) {
        x <- df[df[, subgroup_variable] == i, 1]
        y <- df[df[, subgroup_variable] == i, 2]
        if(length(x) > 0 || length(y) > 0) return(lm(y ~ x))
      })
    }
    if(regression_method == "rlm"){
      robustfit <- rlm(df[,2] ~ df[,1] * df[,3])
      # i.e.: gene_A ~ gene_B * subgroup
      p_interaction <- f.robftest(robustfit, var = 3)$p.value
      fit <- lapply(subgroups, function(i) {
        x <- df[df[, subgroup_variable] == i, 1]
        y <- df[df[, subgroup_variable] == i, 2]
        if(length(x) > 0 || length(y) > 0) return(rlm(y ~ x))
      })
    }
  }
  if(subgroup_length >= 3){   
    # one-way ANOVA
    # i.e.: gene_A ~ gene_B * subgroup
    res_aov <- aov(df[,2] ~ df[,1] * df[,3], data = df)
    p_interaction <- summary(res_aov)[[1]][["Pr(>F)"]][3]
    fit <- lapply(subgroups, function(i) {
      x <- df[df[, subgroup_variable] == i, 1]
      y <- df[df[, subgroup_variable] == i, 2]
      if(regression_method == "lm"){
        if(length(x) > 0 || length(y) > 0) return(lm(y ~ x));
      }
      if(regression_method == "rlm"){
        if(length(x) > 0 || length(y) > 0) return(rlm(y ~ x))
      }
      
    })
  }
  
  # Plot 
  col <- viridis::viridis(n = subgroup_length)
  x_adj <- (max(df[,1], na.rm = T) - min(df[,1], na.rm = T)) * 0.05
  new_x <- seq(min(df[,1], na.rm = T) - x_adj, 
               max(df[,1], na.rm=T) + x_adj, 
               length.out=100)
  
  # prediction of the fitted model
  preds <- lapply(fit,
                  function(i) stats::predict(i, 
                                             newdata = data.frame(x = new_x), 
                                             interval = 'confidence'))
  
  par(mar = c(5.2, 6, 3.3, 10.7), xpd = TRUE)
  plot(df[,1], df[,2], type = 'n', bty = 'l', las = 1, cex.axis = 1.1, 
       font.main = 1, cex.lab = 1.3, xlab = colnames(df)[1], 
       ylab = colnames(df)[2])
  par(xpd = FALSE)
  
  
  cols <- col[df[, subgroup_variable]]
  pch <- c(16 : (16 + subgroup_length -1))[df[, subgroup_variable]]
  for(i in 1:subgroup_length) {
    
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
  
  legend("topright", bty = "n", inset=c(legend_offset, 0), legend = subgroups, 
         xpd=TRUE, col = col, lty = 1, cex = 1.2)
  
  
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
  edges$title <- formatC(edges$p.value, format = "e", digits = 2)
  # normalise p value between 0 and 1
  edges$width <- edges$p.value - min(edges$p.value) /
    (max(edges$p.value) - min(edges$p.value))
  # invert values (and multiply by 4 to increase width)
  edges$width <- (1 - edges$width) * 4
  
  nodes <- data.frame("id" = unique(c(edges$from, edges$to)),
                      "label" = unique(c(edges$from, edges$to)),
                      "title" = unique(c(edges$from, edges$to)))
  
  # sever
  server <- shiny::shinyServer(function(input, output) {
    
    output$network <- visNetwork::renderVisNetwork({
      visNetwork::visNetwork(nodes, edges) %>%
        visNetwork::visIgraphLayout(physics = TRUE, smooth = TRUE, type = "full") %>%
        visNetwork::visNodes(color = list("border" = 'white'),
                             font = list("size" = 16)) %>%
        visNetwork::visEdges(arrows ="to",
                             color = list("inherit" = FALSE)) %>%
        
        visNetwork::visLayout(randomSeed = 12) %>% # to have always the same network
        visNetwork::visOptions(highlightNearest = TRUE,
                               nodesIdSelection = TRUE)  %>%
        visNetwork::visInteraction(hover = TRUE, tooltipDelay = 20)  %>%
        visNetwork::visEvents(select = "function(data) {
                Shiny.onInputChange('current_nodes_selection', data.nodes);
                Shiny.onInputChange('current_edges_selection', data.edges);
                ;}")
    })
    
    # render data table restricted to selected nodes
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
                                                             height = "60vh"), # column widths in a fluidRow should sum to 12
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
