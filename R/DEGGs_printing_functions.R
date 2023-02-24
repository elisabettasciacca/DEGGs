#' Plot regressions for a specific gene-gene interaction
#'
#' @param deggs_object an object of class `deggs` generated from
#' `generate_subnetworks`
#' @param gene_A character. Name of the first gene  
#' @param gene_B character. Name of the second gene
#' @param legend_position position of the legend in the plot. It can be
#' specified by keyword or in any way which is accepted by `xy.coords` (defalut
#' "topright")
#' @importFrom methods is
#' @importFrom grDevices adjustcolor
#' @importFrom graphics abline axis boxplot legend mtext par points polygon text
#' @return base graphics plot showing gene-gene robust regressions for each 
#' subgroup. The p value of the interaction term of
#' gene A ~ gene B \* subgroup is reported on top.
#' @export
#' @seealso [`deggs-class`].
print_regressions <- function (deggs_object,
                               gene_A,
                               gene_B,
                               legend_position = "topright"){

  use_qvalues <- deggs_object@use_qvalues
  sig_var <- ifelse(use_qvalues, "q.value", "p.value")
  metadata <- deggs_object@metadata
  normalised_counts <- deggs_object@normalised_counts
  subgroups = deggs_object@subgroups
  subgroup_variable <- deggs_object@subgroup_variable
  regression_method <- deggs_object@regression_method

  # integrity checks
  if(!is(deggs_object, "deggs")){
    stop("deggs_object must be of class deggs")
  }

  if(!gene_A %in% rownames(deggs_object@normalised_counts))(
    stop("gene_A is not in rownames(normalised_counts)")
  )
  if(!gene_B %in% rownames(deggs_object@normalised_counts))(
    stop("gene_B is not in rownames(normalised_counts)")
  )

  if(is.null(subgroups)) {
    subgroups <- levels(metadata[, subgroup_variable])
  }

  subgroup_length <- length(subgroups)

  # prepare data frame
  df <- data.frame(t(normalised_counts[gene_A, ]),
                   t(normalised_counts[gene_B, ]),
                   metadata[, subgroup_variable],
                   check.names = FALSE)
  colnames(df)[3] <- subgroup_variable

  # compute gene-gene regression
  if(subgroup_length == 2){
    if(regression_method == "lm"){
      lmfit <- stats::lm(df[,2] ~ df[,1] * df[,3])
      # i.e.: gene_A ~ gene_B * subgroup
      p_interaction <- stats::coef(summary(lmfit))[4,4]
      fit <- lapply(subgroups, function(i) {
        x <- df[df[, subgroup_variable] == i, 1]
        y <- df[df[, subgroup_variable] == i, 2]
        if(length(x) > 0 || length(y) > 0) return(stats::lm(y ~ x))
      })
    }
    if(regression_method == "rlm"){
      robustfit <- MASS::rlm(df[,2] ~ df[,1] * df[,3])
      # i.e.: gene_A ~ gene_B * subgroup
      p_interaction <- sfsmisc::f.robftest(robustfit, var = 3)$p.value
      fit <- lapply(subgroups, function(i) {
        x <- df[df[, subgroup_variable] == i, 1]
        y <- df[df[, subgroup_variable] == i, 2]
        if(length(x) > 0 || length(y) > 0) return(MASS::rlm(y ~ x))
      })
    }
  }
  if(subgroup_length >= 3){
    # one-way ANOVA
    # i.e.: gene_A ~ gene_B * subgroup
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
  prefix <- ifelse(use_qvalues, "Padj", "P")
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

  all_interactions <- do.call(rbind, deggs_object@subnetworks)
  Padj_interaction <- all_interactions[which(all_interactions$from == gene_A &
                                            all_interactions$to == gene_B),
                                       sig_var]

  sig_interaction <- ifelse(deggs_object@use_qvalues,
                            Padj_interaction,
                            p_interaction)

  op <- par(mar = c(5.2, 6, 3.3, 5), xpd = TRUE)
  plot(df[,1], df[,2], type = 'n', bty = 'l', las = 1, cex.axis = 1.1,
       font.main = 1, cex.lab = 1.3, xlab = colnames(df)[1],
       ylab = colnames(df)[2])
  op <- par(xpd = FALSE)

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
  mtext(bquote(paste(.(prefix)["interaction"]*"=",
                     .(format(sig_interaction, digits = 2)))),
        cex = 1.2, side = 3, adj = 0.04)
  legend(x = legend_position, legend = subgroups,
         col = col, lty = 1,
         bty = "o", box.lty = 0,
         cex = 0.8)
  on.exit(par(op))
}

#' Boxplots of expression levels of single nodes (genes) per subgroups
#'
#' This function is used only by `View_interactive_subnetwork`
#'
#' @param gene gene name  (must be contained in `rownames(normalised_counts)`)
#' @param deggs_object  an object of class `deggs` generated from
#' `generate_subnetworks`
#' @return the gene boxplot
node_boxplot <- function(gene,
                         deggs_object) {

  metadata <- deggs_object@metadata
  normalised_counts <- deggs_object@normalised_counts
  subgroup_variable <- deggs_object@subgroup_variable
  x <- metadata[, subgroup_variable]
  y <- as.numeric(normalised_counts[gene, ])
  col  <- viridis::viridis(n = nlevels(metadata[, subgroup_variable]))
  cols <- col[metadata[, subgroup_variable]]

  op <- par(bty='l', mar = c(5.2, 6, 3.3, 5))
  boxplot(y ~ x, outline = FALSE, whisklty = 1, medlwd = 2, cex.axis = 1.1,
          col = NA, cex.lab = 1.3, ylab = gene, las = 2, boxwex=.5,
          ylim = range(y), xaxt = 'n', xlab = '')
  points(jitter(as.numeric(x), amount = 0.15), y, pch = 20,
         col = adjustcolor(cols, alpha.f = 0.7), cex = 2)
  boxplot(y ~ x, add = TRUE, col = NA, outline = FALSE, xaxt = 'n', yaxt = 'n',
          whisklty = 1, medlwd = 2, las = 2, xlab = '', boxwex=.5)
  xtick <- levels(x)
  axis(1, 1:length(xtick), labels = FALSE)
  text(x = 1:length(xtick), y = par()$usr[3]-0.04*(par()$usr[4]-par()$usr[3]),
       labels = xtick, srt = 30, adj = 1, xpd = NA, cex = 1.2)
  par(op)
}


#' Interactive shiny app to visualise subnetworks
#'
#' Explore subnetworks and interactively select regression and box plots
#'
#' @param deggs_object an object of class `deggs` generated from
#' `generate_subnetworks`
#' @import igraph
#' @import knitr
#' @import rmarkdown
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @return a shiny interface showing networks with selectable nodes and links
#' @export
View_interactive_subnetwork <- function(deggs_object, host = NULL, port = NULL){
  
  use_qvalues <- deggs_object@use_qvalues
  sig_var <- ifelse(use_qvalues, "q.value", "p.value")

  ################# server ########################
  server <- function(input, output, session) {

    outVar = shiny::reactive({
      if(is.data.frame(deggs_object@subnetworks[[input$subgroup]]))(
        deggs_object@subnetworks[[input$subgroup]]
      ) else (
        data.frame("from" = NA, "to" = NA, "p.value" = 1, "q.value" = 1)
      )
    })
    nodes_selection <- shiny::reactiveValues(current_node = NULL)

    # Set up reactive table
    edges <- shiny::reactive({
      if(is.data.frame(deggs_object@subnetworks[[input$subgroup]])){
        edges <- deggs_object@subnetworks[[input$subgroup]]
        edges <- edges[edges[, sig_var] < input$slider, ]
        edges$id <- rownames(edges)

        if(use_qvalues) {
          edges$`q value` <- formatC(edges$q.value, format = "e", digits = 3)
          edges <- edges[order(edges$q.value),]
        } else {
          edges$`p value` <- formatC(edges$p.value, format = "e", digits = 3)
          edges <- edges[order(edges$p.value),]
          }
        
        if (length(input$current_edges_selection) == 0)(
          DT::datatable(edges,
                        options = list(lengthChange = FALSE, scrollX = T,
                                       columnDefs = list(list(visible=FALSE,
                                                              targets=c(3:4)))),
                        rownames = TRUE)
        ) else (
          DT::datatable(edges %>%
                          dplyr::filter(.data$id %in% input$current_edges_selection),
                        options = list(lengthChange = FALSE, scrollX = T,
                                       columnDefs = list(list(visible=FALSE,
                                                              targets=c(3:4)))),
                        rownames = TRUE)
        )
        }
      })

    # Network
    output$network <- visNetwork::renderVisNetwork({
      if(is.data.frame(deggs_object@subnetworks[[input$subgroup]])){
        edges <- deggs_object@subnetworks[[input$subgroup]]
        edges$id <- rownames(edges)

        # Set up tooltip
        prefix <- ifelse(use_qvalues == TRUE, "Padj=", "P=")
        edges$title <- paste0(prefix,
                              formatC(edges[, sig_var], format = "e", digits = 2) )

        # Set up edges width
        # normalise p value between 0 and 1
        edges$width <- edges[, sig_var] - min(edges[, sig_var]) /
          (max(edges[, sig_var]) - min(edges[, sig_var]))
        # invert values (and multiply by 4 to increase width)
        edges$width <- (1 - edges$width) * 4

        # Set up edges color
        edges$color <- ifelse(edges[, sig_var] < 0.05, "royalblue", "gray")

        # Slider
        edges <- edges[edges[, sig_var] < input$slider, ]

        nodes <- data.frame("id"    = unique(c(edges$from, edges$to)),
                            "label" = unique(c(edges$from, edges$to)),
                            "title" = unique(c(edges$from, edges$to)))

        visNetwork::visNetwork(nodes, edges) %>%
          visNetwork::visIgraphLayout(physics = TRUE, smooth = TRUE, type = "full") %>%
          visNetwork::visNodes(color = list("border" = 'white'),
                               font = list("size" = 16)) %>%
          visNetwork::visEdges(arrows ="to",
                               color = list("inherit" = FALSE)) %>%

          visNetwork::visLayout(randomSeed = 12) %>% # to have always the same network
          visNetwork::visLegend(useGroups = FALSE,
                                position = 'right',
                                addEdges = data.frame(color = c("gray", "royalblue"),
                                                      label = c("not significant", "significant"),
                                                      font.align = "top",
                                                      font.size = 9),
                                stepY = 50, width = 0.15,
                                zoom = FALSE) %>%
          visNetwork::visOptions(highlightNearest = TRUE)  %>%
          visNetwork::visInteraction(hover = TRUE, tooltipDelay = 20)  %>%
          visNetwork::visEvents(select = "function(data) {
                Shiny.onInputChange('current_nodes_selection', data.nodes);
                Shiny.onInputChange('current_edges_selection', data.edges);
                ;}")
      } else {
        nodes <- data.frame()
        edges <- data.frame()
        network.title = "No specific gene-gene interaction for this subgroup"
        visNetwork::visNetwork(nodes, edges, main = network.title)
      }
    })


    # Table
    output$tbl <- DT::renderDT({
      edges()})


    # Plot
    output$edge_or_node_plot <-  shiny::renderPlot({
      edges <- deggs_object@subnetworks[[input$subgroup]]
    # try(
      if(is.null(input$current_nodes_selection) &
         length(input$current_edges_selection) == 1) {

          print_regressions(gene_A = edges[input$current_edges_selection, "from"],
                            gene_B = edges[input$current_edges_selection, "to"],
                            deggs_object = deggs_object)
            } else {
              shiny::req(input$current_nodes_selection != "")
              node_boxplot(input$current_nodes_selection, 
                           deggs_object = deggs_object)
            }
      # ,silent = TRUE)
    })

    # Highligh the searched node in the network
    shiny::observe({
      if(is.data.frame(deggs_object@subnetworks[[input$subgroup]])){
        if(input$searchButton > 0){
          shiny::isolate({
            edges <- deggs_object@subnetworks[[input$subgroup]]
            nodes <- data.frame("id"    = unique(c(edges$from, edges$to)),
                                "label" = unique(c(edges$from, edges$to)),
                                "title" = unique(c(edges$from, edges$to)))

            nodes_selection$current_node <- nodes[grep(input$searchText, nodes$label,
                                                       ignore.case = TRUE), "id"]
            visNetwork::visNetworkProxy("network") %>%
              visNetwork::visSelectNodes(id  = nodes_selection$current_node )
          })
        }
      }
    })

    shiny::observe({
      shiny::updateSliderInput(session, inputId = "slider",
                        max = max(round(outVar()[, sig_var], digits = 3))
      )})

  }

  ################# UI ########################
  ui <- shiny::fluidPage(
    shiny::titlePanel("DEGGs subnetworks"),

    shiny::sidebarLayout(
      shiny::sidebarPanel(
        width = 4, # this will leave more space for the network

        # dropdown menu for subgroup
        shiny::selectInput(inputId = "subgroup", label = "Subgroup",
                           choices = levels(deggs_object@metadata[, deggs_object@subgroup_variable]),
                           selected = levels(deggs_object@metadata[, deggs_object@subgroup_variable])[1]),

       # Slider
        shiny::sliderInput("slider",
                    label = ifelse(use_qvalues, "q values", "p values"),
                    min = 0.01,
                    max = 10,   # fake, it will be updated by updateSliderInput
                    value = 0.05, step = 0.01),

        # Table
        shiny::tags$div(DT::DTOutput('tbl'), style = "font-size: 75%"),

        # Searchbox
        shinydashboard::sidebarSearchForm(textId = "searchText",
                                          buttonId = "searchButton",
                                          label = "Search node..."),

        # Plots
        shiny::tags$div(shiny::plotOutput('edge_or_node_plot'))
      ),

      shiny::mainPanel(
        visNetwork::visNetworkOutput("network",
                                     height = "700px", width = "800px")
      )
    ),
    shiny::tags$script(shiny::HTML("$(function() {
              $('#searchText').keypress(function(e) {
              if (e.which == 13) {
                $('#searchButton').click();
              }
              });
    });"))
   )
  shiny::shinyApp(ui = ui, server = server, options = list(host = host, port = port))
}
