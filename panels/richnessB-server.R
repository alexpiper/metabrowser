betaDist <- metaReactive2({
  data <- physeq()
  metaExpr({
    distance(data, method = ..(input$betaDistance))
    })
  }, varname = "beta.dist")

## MDS
output$betaMdsUI <- renderUI({
  validate(need(physeq(), ""))
  box(
    title = "Setting : " ,
    width = NULL,
    status = "primary",
    checkboxGroupInput(
      "betaMdsAxes",
      label = "Axes : ",
      choices = seq(10),
      selected = c(1, 2),
      inline = TRUE
    ),
    selectInput(
      "betaMdsMethod",
      label = "Method : ",
      selected = "MDS",
      choices = list("DCA (Detrended Correspondence Analysis)" = "DCA", 
                     "CCA (Constrained Correspondence Analysis)" = "CCA", 
                     "RDA (Redundancy Analysis)" = "RDA", 
                     "CAP (Constrained Analysis of Principal Coordinates)" = "CAP", 
                     "DPCoA (Double Principle Coordinate Analysis)" = "DPCoA", 
                     "NMDS (Non-metric MultiDimenstional Scaling)" = "NMDS", 
                     "MDS / PCoA (Principal Coordinate Analysis)" = "MDS")
    ),
    textInput("betaMdsTitle",
              label = "Title : ",
              value = "Samples ordination graphic"),
    selectInput(
      "betaMdsLabel",
      label = "Label : ",
      choices = c("..." = 0, sample_variables(physeq()))
    ),
    selectInput(
      "betaMdsCol",
      label = "Color : ",
      choices = c("..." = 0, sample_variables(physeq()))
    ),
    selectInput(
      "betaMdsShape",
      label = "Shape : ",
      choices = c("..." = 0, sample_variables(physeq()))
    ),
    selectInput(
      "betaMdsEllipse",
      label = "Ellipses : ",
      choices = c("..." = 0, sample_variables(physeq()))
    )
  )
})

output$betaMds <- metaRender2(renderPlot, {
  validate(need(physeq(), "Requires an abundance dataset"),
           need(length(input$betaMdsAxes) == 2, "Requires two projections axes"))
  data <- physeq()
  
  betaMdsEllipse <- if (!is.null(checkNull(input$betaMdsEllipse))) {
    metaExpr({
      stat_ellipse(aes_string(group = ..(input$betaMdsEllipse)))
    })
  }
  
  metaExpr({
    ord <- ordinate(data,
                    method = ..(input$betaMdsMethod),
                    distance = ..(betaDist())
                    )
    p <- plot_ordination(
      physeq = data,
      ordination = ord,
      type = "samples",
      axes = ..(as.numeric(input$betaMdsAxes)),
      color = ..(checkNull(input$betaMdsCol)),
      shape = ..(checkNull(input$betaMdsShape)),
      label = ..(checkNull(input$betaMdsLabel)),
      title = ..(checkNull(input$betaMdsTitle))
    )
    p <- p + ..(betaMdsEllipse)
    p + theme_bw()
  })
})

observeEvent(input$betaMds_output_code,
             {
               displayCodeModal(
                 expandChain(
                   quote(library(phyloseq)),
                   quote(library(phyloseq.extended)),
                   "# Replace `data` with you own data.",
                   output$betaMds()
                 ), clip = NULL
               )
             }
)

## Cluster
output$betaClusterUI <- renderUI({
  validate(need(physeq(), ""))
  box(
    title = "Setting : " ,
    width = NULL,
    status = "primary",
    selectInput("betaClusterMethod",
                label = "Method : ",
                choices = list("ward.D2", "ward.D", "single", "complete", "average", "mcquitty", "median", "centroid")),
    selectInput("betaClusterCol",
                label = "Color : ",
                choices = c("..." = 0, sample_variables(physeq())))
  )
})

output$betaCluster <- metaRender2(renderPlot, {
  validate(need(physeq(), "Requires an abundance dataset"))
  data <- physeq()
  
  metaExpr({
    p <- plot_clust(physeq = data,
                    dist = ..(betaDist()),
                    method = ..(input$betaClusterMethod),
                    color = ..(checkNull(input$betaClusterCol))
                    )
    p
  })
})

observeEvent(input$betaCluster_output_code,
             {
               displayCodeModal(
                 expandChain(
                   quote(library(phyloseq)),
                   quote(library(phyloseq.extended)),
                   "# Replace `data` with you own data.",
                   output$betaCluster()
                 ), clip = NULL
               )
             }
)

## MANOVA
output$betaManovaUI <- renderUI({
  validate(need(physeq(), ""))
  box(
    title = "Setting : " ,
    width = NULL,
    status = "primary",
    selectizeInput("betaManovaVar",
                   label = "Covariates (ordered list ; max=3) :",
                   multiple = TRUE,
                   choices = c(sample_variables(physeq())),
                   selected = NULL,
                   options = list(maxItems = 3)
                   ),
    checkboxGroupInput("betaManovaInteraction",
                       label = "Add interaction :",
                       choices = NULL,
                       inline = TRUE)#,
    # actionButton("betaManovaButton",
    #              "Compute"),
  )
})

observe({
  if (is.null(input$betaManovaVar)) {
    interactions <- NULL
  } else {
    interactions <- list(NULL,
                         paste(input$betaManovaVar, collapse = ":"),
                         c(paste(input$betaManovaVar[c(1, 2)], collapse = ":"),
                           paste(input$betaManovaVar[c(1, 3)], collapse = ":"),
                           paste(input$betaManovaVar[c(2, 3)], collapse = ":"),
                           paste(input$betaManovaVar, collapse = ":"))
                         )[[length(input$betaManovaVar)]]
  }
  
  updateCheckboxGroupInput(session, 
                           inputId = "betaManovaInteraction",
                           choices = interactions,
                           inline = TRUE
  )
})

formula <- reactive({as.formula(paste("betaDist()", "~", paste(c(input$betaManovaVar, input$betaManovaInteraction), collapse = "+")))})

output$betaManova <- renderPrint({
  validate(need(physeq(), "Requires an abundance dataset"))
  validate(need(input$betaManovaVar, "Requires at least one covariate"))

  sample_data <- as(sample_data(physeq()), "data.frame")
  
  manova <- vegan::adonis(
    formula = formula(),
    data = sample_data,
    perm = 9999)
  
  manova$call$formula <- formula()
  manova
})

## Heatmap
output$betaHeatmapUI <- renderUI({
  validate(need(physeq(), ""))
  box(
    title = "Setting : " ,
    width = NULL,
    status = "primary",
    selectInput("betaHeatmapOrder",
                label = "Sorting sample : ",
                choices = c("..." = 0, sample_variables(physeq()))
    ),
    textInput("betaHeatmapTitle",
              label = "Title : ",
              value = "Beta diversity heatmap")
  )
})

output$betaHeatmap <- metaRender2(renderPlot, {
  validate(need(physeq(), "Requires an abundance dataset"))
  data <- physeq()

  if (!is.null(checkNull(input$betaHeatmapOrder)))
  {
    metaExpr({
      beta <- reshape2::melt(as.matrix(..(betaDist())))
      colnames(beta) <- c("x", "y", "distance")
      new_factor = as.factor(get_variable(data, ..(input$betaHeatmapOrder)))
      variable_sort <- as.factor(get_variable(data, ..(input$betaHeatmapOrder))[order(new_factor)])
      L = levels(reorder(sample_names(data), as.numeric(new_factor)))
      beta$x <- factor(beta$x, levels = L)
      beta$y <- factor(beta$y, levels = L)
      palette <- scales::hue_pal()(length(levels(new_factor)))
      tipColor <- scales::col_factor(palette, levels = levels(new_factor))(variable_sort)
      p <- ggplot(beta, aes(x = x, y = y, fill = distance)) + geom_tile()
      p <- p + ggtitle(..(input$betaHeatmapTitle))
      p <- p + theme(
        axis.text.x = element_text(angle = 90, hjust = 1, color = tipColor),
        axis.text.y = element_text(color = tipColor),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
      p + scale_fill_gradient2()
      })
  } else {
    metaExpr({
      beta <- reshape2::melt(as.matrix(..(betaDist())))
      colnames(beta) <- c("x", "y", "distance")
      p <- ggplot(beta, aes(x = x, y = y, fill = distance)) + geom_tile()
      p <- p + ggtitle(..(input$betaHeatmapTitle))
      p <- p + theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
      p + scale_fill_gradient2()
    })
  }
})

observeEvent(input$betaHeatmap_output_code,
             {
               displayCodeModal(
                 expandChain(
                   quote(library(phyloseq)),
                   quote(library(phyloseq.extended)),
                   "# Replace `data` with you own data.",
                   output$betaHeatmap()
                 ), clip = NULL
               )
             }
)

## Network
output$betaNetworkUI <- renderUI({
  validate(need(physeq(), ""))
  box(
    title = "Setting : " ,
    width = NULL,
    status = "primary",
    sliderInput("betaNetworkMax",
                label = "Threshold : ",
                min = 0,
                max = 1,
                value = 0.7),
    checkboxInput("betaNetworkOrphan",
                  label = "Keep orphans",
                  value = TRUE),
    textInput("betaNetworkTitle",
              label = "Title : ",
              value = "Beta diversity network"),
    selectInput("betaNetworkCol",
                label = "Color : ",
                choices = c("..." = 0, sample_variables(physeq()))),
    selectInput("betaNetworkShape",
                label = "Shape : ",
                choices = c("..." = 0, sample_variables(physeq()))),
    selectInput("betaNetworkLabel",
                label = "Label : ",
                choices = c("..." = 0, "Sample name" = "value", sample_variables(physeq())))
    )
})

output$betaNetwork <- metaRender2(renderPlot, {
  validate(need(physeq(), "Requires an abundance dataset"))
  data <- physeq()
  
  metaExpr({
    g <- make_network(data,
                      distance = ..(betaDist()),
                      max.dist = ..(input$betaNetworkMax),
                      keep.isolates = ..(input$betaNetworkOrphan)
                      )
    p <- plot_network(g,
                      physeq = data,
                      color = ..(checkNull(input$betaNetworkCol)),
                      shape = ..(checkNull(input$betaNetworkShape)),
                      label = ..(checkNull(input$betaNetworkLabel)),
                      hjust = 2,
                      title = ..(checkNull(input$betaNetworkTitle))
                      )
    p
  })
})

observeEvent(input$betaNetwork_output_code,
             {
               displayCodeModal(
                 expandChain(
                   quote(library(phyloseq)),
                   quote(library(phyloseq.extended)),
                   "# Replace `data` with you own data.",
                   output$betaNetwork()
                 ), clip = NULL
               )
             }
)

## Table
output$betaTable <- renderUI({
  validate(need(physeq(), "Requires an abundance dataset"))
  box(
    title = "Distance table",
    width = NULL,
    status = "primary",
    beautifulTable(tibble::rownames_to_column(as.data.frame(round(as.matrix(betaDist()), digits = 2)), var = "SAMPLE"))
  )
})
