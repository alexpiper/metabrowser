output$heatmapUI <- renderUI({
  validate(need(physeq(), ""))
  box(
    title = "Setting : " ,
    width = NULL,
    status = "primary",
    textInput("heatmapTitle",
              label = "Title : ",
              value = "Taxa heatmap by samples"),
    radioButtons(
      inputId = "plot_rank",
      label = "Taxonomic rank for y axis label : ",
      inline = TRUE,
      choices = c(
        phyloseq::rank_names(subset_physeq(), errorIfNULL = FALSE),
        "OTU"
      ),
      selected = phyloseq::rank_names(subset_physeq(), errorIfNULL = FALSE)[length(phyloseq::rank_names(subset_physeq(), errorIfNULL = FALSE))]
    ),
    selectInput(
      "heatmapGrid",
      label = "Subplot : ",
      choices = c("..." = 0, sample_variables(physeq()))
    ),
    selectInput(
      "heatmapX",
      label = "X : ",
      choices = c("..." = 0, sample_variables(physeq()))
    ),
    sliderInput(
      "heatmapTopOtu",
      label = "Show the n most abundant OTU : ",
      min = 1,
      max = ntaxa(physeq()),
      value = 250
    ),
    selectInput(
      "heatmapDist",
      label = "Distance : ",
      selected = "bray",
      choices = list(
        "bray",
        "jaccard",
        "unifrac",
        "wunifrac",
        "dpcoa",
        "jsd",
        "euclidean"
      )
    ),
    selectInput(
      "heatmapMethod",
      label = "Method : ",
      selected = "NMDS",
      choices = list(
        "NMDS",
        "ward.D2",
        "ward.D",
        "single",
        "complete",
        "average",
        "mcquitty",
        "median",
        "centroid"
      )
    )
  )
})

output$heatmap <- metaRender2(renderPlot, {
  validate(need(physeq(), "Requires an abundance dataset"))
  data <- physeq()
  
  heatmapGrid <- if (!is.null(checkNull(input$heatmapGrid))) {
    metaExpr({
      facet_grid(..(paste(".", "~", input$heatmapGrid)), scales = "free_x", space = "free")
    })
  }
  
  taxa_label <- NULL
   if(!input$plot_rank == "OTU"){
     taxa_label <- input$plot_rank
   } 
  
  metaExpr({
    data_select <- prune_taxa(names(sort(taxa_sums(data), decreasing = TRUE)[1:..(input$heatmapTopOtu)]), data)
    p <- plot_heatmap(
      physeq = data_select,
      distance = ..(input$heatmapDist),
      method = ..(input$heatmapMethod),
      title = ..(checkNull(input$heatmapTitle)),
      taxa.label = ..(checkNull(taxa_label)),
      sample.order = ..(checkNull(input$heatmapX)),
      low = "yellow",
      high = "red",
      na.value = "white"
    ) #+
     # scale_fill_viridis_c(na.value = "white")
    p + ..(heatmapGrid)
  })
})


observeEvent(input$heatmap_output_code,
             {
               displayCodeModal(
                 expandChain(
                   quote(library(phyloseq)),
                   "# Replace `data` with you own data.",
                   output$heatmap()
                 ), clip = NULL
               )
             }
)
