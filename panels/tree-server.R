output$treeUI <- renderUI({
  validate(need(phy_tree(physeq(), errorIfNULL = FALSE), ""))
  box(
    title = "Setting : " ,
    width = NULL,
    status = "primary",
    radioButtons(
      "treeRank",
      label = "Taxonomic rank captioned : ",
      choices = c(None = "",
                  rank_names(physeq()),
                  OTU = "taxa_names"),
      inline = TRUE
    ),
    sliderInput(
      "treeTopOtu",
      label = "Show the n most abundant OTU : ",
      min = 1,
      max = ntaxa(physeq()),
      value = 20
    ),
    sliderInput(
      "treeSize",
      label = "Abundance size based on a logarithm basis : ",
      min = 2,
      max = 10,
      value = 5
    ),
    checkboxInput("treeRadial", label = "Radial tree", value = FALSE),
    checkboxInput("treeSample", label = "Show samples", value = TRUE),
    textInput("treeTitle",
              label = "Title : ",
              value = "Phylogenetic tree"),
    selectInput(
      "treeCol",
      label = "Color : ",
      choices = c("..." = 0, sample_variables(physeq()), rank_names(physeq()))
    ),
    selectInput(
      "treeShape",
      label = "Shape : ",
      choices = c("..." = 0, sample_variables(physeq()), rank_names(physeq()))
    )
  )
})

output$tree <- metaRender2(renderPlot, {
  validate(need(physeq(), "Requires an abundance dataset"),
           need(phy_tree(physeq(), errorIfNULL = FALSE), "Requires a phylogenetic tree"))
  data <- physeq()
  
  treeRadial <- if (input$treeRadial) {
    metaExpr({
      coord_polar(theta = "y")
    })
  }
  
  metaExpr({
    data_select <- prune_taxa(names(sort(taxa_sums(data), decreasing = TRUE)[1:..(input$treeTopOtu)]), data)
    p <- plot_tree(
      physeq = data_select,
      method = ..(ifelse(input$treeSample, "sampledodge", "treeonly")),
      color = ..(checkNull(input$treeCol)),
      shape = ..(checkNull(input$treeShape)),
      size = "abundance",
      label.tips = ..(checkNull(input$treeRank)),
      sizebase = ..(checkNull(input$treeSize)),
      ladderize = "left",
      plot.margin = 0.1,
      title = ..(checkNull(input$treeTitle))
    )
    p + ..(treeRadial)
  })
})

observeEvent(input$tree_output_code,
             {
               displayCodeModal(
                 expandChain(
                   quote(library(phyloseq)),
                   "# Replace `data` with you own data.",
                   output$tree()
                 ), clip = NULL
               )
             }
)
