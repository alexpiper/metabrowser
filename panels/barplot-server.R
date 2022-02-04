output$barplotShowRankUI <- renderUI({
  validate(need(physeq(), ""))
  radioButtons(
    "barplotShowRank",
    label = "Taxonomic rank used for coloring : ",
    choices = c(rank_names(physeq()), "OTU"),
    selected = "Genus",
    inline = TRUE
  )
})

output$barplotNbTaxaUI <- renderUI({
  validate(need(physeq(), ""))
  sliderInput(
    "barplotNbTaxa",
    label = "Show the n most abundant OTU  : ",
    min = 1,
    max = get_max_taxa(ps = physeq(), filt_rank = input$barplotFilterRank, filt_name =  input$barplotTaxa),
    value = get_max_taxa(ps = physeq(), filt_rank = input$barplotFilterRank, filt_name =  input$barplotTaxa)
  )
})

output$barplotGridUI <- renderUI({
  validate(need(physeq(), ""))
  selectInput("barplotGrid",
              label = "Subplot : ",
              choices = c("..." = 0, sample_variables(physeq())))
})

output$barplotXUI <- renderUI({
  validate(need(physeq(), ""))
  selectInput("barplotX",
              label = "X : ",
              choices = c("..." = "Sample", sample_variables(physeq())))
})

output$barplotUI <- renderUI({
  validate(need(physeq(), ""))
  box(
    title = "Setting : ",
    width = NULL,
    status = "primary",
    uiOutput("barplotShowRankUI"),
    uiOutput("barplotNbTaxaUI"),
    uiOutput("barplotGridUI"),
    uiOutput("barplotXUI")
  )
})


# Display plot ---------------------------------------------------------
output$barplot <- metaRender2(renderPlot, {
  validate(need(physeq(), "Requires an abundance dataset"),
           need(input$barplotShowRank, ""))
  data <- physeq()

  # Set up facetting
  barplotGrid <- if (!is.null(checkNull(input$barplotGrid))) {
    metaExpr({
      facet_grid(..(paste(".", "~", input$barplotGrid)), scales = "free_x", space = "free")
    })
  }
  # Create plot
  metaExpr({
    p <- plot_composition(
      physeq = data,
      taxaRank1 = NULL,
      taxaSet1 = NULL,
      taxaRank2 = ..(input$barplotShowRank),
      fill = ..(input$barplotShowRank),
      numberOfTaxa = ..(input$barplotNbTaxa),
      x = ..(input$barplotX)
    ) +
      scale_x_discrete(expand=c(0,0))+
      scale_y_continuous(expand=c(0,0))
    p + ..(barplotGrid)
  })
})


# Display code ------------------------------------------------------------
observeEvent(input$barplot_output_code,
             {
               displayCodeModal(
                 expandChain(
                   quote(library(phyloseq)),
                   "# Replace `data` with you own data.",
                   output$barplot()
                 ), clip = NULL
               )
             }
)
