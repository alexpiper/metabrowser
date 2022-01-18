output$barplotShowRankUI <- renderUI({
  validate(need(physeq(), ""))
  radioButtons(
    "barplotShowRank",
    label = "Taxonomic rank used for coloring : ",
    choices = c(rank_names(physeq()), "OTU"),
    selected = "Phylum",
    inline = TRUE
  )
})

output$barplotFilterRankUI <- renderUI({
  validate(need(physeq(), ""))
  radioButtons(
    "barplotFilterRank",
    label = "Taxonomic rank used for filtering : ",
    choices = c("NULL" = 0, rank_names(physeq())),
    inline = TRUE
  )
})

output$barplotTaxaUI <- renderUI({
  validate(need(physeq(), ""),
           need(input$barplotFilterRank, ""),
           need(input$barplotFilterRank!=0, ""))
  selectInput(
    "barplotTaxa",
    label = "Selected filter taxa : ",
    choices = unique(as.vector(tax_table(physeq())[, input$barplotFilterRank])),
    selected = TRUE
  )
})

output$barplotNbTaxaUI <- renderUI({
  validate(need(physeq(), ""))
  sliderInput(
    "barplotNbTaxa",
    label = "Number of sub-taxa : ",
    min = 1,
    #max = sum(tax_table(tax_glom(physeq(), rank_names(physeq())[1+as.integer(input$barplotFilterRank)]))[, as.integer(input$barplotFilterRank)]==input$barplotTaxa)
    max = 30,
    value = 10
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
    uiOutput("barplotFilterRankUI"),
    uiOutput("barplotTaxaUI"),
    uiOutput("barplotNbTaxaUI"),
    uiOutput("barplotShowRankUI"),
    uiOutput("barplotGridUI"),
    uiOutput("barplotXUI")
  )
})

output$barplot <- metaRender2(renderPlot, {
  validate(need(physeq(), "Requires an abundance dataset"),
           need(input$barplotShowRank, ""))
  data <- physeq()
  
  barplotGrid <- if (!is.null(checkNull(input$barplotGrid))) {
    metaExpr({
      facet_grid(..(paste(".", "~", input$barplotGrid)), scales = "free_x", space = "free")
    })
  }
  
  metaExpr({
    p <- plot_composition(
      physeq = data,
      taxaRank1 = ..(checkNull(input$barplotFilterRank)),
      taxaSet1 = ..(input$barplotTaxa),
      taxaRank2 = ..(input$barplotShowRank),
      numberOfTaxa = ..(input$barplotNbTaxa),
      x = ..(input$barplotX)
    )
    p + ..(barplotGrid)
  })
})

observeEvent(input$barplot_output_code,
             {
               displayCodeModal(
                 expandChain(
                   quote(library(phyloseq)),
                   quote(library(phyloseq.extended)),
                   "# Replace `data` with you own data.",
                   output$barplot()
                 ), clip = NULL
               )
             }
)
