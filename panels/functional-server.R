# Inputs

# Spreadsheet

# Annotation rank
output$funcAnnotInput <- renderUI({
  validate(need(physeq(), ""))
  fileInput(
    inputId = "fileAnnot",
    label = "Annotation file : ",
    placeholder = "annot.csv"
  )
})

# Annotation rank
output$funcYui <- renderUI({
  validate(need(physeq(), ""))
  radioButtons(
    inputId = "funcY",
    label = "What to put on Y axis : ",
    inline = TRUE,
    choices = list(
      "Relative abundance" = "ab",
      "Unique OTUs" = "otu"
    ),
    selected = "ab"
  )
})


output$funcAnnotRankUI <- renderUI({
  validate(need(physeq(), ""))
  radioButtons(
    "funcFilterRank",
    label = "Taxonomic rank used for annotation : ",
    choices = c("NULL" = 0, rank_names(physeq())),
    inline = TRUE
  )
})

output$funcGridUI <- renderUI({
  validate(need(physeq(), ""))
  selectInput("funcGrid",
              label = "Subplot : ",
              choices = c("..." = 0, sample_variables(physeq())))
})

output$funcXUI <- renderUI({
  validate(need(physeq(), ""))
  selectInput("funcX",
              label = "X : ",
              choices = c("..." = "Sample", sample_variables(physeq())))
})

output$funcUI <- renderUI({
  validate(need(physeq(), ""))
  box(
    title = "Setting : ",
    width = NULL,
    status = "primary",
    uiOutput("funcAnnotInput"),
    uiOutput("funcAnnotRankUI"),
    uiOutput("funcYUI"),
    uiOutput("funcGridUI"),
    uiOutput("funcXUI")
  )
})


# Display plot ---------------------------------------------------------
output$funcplot <- metaRender2(renderPlot, {
  validate(need(physeq(), "Requires an abundance dataset"))
  data <- physeq()

  # Import
  annot <- NULL
  observeEvent(input$fileAnnot, {
    annot <<- read_csv(input$fileAnnot$datapath)
  })
  
  # Set up facetting
  funcGrid <- if (!is.null(checkNull(input$funcGrid))) {
    metaExpr({
      facet_grid(..(paste(".", "~", input$funcGrid)), scales = "free_x", space = "free")
    })
  }
  # Create plot
  metaExpr({
    p <- plot_function(
      physeq = data,
      annot_table = ..(checkNull(annot)),
      annot_rank = ..(input$funcAnnotRankUI),
      x = ..(input$funcX),
      y = ..(input$funcY)
    )
    p + ..(funcGrid)
  })
})



# Display code ------------------------------------------------------------
observeEvent(input$funcplot_output_code,
             {
               displayCodeModal(
                 expandChain(
                   quote(library(phyloseq)),
                   "# Replace `data` with you own data.",
                   output$funcplot()
                 ), clip = NULL
               )
             }
)
