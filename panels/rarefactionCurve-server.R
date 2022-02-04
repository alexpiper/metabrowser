output$rarefactionCurveUI <- renderUI({
  validate(need(physeq(), ""))
  box(
    title = "Setting : " ,
    width = NULL,
    status = "primary",
    checkboxInput("rarefactionMin", label = "Show min sample threshold", value = FALSE),
    textInput("rarefactionTitle",
              label = "Title : ",
              value = "Rarefaction curves"),
    selectInput(
      "rarefactionColor",
      label = "Color : ",
      choices = c("..." = 0, sample_variables(physeq()))
    ),
    selectInput(
      "rarefactionLabel",
      label = "Label : ",
      choices = c("..." = 0, sample_variables(physeq()))
    ),
    selectInput(
      "rarefactionGrid",
      label = "Subplot : ",
      choices = c("..." = 0, sample_variables(physeq()))
    )
  )
})

output$rarefactionCurve <- metaRender2(renderPlot, {
  validate(need(physeq(), "Requires an abundance dataset"))
  data <- physeq()

  rarefactionMin <- if (input$rarefactionMin) {
    metaExpr({
      geom_vline(xintercept = min(sample_sums(data)), color = "gray60")
    })
  }
  
  rarefactionGrid <- if (!is.null(checkNull(input$rarefactionGrid))) {
    metaExpr({
      facet_grid(..(paste(".", "~", input$rarefactionGrid)))
    })
  }
  
  metaExpr({
    p <- ggrare(
      physeq = data,
      step = 100,
      color = ..(checkNull(input$rarefactionColor)),
      label = ..(checkNull(input$rarefactionLabel)),
      se = FALSE
      )
    p <- p + ..(rarefactionMin)
    p <- p + ..(rarefactionGrid)
    p + ggtitle(..(input$rarefactionTitle))
  })
})

observeEvent(input$rarefactionCurve_output_code,
             {
               displayCodeModal(
                 expandChain(
                   quote(library(phyloseq)),
                   "# Replace `data` with you own data.",
                   output$rarefactionCurve()
                 ), clip = NULL
               )
             }
)
