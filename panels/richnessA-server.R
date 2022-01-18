output$alphaPlotUI <- renderUI({
  validate(need(physeq(), ""))
  box(
    title = "Setting : " ,
    width = NULL,
    status = "primary",
    checkboxGroupInput(
      "richnessMeasures",
      label = "Measures : ",
      choices = c(
        "Observed",
        "Chao1",
        "ACE",
        "Shannon",
        "Simpson",
        "InvSimpson",
        "Fisher"
      ),
      selected = c(
        "Observed",
        "Chao1",
        "ACE",
        "Shannon",
        "Simpson",
        "InvSimpson",
        "Fisher"
      ),
      inline = TRUE
    ),
    radioButtons(
      "richnessBoxplot",
      label = "Representation : ",
      choices = list(
        "Dots only" = 1,
        "Dots and boxplot" = 2,
        "Boxplot only" = 3
      ),
      selected = 2,
      inline = TRUE
    ),
    textInput("richnessTitle",
              label = "Title : ",
              value = "Alpha diversity graphics"),
    selectInput(
      "richnessX",
      label = "X : ",
      choices = c("..." = "samples", sample_variables(physeq()))
    ),
    selectInput(
      "richnessColor",
      label = "Color : ",
      choices = c("..." = 0, sample_variables(physeq()))
    ),
    selectInput(
      "richnessShape",
      label = "Shape : ",
      choices = c("..." = 0, sample_variables(physeq()))
    )
  )
})

output$alphaPlot <- metaRender2(renderPlot, {
  validate(need(physeq(), "Requires an abundance dataset"))
  data <- physeq()

  alphaPlot_boxplot <- if (input$richnessBoxplot >= 2) {
    metaExpr({geom_boxplot()})
  }

  alphaPlot_point <- if (input$richnessBoxplot <= 2) {
    metaExpr({geom_point()})
  }

  metaExpr({
    p <- plot_richness(
      physeq = data,
      x = ..(checkNull(input$richnessX)),
      color = ..(checkNull(input$richnessColor)),
      shape = ..(checkNull(input$richnessShape)),
      title = ..(checkNull(input$richnessTitle)),
      measures = ..(checkNull(input$richnessMeasures))
    )
    p + ..(alphaPlot_boxplot) + ..(alphaPlot_point)
  })
})

observeEvent(input$alphaPlot_output_code,
             {
               displayCodeModal(
                 expandChain(
                   quote(library(phyloseq)),
                   quote(library(phyloseq.extended)),
                   "# Replace `data` with you own data.",
                   output$alphaPlot()
                 ), clip = NULL
               )
             }
)

output$alphaTable <- renderUI({
  validate(need(physeq(), "Requires an abundance dataset"))
  box(
    title = "Alpha diversity estimation",
    width = NULL,
    status = "primary",
    beautifulTable(tibble::rownames_to_column(round(estimate_richness(physeq()), digits = 2), var = "SAMPLE"))
  )
})

output$alphaAnovaUI <- renderUI({
  validate(need(physeq(), "Requires an abundance dataset"))
  box(
    title = "Alpha diversity Anova analysis",
    width = NULL,
    status = "primary",
    selectInput(
      inputId = "anovaRichnessMeasure",
      label = "Richness measure : ",
      choices = c(
        "Observed",
        "Chao1",
        "ACE",
        "Shannon",
        "Simpson",
        "InvSimpson",
        "Fisher"
      ),
      selected = "Observed",
      multiple = FALSE
    ),
    selectInput(
      inputId = "anovaCovariate",
      label = "Covariate : ",
      choices = c(sample_variables(physeq()))
    )
  )
})

model <- reactive({
  validate(
    need(physeq(), "Requires an abundance dataset"),
    need(input$anovaRichnessMeasure, ""),
    need(input$anovaCovariate, "Select a covariate")
  )
  anova_data <- cbind(estimate_richness(physeq(), measures = input$anovaRichnessMeasure),
                      as(sample_data(physeq()), "data.frame"))
  formula <-  paste0(input$anovaRichnessMeasure, " ~ ", input$anovaCovariate)
  do.call("lm", list(formula = as.formula(formula), data = as.name("anova_data")))
})

output$alphaAnova <- renderPrint({
  validate(need(model(), ""))
  print(model())
  anova(model())
})

output$anovaHSD <- DT::renderDT({
  validate(
    need(model(), ""),
    need(physeq(), "Requires an abundance dataset"),
    need(input$anovaRichnessMeasure, ""),
    need(input$anovaCovariate, "Select a covariate")
  )
  var <- get_variable(physeq(), input$anovaCovariate)
  if (! (is.factor(var) || is.character(var))) {
    tibble::tibble(.rows = 0)
  } else {
    TukeyHSD(aov(model()))[[1]] %>% as_tibble(rownames = "Pair") %>%
      beautifulTable() %>%
      DT::formatSignif(columns = c("diff", "lwr", "upr", "p adj"), digits = 4)
  }
})
