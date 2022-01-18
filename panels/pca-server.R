output$pcaUI <- renderUI({
  validate(need(physeq(), ""))
  box(
    title = "Setting : " ,
    width = NULL,
    status = "primary",
    checkboxGroupInput(
      "pcaSetting",
      label = "PCA setting",
      choices = list(
        "Center" = "center",
        "Scale" = "scale"
      ),
      selected = c("norm", "sqrt", "center", "scale"),
      inline = TRUE
    ),
    radioButtons(
      "pcaType",
      label = "Type of graph : ",
      choices = list(
        "Biplot of individuals and variables" = "biplot",
        "Graph of individuals" = "ind",
        "Graph of variables" = "var"
      ),
      selected = "biplot",
      inline = TRUE
    ),
    checkboxGroupInput(
      "pcaAxes",
      label = "Axes : ",
      choices = seq(10),
      selected = c(1, 2),
      inline = TRUE
    ),
    textInput("pcaTitle",
              label = "Title : ",
              value = "Principal Component Analysis"),
    h4(strong("Individuals ( = Samples)")),
    checkboxGroupInput(
      "pcaGeomInd",
      label = "Geometry for individuals  : ",
      choices = c("point", "text"),
      selected = c("point", "text"),
      inline = TRUE
    ),
    selectInput(
      "pcaHabillage",
      label = "Group : ",
      choices = c("..." = 0, sample_variables(physeq()))
    ),
    checkboxInput("pcaEllipse",
                  label = "Add ellipses",
                  value = FALSE),
    h4(strong("Variables ( = OTU)")),
    checkboxGroupInput(
      "pcaGeomVar",
      label = "Geometry for variables : ",
      choices = c("arrow", "text"),
      selected = c("arrow", "text"),
      inline = TRUE
    ),
    sliderInput(
      "pcaSelect",
      label = "Represent top contrib variables : ",
      min = 1,
      max = ntaxa(physeq()),
      value = 50,
      step = 1
    )
  )
})

output$pca <- metaRender2(renderPlot, {
  validate(
    need(physeq(), "Requires an abundance dataset"),
    need(length(input$pcaAxes) == 2, "Requires two projections axes"))
  data <- physeq()
  
  pca <- metaExpr({
    data_matrix <- as.data.frame(t(otu_table(data)))
    pca <- prcomp(data_matrix[colSums(data_matrix) != 0],
                  center = ..("center" %in% input$pcaSetting),
                  scale = ..("scale" %in% input$pcaSetting)
                  )
    })
  
  habillage <- if (!is.null(checkNull(input$pcaHabillage))) {
    metaExpr({
      get_variable(data, ..(input$pcaHabillage))
    })
  } else metaExpr({"none"})
  
  pcaType <- 
    switch(input$pcaType,
           "biplot" = metaExpr({
             fviz_pca_biplot(pca,
                             axes = ..(as.numeric(input$pcaAxes)),
                             geom.ind = ..(c(input$pcaGeomInd, "")),
                             geom.var = ..(c(input$pcaGeomVar, "")),
                             habillage = habillage,
                             invisible = "quali",
                             addEllipses = ..(input$pcaEllipse),
                             title = ..(input$pcaTitle),
                             select.var = ..(list(contrib = input$pcaSelect))
                             )
             }),
           "ind" = metaExpr({
             fviz_pca_ind(pca,
                          axes = ..(as.numeric(input$pcaAxes)),
                          geom.ind = ..(c(input$pcaGeomInd, "")),
                          habillage = habillage,
                          invisible = "quali",
                          addEllipses = ..(input$pcaEllipse),
                          title = ..(input$pcaTitle)
                          )
             }),
           "var" = metaExpr({
             fviz_pca_var(pca,
                          axes = ..(as.numeric(input$pcaAxes)),
                          geom.var = ..(c(input$pcaGeomVar, "")),
                          invisible = "quali",
                          title = ..(input$pcaTitle),
                          select.var = ..(list(contrib = input$pcaSelect))
                          )
             })
    )

  metaExpr({
    ..(pca)
    habillage <- ..(habillage)
    p <- ..(pcaType)
    p + theme_bw()
  })
})

observeEvent(input$pca_output_code,
             {
               displayCodeModal(
                 expandChain(
                   quote(library(phyloseq)),
                   quote(library(phyloseq.extended)),
                   quote(library(factoextra)),
                   "# Replace `data` with you own data.",
                   output$pca()
                 ), clip = NULL
               )
             }
)
