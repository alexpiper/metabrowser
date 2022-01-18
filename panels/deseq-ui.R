deseq <- fluidPage(outputCodeButton(withLoader(plotOutput("deseq", height = 700, brush = brushOpts(id = "deseq_brush")))),
                   box(
                     title = "Table of Brushed OTUs",
                     width = NULL,
                     status = "primary",
                     collapsible = TRUE,
                     collapsed = FALSE,
                     DT::DTOutput("deseqBrushed")),
                   box(
                     title = "Table of OTUs with significant effect",
                     width = NULL,
                     status = "primary",
                     collapsible = TRUE,
                     collapsed = TRUE,
                     DT::DTOutput("deseqTable")),
                   uiOutput("deseqUI"))
