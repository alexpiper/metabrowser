alphaPlot <- fluidPage(outputCodeButton(withLoader(plotOutput("alphaPlot", height = 700))),
                       uiOutput("alphaPlotUI"))
alphaTable <- fluidPage(uiOutput("alphaTable"))
alphaAnova <- fluidPage(uiOutput("alphaAnovaUI"),
                        verbatimTextOutput("alphaAnova"),
                        box(
                          title = "Table of pairwise comparisons",
                          width = NULL,
                          status = "primary",
                          collapsible = TRUE,
                          collapsed = TRUE,
                          DT::DTOutput("anovaHSD"))
                        )
