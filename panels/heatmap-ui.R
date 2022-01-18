heatmap <- fluidPage(outputCodeButton(withLoader(plotOutput("heatmap", height = 700))),
                     uiOutput("heatmapUI"))
