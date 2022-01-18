barplot <- fluidPage(outputCodeButton(withLoader(plotOutput("barplot", height = 700))),
                     uiOutput("barplotUI"))
