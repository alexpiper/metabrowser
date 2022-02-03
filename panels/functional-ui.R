funcplot <- fluidPage(outputCodeButton(withLoader(plotOutput("funcplot", height = 700))),
                     uiOutput("funcUI"))
