pca <- fluidPage(outputCodeButton(withLoader(plotOutput("pca", height = 700))),
                 uiOutput("pcaUI"))
