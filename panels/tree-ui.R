tree <- fluidPage(outputCodeButton(withLoader(plotOutput("tree", height = 700))),
                     uiOutput("treeUI"))
