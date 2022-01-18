betaMds <- fluidPage(outputCodeButton(withLoader(plotOutput("betaMds", height = 700))),
                     uiOutput("betaMdsUI"))
betaCluster <- fluidPage(outputCodeButton(withLoader(plotOutput("betaCluster", height = 700))),
                        uiOutput("betaClusterUI"))
betaManova <- fluidPage(withLoader(verbatimTextOutput("betaManova")),
                           uiOutput("betaManovaUI"),
                           verbatimTextOutput("betaManova2"))
betaNetwork <- fluidPage(outputCodeButton(withLoader(plotOutput("betaNetwork", height = 700))),
                          uiOutput("betaNetworkUI"))
betaHeatmap <- fluidPage(outputCodeButton(withLoader(plotOutput("betaHeatmap", height = 700))),
                         uiOutput("betaHeatmapUI"))
betaTable <- fluidPage(uiOutput("betaTable"))
