rarefactionCurve <- fluidPage(outputCodeButton(withLoader(plotOutput("rarefactionCurve", height = 700))),
                              uiOutput("rarefactionCurveUI"))
