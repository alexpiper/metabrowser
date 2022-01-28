Summary <- fluidPage(
  withLoader(plotOutput("filtplot", height = 700)),
  uiOutput("rawphyloseqPrint"),
  uiOutput("phyloseqPrint"),
  tags$footer(
    "Questions regarding this application should be sent to ",
    a(href = "mailto:alexander.piper@agriculture.vic.gov.au?subject=[metabrowser]", "alexander.piper@agriculture.vic.gov.au"),
    align = "center",
    style = "position:absolute;bottom: 0;width: 100%;color: grey;padding: 10px;# background-color: white;z-index: 1000;"
  )
)