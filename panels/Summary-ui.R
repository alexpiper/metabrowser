Summary <- fluidPage(
  withLoader(verbatimTextOutput("phyloseqPrint")),
  tags$footer(
    "Questions, problems or comments regarding this application should be sent to ",
    a(href = "mailto:cedric.midoux@inrae.fr?subject=[Easy16S]", "cedric.midoux@inrae.fr"),
    align = "center",
    style = "position:absolute;bottom: 0;width: 100%;color: grey;padding: 10px;# background-color: white;z-index: 1000;"
  )
)
