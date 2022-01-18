output$otuTable <- renderUI({
  validate(need(physeq(), "Requires an abundance dataset"))
  box(
    title = "OTU table",
    width = NULL,
    status = "primary",
    beautifulTable(tibble::as_tibble(as.data.frame(otu_table(physeq())), rownames = "OTU"))
  )
})

output$taxTable <- renderUI({
  validate(need(physeq(), "Requires an abundance dataset"))
  box(
    title = "Taxonomy table",
    width = NULL,
    status = "primary",
    beautifulTable(tibble::as_tibble(as.data.frame(tax_table(physeq())), rownames = "OTU"))
  )
})

output$glomTable <- renderUI({
  validate(need(physeq(), "Requires an abundance dataset"))
  box(
    title = "Agglomerate OTU table",
    width = NULL,
    status = "primary",
    radioButtons("glomRank",
                 label = "Taxonomic rank : ",
                 choices = rank_names(physeq()),
                 inline = TRUE),
    withLoader(dataTableOutput("tableGlom"))
  )
})

output$tableGlom <- renderDataTable({
  Glom <- tax_glom(physeq(), input$glomRank, NArm = FALSE)
  taxTableGlom <- Glom %>%
    tax_table() %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    dplyr::select(1:input$glomRank) %>%
    tibble::rownames_to_column()
  otuTableGlom <- Glom %>%
    otu_table() %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    tibble::rownames_to_column()
  joinGlom <- dplyr::left_join(taxTableGlom, otuTableGlom, by = "rowname") %>%
    dplyr::select(-rowname)
  beautifulTable(joinGlom)
})

output$sampleTable <- renderUI({
  validate(need(physeq(), "Requires an abundance dataset"))
  box(
    title = "Sample data table",
    width = NULL,
    status = "primary",
    beautifulTable(tibble::as_tibble(sample_data(physeq()), rownames = "SAMPLE")),
    box(
      title = "Class of sample data",
      width = NULL,
      status = "primary",
      collapsible = TRUE,
      collapsed = TRUE,
      renderTable({(sapply(sample_data(physeq()), class))}, rownames = TRUE, colnames = FALSE)
      )
    )
})
