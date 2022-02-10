cat(file=stderr(), paste0(R.version))
options(shiny.maxRequestSize = 30 * 1024 ^ 2)

suppressMessages(suppressWarnings(library(shinydashboard)))
suppressMessages(suppressWarnings(library(shinymeta)))
suppressMessages(suppressWarnings(library(shinyWidgets)))
suppressMessages(suppressWarnings(library(shinycustomloader)))
suppressMessages(suppressWarnings(library(phyloseq)))
suppressMessages(suppressWarnings(library(ape)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(DT)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(magrittr)))
suppressMessages(suppressWarnings(library(factoextra)))
suppressMessages(suppressWarnings(library(vegan)))

shinyServer
(function(input, output, session)
{  
  source("internals.R")
  source("graphical_methods.R")
  source("filters.R")
  source("panels/dataInput.R", local = TRUE)
  source("panels/Summary-server.R", local = TRUE)
  source("panels/table-server.R", local = TRUE)
  source("panels/barplot-server.R", local = TRUE)
  source("panels/heatmap-server.R", local = TRUE)
  source("panels/rarefactionCurve-server.R", local = TRUE)
  source("panels/richnessA-server.R", local = TRUE)
  source("panels/richnessB-server.R", local = TRUE)
  source("panels/pca-server.R", local = TRUE)
  source("panels/deseq-server.R", local = TRUE)
  source("panels/tree-server.R", local = TRUE)
  source("panels/functional-server.R", local = TRUE)

  physeq <- reactiveVal()
  raw_physeq <- reactiveVal()
  select_physeq <- reactiveVal()
  subset_physeq <- reactiveVal()
  filter_physeq <- reactiveVal()
  transform_physeq <- reactiveVal()
  
  showModal(dataInput())

  observeEvent(input$dataButton, {
    showModal(dataInput())
  })
  
  observeEvent(input$selectButton, {
    showModal(subsetSample())
  })
  observeEvent(input$filterButton, {
    showModal(filterSample())
  })
  observeEvent(input$taxaButton, {
    showModal(subsetTaxa())
  })
  observeEvent(input$transformButton, {
    showModal(transformSample())
  })
  
  observeEvent(input$downloadButton, {
    showModal(dataDownload())
  })
  
  observeEvent(input$plotButton, {
    showModal(plotDownload())
  })
  
})
