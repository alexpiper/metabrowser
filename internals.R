## High level Helper functions

.import_biom <- function(input) {
  ## Select appropriate import_function
  import_function <- switch(input$biomFormat,
                            "std"   = import_biom,
                            "frogs" = import_frogs)
  ## Import
  return(import_function(input$fileBiom$datapath, ## biom file
                         input$fileTree$datapath  ## tree file
  ))
}

.import_sample_data <- function(input, physeq) {
  ## Unhappy path
  if (is.null(input$fileMeta)) {
    return(data.frame(SampleID = sample_names(physeq) , row.names = sample_names(physeq)))
  }

  ## Happy path: excel version
  if (input$CSVsep == "excel") {
    sdf <- as.data.frame(readxl::read_excel(input$fileMeta$datapath))
    row.names(sdf) <- sdf[, 1]
    sdf <- sdf[, -1]
    sdf$SampleID <- rownames(sdf)
    return(sdf)
  }

  ## Happy path: csv version
  sdf <- read.csv(
    input$fileMeta$datapath,
    header = TRUE,
    sep = input$CSVsep,
    row.names = 1,
    na.strings = NA
  )
  sdf$SampleID <- rownames(sdf)
  return(sdf)
}

.format_tax_table <- function(tdf) {
  ## explicit rank names
  colnames(tdf) <- c("Kingdom", "Phylum", "Class", "Order",
                     "Family", "Genus", "Species", "Strain")[1:ncol(tdf)]
  ## Replace unknown by NA
  tdf[grep("unknown ", tdf)] <- NA
  #tdf[grep("Unclassified", tdf)] <- NA
  return(tdf)
}

.import_from_rdata <- function(input) {
  ## Happy path
  ne <- new.env() ## new env to store RData content and avoid border effects
  if (!is.null(input$fileRData))
    load(input$fileRData$datapath, envir = ne)
  if (class(ne$data) == "phyloseq")
    return(ne$data)

  ## Unhappy paths: everything else
  return()
}

###

checkNull <- function(x) {
  print(as.character(substitute(x)))
  #if (!exists(as.character(substitute(x)))) { #NOTE: This was breaking
  #if (!exists(as.character(x))) { 
  if (!any(sapply(as.character(substitute(x)), exists))) { #NOTE: This was breaking
    return(NULL)
  } else if (is.null(x)) {
    return(NULL)
  } else if (length(x) > 1) {
    return(x)
  }
  else if (x %in% c(0, "", NA, "NULL")) {
    return(NULL)
  } else {
    return(x)
  }
}

beautifulTable <- function(data)  {
  DT::datatable(
    data = data,
    rownames = FALSE,
    filter = "top",
    extensions = c("Buttons", "ColReorder", "FixedColumns"),
    options = list(
      dom = "Btlip",
      pageLength = 10,
      lengthMenu = list(c(10, 25, 50, 100,-1), list('10', '25', '50', '100', 'All')),
      buttons = list(
        'colvis',
        list(
          extend = 'collection',
          buttons = c('copy', 'csv', 'excel', 'pdf'),
          text = 'Download'
        )
      ),
      colReorder = TRUE,
      scrollX = TRUE,
      fixedColumns = list(leftColumns = 1, rightColumns = 0)
    ),
    width = "auto",
    height = "auto"
  )
}
