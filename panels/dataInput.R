### Upload Data ###
dataInput <- function(failed = FALSE) {
  modalDialog(
    title = "Select your data",
    
    if (failed)
    {
      div(
        "Invalide dataset. \n Please try again \n", 
        style = "font-weight: bold; color: red; text-align: center; white-space: pre-line"
        )
    }, 
    
    "Firstly, you should select a demo dataset or upload an abundance BIOM file. For example, with Galaxy, a BIOM file can be obtained at the end of FROGS workflow with the `FROGS BIOM to std BIOM` tool. Make sure that the phyloseq object in the RData file is called `data`.",
    
    radioButtons(
      inputId = "dataset",
      label = "Select dataset : ",
      inline = TRUE,
      choices = list(
        "Demo" = "demo",
        "Input data" = "input",
        "RData" = "rdata",
        "RDS" = "rds"
      ),
      selected = "demo"
    ),
    
    wellPanel(uiOutput("dataUI")),

    footer = tagList(modalButton("Cancel"),
                     actionButton(inputId = "okData", label = "OK"))
  )
}

output$dataUI <- renderUI({
  if (is.null(input$dataset))
    return()
  
  switch(
    input$dataset,
    "demo" = selectInput(
      inputId = "demo",
      label = "Select a demo dataset",
      choices = c("Chaillou et al., 2015" = "food",
      # "Mach et al., 2015" = "kinetic",
      # "Morton et al., 2017" = "soil",
      # "Ravel et al., 2011" = "ravel",
      # "biorare" = "biorare",
      "GlobalPatterns" = "GlobalPatterns")
    ),
    "input" = tags$div(
      tags$div(
        title = "Abundance BIOM file come from FROGS with 'FROGS BIOM to std BIOM', Qiime or another metagenomic tool.",
        fileInput(
          inputId = "fileBiom",
          label = "Abundance BIOM file : ",
          placeholder = "data.biom"
        )
      ),
      radioButtons(
        inputId = "biomFormat",
        label = NULL,
        inline = TRUE,
        choices = c("STD BIOM" = "std",
                    "FROGS BIOM" = "frogs"),
        selected = "std"
      ),
      # tags$div(
      # style = "text-align:center",
      # title = "Resample dataset such that all samples have the same library size. \nIt's using an random sampling without replacement.",
      checkboxInput(
        inputId = "rareData",
        label = "Rarefy dataset",
        value = FALSE
      ),
      # textOutput("rarefactionMin")
      # ),
      tags$div(
        title = "Metadata table with variables (in columns) and samples (in rows). \nMake sure you follow the exact spelling of the sample names (1st column). \nThe import of an excel table is possible but not recommended.",
        fileInput(
          inputId = "fileMeta",
          label = "Metadata table : ",
          placeholder = "data.csv"
        )
      ),
      radioButtons(
        inputId = "CSVsep",
        label = "CSV separator : ",
        inline = TRUE,
        choices = c("<tab>" = "\t",
                    "," = ",",
                    ";" = ";",
                    excel = "excel"
                    )
      ),
      fileInput(
        inputId = "fileTree",
        label = "Phylogenetic tree : ",
        placeholder = "data.nwk"),
      fileInput(
        inputId = "fileSeq",
        label = "Representative FASTA sequences of OTU : ",
        placeholder = "data.fasta"
      )
    ),
    "rdata" = fileInput(
      inputId = "fileRData",
      label = "RData where 'data' is a phyloseq object : ",
      placeholder = "data.RData"
    ),
    "rds" = fileInput(
      inputId = "fileRDS",
      label = "RDS with a phyloseq object : ",
      placeholder = "phyloseq.RDS"
    )
    
  )
})

observeEvent(input$okData, {
  raw_physeq(NULL)
  try(
    raw_physeq(
      switch(
        input$dataset,
        "demo" =
          {
            message <- as.character(input$demo)
            readRDS(paste0("demo/", input$demo,".RDS"))
          },
        "input" =
          {
            message <- as.character(input$fileBiom$name)
            d <- .import_biom(input)
            ## Format tax table
            tax_table(d) <- .format_tax_table(tax_table(d))
      
            ## import metadata and store it in phyloseq object
            sample_data(d) <- .import_sample_data(input, d)
      
            ## Rarefy data
            if (input$rareData)
            {
              d <- rarefy_even_depth(
                d,
                replace = FALSE,
                rngseed = 314,
                verbose = FALSE
              )
            }
            d
          },
        "rdata" =
          {
            message <- as.character(input$fileRData$name)
            ne <- new.env()
            if (!is.null(input$fileRData))
              {load(input$fileRData$datapath, envir = ne)}
            ne$data
          },
        "rds" = 
          {
            message <- as.character(input$fileRDS$name)
            if (!is.null(input$fileRDS))
              {readRDS(input$fileRDS$datapath)}
          }
        )
      ),
    silent = TRUE, 
    outFile = showModal(dataInput(failed = TRUE))
    )
  
  if (inherits(raw_physeq(), "phyloseq")) {
    select_physeq(raw_physeq())
    physeq(raw_physeq())
    transform_physeq(NULL)
    shinyWidgets::updateSwitchInput(session = session, inputId = "useTransf", disabled = TRUE, value = FALSE)
    message(paste("[Easy16S] Correct upload with", input$dataset, "mode :", message))
    removeModal()
  } else {
    showModal(dataInput(failed = TRUE))
  }
})

### Filter ###
filterSample <- function() {
  modalDialog(
    title = "Select some sample",
    
    radioButtons(
      inputId = "filterCriteria",
      label = "Filter criteria : ",
      inline = TRUE,
      choices = c(
        "Sample" = "sample",
        sample_variables(select_physeq(), errorIfNULL = FALSE)
      ),
      selected = "sample"
    ),
    
    wellPanel(uiOutput("filterUI"), 
              actionLink(inputId = "selectAll", label = "(Un)select All")),
    
    footer = tagList(modalButton("Cancel"),
                     actionButton(inputId = "okData", label = "Refresh filter"),
                     actionButton(inputId = "selectData", label = "Select")
    )
  )
}

output$filterUI <- renderUI({
  if (is.null(input$filterCriteria))
    return()
  
  if (input$filterCriteria == "sample") {
    label <- "Sample to keep :"
    choices <- sample_names(select_physeq())
  } else {
    label <- "Variable to keep :"
    choices <- levels(get_variable(select_physeq(), input$filterCriteria))
  }
  
  observe({
    if(input$selectAll == 0) return(NULL)
    else if (input$selectAll%%2 == 0)
    {
      updateCheckboxGroupInput(session = session,
                               inputId = "filterCheck",
                               label = label,
                               choices = choices,
                               selected = choices,
                               inline = TRUE)
    } else {
      updateCheckboxGroupInput(session = session,
                               inputId = "filterCheck",
                               label = label,
                               choices = choices,
                               inline = TRUE)
    }
  })
  
  checkboxGroupInput(inputId = "filterCheck",
                     label = label,
                     choices = choices,
                     selected = choices,
                     inline = TRUE)
})

observeEvent(input$selectData, {
  if (is.null(input$filterCheck)) 
  {
    showModal(dataInput(failed = TRUE))
  } else {
    try(
      if (input$filterCriteria == "sample") {
        select_physeq(prune_samples(samples = input$filterCheck, select_physeq()))
      } else {
        criteria <<- input$filterCriteria
        check <<- input$filterCheck
        select_physeq(subset_samples(select_physeq(), eval(parse(text = criteria)) %in% check))
        },
      silent = TRUE,
      outFile = showModal(dataInput(failed = TRUE)))
    
    if (inherits(select_physeq(), "phyloseq")) {
      message <- paste(input$filterCheck, collapse = ", ")
      message(paste("[Easy16S] Select some samples :", message))
      physeq(select_physeq())
      shinyWidgets::updateSwitchInput(session = session, inputId = "useTransf", value = FALSE, disabled = TRUE)
      removeModal()
    } else {
      showModal(dataInput(failed = TRUE))
    }    
  }
})

### Transformation ###
transformSample <- function() {
  modalDialog(
    title = "Transform abundance data",
    
    radioButtons(
      inputId = "transformFun",
      label = "Transform function : ",
      selected = NULL,
      choices = c("Rarefaction" = "rare",
                  "Proportional Transformation" = "prop", 
                  "Square Root Transformation" = "sqrt", 
                  "Square Root Proportional Transformation" = "sqrtprop", 
                  "Centered Log-Ratio (CLR) Transformation" = "clr")
    ),
    
    wellPanel(verbatimTextOutput("transformFun")),
    
    footer = tagList(modalButton("Cancel"),
                     actionButton(inputId = "transformData", label = "Transform")
    )
  )
}

output$transformFun <- renderText({
  validate(need(input$transformFun, ""))
  switch (input$transformFun,
          "rare" = "data_rare <- rarefy_even_depth(data, rngseed = 314, replace = TRUE)",
          "prop" = paste("count_to_prop <- function(x) {return( x / sum(x) )}", 
                         "data_prop <- transform_sample_counts(data, count_to_prop)",
                         sep = "\n"),
          "sqrt" = "data_sqrt <- transform_sample_counts(data, sqrt)",
          "sqrtprop" = paste("count_to_sqrtprop <- function(x) {return(sqrt(x / sum(x)))}", 
                         "data_sqrtprop <- transform_sample_counts(data, count_to_sqrtprop)",
                         sep = "\n"),
          "clr" = paste("gm_mean <- function(x, na.rm=TRUE) {",
                        "  return(exp(mean(log(x), na.rm=na.rm)))",
                        "}",
                        "clr <- function(x, base=exp(1)) {",
                        "  x <- x+1",
                        "  x <- log(x/gm_mean(x), base)",
                        "  return(x)",
                        "}", 
                        "data_clr <- transform_sample_counts(data, clr)",
                        sep = "\n")
          )
})

observeEvent(input$transformData, {
  if (is.null(input$transformData)) 
  {
    removeModal()
  } else {
    shinyWidgets::updateSwitchInput(session = session, inputId = "useTransf", value = FALSE)
    try(
      switch (input$transformFun,
              "rare" = {
                transform_physeq(rarefy_even_depth(select_physeq(), rngseed = 314, replace = TRUE))
                },
              "prop" = {
                count_to_prop <- function(x) {return( x / sum(x) )}
                transform_physeq(transform_sample_counts(select_physeq(), count_to_prop))
              },
              "sqrt" = {
                transform_physeq(transform_sample_counts(select_physeq(), sqrt))
              },
              "sqrtprop" = {
                count_to_sqrtprop <- function(x) {return(sqrt(x / sum(x)))}
                transform_physeq(transform_sample_counts(select_physeq(), count_to_sqrtprop))
              },
              "clr" = {
                gm_mean <- function(x, na.rm=TRUE) {
                  return(exp(mean(log(x), na.rm=na.rm)))
                }
                clr <- function(x, base=exp(1)) {
                  x <- x+1
                  x <- log(x/gm_mean(x), base)
                  return(x)
                }
                transform_physeq(transform_sample_counts(select_physeq(), clr))
              }
      ),
      silent = TRUE,
      outFile = showModal(dataInput(failed = TRUE)))
    
    if (class(transform_physeq()) == "phyloseq") {
      message(paste("[Easy16S] Transfom data :", input$transformFun))
      shinyWidgets::updateSwitchInput(session = session, inputId = "useTransf", disabled = FALSE)
      shinyWidgets::updateSwitchInput(session = session, inputId = "useTransf", value = TRUE)
      removeModal()
    } else {
      showModal(dataInput(failed = TRUE))
    }
  }
})

observeEvent(input$useTransf,
             if (input$useTransf) {
               physeq(transform_physeq())
             } else {
               physeq(select_physeq())
             }
             )

### Download Data ###
dataDownload <- function() {
  modalDialog(
    title = "Download data",
    size = "s",
    textInput("dataName", "File name : ", value = paste("Easy16S-data", Sys.Date(), sep = "-")),
    radioButtons("dataFormat", "File format : ", choices = c("RData", "RDS", "biom"), selected = "RData", inline = TRUE),
    footer = tagList(modalButton("Cancel"),
                     downloadButton("okDownload", "Download")
    )
  )
}

output$okDownload <- downloadHandler(
  filename = function() {
    paste(input$dataName, input$dataFormat, sep = ".")
  },
  content = function(file) {
    if (input$dataFormat == "RData") {
      data <- physeq()
      save(data, file = file)
    } else if (input$dataFormat == "RDS") {
      saveRDS(physeq(), file = file)
    } else if (input$dataFormat == "biom") {
      write_phyloseq(physeq = physeq(), biom_file = file, biom_format = "frogs") #"standard"
    }
  }
)

### Download Plot ###
plotDownload <- function() {
  modalDialog(
    title = "Download last modified plot",
    size = "s",
    textInput("plotName", "File name : ", value = "plot"),
    numericInput("plotWidth", "Width (cm) : ", value = 7),
    numericInput("plotHeight", "Height (cm) : ", value = 7),
    numericInput("plotDPI", "DPI : ", value = 300),
    radioButtons("plotFormat", "File format : ", choices = c("png", "pdf", "jpeg", "svg", "wmf"), selected = "png", inline = TRUE),
    footer = tagList(modalButton("Cancel"),
                     downloadButton("okPlot", "Download")
                     )
  )
}

output$okPlot <- downloadHandler(
  filename = function() {paste(input$plotName, input$plotFormat, sep = ".")},
  content = function(file) {ggsave(file, width = input$plotWidth, height = input$plotHeight, dpi = input$plotDPI, units = "cm")}
)

# output$rarefactionMin <- renderText({
#   validate(need(input$fileBiom, ""),
#            need(input$dataset == "input", ""))
#   paste("(min sample =", format(min(sample_sums(data16S(
#     
#   ))), big.mark = " "), "reads)")
# })
