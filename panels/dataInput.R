### Upload Data ###
dataInput <- function(failed = FALSE) {
  modalDialog(
    title = "Select your data",
    
    if (failed)
    {
      div(
        "Invalid dataset. \n Please try again \n", 
        style = "font-weight: bold; color: red; text-align: center; white-space: pre-line"
        )
    }, 
    
    "Firstly, you should select a demo dataset or either upload the necessary files to create a phyloseq object or an already made phyloseq object",
    
    radioButtons(
      inputId = "dataset",
      label = "Select dataset : ",
      inline = TRUE,
      choices = list(
        "Demo" = "demo",
        "Input data" = "input",
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
        title = "Input seqtab file",
        fileInput(
          inputId = "fileSeqtab",
          label = "seqtab file : ",
          placeholder = "seqtab.csv"
        )
      ),
      tags$div(
        title = "Input taxtab file",
        fileInput(
          inputId = "fileTaxtab",
          label = "taxtab file : ",
          placeholder = "taxtab.csv"
        )
      ),
      #radioButtons(
      #  inputId = "biomFormat",
      #  label = NULL,
      #  inline = TRUE,
      #  choices = c("STD BIOM" = "std",
      #              "FROGS BIOM" = "frogs"),
      #  selected = "std"
      #),
      tags$div(
        title = "Metadata table with variables (in columns) and samples (in rows). \nMake sure you follow the exact spelling of the sample names (1st column). \nThe import of an excel table is possible but not recommended.",
        fileInput(
          inputId = "FileSamdf",
          label = "Metadata table : ",
          placeholder = "samdf.csv"
        )
      ),
      fileInput(
        inputId = "fileTree",
        label = "Phylogenetic tree : ",
        placeholder = "data.nwk"),
      fileInput(
        inputId = "fileSeq",
        label = "Representative FASTA sequences of OTUs : ",
        placeholder = "data.fasta"
      ),
      checkboxInput(
        inputId = "renameASV",
        label = "Rename ASV's",
        value = FALSE
      ),
      checkboxInput(
        inputId = "rareData",
        label = "Rarefy dataset",
        value = FALSE
      )
      # tags$div(
      # style = "text-align:center",
      # title = "Resample dataset such that all samples have the same library size. \nIt's using an random sampling without replacement.",
      # textOutput("rarefactionMin")
      # ),
    ),
    "rds" = fileInput(
      inputId = "fileRDS",
      label = "RDS with a phyloseq object : ",
      placeholder = "phyloseq.RDS"
    )
    
  )
})

# Inport data when OK button is pressed
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
            # Process seqtabs
            message <- as.character(input$fileBiom$name)
            d <- .import_data(input)
            
            print(d)
            ## Format tax table
            tax_table(d) <- .format_tax_table(tax_table(d)) #Need to detect double underscores here as well
      
            ## import metadata and store it in phyloseq object
            #sample_data(d) <- .import_sample_data(input, d)
      
            ## Rename ASV's to numbered lowest rank
            # Number asv's then addd
            if (input$renameASV) {
              # Check if taxa names are sequences - if so add them to refseq
              #if(all(grepl("[^ACTG]", colnames(phyloseq::get_taxa(d))))){
              #  new_seqs <- Biostrings::DNAStringSet(colnames(phyloseq::get_taxa(d)))
              #}
              # Get lowest classified rank
              lowest_name <- seqateurs::na_to_unclassified(as(tax_table(d), "matrix")) %>%
                as.data.frame()%>%
                dplyr::pull(Species) %>%
                stringr::str_replace_all(" ", "_")
              taxa_names(d) <-  paste0("SV", 1:phyloseq::ntaxa(d), "_", lowest_name)
            }
            
            ## Rarefy data
            if (input$rareData) {
              d <- rarefy_even_depth(
                d,
                replace = FALSE,
                rngseed = 314,
                verbose = FALSE
              )
            }
            d
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
    subset_physeq(raw_physeq())
    physeq(raw_physeq())
    transform_physeq(NULL)
    shinyWidgets::updateSwitchInput(session = session, inputId = "useTransf", disabled = TRUE, value = FALSE)
    shinyWidgets::updateSwitchInput(session = session, inputId = "useFiltf", disabled = TRUE, value = FALSE)
    message(paste("[metabrowser] Correct upload with", input$dataset, "mode :", message))
    removeModal()
  } else {
    showModal(dataInput(failed = TRUE))
  }
})


# Subset samples ----------------------------------------------------------
subsetSample <- function() {
  modalDialog(
    title = "Subset samples",
    
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
    
    wellPanel(uiOutput("subsetUI"), 
              actionLink(inputId = "selectAll", label = "(Un)select All")),
    
    footer = tagList(modalButton("Cancel"),
                     actionButton(inputId = "okData", label = "Refresh filter"),
                     actionButton(inputId = "selectData", label = "Select")
    )
  )
}

output$subsetUI <- renderUI({
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
        select_physeq(prune_samples(samples = input$filterCheck, select_physeq())%>% 
                        filter_taxa(function(x) mean(x) > 0, TRUE) #Drop missing taxa from table
        )
      } else {
        criteria <<- input$filterCriteria
        check <<- input$filterCheck
        select_physeq(subset_samples(select_physeq(), eval(parse(text = criteria)) %in% check)%>% 
                        filter_taxa(function(x) mean(x) > 0, TRUE) #Drop missing taxa from table
        )
        },
      silent = TRUE,
      outFile = showModal(dataInput(failed = TRUE)))
    
    if (inherits(select_physeq(), "phyloseq")) {
      message <- paste(input$filterCheck, collapse = ", ")
      message(paste("[metabrowser] Select some samples :", message))
      physeq(select_physeq())
      subset_physeq(select_physeq())
      shinyWidgets::updateSwitchInput(session = session, inputId = "useTransf", value = FALSE, disabled = TRUE)
      removeModal()
    } else {
      showModal(dataInput(failed = TRUE))
    }    
  }
})

# Subset taxa ----------------------------------------------------------
#TODO Make this select desired taxa, either manually or using a list
subsetTaxa <- function() {
  modalDialog(
    title = "Subset taxa",

    radioButtons(
      inputId = "subsetCriteria",
      label = "Taxa subset rank : ",
      inline = TRUE,
      choices = c(
        phyloseq::rank_names(subset_physeq(), errorIfNULL = FALSE),
        "OTU"
      ),
      selected = phyloseq::rank_names(subset_physeq(), errorIfNULL = FALSE)[length(phyloseq::rank_names(subset_physeq(), errorIfNULL = FALSE))]
    ),
    
    wellPanel(uiOutput("taxaUI"), 
              fileInput(
                inputId = "fileTaxlist",
                label = "Taxon list : ",
                placeholder = "taxlist.txt"
              ),
              radioButtons(
                inputId = "select_type", 
                label = "Choose selection:",
                choices = c("Select all" = "sel_all",
                  "Unselect all" = "exc_all",
                  "Only include list" = "sel_list",
                  "All except list" = "exc_list"),
                selected = "sel_all")
              ),
    
    footer = tagList(modalButton("Cancel"),
                     actionButton(inputId = "okData", label = "Refresh filter"),
                     actionButton(inputId = "selectTaxa", label = "Select")
    )
  )
}

output$taxaUI <- renderUI({
  # TODO This needs to change depending on what rank you want to look at
  if (is.null(input$subsetCriteria))
    return()
  label <- "Taxa to keep :"
  if (input$subsetCriteria == "OTU") {
    choices <- taxa_names(subset_physeq())
  } else {
    choices <- unique(as.data.frame(as(tax_table(subset_physeq()), "matrix"))[,input$subsetCriteria])
    choices <- choices[!is.na(choices)]
  }
  
  # Import
  taxlist <- NULL
  observeEvent(input$fileTaxlist, {
    taxlist <<- readLines(input$fileTaxlist$datapath)
  })

  # Handle the select all checkbox
  observe({
    curr_select <- choices
  if (input$select_type == "sel_all") {
    curr_select <- choices
    print("selecting all")
  } else if (input$select_type == "exc_all") {
      curr_select <- NULL
      print("unselecting all")
  } else if (input$select_type == "sel_list") {
    validate(
      need(!is.null(taxlist), "Must load taxlist first")
    )
    curr_select <- taxlist
    print("including taxlist only")
  }else if (input$select_type == "exc_list") {
    validate(
      need(!is.null(taxlist), "Must load taxlist first")
    )
    curr_select <- choices[!choices %in% taxlist]
    print("excluding taxlist")
  }
    updateCheckboxGroupInput(session = session,
                             inputId = "subsetCheck",
                             label = label,
                             choices = choices,
                             selected = curr_select,
                             inline = TRUE)
  })
  
  # Handle manual selecting
  checkboxGroupInput(inputId = "subsetCheck",
                     label = label,
                     choices = choices,
                     selected = choices,
                     inline = TRUE)
})

observeEvent(input$selectTaxa, {
  if (is.null(input$subsetCheck)) {
    showModal(dataInput(failed = TRUE))
  } else {
    taxtab <- as.data.frame(as(tax_table(subset_physeq()), "matrix"))
    filt_pos <- which(rank_names(subset_physeq()) == input$subsetCriteria)
    to_subset <- rownames(taxtab)[taxtab[,filt_pos] %in% input$subsetCheck]
    
    try(subset_physeq(prune_taxa(taxa = to_subset, subset_physeq())),
      silent = TRUE,
      outFile = showModal(dataInput(failed = TRUE)))
    
    if (inherits(subset_physeq(), "phyloseq")) {
      message <- paste(input$subsetCheck, collapse = ", ")
      message(paste("[metabrowser] Select some taxa :", message))
      physeq(subset_physeq())
      select_physeq(subset_physeq())
      shinyWidgets::updateSwitchInput(session = session, inputId = "useTransf", value = FALSE, disabled = TRUE)
      removeModal()
    } else {
      showModal(dataInput(failed = TRUE))
    }    
  }
})

# filter taxa by abundance  ---------------------------------------------------------------
filterSample <- function() {
  modalDialog(
    title = "Filter abundance data",
    h4('OTU filtering'),
    fluidRow(column(width=12,
                    div(class="col-md-6",
                        numericInput("filter_taxa_sums_threshold", "Taxa Min reads per sample.",
                                        value=0, min=0, step=1))
    )),
    fluidRow(column(width=12,
                    div(class="col-md-6",
                        numericInputIcon("filter_taxa_ra_threshold", "Taxa Min relative abundance per sample.",
                                     value=0, min=0, max=100, step=NA, icon = list(NULL, icon("percent"))))
    )),
    fluidRow(column(width=12,
                    div(class="col-md-6",
                        numericInputIcon("filter_taxa_mean_threshold", "Taxa mean relative abundance over whole dataset.",
                                     value=0, min=0, max=100, step=NA, icon = list(NULL, icon("percent"))))
    )),
    fluidRow(column(width=12,
                    div(class="col-md-6",
                        numericInputIcon("filter_taxa_total_threshold", "Taxa relative abundance above fraction over whole dataset",
                                     value=0, min=0, max=100, step=NA, icon = list(NULL, icon("percent"))))
    )),
    h4('Include only taxa with more than A reads (on average) in at least k% samples'),
    fluidRow(column(width=12,
                    div(class="col-md-6",
                        numericInputIcon("filter_kOverA_count_threshold", "A -- The relative abundance value minmum threshold",
                                     value=0, min=0, max=100, step=NA, icon = list(NULL, icon("percent"))), # TODO Change the step to 1 when reads is selected
                    div(class="col-md-6",
                        numericInputIcon("filter_kOverA_sample_threshold", "k - The proportion of samples in which a taxa exceeded A",
                                     min=0, max=100, value=0, step=1, icon = list(NULL, icon("percent"))))
    ))),
    h4('Keep top N taxa'),
    fluidRow(column(width=12,
                    div(class="col-md-6",
                        numericInput("filter_top_taxa", "Number of most abundant taxa to keep",
                                     value=Inf, min=0, step=1))
    )),
    h4('Sample filtering'),
    fluidRow(column(width=12,
                    div(class="col-md-6",
                        numericInput("filter_sample_sums_threshold", "Sample Min reads.",
                                     value=0, min=0, step=1))
    )),
    fluidRow(column(width=12,
                    div(class="col-md-6",
                        numericInput("filter_sample_taxa_threshold", "Sample Min taxa",
                                     value=0, min=0, step=1))
    )),
    footer = tagList(modalButton("Cancel"),
                     actionButton(inputId = "filterData", label = "Filter")
    )
  )
}

#Update the use filtered data option on left
observeEvent(input$filterData, {
  if (is.null(input$filterData)) 
  {
    removeModal()
  } else {
    shinyWidgets::updateSwitchInput(session = session, inputId = "useFiltf", value = FALSE)
    
    # Define filters
    filter_all <- function(input, physeq){
      
      # Get initial phyloseq object
      ps0 <- physeq
      
      # Per sample counts
      if( input$filter_taxa_sums_threshold > 0 ){
        ps0 <- phyloseq_filter_sample_wise_abund_trim(
          ps0, 
          minabund = input$filter_taxa_sums_threshold,
          relabund = FALSE,
          rm_zero_OTUs = TRUE)
        
      }
      
      # Per sample RA
      if( input$filter_taxa_ra_threshold > 0 ){
        ps0 <- phyloseq_filter_sample_wise_abund_trim(
          ps0, 
          minabund = 1/input$filter_taxa_ra_threshold,
          relabund = TRUE,
          rm_zero_OTUs = TRUE)
        
      }
      
      # Total dataset mean relative abundance.
      if( input$filter_taxa_mean_threshold > 0 ){
        ps0 <- phyloseq_filter_taxa_rel_abund(
          ps0, 
          frac = 1/input$filter_taxa_mean_threshold)
      }
      
      #Remove taxa with abundance less then a certain fraction of total abundance.
      if( input$filter_taxa_total_threshold > 0 ){
        ps0 <- phyloseq_filter_taxa_tot_fraction(
          ps0, 
          frac = 1/input$filter_taxa_total_threshold)
      }
      
      # K over A filtering 
      if( input$filter_kOverA_sample_threshold > 1){
        ps0 <- phyloseq_filter_prevalence(
          ps0, 
          prev.trh = 1/input$filter_kOverA_sample_threshold,
          abund.trh = 1/input$filter_kOverA_count_threshold,
          threshold_condition = "AND",
          abund.type = "mean")
      }
      
      # Extract most abundant taxa
      if( input$filter_top_taxa > 0){
        ps0 <- phyloseq_filter_top_taxa(
          ps0, 
          n = input$filter_top_taxa
        )
      }
      
      # then filter sample sums
      if( input$filter_sample_sums_threshold > 0 ){
        # Sample sums filtering
        ps0 <- prune_samples({sample_sums(ps0) > input$filter_sample_sums_threshold}, ps0)
      }
      phyloseq_richness_filter <- function(physeq, mintaxa = 10){
        
        ## Estimate number of OTUs per sample
        sp <- phyloseq::estimate_richness(physeq, measures = "Observed")
        samples_to_keep <- rownames(sp)[ which(sp$Observed >= mintaxa) ]
        
        if(length(samples_to_keep) == 0){
          stop("All samples will be removed.\n")
        }
        
        if(length(samples_to_keep) == phyloseq::nsamples(physeq)){
          cat("All samples will be preserved\n")
          res <- physeq
        }
        
        if(length(samples_to_keep) < phyloseq::nsamples(physeq)){
          res <- phyloseq::prune_samples(samples = samples_to_keep, x = physeq)
        }
        
        return(res)
      }
      # Sample richness filtering
      if( input$filter_sample_taxa_threshold < Inf ){
        ps0 <- phyloseq_richness_filter(ps0, mintaxa = input$filter_sample_taxa_threshold)
      }
      return(ps0)
    }
      
    #apply_filters 
    try(
      filter_physeq(filter_all(input, select_physeq())),
      silent = TRUE,
      outFile = showModal(dataInput(failed = TRUE)))
    
    if (class(filter_physeq()) == "phyloseq") {
      message("[metabrowser] Filtered data")
      shinyWidgets::updateSwitchInput(session = session, inputId = "useFiltf", disabled = FALSE)
      shinyWidgets::updateSwitchInput(session = session, inputId = "useFiltf", value = TRUE)
      removeModal()
    } else {
      showModal(dataInput(failed = TRUE))
    }
  }
})

# Use filtered physeq if selected
observeEvent(input$useFiltf,
             if (input$useFiltf) {
               physeq(filter_physeq())
             } else {
               physeq(select_physeq())
             }
)


# Transform ---------------------------------------------------------------
# TODO: Need to make this pull from filter_physeq rather than select_physeq
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
      message(paste("[metabrowser] Transfom data :", input$transformFun))
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


# Download ----------------------------------------------------------------
dataDownload <- function() {
  modalDialog(
    title = "Download data",
    size = "s",
    textInput("dataName", "File name : ", value = paste("metabrowser-data", Sys.Date(), sep = "-")),
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
