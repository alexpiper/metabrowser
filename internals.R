## High level Helper functions

# Create phyloseq ---------------------------------------------------------
make_ps <- function(seqtabPath, taxtabPath, samdfPath, treePath=NULL, refseqpath=NULL) {
  # Import seqtab
  if (grepl(tolower(seqtabPath), pattern = ".txt|.tsv")) {
    seqtab <- read.table(seqtabPath, header = 1,
               sep = "\t", stringsAsFactors = F,
               quote = "", comment.char = "")
  } else if (grepl(tolower(seqtabPath), pattern = ".csv")) {
    seqtab <- read.csv(seqtabPath, header = 1,
             stringsAsFactors = F,
             quote = "", comment.char = "")%>%
      tibble::column_to_rownames("sample_id") %>%
      as.matrix()
  } else if (grepl(tolower(seqtabPath), pattern = ".rds")) {
    seqtab <- readRDS(seqtab)
  } else {
    seqtab <- NULL
  }
  
  # Import taxtab
  if (grepl(tolower(taxtabPath), pattern = ".txt|.tsv")) {
    taxtab <- read.table(taxtabPath, header = 1,
               sep = "\t", stringsAsFactors = F,
               quote = "", comment.char = "")
  } else if (grepl(tolower(taxtabPath), pattern = ".csv")) {
    taxtab <- read.csv(taxtabPath, header = 1,
             stringsAsFactors = F,
             quote = "", comment.char = "")%>%
      tibble::column_to_rownames("OTU") %>%
      as.matrix()
    
  }else if (grepl(tolower(taxtabPath), pattern = ".rds")) {
    taxtab <- readRDS(taxtabPath)
  } else {
    taxtab <- NULL
  }
  
  # Import samdf
  if (grepl(tolower(samdfPath), pattern = ".txt|.tsv")) {
    samdf <- read.table(samdfPath,
               header = 1, sep = "\t", stringsAsFactors = F,
               quote = "", comment.char = "")
  } else if (grepl(tolower(samdfPath), pattern = ".csv")) {
    samdf <- read.csv(samdfPath, header = 1,
             stringsAsFactors = F,
             quote = "", comment.char = "")
  }else {
    samdf <- NULL
  }
  
  # Import tree
  if(!is.null(treePath)){
    if (grepl(tolower(treePath), pattern = ".nwk")) {
      tree <- phyloseq::read_tree(treePath)
    }else if (grepl(tolower(treePath), pattern = ".rds")) {
      tree <- readRDS(tree)
      if (!class(tree) == "phylo"){
        tree <- NULL
        warning("Tree .rds file is not a phylo object")
      }
    } else {
      warning("Tree file is not .nwk or .rds file")
      tree <- NULL
    }
  }else {
    tree <- NULL
  }
  
  # refseq
  if(!is.null(refseqpath)){
    if (grepl(tolower(refseqpath), pattern = ".fasta|.fasta.gz|.fa|.fa.gz")) {
      seqs <- Biostrings::readDNAStringSet(refseqpath)
    } else {
      warning("refseq file is not a fasta file")
      seqs <- NULL
    }
  }else {
    seqs <- NULL
  }
  
  # Create phyloseq object
  if(!is.null(tree) & !is.null(seqs)){
    ps <- phyloseq(tax_table(taxtab), 
                   sample_data(samdf),
                   otu_table(seqtab, taxa_are_rows = FALSE),
                   phy_tree(tree),
                   refseq(seqs)
                   ) 
  } else if (!is.null(tree) & is.null(seqs)){
    ps <- phyloseq(tax_table(taxtab), 
                   sample_data(samdf),
                   otu_table(seqtab, taxa_are_rows = FALSE),
                   phy_tree(tree)
    ) 
  } else if (is.null(tree) & is.null(seqs)){
    ps <- phyloseq(tax_table(taxtab), 
                   sample_data(samdf),
                   otu_table(seqtab, taxa_are_rows = FALSE)
    ) 
  } else {
    stop("An error occured creating phyloseq object")
  }
  return(ps)
}

# Import data file --------------------------------------------------------
.import_data <- function(input) {
  ## Import
  return(
    make_ps(
      seqtabPath = input$fileSeqtab$datapath,
      taxtabPath = input$fileTaxtab$datapath,
      samdfPath = input$FileSamdf$datapath,
      treePath = input$fileTree$datapath,
      refseqpath = input$fileSeq$datapath) 
  )
}

# Import biom file --------------------------------------------------------
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


# Import sample metadata file ---------------------------------------------
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


# Reformat tax table -----------------------------------------------------
.format_tax_table <- function(tdf) {
  ## explicit rank names
  colnames(tdf) <- c("Root", "Kingdom", "Phylum", "Class", "Order",
                     "Family", "Genus", "Species")[1:ncol(tdf)]
  ## Replace unknown by NA
  tdf[grep("unknown ", tdf)] <- NA
  tdf[grep("__", tdf)] <- NA
  #tdf[grep("Unclassified", tdf)] <- NA
  return(tdf)
}


# import from r data file -------------------------------------------------
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
