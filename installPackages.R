
# List of cran packages
cran_packages <-c(
"shinydashboard",
"shinymeta",
"shinycustomloader",
"shinyWidgets",
"shinyAce",
"phyloseq",
"ape",
"tidyverse",
"DT",
"reshape2",
"data.table",
"magrittr",
"factoextra",
"gridExtra",
"readxl",
"vegan",
"plotly",
"devtools"
)

new_cran_packages <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
if(length(new_cran_packages)) install.packages(new_cran_packages)

# Bioconductor packages
if(!"phyloseq" %in% installed.packages()[,"Package"]) {BiocManager::install("phyloseq")}
if(!"DESeq2" %in% installed.packages()[,"Package"]) {BiocManager::install("DESeq2")}
