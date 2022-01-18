list.of.packages <- c("shiny", "shinydashboard", "ggplot2", "reshape2", "ape", "gridExtra", "readxl", "shinycustomloader", "dplyr", "glue", "DT", "tibble", "plotly", "BiocManager", "devtools", "shinyAce")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

if(!"factoextra" %in% installed.packages()[,"Package"]) {devtools::install_github("kassambara/factoextra")}
if(!"phyloseq" %in% installed.packages()[,"Package"]) {BiocManager::install("phyloseq")}
if(!"DESeq2" %in% installed.packages()[,"Package"]) {BiocManager::install("DESeq2")}

devtools::install_github("mahendra-mariadassou/phyloseq-extended")
