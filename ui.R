suppressMessages(suppressWarnings(library(shinydashboard)))
suppressMessages(suppressWarnings(library(shinymeta)))
suppressMessages(suppressWarnings(library(shinycustomloader)))

source("panels/Summary-ui.R", local = TRUE)
source("panels/table-ui.R", local = TRUE)
source("panels/barplot-ui.R", local = TRUE)
source("panels/heatmap-ui.R", local = TRUE)
source("panels/rarefactionCurve-ui.R", local = TRUE)
source("panels/richnessA-ui.R", local = TRUE)
source("panels/richnessB-ui.R", local = TRUE)
source("panels/pca-ui.R", local = TRUE)
source("panels/deseq-ui.R", local = TRUE)
source("panels/tree-ui.R", local = TRUE)
source("panels/Help-ui.R", local = TRUE)

shinyUI(dashboardPage(
dashboardHeader(title = "Easy16S"),
  dashboardSidebar(
    actionButton("dataButton",
                 "Select your data",
                 icon = icon("upload"),
                 width = 200
                 ),
    actionButton("filterButton",
                 "Select some samples",
                 icon = icon("filter"),
                 width = 200
                 ),
    actionButton("transformButton",
                 "Transform abundance",
                 icon = icon("square-root-alt"),
                 width = 200
                 ),
    shinyWidgets::switchInput("useTransf",
                              value = FALSE,
                              label = HTML("Transformed&nbsp;data"),
                              disabled = TRUE,
                              size = "mini",
                              width = 200
                              ),
    actionButton("downloadButton",
                 "Download data",
                 icon = icon("download"),
                 width = 200
                 ),
    actionButton("plotButton",
                 "Download last plot",
                 icon = icon("file-image"),
                 width = 200
                 ),
    sidebarMenu(
      menuItem("Summary", tabName = "Summary", icon = icon("dna")),
      menuItem("Tables", icon = icon("table"),
               menuSubItem("OTU table", tabName = "otuTable"),
               menuSubItem("Taxonomy table", tabName = "taxtable"),
               menuSubItem("Agglomerate OTU table", tabName = "glomTable"),
               menuSubItem("Sample data table", tabName = "sampleTable")
      ),
      menuItem("Barplot", tabName = "barplot", icon = icon("chart-bar")),
      menuItem("Heatmap", tabName = "heatmap", icon = icon("chess-board")),
      menuItem("Rarefaction curves", tabName = "rarefactionCurve", icon = icon("chart-line")),
      menuItem(HTML("&alpha;-diversity"), icon = icon("th"),
               menuSubItem("Plots", tabName = "alphaPlot"),
               menuSubItem("Table", tabName = "alphaTable"),
               menuSubItem("ANOVA", tabName = "alphaAnova")
      ),
      menuItem(HTML("&beta;-diversity"), icon = icon("th"),
               selectInput("betaDistance",
                           label = "Distance : ",
                           choices = list("bray", "jaccard" = 'cc', "unifrac", "wunifrac", "dpcoa", "jsd", "euclidean")),
               menuSubItem("MultiDimensional Scaling", tabName = "betaMds"),
               menuSubItem("Samples clustering", tabName = "betaCluster"),
               menuSubItem("Multivariate ANOVA", tabName = "betaManova"),
               menuSubItem("Network", tabName = "betaNetwork"),
               menuSubItem("Samples heatmap", tabName = "betaHeatmap"),
               menuSubItem("Table", tabName = "betaTable")
      ),
      menuItem("PCA", tabName = "pca", icon = icon("bullseye")),
      menuItem("Differential expression analysis", tabName = "deseq", icon = icon("balance-scale-left")),
      menuItem("Phylogenetic tree", tabName = "tree", icon = icon("tree")),
      menuItem("Help", tabName = "Help", icon = icon("info-circle"))
  )),
  dashboardBody(
    tabItems(
      tabItem(tabName = "Summary", Summary),
      tabItem(tabName = "otuTable", otuTable),
      tabItem(tabName = "taxtable", taxtable),
      tabItem(tabName = "glomTable", glomTable),
      tabItem(tabName = "sampleTable", sampleTable),
      tabItem(tabName = "barplot", barplot),
      tabItem(tabName = "heatmap", heatmap),
      tabItem(tabName = "rarefactionCurve", rarefactionCurve),
      tabItem(tabName = "alphaPlot", alphaPlot),
      tabItem(tabName = "alphaTable", alphaTable),
      tabItem(tabName = "alphaAnova", alphaAnova),
      tabItem(tabName = "betaMds", betaMds),
      tabItem(tabName = "betaCluster", betaCluster),
      tabItem(tabName = "betaManova", betaManova),
      tabItem(tabName = "betaNetwork", betaNetwork),
      tabItem(tabName = "betaHeatmap", betaHeatmap),
      tabItem(tabName = "betaTable", betaTable),
      tabItem(tabName = "pca", pca),
      tabItem(tabName = "deseq", deseq),
      tabItem(tabName = "tree", tree),
      tabItem(tabName = "Help", Help)
      )
    )
  )
)
