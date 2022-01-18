output$deseqContrastVarUI <- renderUI({
  validate(need(physeq(), ""))
  selectInput(
    "deseqContrastVar",
    label = "Experimental design : ",
    choices = c(sample_variables(physeq()))
  )
})

output$deseqContrastModUI <- renderUI({
  validate(need(physeq(), ""), 
           need(input$deseqContrastVar, ""),
           need(class(get_variable(physeq(), input$deseqContrastVar)) != "numeric", "")
           )
  checkboxGroupInput(
    "deseqContrastMod",
    label = "Contrast (exactly two required) : ",
    choices = NULL,
    inline = TRUE
  )
})

observe({
  validate(need(physeq(), ""), need(input$deseqContrastVar, ""))
  var <- levels(as.factor(get_variable(physeq(), input$deseqContrastVar)))
  updateCheckboxGroupInput(session,
                           inputId = "deseqContrastMod",
                           choices = var,
                           selected = var[c(1, 2)],
                           inline = TRUE
  )
})

output$deseqTitleUI <- renderUI({
  validate(need(physeq(), ""))
  textInput("deseqTitle",
            label = "Title : ",
            value = "Volcano Plot")
})

output$deseqUI <- renderUI({
  validate(need(physeq(), ""))
  box(
    title = "Setting : " ,
    width = NULL,
    status = "primary",
    uiOutput("deseqContrastVarUI"),
    uiOutput("deseqContrastModUI"),
    # actionButton("deseqButton", label = "Execute", icon = icon("check"), class = "btn-primary"),
    uiOutput("deseqTitleUI")
  )
})

design <- metaReactive2({
  req(input$deseqContrastVar)
  metaExpr({as.formula(..(paste("~", input$deseqContrastVar)))})
})

cds <- metaReactive2({
  req(design())
  req(physeq())
  
  data <- physeq()
  metaExpr({phyloseq_to_deseq2(data, design = ..(design()))})
  })

dds <- metaReactive2({
  req(cds())
  
  metaExpr({DESeq2::DESeq(..(cds()), sfType = "poscounts")})
})

results <- metaReactive2({
  req(input$deseqContrastVar)
  req(dds())
  req(physeq())
  req({class(get_variable(physeq(), input$deseqContrastVar)) == "numeric" || length(input$deseqContrastMod) == 2})
  
  data <- physeq()

  if (class(get_variable(data, input$deseqContrastVar)) == "numeric") {
    # First case: regression against a continuous variable
    metaExpr({
      r <- DESeq2::results(object = ..(dds()),
                           name = ..(input$deseqContrastVar),
                           tidy = TRUE) %>%
        as_tibble() %>%
        rename(OTU = row) %>%
        inner_join(tax_table(data) %>% as("matrix") %>% as_tibble(rownames = "OTU"), by = "OTU")
      comment(r) <- ..(paste0("You compare low and high values of the continuous variable ", input$deseqContrastVar, ".\nA positive log2FoldChange means more abundant for high values of ", input$deseqContrastVar, "."))
      r
      }, bindToReturn = TRUE)
    } else {
      validate(need(length(input$deseqContrastMod) == 2, "Invalid input."))
      
      if (length(levels(as.factor(get_variable(data, input$deseqContrastVar)))) == 2) {
        # Second case: regression against a binary variable
        metaExpr({
          r <- DESeq2::results(object = ..(dds()),
                               name = DESeq2::resultsNames(..(dds()))[-1],
                               tidy = TRUE) %>%
            as_tibble() %>% rename(OTU = row) %>%
            inner_join(tax_table(data) %>% as("matrix") %>% as_tibble(rownames = "OTU"), by = "OTU")
          comment(r) <- ..(paste0("You compare ", input$deseqContrastMod[1], " to ", input$deseqContrastMod[2], " for the variable ", input$deseqContrastVar, ".\nA positive log2FoldChange means more abundant in ", input$deseqContrastMod[2], " than in ", input$deseqContrastMod[1], "."))
          r
          }, bindToReturn = TRUE)
        } else {
          # Third case: regression against a qualiative variable with three or more levels
          metaExpr({
            r <- DESeq2::results(object = ..(dds()),
                                 contrast = ..(c(input$deseqContrastVar, input$deseqContrastMod[1], input$deseqContrastMod[2])),
                                 tidy = TRUE) %>%
              as_tibble() %>% rename(OTU = row) %>%
              inner_join(tax_table(data) %>% as("matrix") %>% as_tibble(rownames = "OTU"), by = "OTU")
            comment(r) <- ..(paste0("You choose to compare ", input$deseqContrastMod[1], " to ", input$deseqContrastMod[2], " for the variable ", input$deseqContrastVar, ".\nA positive log2FoldChange means more abundant in ", input$deseqContrastMod[1], " than in ", input$deseqContrastMod[2], "."))
            r
            }, bindToReturn = TRUE)
        }
    }
})


output$deseq <- metaRender2(renderPlot, {
  validate(
    need(physeq(), "Requires an abundance dataset"),
    need(class(get_variable(physeq(), input$deseqContrastVar)) == "numeric" || 
           length(input$deseqContrastMod) == 2, "Requires a continuous design or a selection of two modalities for a discrete design.")#,
    #need(class(results()) == "DESeqResults", "Invalid input.")
    )
  data <- physeq()
  
  deseqPlot <- metaExpr({
    ggplot(..(results()) %>% mutate(evidence = -log10(padj), 
                              evolution = case_when(
                                padj <= 0.05 & log2FoldChange < 0 ~ "Down", 
                                padj <= 0.05 & log2FoldChange > 0 ~ "Up", 
                                TRUE                              ~ "Not DA"
                              )),
           aes(x = log2FoldChange, y = evidence)) + 
      geom_point(aes(color = evolution), size = 1.75, alpha = 0.8, na.rm = T) +     # base layer 
      theme_bw(base_size = 16) +                                                    # clean up theme
      theme(legend.position = "none",                                               # remove legend 
            plot.subtitle = element_text(size = 12)) +                              # add subtitle
      ggtitle(label = ..(input$deseqTitle), subtitle = comment(..(results()))) +              # add informative title
      xlab(expression(log[2]("FoldChange"))) +                                      # x-axis label
      ylab(expression(-log[10]("adjusted p-value"))) +                              # y-axis label
      geom_vline(xintercept = 0, colour = "grey80", linetype = 2) +                 # add line at 0
      geom_hline(yintercept = -log10(0.05), colour = "grey80", linetype = 2) +
      scale_color_manual(values = c("Down" = "red", "Not DA" = "grey20", "Up" = "green")) # change colors
  })
  
  metaExpr({
    p <- ..(deseqPlot)
    p
  })
})

observeEvent(input$deseq_output_code,
             {
               displayCodeModal(
                 expandChain(
                   quote(library(phyloseq)),
                   quote(library(phyloseq.extended)),
                   quote(library(DESeq2)),
                   quote(library(ggplot2)),
                   quote(library(magrittr)),
                   quote(library(dplyr)),
                   
                   "# Replace `data` with you own data.",
                   output$deseq()
                 ), clip = NULL
               )
             }
)

output$deseqTable <- renderDT({
  validate(
    need(physeq(), "Requires an abundance dataset"),
    need(class(get_variable(physeq(), input$deseqContrastVar)) == "numeric" || 
           length(input$deseqContrastMod) == 2, "Requires a continuous design or a selection of two modalities for a discrete design."),
    need(results(), "Invalid input.")
  )
  
  results() %>%
    filter(padj <= 0.05) %>% ## Only significant OTUs
    datatable(rownames = FALSE,
              filter = "top",
              extensions = c("Buttons", "ColReorder", "FixedColumns"),
              options = list(dom = "Btlip",
                             pageLength = 10,
                             lengthMenu = list(c(10, 25, 50, 100, -1), list('10', '25', '50', '100', 'All')),
                             buttons = list('colvis',
                                            list(extend = 'collection',
                                                 buttons = c('copy', 'csv', 'excel', 'pdf'),
                                                 text = 'Download')),
                             colReorder = TRUE,
                             scrollX = TRUE,
                             fixedColumns = list(leftColumns = 1, rightColumns = 0)
                             ),
              width = "auto",
              height = "auto") %>%
    formatSignif(columns = c("baseMean", "log2FoldChange", "lfcSE", "stat", "padj", "pvalue"), digits = 4) %>%
    formatStyle(columns = "log2FoldChange", color = DT::styleInterval(0, c('red', 'green')))
  })

output$deseqBrushed <- renderDT({
  validate(
    need(physeq(), "Requires an abundance dataset"),
    need(class(get_variable(physeq(), input$deseqContrastVar)) == "numeric" || 
           length(input$deseqContrastMod) == 2, "Requires a continuous design or a selection of two modalities for a discrete design."),
    need(results(), "Invalid input."),
    need(input$deseq_brush, "")
  )
  
  results() %>% mutate(evidence = -log10(padj),
                       evolution = case_when(
                         padj <= 0.05 & log2FoldChange < 0 ~ "Down",
                         padj <= 0.05 & log2FoldChange > 0 ~ "Up", 
                         TRUE                              ~ "Not DA"
                         )) %>%
    brushedPoints(brush = input$deseq_brush) %>%
    #filter(padj <= 0.05) %>% ## Only significant OTUs
    datatable(rownames = FALSE,
              filter = "top",
              extensions = c("Buttons", "ColReorder", "FixedColumns"),
              options = list(dom = "Btlip",
                             pageLength = 10,
                             lengthMenu = list(c(10, 25, 50, 100, -1), list('10', '25', '50', '100', 'All')),
                             buttons = list('colvis',
                                            list(extend = 'collection',
                                                 buttons = c('copy', 'csv', 'excel', 'pdf'),
                                                 text = 'Download')),
                             colReorder = TRUE,
                             scrollX = TRUE,
                             fixedColumns = list(leftColumns = 1, rightColumns = 0)
              ),
              width = "auto",
              height = "auto") %>%
    formatSignif(columns = c("baseMean", "log2FoldChange", "lfcSE", "stat", "padj", "pvalue"), digits = 4) %>%
    formatStyle(columns = "log2FoldChange", color = DT::styleInterval(c(-0.5, +0.5), c('red', 'black', 'green')))
})
