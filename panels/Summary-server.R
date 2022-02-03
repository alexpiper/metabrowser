output$rawphyloseqPrint <- renderPrint({
  validate(
    need(
      raw_physeq(),
      "Firstly, you should select a demo dataset or upload an abundance BIOM file.\nFor example, with Galaxy, a BIOM file can be obtained at the end of FROGS workflow with the 'FROGS BIOM to std BIOM' tool. \nMake sure that the phyloseq object in the RData file is called 'data'."
    )
  )
  print(paste0("Original data: ", phyloseq::nsamples(raw_physeq()), " Samples, ", phyloseq::ntaxa(raw_physeq()), " OTUs."))
})

output$phyloseqPrint <- renderPrint({
  validate(
    need(
      physeq(),
      ""
    )
  )
  print(paste0("Filtered data: ", phyloseq::nsamples(physeq()), " Samples, ", phyloseq::ntaxa(physeq()), " OTUs."))
})

# Display plot ---------------------------------------------------------
output$filtplot <- metaRender2(renderPlot, {
  validate(need(physeq(), "Requires an abundance dataset"))
  validate(need(raw_physeq(), "Requires an abundance dataset"))
  raw <- raw_physeq()
  filt <- physeq()
  
  # Create plot
  metaExpr({
    # Unfiltered data
    gg.lib_hist <- tibble::enframe(sample_sums(raw)) %>%
      ggplot(aes(x=value)) + 
      geom_histogram() +
      labs(x = "Library size", y = "Sample count") +
      ggtitle("Original library sizes")
    
    gg.otu_hist <- tibble::enframe(taxa_sums(raw)) %>%
      ggplot(aes(x=value)) + 
      geom_histogram() +
      labs(x = "Abundance", y = "OTU count") +
      ggtitle("Original OTU abundance")
    
    # Filtered data
    gg.lib_hist_filtered <- tibble::enframe(sample_sums(filt)) %>%
      ggplot(aes(x=value)) + 
      geom_histogram() +
      labs(x = "Library size", y = "Sample count") +
      ggtitle("Filtered library sizes")
    
    gg.otu_hist_filtered  <- tibble::enframe(taxa_sums(filt)) %>%
      ggplot(aes(x=value)) + 
      geom_histogram() +
      labs(x = "Abundance", y = "OTU count")+
      ggtitle("Filtered OTU abundance")
    
    gridExtra::grid.arrange(gg.lib_hist, gg.otu_hist, gg.lib_hist_filtered, gg.otu_hist_filtered, ncol=2) #, main="Histograms: Before and After Filtering")
  })
})
