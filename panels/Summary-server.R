output$phyloseqPrint <- renderPrint({
  validate(
    need(
      physeq(),
      "Firstly, you should select a demo dataset or upload an abundance BIOM file.\nFor example, with Galaxy, a BIOM file can be obtained at the end of FROGS workflow with the 'FROGS BIOM to std BIOM' tool. \nMake sure that the phyloseq object in the RData file is called 'data'."
    )
  )
  physeq()
})
