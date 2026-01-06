library(shiny)
library(CellTFusion)

server <- function(input, output, session) {

  results <- reactiveValues(
    raw = NULL,
    res = NULL
  )

  ## Load data
  observeEvent(input$counts, {
    results$raw <- read.csv(
      input$counts$datapath,
      row.names = 1,
      check.names = FALSE
    )
  })

  ## Run pipeline ONCE
  observeEvent(input$run_all, {

    req(results$raw)

    results$res <- CellTFusion(
      raw.counts = results$raw,

      ## normalization
      normalized = input$normalized,

      ## deconvolution
      deconv_methods = input$deconv_methods,

      ## TF activity
      TF.collection = input$TF.collection,
      min_targets_size = input$min_targets_size,
      tfs.pruned = input$tfs.pruned,

      ## network
      minMod = input$minMod,
      corr_mod = input$corr_mod,

      ## general
      verbose = TRUE,
      return = TRUE
    )
  })

  ## Outputs
  output$status <- renderPrint({
    list(
      data_loaded = !is.null(results$raw),
      results_ready = !is.null(results$res)
    )
  })

  output$deconv_out <- renderTable({
    req(results$res)
    head(results$res$Deconvolution)
  })

  output$tfs_out <- renderTable({
    req(results$res)
    head(results$res$TFs_matrix)
  })

  output$network_out <- renderPrint({
    req(results$res)
    str(results$res$TF_network)
  })
}
