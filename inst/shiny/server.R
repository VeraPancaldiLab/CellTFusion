library(shiny)
library(CellTFusion)   # <-- CHANGE if your package name differs

server <- function(input, output, session) {

  # --------------------------------------------------
  # Store results between steps
  # --------------------------------------------------
  results <- reactiveValues(
    raw = NULL,
    counts = NULL,
    deconv = NULL,
    tfs = NULL,
    network = NULL,
    pathways = NULL,
    cell_groups = NULL,
    latent = NULL
  )

  # --------------------------------------------------
  # Load input data
  # --------------------------------------------------
  observeEvent(input$counts, {
    results$raw <- read.csv(
      input$counts$datapath,
      row.names = 1,
      check.names = FALSE
    )
  })

  # --------------------------------------------------
  # 1. Normalize counts
  # --------------------------------------------------
  observeEvent(input$run_norm, {
    req(results$raw)

    results$counts <- normalize_counts(
      raw.counts = results$raw,
      normalized = input$normalized
    )
  })

  # --------------------------------------------------
  # 2. Deconvolution
  # --------------------------------------------------
  observeEvent(input$run_deconv, {
    req(results$counts)

    results$deconv <- run_deconvolution(
      raw.counts = results$counts,
      normalized = input$normalized,
      deconv_methods = c("Quantiseq", "DWLS")
    )
  })

  # --------------------------------------------------
  # 3. TF activity
  # --------------------------------------------------
  observeEvent(input$run_tfs, {
    req(results$counts)

    results$tfs <- run_tf_activity(
      counts.norm = results$counts,
      TF.collection = "CollecTRI",
      min_targets_size = 10,
      tfs.pruned = FALSE,
      universe = NULL
    )
  })

  # --------------------------------------------------
  # 4. TF network
  # --------------------------------------------------
  observeEvent(input$run_network, {
    req(results$tfs)

    results$network <- run_tf_network(
      tfs = results$tfs,
      batch = FALSE,
      minMod = 10,
      corr_mod = 0.9,
      return = TRUE
    )
  })

  # --------------------------------------------------
  # Run EVERYTHING (wrapper)
  # --------------------------------------------------
  observeEvent(input$run_all, {
    req(results$raw)

    res <- CellTFusion(
      raw.counts = results$raw,
      normalized = input$normalized,
      verbose = TRUE
    )

    results$deconv <- res$Deconvolution
    results$tfs <- res$TFs_matrix
    results$network <- res$TF_network
    results$pathways <- res$Pathways_scores
    results$cell_groups <- res$Cell_groups
    results$latent <- res$Latent_spaces
  })

  # --------------------------------------------------
  # Outputs
  # --------------------------------------------------
  output$status <- renderPrint({
    list(
      raw_loaded = !is.null(results$raw),
      normalized = !is.null(results$counts),
      deconvolution = !is.null(results$deconv),
      tfs = !is.null(results$tfs),
      network = !is.null(results$network)
    )
  })

  output$deconv_out <- renderTable({
    req(results$deconv)
    head(results$deconv)
  })

  output$tfs_out <- renderTable({
    req(results$tfs)
    head(results$tfs)
  })

  output$network_out <- renderPrint({
    req(results$network)
    str(results$network)
  })

}
