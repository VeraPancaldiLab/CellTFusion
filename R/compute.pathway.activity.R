compute.pathway.activity <- function(RNA.tpm, gene_sets = NULL, paths = NULL){

  RNA.tpm = as.matrix(RNA.tpm)
  #Get universe
  if(is.null(paths)){
    paths <- get_progeny(organism = 'human', top = 500)
  }

  # Run mlm
  progeny <- run_mlm(mat=RNA.tpm, net=paths, .source='source', .target='target', .mor='weight', minsize = 5)

  #Remove variable
  rm(paths)
  gc()

  # Transform to wide matrix
  sample_acts_progeny <- progeny %>%
    pivot_wider(id_cols = 'condition', names_from = 'source',
                values_from = 'score') %>%
    column_to_rownames('condition') %>%
    as.matrix()

  if(is.null(gene_sets)==F){

    cat("Computing GSVA analysis using provided gene sets.....................................................\n")

    gsva_results <- gsva(
      RNA.tpm,
      gene_sets,
      method = "gsva",
      kcdf = "Gaussian",
      min.sz = 1,
      mx.diff = TRUE,
      verbose = TRUE
    )
    sample_acts_hallmarks <- data.frame(scale(t(gsva_results)))
    return(list(sample_acts_progeny, sample_acts_hallmarks))
  }

  # Scale per feature
  sample_acts_progeny <- data.frame(scale(sample_acts_progeny))

  message("Pathways scores computed")
  return(sample_acts_progeny)

}
