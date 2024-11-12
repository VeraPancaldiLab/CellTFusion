#' Computes TF-modules pathway activities
#'
#'Pathways activity computation is done by using a multivariate linear model (mlm) and a resource that leverages a large compendium of publicly available signaling perturbation experiments to yield a common core of pathway responsive genes for human PROGENy. Users can also select to do an additional Gene Set Variation Analysis (GSVA) using hallmark signatures or predefined gene sets.
#'
#' @param RNA.tpm gene expression matrix normalized with genes as rows and samples as columns.
#' @param gene_sets either hallmark signatures extracted from the Molecular Signature Database (MSigDB) [2] or predefined gene sets signatures. Default is NULL.
#' @param paths whether external database wants to be provided for computing the pathway activities
#'
#' @return A matrix with PRGENy pathway activities (samples as rows and pathways as columns). If gene_sets is not NULL it will also return a matrix containing the signatures scores (samples as rows and signatures as columns).
#'
#' @export
#'
#' @examples
#'
#' pathways = compute.pathway.activity(counts.norm)
#'
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
