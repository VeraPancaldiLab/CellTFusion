.onLoad <- function(libname, pkgname) {
  results_dir <- file.path("Results")
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("DeconRNASeq", force = TRUE)
}
