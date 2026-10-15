# Tutorial data shared by the tests (no download needed)
data("tfs.tuto", package = "CellTFusion")
data("network.tuto", package = "CellTFusion")
data("deconv_subgroups.tuto", package = "CellTFusion")
data("deconv.tuto", package = "CellTFusion")
data("traitdata.tuto", package = "CellTFusion")
data("counts.norm.tuto", package = "CellTFusion")

# Two pseudo-cohorts of the tutorial samples, with TF activities scaled within each cohort (as CellTFusion() does)
tuto_batch <- stats::setNames(rep(c("A", "B"), length.out = nrow(tfs.tuto)), rownames(tfs.tuto))
tfs_batch <- tfs.tuto
for (cohort in unique(tuto_batch)) {
  tfs_batch[tuto_batch == cohort, ] <- scale(tfs.tuto[tuto_batch == cohort, ])
}

# Cell groups and latent factors of the tutorial data, computed once for the tests
cg_tuto <- suppressMessages(construct_cell_groups(network.tuto, deconv_subgroups.tuto, return = FALSE))
lf_tuto <- suppressMessages(compute.latent_factors(cg_tuto$Cell_groups, return = FALSE))

# Small offline networks for TF and pathway scoring (genes of the tutorial expression data)
make_toy_network <- function(n_sources, n_targets, weight_col) {
  set.seed(1)
  genes <- rownames(counts.norm.tuto)[1:2000]
  net <- do.call(rbind, lapply(seq_len(n_sources), function(i) {
    data.frame(source = paste0("S", i), target = sample(genes, n_targets), stringsAsFactors = FALSE)
  }))
  net[[weight_col]] <- sample(c(-1, 1), nrow(net), replace = TRUE)
  net
}
