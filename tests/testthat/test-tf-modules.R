test_that("compute.WTCNA() builds modules with the documented defaults and stores its parameters", {
  net <- suppressMessages(compute.WTCNA(tfs.tuto, return = FALSE, verbose = FALSE))
  expect_named(net, c("TFs module matrix", "TFs colors", "TFs per module", "Proportion of variance", "TFs_matrix", "Network parameters"))
  expect_identical(rownames(net[["TFs module matrix"]]), rownames(tfs.tuto))
  expect_identical(colnames(net[["TFs module matrix"]]), names(net[["TFs per module"]]))
  expect_false("grey" %in% names(net[["TFs per module"]]))
  expect_equal(net[["Network parameters"]]$corr_type, "spearman")
  expect_true(all(lengths(net[["TFs per module"]]) >= 5))
  # module scores are z-scored
  expect_equal(unname(colMeans(net[["TFs module matrix"]])), rep(0, ncol(net[["TFs module matrix"]])), tolerance = 1e-8)
})

test_that("compute.WTCNA() results do not depend on the scaling of the TF activities", {
  raw <- as.data.frame(sweep(sweep(as.matrix(tfs.tuto), 2, seq_len(ncol(tfs.tuto)) / 100, "*"), 2, 3, "+"))
  n_scaled <- suppressMessages(compute.WTCNA(tfs.tuto, return = FALSE, verbose = FALSE))
  n_raw <- suppressMessages(compute.WTCNA(raw, return = FALSE, verbose = FALSE))
  expect_identical(n_raw[["TFs per module"]], n_scaled[["TFs per module"]])
  expect_equal(n_raw[["TFs module matrix"]], n_scaled[["TFs module matrix"]], tolerance = 1e-6)
})

test_that("compute.WTCNA() batch mode takes a cohort vector and returns scores in the input order", {
  net <- suppressMessages(compute.WTCNA(tfs_batch, batch = tuto_batch, return = FALSE, verbose = FALSE))
  M <- net[["TFs module matrix"]]
  expect_identical(rownames(M), rownames(tfs_batch))
  expect_true(is.list(net$TFs_matrix) && !is.data.frame(net$TFs_matrix))
  for (cohort in unique(tuto_batch)) {
    expect_equal(unname(colMeans(M[tuto_batch == cohort, , drop = FALSE])), rep(0, ncol(M)), tolerance = 1e-8)
  }
  # a single cohort is the same as no batch
  n1 <- suppressMessages(compute.WTCNA(tfs.tuto, batch = rep("A", nrow(tfs.tuto)), return = FALSE, verbose = FALSE))
  n0 <- suppressMessages(compute.WTCNA(tfs.tuto, return = FALSE, verbose = FALSE))
  expect_identical(n1, n0)
})

test_that("compute.WTCNA() validates its arguments", {
  expect_error(compute.WTCNA(tfs.tuto, corr_type = "s", return = FALSE, verbose = FALSE), "corr_type")
  expect_error(compute.WTCNA(tfs_batch, batch = tuto_batch, consensus_quantile = 2, return = FALSE, verbose = FALSE), "consensus_quantile")
  expect_error(compute.WTCNA(tfs.tuto, batch = c("A", "B"), return = FALSE, verbose = FALSE), "one value per sample")
})

test_that("identify_hub_TFs() returns hubs that belong to their module", {
  hubs <- suppressMessages(identify_hub_TFs(network.tuto))
  expect_named(hubs, c("hubGenes", "detailedData"))
  for (m in names(hubs$hubGenes)) {
    expect_true(all(hubs$hubGenes[[m]] %in% network.tuto[["TFs per module"]][[m]]))
  }
})
