test_that("construct_cell_groups() returns scores, compositions and weights", {
  cg <- cg_tuto
  expect_named(cg, c("Cell_groups", "Composition", "Weights"))
  expect_identical(colnames(cg$Cell_groups), names(cg$Composition))
  expect_identical(names(cg$Weights), names(cg$Composition))
  expect_true(all(grepl("^Dendrogram_.+_group_[0-9]+$", names(cg$Composition))))
  expect_named(cg$Weights[[1]], c("xcoef", "train_means", "train_sds"))
})

test_that("each cell group score is positively related to its TF module", {
  cg <- cg_tuto
  modules <- sub("^Dendrogram_(.*)_group_.*$", "\\1", colnames(cg$Cell_groups))
  M <- network.tuto[["TFs module matrix"]][rownames(cg$Cell_groups), ]
  r <- mapply(function(g, m) stats::cor(cg$Cell_groups[[g]], M[, m]), colnames(cg$Cell_groups), modules)
  expect_true(all(r > 0))
})

test_that("Pearson correlation reproduces the same kind of output and Spearman is the default", {
  cg_p <- suppressMessages(construct_cell_groups(network.tuto, deconv_subgroups.tuto, corr_type = "pearson", return = FALSE))
  expect_gt(ncol(cg_p$Cell_groups), 0)
  expect_error(construct_cell_groups(network.tuto, deconv_subgroups.tuto, corr_type = "p", return = FALSE), "corr_type")
})

test_that("construct_cell_groups() requires the same samples in the network and the deconvolution", {
  dsub <- deconv_subgroups.tuto
  dsub[[1]] <- dsub[[1]][1:150, ]
  expect_error(construct_cell_groups(network.tuto, dsub, return = FALSE), "same samples")
})

test_that("compute.test.set() reproduces the training scores of the training samples (single cohort)", {
  cg <- cg_tuto
  groups <- names(cg$Composition)[1:10]
  test <- suppressWarnings(compute.test.set(deconv_subgroups.tuto, cg, group_names = groups, deconv_test = deconv.tuto))
  expect_equal(as.matrix(test), as.matrix(cg$Cell_groups[rownames(test), colnames(test)]), tolerance = 1e-8, ignore_attr = TRUE)
  # a single sample can be scored with the training scaling
  one <- suppressWarnings(compute.test.set(deconv_subgroups.tuto, cg, group_names = groups, deconv_test = deconv.tuto[1, , drop = FALSE]))
  expect_equal(unlist(one), unlist(cg$Cell_groups[rownames(one), colnames(one)]), tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("batch mode: features z-scored per cohort, only xcoef stored, cohorts re-projected exactly", {
  nb <- suppressWarnings(suppressMessages(compute.WTCNA(tfs_batch, batch = tuto_batch, return = FALSE, verbose = FALSE))) # foreach: "executing %dopar% sequentially"
  cb <- suppressMessages(construct_cell_groups(nb, deconv_subgroups.tuto, batch = tuto_batch, return = FALSE))
  expect_identical(unique(unlist(lapply(cb$Weights, names))), "xcoef")
  S <- as.matrix(cb$Cell_groups)
  b <- tuto_batch[rownames(S)]
  for (cohort in unique(b)) {
    expect_equal(unname(colMeans(S[b == cohort, , drop = FALSE])), rep(0, ncol(S)), tolerance = 1e-8)
  }
  iA <- names(b)[b == "A"]
  pA <- suppressWarnings(compute.test.set(deconv_subgroups.tuto, cb, deconv_test = deconv.tuto[iA, ]))
  expect_equal(as.matrix(pA), S[iA, colnames(pA)], tolerance = 1e-8, ignore_attr = TRUE)
  # few test samples: warning
  expect_warning(compute.test.set(deconv_subgroups.tuto, cb, deconv_test = deconv.tuto[iA[1:5], ]), "fewer than 10")
})
