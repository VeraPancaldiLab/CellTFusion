test_that("compute.latent_factors() splits the scores and returns non-negative factors", {
  cg <- cg_tuto
  lf <- suppressMessages(compute.latent_factors(cg$Cell_groups, return = FALSE))
  expect_true(all(c("Z", "W", "nmf_input", "nmf_model", "patient_mixture", "stability", "rank_selection") %in% names(lf)))
  # rank selection by bootstrap NMF: the selected rank is the largest with all factors robust and none redundant
  rs <- lf$rank_selection
  expect_identical(ncol(lf$Z), as.integer(if (any(rs$stable_distinct)) max(rs$rank[rs$stable_distinct]) else rs$rank[which.max(rs$robust_factors)]))
  expect_identical(nrow(lf$stability), ncol(lf$Z))
  expect_true(all(lf$stability$recovery >= 0 & lf$stability$recovery <= 1))
  expect_identical(rownames(lf$Z), rownames(cg$Cell_groups))
  expect_true(all(lf$Z >= 0) && all(lf$W >= 0))
  # positive part of the scores only (high-activity side), one column per cell group
  expect_identical(colnames(lf$nmf_input), colnames(cg$Cell_groups))
  expect_equal(unname(lf$nmf_input), unname(pmax(as.matrix(cg$Cell_groups), 0)))
  # reproducible with the same seed
  lf2 <- suppressMessages(compute.latent_factors(cg$Cell_groups, return = FALSE))
  # a given rank: no rank selection, stability of its factors only
  lf3 <- suppressMessages(compute.latent_factors(cg$Cell_groups, rank = 3, n_boot = 10, return = FALSE))
  expect_null(lf3$rank_selection)
  expect_identical(nrow(lf3$stability), 3L)
  expect_identical(lf$Z, lf2$Z)
  expect_error(compute.latent_factors(cbind(cg$Cell_groups, bad = NA), return = FALSE), "missing values")
})

test_that("project_factors() places the training samples close to their factors", {
  cg <- cg_tuto
  lf <- lf_tuto
  pz <- project_factors(lf, cg$Cell_groups)
  expect_identical(dim(pz), dim(lf$Z))
  Zd <- sweep(lf$Z, 2, lf$nmf_model$d, "*")
  expect_true(all(diag(stats::cor(pz, Zd)) > 0.9))
})

test_that("compute_cells_niches() returns the top cell groups of each factor with their CCA coefficients", {
  withr::local_dir(withr::local_tempdir()) # compute_cells_niches() creates Results/
  cg <- cg_tuto
  lf <- lf_tuto
  niches <- compute_cells_niches(lf, cg, top_n = 3, file_name = "test")
  expect_named(niches, c("factor", "cell_group", "module", "weight", "feature", "coefficient"))
  expect_setequal(unique(niches$factor), colnames(lf$W))
  for (f in colnames(lf$W)) {
    top <- names(sort(lf$W[, f], decreasing = TRUE))[1:3]
    expect_identical(unique(niches$cell_group[niches$factor == f]), top)
  }
  g <- niches$cell_group[1]
  expect_equal(niches$coefficient[niches$factor == niches$factor[1] & niches$cell_group == g], unname(cg$Weights[[g]]$xcoef[, 1]))
  expect_true(file.exists("Results/Cells_niches_test.pdf"))
  expect_true(file.exists("Results/Cells_niches_heatmap_test.pdf"))
})
