test_that("compute.latent_factors() splits the scores and returns non-negative factors", {
  cg <- cg_tuto
  lf <- suppressMessages(compute.latent_factors(cg$Cell_groups, return = FALSE))
  expect_true(all(c("Z", "W", "nmf_input", "nmf_model", "patient_mixture") %in% names(lf)))
  expect_identical(rownames(lf$Z), rownames(cg$Cell_groups))
  expect_true(all(lf$Z >= 0) && all(lf$W >= 0))
  expect_identical(colnames(lf$nmf_input), c(paste0(colnames(cg$Cell_groups), "_pos"), paste0(colnames(cg$Cell_groups), "_neg")))
  # reproducible with the same seed
  lf2 <- suppressMessages(compute.latent_factors(cg$Cell_groups, return = FALSE))
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

test_that("compute_cells_niches() returns cell types per factor", {
  withr::local_dir(withr::local_tempdir()) # compute_cells_niches() creates Results/
  cg <- cg_tuto
  lf <- lf_tuto
  niches <- suppressMessages(compute_cells_niches(lf, deconv_subgroups.tuto, cg, return = FALSE))
  expect_type(niches, "list")
  expect_true(all(names(niches) %in% colnames(lf$Z)))
  expect_true(all(vapply(niches, function(x) is.numeric(x) && !is.null(names(x)), logical(1))))
})
