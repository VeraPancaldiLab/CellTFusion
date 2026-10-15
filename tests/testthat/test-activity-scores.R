test_that("compute.TFs.activity() raw scores do not depend on the other samples", {
  universe <- make_toy_network(20, 40, "mor")
  X <- counts.norm.tuto[, 1:40]
  all_raw <- suppressMessages(compute.TFs.activity(X, universe = universe, scale = FALSE, return = FALSE))
  part_raw <- suppressMessages(compute.TFs.activity(X[, 1:10], universe = universe, scale = FALSE, return = FALSE))
  expect_equal(as.matrix(part_raw), as.matrix(all_raw[1:10, colnames(part_raw)]), tolerance = 1e-8, ignore_attr = TRUE)
  scaled <- suppressMessages(compute.TFs.activity(X, universe = universe, return = FALSE))
  expect_equal(unname(apply(scaled, 2, stats::sd)), rep(1, ncol(scaled)), tolerance = 1e-8)
})

test_that("compute.pathway.activity() with batch scales within each cohort and keeps the sample order", {
  paths <- make_toy_network(6, 60, "weight")
  X <- counts.norm.tuto[, 1:60]
  b <- rep(c("A", "B", "C"), length.out = ncol(X))
  pb <- suppressMessages(compute.pathway.activity(X, paths = paths, batch = b, return = FALSE))
  expect_identical(rownames(pb), colnames(X))
  for (cohort in unique(b)) {
    expect_equal(unname(colMeans(pb[b == cohort, ])), rep(0, ncol(pb)), tolerance = 1e-8)
  }
  per_cohort <- do.call(rbind, unname(lapply(split(seq_len(ncol(X)), b), function(i)
    suppressMessages(compute.pathway.activity(X[, i], paths = paths, return = FALSE)))))
  expect_equal(pb, per_cohort[colnames(X), ], tolerance = 1e-8)
  # a single cohort is the same as no batch
  expect_identical(suppressMessages(compute.pathway.activity(X, paths = paths, batch = rep("A", ncol(X)), return = FALSE)),
                   suppressMessages(compute.pathway.activity(X, paths = paths, return = FALSE)))
})
