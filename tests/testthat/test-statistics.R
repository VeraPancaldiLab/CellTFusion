test_that("compute.modules.relationship() returns correlations and BH adjusted p-values", {
  rel <- compute.modules.relationship(network.tuto[[1]], deconv_subgroups.tuto[[1]], return = FALSE)
  expect_named(rel, c("correlation", "significant", "p_value", "p_adjusted"))
  expect_identical(rel$p_value, rel$p_adjusted)
  rel_bh <- compute.modules.relationship(network.tuto[[1]], deconv_subgroups.tuto[[1]], return = FALSE, padj = TRUE)
  expect_true(all(rel_bh$p_adjusted >= rel_bh$p_value, na.rm = TRUE))
  expect_true(sum(lengths(rel_bh$significant)) <= sum(lengths(rel$significant)))
  expect_error(compute.modules.relationship(network.tuto[[1]], deconv_subgroups.tuto[[1]], return = FALSE, padj = "BH"), "padj")
})

test_that("scores.stat.analysis() runs the chosen test and validates padj", {
  withr::local_dir(withr::local_tempdir()) # the scores.* tests save their plots in Results/
  dir.create("Results")
  scores <- network.tuto[[1]]
  res <- suppressMessages(scores.stat.analysis(scores, traitdata.tuto, trait = "Best.Confirmed.Overall.Response", test = "kruskal"))
  expect_true("Statistics" %in% names(res))
  expect_identical(res$Statistics$feature, colnames(scores))
  res_bh <- suppressMessages(scores.stat.analysis(scores, traitdata.tuto, trait = "binaryResponse", test = "wilcox", padj = TRUE))
  expect_true(all(res_bh$Statistics$p_adjusted >= res_bh$Statistics$p_value, na.rm = TRUE))
  expect_error(scores.stat.analysis(scores, traitdata.tuto, trait = "binaryResponse", test = "wilcox", padj = "BH"), "padj")
})

test_that("compute.metadata.association() tests categorical traits without saving when no trait is numeric", {
  withr::local_dir(withr::local_tempdir())
  dir.create("Results")
  res <- suppressMessages(compute.metadata.association(network.tuto[[1]], traitdata.tuto[, 1:3], plot_grid = TRUE, padj = TRUE))
  expect_true(is.data.frame(res$anova))
  expect_true(all(res$anova$p_adjusted >= res$anova$p, na.rm = TRUE))
})

test_that("compute.survival.analysis() works on the tutorial overall survival", {
  skip_if_not_installed("survival")
  skip_if_not_installed("survminer")
  skip_if_not_installed("gridExtra")
  withr::local_dir(withr::local_tempdir())
  dir.create("Results")
  res <- suppressWarnings(suppressMessages(compute.survival.analysis(traitdata.tuto, time_column = "os", event_column = "censOS",
                                                                     group_column = "binaryResponse")))
  expect_true(is.numeric(res$p_value) && length(res$p_value) == 1)
  expect_error(compute.survival.analysis(traitdata.tuto, time_column = "os", event_column = "censOS", group_column = "binaryResponse", padj = "BH"), "padj")
})
