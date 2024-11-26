dir.create(file.path(getwd(), "Results"))

test_that("TFs activity inference works", {
  tfs_test <- read.csv("../TFs_test.csv", row.names = 1)

  raw.counts = read.csv("../Counts_test.csv", row.names = 1)
  counts.norm = data.frame(ADImpute::NormalizeTPM(raw.counts, log = T))
  tfs = compute.TFs.activity(counts.norm)

  expect_equal(
    info = "rows of TFs equal to columns of counts (same samples in same order)",
    object = sort(rownames(tfs)), expected = sort(rownames(tfs_test))
  )
  expect_equal(
    info = "columns contains same TFs as in test",
    object = colnames(tfs), expected = colnames(tfs_test)
  )
  expect_equal(
    info = "TFs result is correct", object = tfs,
    expected = tfs_test, tolerance = 1e-1
  )

})

test_that("Network construction works", {
  network_test <- readRDS("../network.rds")

  tfs = read.csv("../TFs_test.csv", row.names = 1)
  network = compute.WTCNA(tfs, corr_mod = 0.9, clustering.method = "ward.D2", return = T)

  expect_equal(
    info = "number of elements from network is the same as in test",
    object = length(network), expected = length(network_test)
  )

  expect_equal(
    info = "number of samples in ME module matrix is the same",
    object = nrow(network[[1]]), expected = nrow(network_test[[1]])
  )

  expect_equal(
    info = "number of colors is the same as the total number of TFs",
    object = length(network[[2]]), expected = ncol(tfs)
  )

  expect_equal(
    info = "TFs module matrix is correct", object = network[[1]],
    expected = network_test[[1]], tolerance = 1e-1
  )

})

test_that("HubTFs identification works", {
  hubTFs_test <- readRDS("../hubTFs.rds")

  tfs_test <- read.csv("../TFs_test.csv", row.names = 1)
  network_test <- readRDS("../network.rds")

  hub_tfs = identify_hub_TFs(t(tfs_test), network_test, MM_thresh = 0.8, degree_thresh = 0.9)

  expect_equal(
    info = "number of elements from hubTFs is the same as in test",
    object = length(hub_tfs), expected = length(hubTFs_test)
  )

  expect_equal(
    info = "number of elements from hubTFs is the same as number of modules",
    object = length(hub_tfs[[1]]), expected = length(network_test[[3]])
  )

  expect_equal(
    info = "Detailed data of hubTFs is correct", object = hub_tfs[[2]],
    expected = hubTFs_test[[2]], tolerance = 1e-1
  )

})

test_that("Pathways activity works", {
  pathways_test <- read.csv("../Pathways_test.csv", row.names = 1)

  raw.counts = read.csv("../Counts_test.csv", row.names = 1)
  counts.norm = data.frame(ADImpute::NormalizeTPM(raw.counts, log = T))
  pathways = compute.pathway.activity(counts.norm)

  expect_equal(
    info = "rows of pathways equal (same samples as in test)",
    object = rownames(pathways), expected = rownames(pathways_test)
  )

  expect_equal(
    info = "Pathways result is correct", object = pathways,
    expected = pathways_test, tolerance = 1e-1
  )

})

test_that("TF network classification works", {
  tfs.modules.clusters_test <- readRDS("../TF_clusters.rds")

  pathways_test <- read.csv("../Pathways_test.csv", row.names = 1)
  network_test <- readRDS("../network.rds")

  tfs.modules.clusters = compute.TF.network.classification(network_test, pathways_test, return = T)

  expect_equal(
    info = "number of elements is the same as in test",
    object = length(tfs.modules.clusters), expected = length(tfs.modules.clusters_test)
  )

  expect_equal(
    info = "names of elements are the same",
    object = names(tfs.modules.clusters), expected = names(tfs.modules.clusters_test)
  )

  expect_equal(
    info = "cluster 1 is the same",
    object = tfs.modules.clusters[[1]], expected = tfs.modules.clusters_test[[1]]
  )

  expect_equal(
    info = "cluster 2 is the same",
    object = tfs.modules.clusters[[2]], expected = tfs.modules.clusters_test[[2]]
  )

})
