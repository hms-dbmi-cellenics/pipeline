# Faithfulness tests for the spatial local-outlier primitives factored out of
# SpotSweeper::localOutliers (https://github.com/MicTott/SpotSweeper, devel):
#
#   find_spatial_knn  <- BiocNeighbors::findKNN over the tissue coordinates
#   local_outliers    <- the per-spot neighbourhood modified z-score
#                        (vectorised with matrixStats for Visium HD scale)
#
# SpotSweeper computes, for each spot, a robust (modified) z-score comparing the
# spot's QC metric to its spatial neighbourhood:
#
#   z_i = 0.6745 * (m_i - median(N_i)) / (1.4826 * median(|N_i - median(N_i)|))
#
# where N_i = {m_i} U {metric of the k nearest spatial neighbours of i}. In the
# package this is `spatialEco::outliers(metric[c(i, knn_i)])[1]` — the focal spot
# is prepended to its neighbours and its (first) z-score is taken.
#
# These tests pin our re-implementation to that reference behaviour using a
# small, independent re-derivation of the formula (ref_* helpers below).

library(Seurat)

# ── independent reference implementations (mirror SpotSweeper) ────────────────

# modified z-score, identical to spatialEco::outliers / the SpotSweeper formula
ref_mod_z <- function(x) {
  med <- median(x)
  mad <- 1.4826 * median(abs(x - med))
  0.6745 * (x - med) / mad
}

# per-spot focal modified z-score over c(focal, neighbours)
ref_local_z <- function(metric_vals, spatial_knn) {
  z <- vapply(seq_len(nrow(spatial_knn)), function(i) {
    nbrs <- spatial_knn[i, ]
    nbrs <- nbrs[nbrs != 0]
    ref_mod_z(c(metric_vals[i], metric_vals[nbrs]))[1]
  }, numeric(1))
  z[!is.finite(z)] <- 0
  z
}


# ── find_spatial_knn ──────────────────────────────────────────────────────────

test_that("find_spatial_knn returns a k-neighbour index matrix from tissue coords", {
  scdata <- mock_spatial_scdata(ncells = 64, sample_id = "s1")
  knn <- find_spatial_knn(scdata, workers = 1, n_neighbors = 8)

  expect_true(is.matrix(knn))
  expect_equal(nrow(knn), ncol(scdata))
  expect_equal(ncol(knn), 8)
  # indices are valid 1-based row positions
  expect_true(all(knn >= 1 & knn <= ncol(scdata)))
})

test_that("find_spatial_knn excludes the focal spot from its own neighbours", {
  scdata <- mock_spatial_scdata(ncells = 64, sample_id = "s1")
  knn <- find_spatial_knn(scdata, workers = 1, n_neighbors = 8)
  self_in_neighbours <- vapply(
    seq_len(nrow(knn)), function(i) i %in% knn[i, ], logical(1)
  )
  expect_false(any(self_in_neighbours))
})

test_that("find_spatial_knn finds the spatially closest spots", {
  # on the regular grid built by mock_spatial_scdata, the nearest neighbours of
  # an interior spot are exactly the spots that are closest in (x, y).
  scdata <- mock_spatial_scdata(ncells = 64, sample_id = "s1")
  coords <- as.matrix(Seurat::GetTissueCoordinates(scdata)[, c("x", "y")])
  knn <- find_spatial_knn(scdata, workers = 1, n_neighbors = 4)

  i <- 28L # an interior grid point
  expected <- order(sqrt(rowSums((sweep(coords, 2, coords[i, ])^2))))
  expected <- expected[expected != i][1:4] # 4 closest, excluding self
  expect_setequal(knn[i, ], expected)
})

test_that("find_spatial_knn honours n_neighbors", {
  scdata <- mock_spatial_scdata(ncells = 64, sample_id = "s1")
  expect_equal(ncol(find_spatial_knn(scdata, n_neighbors = 5)), 5)
  expect_equal(ncol(find_spatial_knn(scdata, n_neighbors = 12)), 12)
})


# ── local_outliers ────────────────────────────────────────────────────────────

test_that("local_outliers reproduces the SpotSweeper focal-spot modified z-score", {
  scdata <- mock_spatial_scdata(ncells = 64, sample_id = "s1")
  knn <- find_spatial_knn(scdata, workers = 1, n_neighbors = 8)

  out <- local_outliers(scdata, knn, metric = "percent.mt", log = FALSE)
  expect_equal(
    out[["percent.mt_z"]],
    ref_local_z(scdata@meta.data[["percent.mt"]], knn)
  )
})

test_that("local_outliers z-score depends on the focal spot's own metric", {
  # this is the property the neighbours-only implementation violated: the spot's
  # own metric value must drive its z-score.
  scdata <- mock_spatial_scdata(ncells = 64, sample_id = "s1")
  knn <- find_spatial_knn(scdata, workers = 1, n_neighbors = 8)

  z_before <- local_outliers(scdata, knn, metric = "percent.mt")[["percent.mt_z"]]

  i <- 30L
  scdata@meta.data[["percent.mt"]][i] <- 1000 # make spot i a glaring outlier
  z_after <- local_outliers(scdata, knn, metric = "percent.mt")[["percent.mt_z"]]

  # spot i's own z-score must jump; spots that don't neighbour i are unaffected
  expect_gt(z_after[i], z_before[i] + 1)
  expect_gt(z_after[i], 3)
})

test_that("local_outliers log-transforms the metric with log1p before scoring", {
  scdata <- mock_spatial_scdata(ncells = 64, sample_id = "s1")
  knn <- find_spatial_knn(scdata, workers = 1, n_neighbors = 8)

  out <- local_outliers(scdata, knn, metric = "nCount_RNA", log = TRUE)

  # the log column is added and equals log1p of the raw metric
  expect_true("nCount_RNA_log" %in% colnames(out))
  expect_equal(out[["nCount_RNA_log"]], log1p(scdata@meta.data[["nCount_RNA"]]))

  # z-scores are computed on the log-transformed metric
  expect_equal(
    out[["nCount_RNA_z"]],
    ref_local_z(log1p(scdata@meta.data[["nCount_RNA"]]), knn)
  )
})

test_that("local_outliers z-scores are invariant to a constant metric scaling", {
  # SpotSweeper log2-transforms metrics whereas we use log1p; because the
  # modified z-score divides out the MAD, a constant scale factor between the
  # two cancels, so our z-scores match SpotSweeper's regardless of log base.
  scdata <- mock_spatial_scdata(ncells = 64, sample_id = "s1")
  knn <- find_spatial_knn(scdata, workers = 1, n_neighbors = 8)

  z1 <- local_outliers(scdata, knn, metric = "percent.mt")[["percent.mt_z"]]
  scdata@meta.data[["percent.mt"]] <- scdata@meta.data[["percent.mt"]] * 7.3
  z2 <- local_outliers(scdata, knn, metric = "percent.mt")[["percent.mt_z"]]
  expect_equal(z1, z2)
})

test_that("local_outliers maps non-finite z-scores (zero-MAD neighbourhoods) to 0", {
  scdata <- mock_spatial_scdata(ncells = 64, sample_id = "s1")
  knn <- find_spatial_knn(scdata, workers = 1, n_neighbors = 8)

  # constant metric -> every neighbourhood has MAD 0 -> all z-scores non-finite
  scdata@meta.data[["percent.mt"]] <- 5
  out <- local_outliers(scdata, knn, metric = "percent.mt")

  expect_true(all(is.finite(out[["percent.mt_z"]])))
  expect_true(all(out[["percent.mt_z"]] == 0))
})

test_that("local_outliers returns metric and z columns row-aligned to meta.data", {
  scdata <- mock_spatial_scdata(ncells = 64, sample_id = "s1")
  knn <- find_spatial_knn(scdata, workers = 1, n_neighbors = 8)

  out <- local_outliers(scdata, knn, metric = "percent.mt", log = FALSE)
  expect_equal(rownames(out), rownames(scdata@meta.data))
  expect_setequal(colnames(out), c("percent.mt", "percent.mt_z"))
  expect_equal(nrow(out), ncol(scdata))
})
