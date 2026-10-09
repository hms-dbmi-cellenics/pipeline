# mock_spatial_config() and add_zscores() live in helper-spatial.R

test_that("filter_spatial_umi_outlier removes cells below -cutoff (lower direction)", {
  scdata_list <- add_tissue_coords(mock_scdata_list())
  sample1_id <- names(scdata_list)[1]
  ncells <- ncol(scdata_list[[sample1_id]])

  # first 10 cells are strong low outliers, the rest are within range
  z <- rep(0, ncells)
  z[1:10] <- -5
  scdata_list <- add_zscores(scdata_list, "nCount_RNA", setNames(list(z), sample1_id))

  cells_id <- mock_ids(scdata_list)

  config <- mock_spatial_config(cutoff = 3, direction = "lower")
  out <- filter_spatial_umi_outlier(scdata_list, config, sample1_id, cells_id)

  # original object is not subset; new_ids drops the 10 outliers
  expect_equal(ncol(out$data[[sample1_id]]), ncells)
  expect_equal(length(out$new_ids[[sample1_id]]), ncells - 10)
  expect_false(any(0:9 %in% out$new_ids[[sample1_id]]))
})

test_that("filter_spatial_umi_outlier can be disabled", {
  scdata_list <- add_tissue_coords(mock_scdata_list())
  sample1_id <- names(scdata_list)[1]
  ncells <- ncol(scdata_list[[sample1_id]])

  z <- rep(-5, ncells) # all would be outliers if enabled
  scdata_list <- add_zscores(scdata_list, "nCount_RNA", setNames(list(z), sample1_id))

  cells_id <- mock_ids(scdata_list)
  config <- mock_spatial_config(cutoff = 3, direction = "lower")
  config$enabled <- FALSE

  out <- filter_spatial_umi_outlier(scdata_list, config, sample1_id, cells_id)
  expect_equal(length(out$new_ids[[sample1_id]]), ncells)
})

test_that("filter_spatial_umi_outlier auto uses the default cutoff of 3", {
  scdata_list <- add_tissue_coords(mock_scdata_list())
  sample1_id <- names(scdata_list)[1]
  ncells <- ncol(scdata_list[[sample1_id]])

  z <- rep(0, ncells)
  z[1:5] <- -4 # beyond default cutoff
  scdata_list <- add_zscores(scdata_list, "nCount_RNA", setNames(list(z), sample1_id))

  cells_id <- mock_ids(scdata_list)
  # manual cutoff is huge, but auto should override it to 3
  config <- mock_spatial_config(cutoff = 100, direction = "lower", auto = TRUE)

  out <- filter_spatial_umi_outlier(scdata_list, config, sample1_id, cells_id)
  expect_equal(out$config$filterSettings$methodSettings$localOutlier$cutoff, 3)
  expect_equal(length(out$new_ids[[sample1_id]]), ncells - 5)
})

test_that("filter_spatial_umi_outlier emits the expected plotData", {
  scdata_list <- add_tissue_coords(mock_scdata_list())
  sample1_id <- names(scdata_list)[1]
  ncells <- ncol(scdata_list[[sample1_id]])

  z <- rnorm(ncells)
  scdata_list <- add_zscores(scdata_list, "nCount_RNA", setNames(list(z), sample1_id))

  cells_id <- mock_ids(scdata_list)
  config <- mock_spatial_config()
  out <- filter_spatial_umi_outlier(scdata_list, config, sample1_id, cells_id)

  task_name <- "spatialUmiOutlier"
  plot0 <- out$plotData[[generate_gui_uuid(sample1_id, task_name, 0)]]
  plot1 <- out$plotData[[generate_gui_uuid(sample1_id, task_name, 1)]]
  stats <- out$plotData[[generate_gui_uuid(sample1_id, task_name, 2)]]

  # plot0: one full per-cell record with cellId, value, zscore and coords (not downsampled)
  expect_equal(length(plot0), ncells)
  expect_setequal(names(plot0[[1]]), c("cellId", "value", "raw", "zscore", "x", "y"))

  # plot1: z-score values for the histogram
  expect_true(length(plot1) > 0)
  expect_setequal(names(plot1[[1]]), c("raw", "zscore"))

  # plot2: before/after filter statistics
  expect_setequal(names(stats), c("before", "after"))
})

test_that("manual (non-auto) cutoff is persisted to config and used for filtering", {
  scdata_list <- add_tissue_coords(mock_scdata_list())
  sample1_id <- names(scdata_list)[1]
  ncells <- ncol(scdata_list[[sample1_id]])

  z <- rep(0, ncells)
  z[1:6] <- -2.5 # outliers only when cutoff is below 2.5
  scdata_list <- add_zscores(scdata_list, "nCount_RNA", setNames(list(z), sample1_id))

  cells_id <- mock_ids(scdata_list)
  # auto = FALSE -> user-provided cutoff of 2 must be honoured (not reset to 3)
  config <- mock_spatial_config(cutoff = 2, direction = "lower", auto = FALSE)
  out <- filter_spatial_umi_outlier(scdata_list, config, sample1_id, cells_id)

  expect_equal(out$config$filterSettings$methodSettings$localOutlier$cutoff, 2)
  # z = -2.5 is beyond -2, so the 6 outliers are removed
  expect_equal(length(out$new_ids[[sample1_id]]), ncells - 6)
})

test_that("lower keeps z > -cutoff, so a cell exactly at -cutoff is removed", {
  scdata_list <- add_tissue_coords(mock_scdata_list())
  sample1_id <- names(scdata_list)[1]
  ncells <- ncol(scdata_list[[sample1_id]])
  cells_id <- mock_ids(scdata_list)

  z <- rep(0, ncells)
  z[1] <- -3.1 # beyond -cutoff -> removed
  z[2] <- -3 # exactly -cutoff -> removed
  z[3] <- -2.9 # within -> kept

  scdata_list <- add_zscores(scdata_list, "nCount_RNA", setNames(list(z), sample1_id))
  out <- filter_spatial_umi_outlier(
    scdata_list, mock_spatial_config(3, "lower"), sample1_id, cells_id
  )

  expect_false(any(c(0, 1) %in% out$new_ids[[sample1_id]]))
  expect_true(2 %in% out$new_ids[[sample1_id]])
  expect_equal(length(out$new_ids[[sample1_id]]), ncells - 2)
})

test_that("missing z-score column is a no-op with a warning message", {
  scdata_list <- add_tissue_coords(mock_scdata_list())
  sample1_id <- names(scdata_list)[1]
  ncells <- ncol(scdata_list[[sample1_id]])
  cells_id <- mock_ids(scdata_list)
  config <- mock_spatial_config()

  # no nCount_RNA_z column added
  expect_message(
    out <- filter_spatial_umi_outlier(scdata_list, config, sample1_id, cells_id),
    "No spatial z-scores"
  )
  expect_equal(length(out$plotData), 0)
  expect_identical(out$config, config)
})

test_that("plot 0 is the full, non-downsampled per-cell array with x/y coords", {
  sample1_id <- "spatial1"
  scdata_list <- mock_spatial_scdata_list(ncells = 64, sample_id = sample1_id)
  ncells <- ncol(scdata_list[[sample1_id]])

  z <- rnorm(ncells)
  scdata_list <- add_zscores(scdata_list, "nCount_RNA", setNames(list(z), sample1_id))
  # add the *_log column the filter colours by for a log-scaled metric
  scdata_list[[sample1_id]]$nCount_RNA_log <- log1p(scdata_list[[sample1_id]]$nCount_RNA)

  cells_id <- mock_ids(scdata_list)
  # small downsample target proves plot 0 is NOT downsampled while plot 1 is
  out <- filter_spatial_umi_outlier(
    scdata_list, mock_spatial_config(), sample1_id, cells_id,
    num_cells_to_downsample = 10
  )

  task_name <- "spatialUmiOutlier"
  plot0 <- out$plotData[[generate_gui_uuid(sample1_id, task_name, 0)]]
  plot1 <- out$plotData[[generate_gui_uuid(sample1_id, task_name, 1)]]

  # plot0 has one record per cell (not downsampled) with coords
  expect_equal(length(plot0), ncells)
  expect_setequal(names(plot0[[1]]), c("cellId", "value", "raw", "zscore", "x", "y"))
  # plot1 is downsampled to the requested number
  expect_equal(length(plot1), 10)
})

# ── absolute thresholds ──────────────────────────────────────────────────────

# sample whose nCount_RNA is set to known values; no z-scores are added, as the
# absolute method must not need them
mock_counts_sample <- function(counts) {
  scdata_list <- mock_spatial_scdata_list(ncells = length(counts))
  scdata_list$spatial1$nCount_RNA <- counts
  scdata_list
}

test_that("absolute min keeps cells at or above min, inclusive, without z-scores", {
  scdata_list <- mock_counts_sample(c(10, 99, 100, 101, rep(500, 60)))
  cells_id <- mock_ids(scdata_list)

  out <- filter_spatial_umi_outlier(
    scdata_list, mock_spatial_absolute_config(min = 100), "spatial1", cells_id
  )

  expect_equal(out$new_ids$spatial1, cells_id$spatial1[-(1:2)])
  expect_equal(out$config$filterSettings$methodSettings$absoluteThreshold$min, 100)
  expect_null(out$config$filterSettings$methodSettings$absoluteThreshold$max)
})

test_that("absolute max keeps cells at or below max, inclusive", {
  scdata_list <- mock_counts_sample(c(rep(500, 60), 999, 1000, 1001, 5000))
  cells_id <- mock_ids(scdata_list)

  out <- filter_spatial_umi_outlier(
    scdata_list, mock_spatial_absolute_config(max = 1000), "spatial1", cells_id
  )

  expect_equal(out$new_ids$spatial1, cells_id$spatial1[1:62])
})

test_that("absolute min and max together keep the cells in between", {
  scdata_list <- mock_counts_sample(c(5, rep(500, 62), 5000))
  cells_id <- mock_ids(scdata_list)

  out <- filter_spatial_umi_outlier(
    scdata_list, mock_spatial_absolute_config(min = 100, max = 1000), "spatial1", cells_id
  )

  expect_equal(out$new_ids$spatial1, cells_id$spatial1[2:63])
})

test_that("absolute with no bounds keeps every cell", {
  scdata_list <- mock_counts_sample(c(1, rep(500, 62), 1e6))
  cells_id <- mock_ids(scdata_list)

  out <- filter_spatial_umi_outlier(
    scdata_list, mock_spatial_absolute_config(), "spatial1", cells_id
  )

  expect_equal(out$new_ids$spatial1, cells_id$spatial1)
})

test_that("absolute auto uses the sample's 3-MAD lower threshold and no max", {
  set.seed(1)
  counts <- c(round(exp(rnorm(60, log(500), 0.2))), 5, 6, 1e5, 2e5)
  scdata_list <- mock_counts_sample(counts)
  cells_id <- mock_ids(scdata_list)

  # manual values are ignored in auto mode
  config <- mock_spatial_absolute_config(min = 1, max = 2, auto = TRUE)
  out <- filter_spatial_umi_outlier(scdata_list, config, "spatial1", cells_id)

  expected_min <- attr(
    scuttle::isOutlier(counts, log = TRUE, type = "lower", nmads = 3), "thresholds"
  )[["lower"]]
  settings <- out$config$filterSettings$methodSettings$absoluteThreshold
  expect_equal(settings$min, unname(expected_min))
  expect_null(settings$max)
  expect_true("max" %in% names(settings))

  # the two low cells are removed, the high ones are kept (no max)
  expect_equal(length(out$new_ids$spatial1), length(counts) - 2)
  expect_false(any(cells_id$spatial1[61:62] %in% out$new_ids$spatial1))
})

test_that("switching method keeps the other method's settings", {
  scdata_list <- mock_counts_sample(rep(500, 64))
  config <- mock_spatial_absolute_config(min = 100)
  config$filterSettings$methodSettings$localOutlier$cutoff <- 2.5

  out <- filter_spatial_umi_outlier(scdata_list, config, "spatial1", mock_ids(scdata_list))

  expect_equal(out$config$filterSettings$method, "absoluteThreshold")
  expect_equal(out$config$filterSettings$methodSettings$localOutlier$cutoff, 2.5)
})

test_that("absolute plotData carries raw values, and no z-scores when there are none", {
  scdata_list <- mock_counts_sample(seq(100, by = 10, length.out = 64))
  out <- filter_spatial_umi_outlier(
    scdata_list, mock_spatial_absolute_config(min = 200), "spatial1", mock_ids(scdata_list)
  )

  plot0 <- out$plotData[[generate_gui_uuid("spatial1", "spatialUmiOutlier", 0)]]
  plot1 <- out$plotData[[generate_gui_uuid("spatial1", "spatialUmiOutlier", 1)]]

  expect_setequal(names(plot0[[1]]), c("cellId", "value", "raw", "x", "y"))
  expect_equal(vapply(plot0, function(cell) cell$raw, numeric(1)), seq(100, by = 10, length.out = 64))
  expect_equal(names(plot1[[1]]), "raw")
})

test_that("disabled absolute filter keeps every cell but still resolves auto defaults", {
  set.seed(1)
  scdata_list <- mock_counts_sample(c(1, round(exp(rnorm(63, log(500), 0.2)))))
  config <- mock_spatial_absolute_config(auto = TRUE)
  config$enabled <- FALSE

  out <- filter_spatial_umi_outlier(scdata_list, config, "spatial1", mock_ids(scdata_list))

  expect_equal(length(out$new_ids$spatial1), 64)
  expect_false(is.null(out$config$filterSettings$methodSettings$absoluteThreshold$min))
})

test_that("a NULL config (step not configured, e.g. single-cell data) is a no-op", {
  scdata_list <- mock_scdata_list()
  cells_id <- mock_ids(scdata_list)
  sample1_id <- names(scdata_list)[1]

  out <- filter_spatial_umi_outlier(scdata_list, NULL, sample1_id, cells_id)

  expect_identical(out$new_ids, cells_id)
  expect_null(out$config)
  expect_equal(length(out$plotData), 0)
})

test_that("an unknown method is an error", {
  scdata_list <- mock_counts_sample(rep(500, 64))
  config <- mock_spatial_config()
  config$filterSettings$method <- "nope"

  expect_error(
    filter_spatial_umi_outlier(scdata_list, config, "spatial1", mock_ids(scdata_list)),
    "Unknown spatial QC method"
  )
})

test_that("generate_default_values_spatial_absolute sets no bound when the MAD is zero", {
  expect_null(generate_default_values_spatial_absolute(rep(100, 50), "lower")$min)
  expect_null(generate_default_values_spatial_absolute(rep(0, 50), "upper")$max)
})
