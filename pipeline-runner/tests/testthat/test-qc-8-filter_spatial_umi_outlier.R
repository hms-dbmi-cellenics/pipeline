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
  expect_equal(out$config$filterSettings$cutoff, 3)
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
  expect_setequal(names(plot0[[1]]), c("cellId", "value", "zscore", "x", "y"))

  # plot1: z-score values for the histogram
  expect_true(length(plot1) > 0)
  expect_equal(names(plot1[[1]]), "zscore")

  # plot2: before/after filter statistics
  expect_setequal(names(stats), c("before", "after"))
})

test_that("filter_spatial_mito_outlier removes cells above the cutoff (upper direction)", {
  scdata_list <- add_tissue_coords(mock_scdata_list())
  sample1_id <- names(scdata_list)[1]
  ncells <- ncol(scdata_list[[sample1_id]])

  z <- rep(0, ncells)
  z[1:8] <- 5 # high outliers
  scdata_list <- add_zscores(scdata_list, "percent.mt", setNames(list(z), sample1_id))

  cells_id <- mock_ids(scdata_list)
  config <- mock_spatial_config(cutoff = 3, direction = "upper")
  out <- filter_spatial_mito_outlier(scdata_list, config, sample1_id, cells_id)

  expect_equal(length(out$new_ids[[sample1_id]]), ncells - 8)
  expect_false(any(0:7 %in% out$new_ids[[sample1_id]]))
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

  expect_equal(out$config$filterSettings$cutoff, 2)
  expect_equal(out$config$filterSettings$direction, "lower")
  # z = -2.5 is beyond -2, so the 6 outliers are removed
  expect_equal(length(out$new_ids[[sample1_id]]), ncells - 6)
})

test_that("upper keeps z < cutoff and lower keeps z > -cutoff at the boundary", {
  scdata_list <- add_tissue_coords(mock_scdata_list())
  sample1_id <- names(scdata_list)[1]
  ncells <- ncol(scdata_list[[sample1_id]])
  cells_id <- mock_ids(scdata_list)

  # upper keeps z < cutoff (strict); lower keeps z > -cutoff (strict).
  # Cells exactly at the boundary are removed.
  z <- rep(0, ncells)
  z[1] <- 3.1 # > cutoff  -> removed (upper)
  z[2] <- 2.9 # < cutoff  -> kept (upper)
  z[3] <- -3.1 # < -cutoff -> removed (lower)
  z[4] <- -2.9 # > -cutoff -> kept (lower)

  sl_up <- add_zscores(scdata_list, "percent.mt", setNames(list(z), sample1_id))
  out_up <- filter_spatial_mito_outlier(
    sl_up, mock_spatial_config(3, "upper"), sample1_id, cells_id
  )
  # cells_id 0 (z=3.1) removed, cells_id 1 (z=2.9) kept
  expect_false(0 %in% out_up$new_ids[[sample1_id]])
  expect_true(1 %in% out_up$new_ids[[sample1_id]])
  expect_equal(length(out_up$new_ids[[sample1_id]]), ncells - 1)

  sl_lo <- add_zscores(scdata_list, "nCount_RNA", setNames(list(z), sample1_id))
  out_lo <- filter_spatial_umi_outlier(
    sl_lo, mock_spatial_config(3, "lower"), sample1_id, cells_id
  )
  # cells_id 2 (z=-3.1) removed, cells_id 3 (z=-2.9) kept
  expect_false(2 %in% out_lo$new_ids[[sample1_id]])
  expect_true(3 %in% out_lo$new_ids[[sample1_id]])
  expect_equal(length(out_lo$new_ids[[sample1_id]]), ncells - 1)
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
  expect_setequal(names(plot0[[1]]), c("cellId", "value", "zscore", "x", "y"))
  # plot1 is downsampled to the requested number
  expect_equal(length(plot1), 10)
})
