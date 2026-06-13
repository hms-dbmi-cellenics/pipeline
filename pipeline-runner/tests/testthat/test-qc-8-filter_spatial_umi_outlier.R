mock_spatial_config <- function(cutoff = 3, direction = "lower", auto = FALSE) {
  list(
    auto = auto,
    enabled = TRUE,
    filterSettings = list(cutoff = cutoff, direction = direction)
  )
}

# adds the <metric>_z column that add_spatial_local_outliers would have computed
add_zscores <- function(scdata_list, metric, zscores_by_sample) {
  for (sample_id in names(zscores_by_sample)) {
    scdata_list[[sample_id]]@meta.data[[paste0(metric, "_z")]] <-
      zscores_by_sample[[sample_id]]
  }
  scdata_list
}

test_that("filter_spatial_umi_outlier removes cells below -cutoff (lower direction)", {
  scdata_list <- mock_scdata_list()
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
  scdata_list <- mock_scdata_list()
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
  scdata_list <- mock_scdata_list()
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
  scdata_list <- mock_scdata_list()
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

  # plot0: one full per-cell record with cellId, value and zscore (not downsampled)
  expect_equal(length(plot0), ncells)
  expect_setequal(names(plot0[[1]]), c("cellId", "value", "zscore"))

  # plot1: z-score values for the histogram
  expect_true(length(plot1) > 0)
  expect_equal(names(plot1[[1]]), "zscore")

  # plot2: before/after filter statistics
  expect_setequal(names(stats), c("before", "after"))
})

test_that("filter_spatial_mito_outlier removes cells above the cutoff (upper direction)", {
  scdata_list <- mock_scdata_list()
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
