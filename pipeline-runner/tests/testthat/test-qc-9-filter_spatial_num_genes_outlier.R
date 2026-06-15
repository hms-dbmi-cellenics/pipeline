# Helpers mock_spatial_config() / add_zscores() live in helper-spatial.R

test_that("filter_spatial_num_genes_outlier removes low local outliers (lower direction)", {
  scdata_list <- mock_scdata_list()
  sample1_id <- names(scdata_list)[1]
  ncells <- ncol(scdata_list[[sample1_id]])

  z <- rep(0, ncells)
  z[1:7] <- -5 # strong low outliers in number of genes
  scdata_list <- add_zscores(scdata_list, "nFeature_RNA", setNames(list(z), sample1_id))

  cells_id <- mock_ids(scdata_list)
  config <- mock_spatial_config(cutoff = 3, direction = "lower")
  out <- filter_spatial_num_genes_outlier(scdata_list, config, sample1_id, cells_id)

  expect_equal(ncol(out$data[[sample1_id]]), ncells) # object not subset
  expect_equal(length(out$new_ids[[sample1_id]]), ncells - 7)
  expect_false(any(0:6 %in% out$new_ids[[sample1_id]]))
})

test_that("filter_spatial_num_genes_outlier auto uses default cutoff of 3", {
  scdata_list <- mock_scdata_list()
  sample1_id <- names(scdata_list)[1]
  ncells <- ncol(scdata_list[[sample1_id]])

  z <- rep(0, ncells)
  z[1:4] <- -4
  scdata_list <- add_zscores(scdata_list, "nFeature_RNA", setNames(list(z), sample1_id))

  cells_id <- mock_ids(scdata_list)
  # manual cutoff huge, but auto overrides to 3
  config <- mock_spatial_config(cutoff = 100, direction = "lower", auto = TRUE)
  out <- filter_spatial_num_genes_outlier(scdata_list, config, sample1_id, cells_id)

  expect_equal(out$config$filterSettings$cutoff, 3)
  expect_equal(out$config$filterSettings$direction, "lower")
  expect_equal(length(out$new_ids[[sample1_id]]), ncells - 4)
})

test_that("filter_spatial_num_genes_outlier can be disabled (no-op)", {
  scdata_list <- mock_scdata_list()
  sample1_id <- names(scdata_list)[1]
  ncells <- ncol(scdata_list[[sample1_id]])

  z <- rep(-5, ncells) # all would be outliers if enabled
  scdata_list <- add_zscores(scdata_list, "nFeature_RNA", setNames(list(z), sample1_id))

  cells_id <- mock_ids(scdata_list)
  config <- mock_spatial_config(cutoff = 3, direction = "lower")
  config$enabled <- FALSE

  out <- filter_spatial_num_genes_outlier(scdata_list, config, sample1_id, cells_id)
  expect_equal(length(out$new_ids[[sample1_id]]), ncells)
})

test_that("filter_spatial_num_genes_outlier emits the three plotData entries", {
  scdata_list <- mock_scdata_list()
  sample1_id <- names(scdata_list)[1]
  ncells <- ncol(scdata_list[[sample1_id]])

  z <- rnorm(ncells)
  scdata_list <- add_zscores(scdata_list, "nFeature_RNA", setNames(list(z), sample1_id))

  cells_id <- mock_ids(scdata_list)
  out <- filter_spatial_num_genes_outlier(
    scdata_list, mock_spatial_config(), sample1_id, cells_id
  )

  task_name <- "spatialNumGenesOutlier"
  plot0 <- out$plotData[[generate_gui_uuid(sample1_id, task_name, 0)]]
  plot1 <- out$plotData[[generate_gui_uuid(sample1_id, task_name, 1)]]
  stats <- out$plotData[[generate_gui_uuid(sample1_id, task_name, 2)]]

  expect_equal(length(plot0), ncells) # full per-cell, not downsampled
  expect_true(all(c("cellId", "value", "zscore") %in% names(plot0[[1]])))
  expect_equal(names(plot1[[1]]), "zscore")
  expect_setequal(names(stats), c("before", "after"))
})
