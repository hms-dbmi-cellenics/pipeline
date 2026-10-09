# mock_spatial_config() and add_zscores() live in helper-spatial.R

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

test_that("upper keeps z < cutoff, so a cell exactly at the cutoff is removed", {
  scdata_list <- add_tissue_coords(mock_scdata_list())
  sample1_id <- names(scdata_list)[1]
  ncells <- ncol(scdata_list[[sample1_id]])
  cells_id <- mock_ids(scdata_list)

  z <- rep(0, ncells)
  z[1] <- 3.1 # beyond cutoff -> removed
  z[2] <- 3 # exactly cutoff -> removed
  z[3] <- 2.9 # within -> kept

  scdata_list <- add_zscores(scdata_list, "percent.mt", setNames(list(z), sample1_id))
  out <- filter_spatial_mito_outlier(
    scdata_list, mock_spatial_config(3, "upper"), sample1_id, cells_id
  )

  expect_false(any(c(0, 1) %in% out$new_ids[[sample1_id]]))
  expect_true(2 %in% out$new_ids[[sample1_id]])
  expect_equal(length(out$new_ids[[sample1_id]]), ncells - 2)
})

test_that("filter_spatial_mito_outlier auto uses the default cutoff of 3", {
  scdata_list <- add_tissue_coords(mock_scdata_list())
  sample1_id <- names(scdata_list)[1]
  ncells <- ncol(scdata_list[[sample1_id]])

  z <- rep(0, ncells)
  z[1:4] <- 4
  scdata_list <- add_zscores(scdata_list, "percent.mt", setNames(list(z), sample1_id))

  cells_id <- mock_ids(scdata_list)
  config <- mock_spatial_config(cutoff = 100, direction = "upper", auto = TRUE)
  out <- filter_spatial_mito_outlier(scdata_list, config, sample1_id, cells_id)

  expect_equal(out$config$filterSettings$cutoff, 3)
  expect_equal(out$config$filterSettings$direction, "upper")
  expect_equal(length(out$new_ids[[sample1_id]]), ncells - 4)
})

test_that("filter_spatial_mito_outlier plots raw percent.mt, not a log value", {
  sample1_id <- "spatial1"
  scdata_list <- mock_spatial_scdata_list(ncells = 64, sample_id = sample1_id)
  ncells <- ncol(scdata_list[[sample1_id]])
  scdata_list <- add_zscores(scdata_list, "percent.mt", setNames(list(rnorm(ncells)), sample1_id))

  out <- filter_spatial_mito_outlier(
    scdata_list, mock_spatial_config(3, "upper"), sample1_id, mock_ids(scdata_list)
  )

  plot0 <- out$plotData[[generate_gui_uuid(sample1_id, "spatialMitoOutlier", 0)]]
  values <- vapply(plot0, function(cell) cell$value, numeric(1))
  expect_equal(values, unname(scdata_list[[sample1_id]]$percent.mt))
})
