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

  expect_equal(out$config$filterSettings$methodSettings$localOutlier$cutoff, 3)
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

test_that("absolute mito keeps cells at or below max percent and ignores any min", {
  scdata_list <- mock_spatial_scdata_list(ncells = 64)
  scdata_list$spatial1$percent.mt <- c(rep(5, 60), 24.9, 25, 25.1, 80)
  cells_id <- mock_ids(scdata_list)

  out <- filter_spatial_mito_outlier(
    scdata_list, mock_spatial_absolute_config(min = 10, max = 25), "spatial1", cells_id
  )

  expect_equal(out$new_ids$spatial1, cells_id$spatial1[1:62])
})

test_that("absolute mito auto uses the sample's 3-MAD upper threshold on percent", {
  set.seed(1)
  mito <- c(runif(60, 2, 6), 40, 50, 60, 70)
  scdata_list <- mock_spatial_scdata_list(ncells = 64)
  scdata_list$spatial1$percent.mt <- mito
  cells_id <- mock_ids(scdata_list)

  out <- filter_spatial_mito_outlier(
    scdata_list, mock_spatial_absolute_config(auto = TRUE), "spatial1", cells_id
  )

  expected_max <- attr(scuttle::isOutlier(mito, type = "higher", nmads = 3), "thresholds")[["higher"]]
  settings <- out$config$filterSettings$methodSettings$absoluteThreshold
  expect_equal(settings$max, unname(expected_max))
  expect_equal(length(out$new_ids$spatial1), 60)
})
