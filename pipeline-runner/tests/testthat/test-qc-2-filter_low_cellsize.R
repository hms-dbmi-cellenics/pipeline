mock_config <- function(mcs = 100, auto_settings = FALSE, enabled_on = TRUE) {
  config <- list(
    auto = auto_settings,
    enabled = enabled_on,
    filterSettings = list(
      minCellSize = mcs
    )
  )
  return(config)
}

test_that("filter_low_cellsize removes cells and works with bpcells", {

  for (use_bpcells in c(FALSE, TRUE)) {
    scdata_list <- mock_scdata_list()
    sample1_id <- names(scdata_list)[1]
    sample2_id <- names(scdata_list)[2]

    config <- mock_config(mcs = 10000)
    cells_id <- mock_ids(scdata_list)

    out <- filter_low_cellsize(scdata_list, config, sample1_id, cells_id)

    expect_equal(ncol(out$data[[sample1_id]]), ncol(scdata_list[[sample1_id]]))
    expect_lt(length(out$new_ids[[sample1_id]]), length(cells_id[[sample1_id]]))

    expect_equal(
      length(out$new_ids[[sample2_id]]),
      length(cells_id[[sample2_id]])
    )
  }
})

test_that("filter_low_cellsize filters only appropriate cells", {
  mcs <- 100
  scdata_list <- mock_scdata_list()
  sample1_id <- names(scdata_list)[1]
  sample2_id <- names(scdata_list)[2]

  config <- mock_config(mcs)
  cells_id <- mock_ids(scdata_list)

  out <- filter_low_cellsize(scdata_list, config, sample1_id, cells_id)


  data_sample_1 <- out$data[[sample1_id]]

  barcode_names_this_sample <-
    colnames(data_sample_1)[data_sample_1$samples == sample1_id]
  sample_subset <- subset(data_sample_1, cells = barcode_names_this_sample)
  sample_subset <- subset_ids(sample_subset, out$new_ids[[sample1_id]])

  expect_false(all(sample_subset$nCount_RNA <= mcs))
})

test_that("filter_low_cellsize is sample aware", {
  mcs <- 200
  scdata_list <- mock_scdata_list()
  sample1_id <- names(scdata_list)[1]
  sample2_id <- names(scdata_list)[2]

  config <- mock_config(mcs)
  cells_id <- mock_ids(scdata_list)

  out <- filter_low_cellsize(scdata_list, config, sample1_id, cells_id)
  data <- out$data
  barcode_names_this_sample <- colnames(data[[sample1_id]])
  expect_equal(length(barcode_names_this_sample), 40)

  barcode_names_this_sample <- colnames(data[[sample2_id]])
  expect_equal(length(barcode_names_this_sample), 40)

  sample1_new_ids <- out$new_ids[[sample1_id]]
  sample2_new_ids <- out$new_ids[[sample2_id]]

  expect_lt(length(sample1_new_ids), 40)
  expect_equal(length(sample2_new_ids), 40)
})

test_that("filter_low_cellsize works on empty data/wrong sample", {
  scdata_list <- mock_scdata_list()
  sample1_id <- names(scdata_list)[1]
  sample2_id <- names(scdata_list)[2]

  config <- mock_config(1000000)
  cells_id <- mock_ids(scdata_list)

  out <- filter_low_cellsize(scdata_list, config, sample1_id, cells_id)

  out <- filter_low_cellsize(out$data, config, sample2_id, out$new_ids)
  new_ids <- out$new_ids

  expect_equal(0, length(new_ids[[sample1_id]]))
  expect_equal(0, length(new_ids[[sample2_id]]))

  out <- filter_low_cellsize(out$data, config, sample1_id, new_ids)
  expect_equal(0, length(new_ids[[sample1_id]]))
  expect_equal(0, length(new_ids[[sample2_id]]))

  # Seurat object wasn't changed
  expect_equal(ncol(out$data), 40)
})

test_that("filter_low_cellsize works with auto", {
  scdata_list <- mock_scdata_list()
  sample1_id <- names(scdata_list)[1]
  sample2_id <- names(scdata_list)[2]
  config <- mock_config(auto_settings = TRUE)
  cells_id <- mock_ids(scdata_list)

  out <- filter_low_cellsize(scdata_list, config, sample2_id, cells_id)

  expect_false(
    config$filterSettings$minCellSize == out$config$filterSettings$minCellSize
  )

  data <- out$data[[sample2_id]]
  barcode_names_this_sample <- colnames(data)[data$samples == sample2_id]
  sample_subset <- subset(data, cells = barcode_names_this_sample)
  sample_subset <- subset_ids(data, out$new_ids[[sample2_id]])

  expect_equal(ncol(data), 40)
  expect_true(
    all(sample_subset$nCount_RNA >= out$config$filterSettings$minCellSize)
  )
})


test_that("filter_low_cellsize can be disabled", {
  scdata_list <- mock_scdata_list()
  sample1_id <- names(scdata_list)[1]
  sample2_id <- names(scdata_list)[2]
  config <- mock_config(mcs = 10000, enabled = FALSE)
  cells_id <- mock_ids(scdata_list)

  out <- filter_low_cellsize(scdata_list, config, sample2_id, cells_id)

  data <- out$data
  expect_equal(sum(purrr::map_int(data, ncol)), 80)
  expect_equal(length(out$new_ids[[sample1_id]]), 40)
  expect_equal(length(out$new_ids[[sample2_id]]), 40)
})
