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

mock_scdata <- function() {
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )

  sample_1_id <- "123abc"
  sample_2_id <- "123def"

  scdata <- Seurat::CreateSeuratObject(counts = pbmc_raw)
  scdata$cells_id <- 0:(ncol(scdata) - 1)

  # add samples
  scdata$samples <- rep(c(sample_1_id, sample_2_id), each = 40)
  scdata <- Seurat::RenameCells(scdata, paste(scdata$samples, colnames(scdata), sep = ""))
  scdata$emptyDrops_FDR <- NA

  scdata_sample1 <- subset(scdata, samples == sample_1_id)
  scdata_sample2 <- subset(scdata, samples == sample_2_id)

  scdata_list <- list(scdata_sample1, scdata_sample2)
  names(scdata_list) <- c(sample_1_id, sample_2_id)

  return(list(scdata_list, sample_1_id, sample_2_id))
}


test_that("filter_low_cellsize removes cells", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  config <- mock_config(mcs = 10000)
  cells_id <- mock_ids()

  out <- filter_low_cellsize(scdata_list, config, sample_1_id, cells_id)

  expect_equal(ncol(out$data[[sample_1_id]]), ncol(scdata_list[[sample_1_id]]))
  expect_lt(length(out$new_ids[[sample_1_id]]), length(cells_id[[sample_1_id]]))

  expect_equal(length(out$new_ids[[sample_2_id]]), length(cells_id[[sample_2_id]]))
})

test_that("filter_low_cellsize filters only appropiate cells", {
  mcs <- 100
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  config <- mock_config(mcs)
  cells_id <- mock_ids()

  out <- filter_low_cellsize(scdata_list, config, sample_1_id, cells_id)


  data_sample_1 <- out$data[[sample_1_id]]

  barcode_names_this_sample <- colnames(data_sample_1)[data_sample_1$samples == sample_1_id]
  sample_subset <- subset(data_sample_1, cells = barcode_names_this_sample)
  sample_subset <- subset_ids(sample_subset, out$new_ids[[sample_1_id]])

  expect_false(all(sample_subset$nCount_RNA <= mcs))
})

test_that("filter_low_cellsize is sample aware", {
  mcs <- 200
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  config <- mock_config(mcs)
  cells_id <- mock_ids()

  out <- filter_low_cellsize(scdata_list, config, sample_1_id, cells_id)
  data <- out$data
  barcode_names_this_sample <- colnames(data[[sample_1_id]])
  expect_equal(length(barcode_names_this_sample), 40)

  barcode_names_this_sample <- colnames(data[[sample_2_id]])
  expect_equal(length(barcode_names_this_sample), 40)

  sample_1_new_ids <- out$new_ids[[sample_1_id]]
  sample_2_new_ids <- out$new_ids[[sample_2_id]]

  expect_lt(length(sample_1_new_ids), 40)
  expect_equal(length(sample_2_new_ids), 40)
})

test_that("filter_low_cellsize works on empty data/wrong sample", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  config <- mock_config(1000000)
  cells_id <- mock_ids()

  out <- filter_low_cellsize(scdata_list, config, sample_1_id, cells_id)

  out <- filter_low_cellsize(out$data, config, sample_2_id, out$new_ids)
  new_ids <- out$new_ids

  expect_equal(0, length(new_ids[[sample_1_id]]))
  expect_equal(0, length(new_ids[[sample_2_id]]))

  out <- filter_low_cellsize(out$data, config, sample_1_id, new_ids)
  expect_equal(0, length(new_ids[[sample_1_id]]))
  expect_equal(0, length(new_ids[[sample_2_id]]))

  # Seurat object wasn't changed
  expect_equal(ncol(out$data), 40)
})

test_that("filter_low_cellsize works with auto", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  config <- mock_config(auto_settings = TRUE)
  cells_id <- mock_ids()

  out <- filter_low_cellsize(scdata_list, config, sample_2_id, cells_id)

  expect_false(config$filterSettings$minCellSize == out$config$filterSettings$minCellSize)

  data <- out$data[[sample_2_id]]
  barcode_names_this_sample <- colnames(data)[data$samples == sample_2_id]
  sample_subset <- subset(data, cells = barcode_names_this_sample)
  sample_subset <- subset_ids(data, out$new_ids[[sample_2_id]])

  expect_equal(ncol(data), 40)
  expect_true(all(sample_subset$nCount_RNA >= out$config$filterSettings$minCellSize))
})


test_that("filter_low_cellsize can be disabled", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  config <- mock_config(mcs = 10000, enabled = FALSE)
  cells_id <- mock_ids()

  out <- filter_low_cellsize(scdata_list, config, sample_2_id, cells_id)

  data <- out$data
  expect_equal(sum(purrr::map_int(data, ncol)), 80)
  expect_equal(length(out$new_ids[[sample_1_id]]), 40)
  expect_equal(length(out$new_ids[[sample_2_id]]), 40)
})

