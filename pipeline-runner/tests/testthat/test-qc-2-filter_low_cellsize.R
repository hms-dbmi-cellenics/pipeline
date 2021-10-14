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

  scdata <- Seurat::CreateSeuratObject(counts = pbmc_raw)
  scdata$cells_id <- 0:(ncol(scdata) - 1)

  # add samples
  scdata$samples <- rep(c("123abc", "123def"), each = 40)
  scdata <- Seurat::RenameCells(scdata, paste(scdata$samples, colnames(scdata), sep = ""))

  scdata$emptyDrops_FDR <- NA
  return(scdata)
}


test_that("filter_low_cellsize removes cells", {
  scdata <- mock_scdata()
  config <- mock_config(mcs = 10000)
  cells_id <- mock_ids()

  out <- filter_low_cellsize(scdata, config, "123abc",cells_id)

  expect_equal(ncol(out$data), ncol(scdata))
  expect_lt(length(out$new_ids$`123abc`),length(cells_id$`123abc`))
  expect_equal(length(out$new_ids$`123def`),length(cells_id$`123def`))
})

test_that("filter_low_cellsize filters only appropiate cells", {
  mcs <- 100
  scdata <- mock_scdata()
  config <- mock_config(mcs)
  cells_id <- mock_ids()

  out <- filter_low_cellsize(scdata, config, "123abc",cells_id)

  data <- out$data
  barcode_names_this_sample <- colnames(data)[data$samples == "123abc"]
  sample_subset <- subset(data, cells = barcode_names_this_sample)
  sample_subset <- subset_ids(sample_subset,out$new_ids$`123abc`)

  expect_false(all(sample_subset$nCount_RNA <= mcs))
})

test_that("filter_low_cellsize is sample aware", {
  mcs <- 10000
  scdata <- mock_scdata()
  config <- mock_config(mcs)
  cells_id <- mock_ids()

  out <- filter_low_cellsize(scdata, config, "123abc",cells_id)
  data <- out$data
  barcode_names_this_sample <- colnames(data)[data$samples == "123abc"]
  expect_equal(length(barcode_names_this_sample), 40)

  barcode_names_this_sample <- colnames(data)[data$samples == "123def"]
  expect_equal(length(barcode_names_this_sample), 40)

  expect_lt(length(out$new_ids$`123abc`),40)
  expect_equal(length(out$new_ids$`123def`),40)

  out <- filter_low_cellsize(scdata, config, "123def",cells_id)
  data <- out$data

  barcode_names_this_sample <- colnames(data)[data$samples == "123def"]
  expect_equal(length(barcode_names_this_sample), 40)

  expect_equal(length(out$new_ids$`123abc`),40)
  expect_lt(length(out$new_ids$`123def`),40)
})

test_that("filter_low_cellsize works on empty data/wrong sample", {
  scdata <- mock_scdata()
  config <- mock_config(100000)
  cells_id <- mock_ids()

  out <- filter_low_cellsize(scdata, config, "123abc",cells_id)

  out <- filter_low_cellsize(out$data, config, "123def",out$new_ids)
  new_ids <- out$new_ids

  expect_equal(0,length(new_ids$`123abc`))
  expect_equal(0,length(new_ids$`123def`))

  out <- filter_low_cellsize(out$data, config, "123abc",new_ids)
  expect_equal(0,length(new_ids$`123abc`))
  expect_equal(0,length(new_ids$`123def`))

  expect_equal(ncol(out$data),80)
})

test_that("filter_low_cellsize works with auto", {
  scdata <- mock_scdata()
  config <- mock_config(auto_settings = TRUE)
  cells_id <- mock_ids()

  out <- filter_low_cellsize(scdata, config, "123def",cells_id)

  expect_false(config$filterSettings$minCellSize == out$config$filterSettings$minCellSize)

  data <- out$data
  barcode_names_this_sample <- colnames(data)[data$samples == "123def"]
  sample_subset <- subset(data, cells = barcode_names_this_sample)
  sample_subset <- subset_ids(data, out$new_ids$`123def`)

  expect_equal(ncol(out$data),80)
  expect_true(all(sample_subset$nCount_RNA >= out$config$filterSettings$minCellSize))
})


test_that("filter_low_cellsize can be disabled", {
  scdata <- mock_scdata()
  config <- mock_config(mcs = 10000, enabled = FALSE)
  cells_id <- mock_ids()

  out <- filter_low_cellsize(scdata, config, "123def",cells_id)

  expect_equal(ncol(out$data), 80)
  expect_equal(length(out$new_ids$`123abc`),40)
  expect_equal(length(out$new_ids$`123def`),40)
})
