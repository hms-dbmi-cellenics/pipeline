mock_config <- function(max_fraction = 0.1) {
  config <- list(
    auto = TRUE,
    enabled = TRUE,
    filterSettings = list(
      method = "absoluteThreshold",
      methodSettings = list(
        absoluteThreshold = list(maxFraction = max_fraction)
      )
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

  # add mitochondrial percent
  scdata@meta.data$percent.mt <- 5
  scdata@meta.data$percent.mt[1:10] <- 11

  # add samples
  scdata$samples <- rep(c(sample_1_id, sample_2_id), each = 40)
  scdata <- Seurat::RenameCells(scdata, paste(scdata$samples, colnames(scdata), sep = ""))

  scdata_sample1 <- subset(scdata, samples == sample_1_id)
  scdata_sample2 <- subset(scdata, samples == sample_2_id)

  scdata_list <- list(scdata_sample1, scdata_sample2)
  names(scdata_list) <- c(sample_1_id, sample_2_id)

  return(list(scdata_list, sample_1_id, sample_2_id))
}

get_threshold <- function(config) config$filterSettings$methodSettings$absoluteThreshold$maxFraction


test_that("filter_high_mito filters based on threshold", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  cells_id <- mock_ids()

  # should filter first 10 cells
  config <- mock_config(0.1)
  config$auto <- FALSE
  out <- filter_high_mito(scdata_list, config, sample_1_id, cells_id)
  expect_equal(ncol(out$data[[sample_1_id]]), 40)
  expect_equal(length(out$new_ids[[sample_1_id]]), 30)

  # should filter all cells in first sample
  config <- mock_config(0.01)
  config$auto <- FALSE

  out <- filter_high_mito(scdata_list, config, sample_1_id, cells_id)
  expect_equal(ncol(out$data[[sample_1_id]]), 40)
  expect_equal(length(out$new_ids[[sample_1_id]]), 0)
})

test_that("filter_high_mito can be disabled", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  cells_id <- mock_ids()
  nstart <- ncol(scdata_list[[sample_1_id]])

  # should filter first 10 cells
  config <- mock_config(0.1)
  config$enabled <- FALSE
  out <- filter_high_mito(scdata_list, config, sample_1_id, cells_id)
  expect_equal(ncol(out$data[[sample_1_id]]), nstart)
  expect_equal(out$new_ids[[sample_1_id]], 0:39)
  expect_equal(out$new_ids[[sample_2_id]], 40:79)
})

test_that("filter_high_mito is sample specific", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  cells_id <- mock_ids()

  # should filter first 10 cells of 123abc only
  config <- mock_config(0.1)
  config$auto <- FALSE

  out1 <- filter_high_mito(scdata_list, config, sample_1_id, cells_id)
  out2 <- filter_high_mito(scdata_list, config, sample_2_id, cells_id)

  expect_equal(ncol(out1$data[[sample_1_id]]), 40)
  expect_equal(length(out1$new_ids[[sample_1_id]]), 30)
  expect_equal(length(out1$new_ids[[sample_2_id]]), 40)
  expect_equal(ncol(out2$data[[sample_2_id]]), 40)
  expect_equal(length(out2$new_ids[[sample_1_id]]), 40)
  expect_equal(length(out2$new_ids[[sample_2_id]]), 40)
})

test_that("filter_high_mito can be set to auto", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  cells_id <- mock_ids()

  # would filter all 40 cells in first sample unless auto
  config <- mock_config(0.01)
  out <- filter_high_mito(scdata_list, config, sample_1_id, cells_id)

  # check that threshold was updated
  expect_true(get_threshold(out$config) > get_threshold(config))

  # check that updated threshold was used
  nkeep <- length(unlist(out$new_ids[[sample_1_id]]))
  expected <- sum(scdata_list[[sample_1_id]]$percent.mt <= get_threshold(out$config) * 100)
  expect_equal(nkeep, expected)

  # didn't subset original data
  expect_equal(ncol(out$data[[sample_1_id]]), 40)
})


test_that("data without percent.mt outliers uses the max percentage as the threshold", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  set.seed(RANDOM_SEED)

  scdata_sample_1 <- scdata_list[[sample_1_id]]

  scdata_list$percent.mt <- rnorm(ncol(scdata_sample_1), 6)
  cells_id <- mock_ids()

  config <- mock_config()
  out <- filter_high_mito(scdata_list, config, sample_1_id, cells_id)

  # nothing filtered
  expect_length(unlist(out$new_ids[[sample_1_id]]), ncol(scdata_sample_1))

  # threshold equals max for sample
  expect_equal(get_threshold(out$config), max(scdata_sample_1$percent.mt / 100))
})

