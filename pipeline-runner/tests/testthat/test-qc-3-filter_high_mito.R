mock_config <- function(max_fraction = 0.1) {
  config <- list(
    auto = TRUE,
    enabled = TRUE,
    filterSettings = list(
      method = "absolute_threshold",
      methodSettings = list(
        absolute_threshold = list(maxFraction = max_fraction)
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

  scdata <- Seurat::CreateSeuratObject(counts = pbmc_raw)
  scdata$cells_id <- 0:(ncol(scdata) - 1)

  # add mitochondrial percent
  scdata@meta.data$percent.mt <- 5
  scdata@meta.data$percent.mt[1:10] <- 11

  # add samples
  scdata$samples <- rep(c("123abc", "123def"), each = 40)
  scdata <- Seurat::RenameCells(scdata, paste(scdata$samples, colnames(scdata), sep = ""))
  return(scdata)
}

get_threshold <- function(config) config$filterSettings$methodSettings$absolute_threshold$maxFraction


test_that("filter_high_mito filters based on threshold", {
  scdata <- mock_scdata()
  cells_id <- mock_ids()

  # should filter first 10 cells
  config <- mock_config(0.1)
  config$auto <- FALSE
  out <- filter_high_mito(scdata, config, "123abc", cells_id)
  expect_equal(ncol(out$data), 80)
  expect_equal(length(out$new_ids$`123abc`), 30)

  # should filter all cells in first sample
  config <- mock_config(0.01)
  config$auto <- FALSE

  out <- filter_high_mito(scdata, config, "123abc", cells_id)
  expect_equal(ncol(out$data), 80)
  expect_equal(length(out$new_ids$`123abc`), 0)
})

test_that("filter_high_mito can be disabled", {
  scdata <- mock_scdata()
  cells_id <- mock_ids()
  nstart <- ncol(scdata)

  # should filter first 10 cells
  config <- mock_config(0.1)
  config$enabled <- FALSE
  out <- filter_high_mito(scdata, config, "123abc", cells_id)
  expect_equal(ncol(out$data), nstart)
  expect_equal(out$new_ids$`123abc`, 0:39)
  expect_equal(out$new_ids$`123def`, 40:79)
})

test_that("filter_high_mito is sample specific", {
  scdata <- mock_scdata()
  cells_id <- mock_ids()
  nstart <- ncol(scdata)

  # should filter first 10 cells of 123abc only
  config <- mock_config(0.1)
  config$auto <- FALSE

  out1 <- filter_high_mito(scdata, config, "123abc", cells_id)
  out2 <- filter_high_mito(scdata, config, "123def", cells_id)

  expect_equal(ncol(out1$data), 80)
  expect_equal(length(out1$new_ids$`123abc`), 30)
  expect_equal(length(out1$new_ids$`123def`), 40)
  expect_equal(ncol(out2$data), 80)
  expect_equal(length(out2$new_ids$`123abc`), 40)
  expect_equal(length(out2$new_ids$`123def`), 40)
})

test_that("filter_high_mito can be set to auto", {
  scdata <- mock_scdata()
  cells_id <- mock_ids()
  nstart <- ncol(scdata)

  # would filter all 40 cells in first sample unless auto
  config <- mock_config(0.01)
  out <- filter_high_mito(scdata, config, "123abc", cells_id)

  # check that threshold was updated
  expect_true(get_threshold(out$config) > get_threshold(config))

  # check that updated threshold was used
  nkeep <- length(unlist(out$new_ids))
  expected <- sum(scdata$percent.mt <= get_threshold(out$config) * 100)
  expect_equal(nkeep, expected)

  # didn't subset original data
  expect_equal(ncol(out$data), 80)
})


test_that("data without percent.mt outliers uses the max percentage as the threshold", {
  scdata <- mock_scdata()
  set.seed(pipeline::gem2s$random.seed)
  scdata$percent.mt <- rnorm(ncol(scdata), 6)
  cells_id <- mock_ids()

  config <- mock_config()
  out <- filter_high_mito(scdata, config, "123abc", cells_id)

  # nothing filtered
  expect_length(unlist(out$new_ids), ncol(scdata))

  # threshold equals max for sample
  expect_equal(get_threshold(out$config), max(scdata$percent.mt[scdata$samples == "123abc"]) / 100)
})

