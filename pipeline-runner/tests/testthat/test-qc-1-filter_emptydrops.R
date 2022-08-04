source("qc_mock.R")

mock_config <- function() {
  config <- list(
    auto = FALSE,
    enabled = TRUE,
    filterSettings = list(
      FDR = 0.01
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

  # add empty drops stuff
  scdata$emptyDrops_FDR <- NA
  scdata$emptyDrops_FDR[1:70] <- 0.009

  # add samples
  scdata$samples <- rep(c("123abc", "123def"), each = 40)
  scdata <- Seurat::RenameCells(scdata, paste(scdata$samples, colnames(scdata), sep = ""))
  return(scdata)
}

test_that("Initial cells_id are correct", {
  scdata <- mock_scdata()
  cells_id <- generate_first_step_ids(scdata)
  expect_equal(unique(cells_id$`123abc`), 0:39)
  expect_equal(unique(cells_id$`123def`), 40:79)
})

test_that("filter_emptydrops removes NA with threshold < 1", {
  scdata <- mock_scdata()
  config <- mock_config()
  cells_id <- mock_ids()

  out <- filter_emptydrops(scdata, config, "123def", cells_id)
  expect_equal(ncol(out$data), 80)
  expect_equal(length(out$new_ids$`123abc`), 40)
  expect_equal(length(out$new_ids$`123def`), 30)
})

test_that("filter_emptydrops is sample aware", {
  scdata <- mock_scdata()
  config <- mock_config()
  cells_id <- mock_ids()

  # NA (empty) drops are in 123def only
  out <- filter_emptydrops(scdata, config, "123abc", cells_id)
  expect_equal(ncol(out$data), 80)
  expect_equal(length(out$new_ids$`123abc`), 40)
  expect_equal(length(out$new_ids$`123def`), 40)
})


test_that("if FDR=1 filter_emptydrops keeps everything", {
  scdata <- mock_scdata()
  config <- mock_config()
  config$filterSettings$FDR <- 1
  cells_id <- mock_ids()

  out <- filter_emptydrops(scdata, config, "123def", cells_id)
  expect_equal(ncol(out$data), 80)
  expect_equal(length(out$new_ids$`123abc`), 40)
  expect_equal(length(out$new_ids$`123def`), 40)
})

test_that("filter_emptydrops can be disabled", {
  scdata <- mock_scdata()
  config <- mock_config()
  config$enabled <- FALSE
  cells_id <- mock_ids()

  out <- filter_emptydrops(scdata, config, "123def", cells_id)
  expect_equal(ncol(out$data), 80)
  expect_equal(length(out$new_ids$`123abc`), 40)
  expect_equal(length(out$new_ids$`123def`), 40)
})


test_that("filter_emptydrops handles missing emptyDrops_FDR", {
  scdata <- mock_scdata()
  config <- mock_config()
  scdata$emptyDrops_FDR <- NULL
  cells_id <- mock_ids()

  out <- filter_emptydrops(scdata, config, "123def", cells_id)
  expect_equal(ncol(out$data), 80)
  expect_equal(length(out$new_ids$`123abc`), 40)
  expect_equal(length(out$new_ids$`123def`), 40)
})
