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
    as.is = TRUE)

  scdata <- Seurat::CreateSeuratObject(counts = pbmc_raw)

  # add empty drops stuff
  scdata$emptyDrops_FDR <- NA
  scdata$emptyDrops_FDR[1:70] <- 0.009

  # add samples
  scdata$samples <- rep(c('123abc', '123def'), each = 40)
  scdata <- Seurat::RenameCells(scdata, paste(scdata$samples, colnames(scdata), sep = ''))
  return(scdata)
}


test_that("filter_emptydrops removes NA with threshold < 1", {
  scdata <- mock_scdata()
  config <- mock_config()

  out <- filter_emptydrops(scdata, config, '123def')
  expect_equal(ncol(out$data), 70)
})

test_that("filter_emptydrops is sample aware", {
  scdata <- mock_scdata()
  config <- mock_config()

  # NA (empty) drops are in 123def only
  out <- filter_emptydrops(scdata, config, '123abc')
  expect_equal(ncol(out$data), 80)
})


test_that("if FDR=1 filter_emptydrops keeps everything", {
  scdata <- mock_scdata()
  config <- mock_config()
  config$filterSettings$FDR <- 1

  out <- filter_emptydrops(scdata, config, '123def')
  expect_equal(ncol(out$data), 80)
})

test_that("filter_emptydrops can be disabled", {
  scdata <- mock_scdata()
  config <- mock_config()
  config$enabled <- FALSE

  out <- filter_emptydrops(scdata, config, '123def')
  expect_equal(ncol(out$data), 80)
})


test_that("filter_emptydrops handles missing emptyDrops_FDR", {
  scdata <- mock_scdata()
  config <- mock_config()
  scdata$emptyDrops_FDR <- NULL

  out <- filter_emptydrops(scdata, config, '123def')
  expect_equal(ncol(out$data), 80)
})
