mock_scdata <- function() {
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )

  scdata <- Seurat::CreateSeuratObject(counts = pbmc_raw)

  # add samples
  scdata$samples <- rep("123abc", 80)
  scdata <- Seurat::RenameCells(scdata, paste(scdata$samples, colnames(scdata), sep = ""))

  # add doublet scores
  scdata$doublet_scores <- rep(c(0.01, 0.9), each = 40)
  scdata$doublet_class <- rep(c("singlet", "doublet"), each = 40)

  # add mitochondrial percent
  scdata$percent.mt <- rnorm(ncol(scdata), mean = 6)
  return(scdata)
}


test_that("cellsize filter is disabled by default and classifier if pre-filtered", {
  scdata <- mock_scdata()

  qc_config <- construct_qc_config(list(scdata), any_filtered = TRUE)

  expect_false(qc_config$classifier$enabled)
  expect_true(qc_config$classifier$prefiltered)
  expect_false(qc_config$cellSizeDistribution$enabled)
})


test_that("cellsize filter is disabled by default and classifier if not pre-filtered", {
  scdata <- mock_scdata()

  qc_config <- construct_qc_config(list(scdata), any_filtered = FALSE)
  expect_false(qc_config$cellSizeDistribution$enabled)
  expect_true(qc_config$classifier$enabled)
})
