mock_scdata <- function() {
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )

  scdata <- Seurat::CreateSeuratObject(counts = pbmc_raw)

  # add samples
  scdata$samples <- rep(c("123abc", "123def"), each = 40)
  return(scdata)
}


test_that("harmony integration works", {
  scdata <- mock_scdata()
  config <- list(
    dimensionalityReduction = list(numPCs = 2),
    dataIntegration = list(method = "harmony", methodSettings = list(harmony = list(numGenes = 10, normalisation = "logNormalize")))
  )

  scdata <- suppressWarnings(run_dataIntegration(scdata, config))
  expect_s4_class(scdata, "Seurat")
})

test_that("SeuratV4 integration doesnt error out with small dataset", {
  scdata <- mock_scdata()
  config <- list(
    dimensionalityReduction = list(numPCs = 2),
    dataIntegration = list(method = "seuratv4", methodSettings = list(seuratv4 = list(numGenes = 1000, normalisation = "logNormalize")))
  )

  scdata <- suppressWarnings(run_dataIntegration(scdata, config))
  expect_s4_class(scdata, "Seurat")
})

test_that("Unisample integration works", {
  scdata <- mock_scdata()
  config <- list(
    dimensionalityReduction = list(numPCs = 2),
    dataIntegration = list(method = "unisample", methodSettings = list(unisample = list(numGenes = 1000, normalisation = "logNormalize")))
  )

  scdata <- suppressWarnings(run_dataIntegration(scdata, config))
  expect_s4_class(scdata, "Seurat")
})

test_that("FastMNN is not working", {
  scdata <- mock_scdata()
  config <- list(
    dimensionalityReduction = list(numPCs = 2),
    dataIntegration = list(method = "fastmnn", methodSettings = list(fastmnn = list(numGenes = 1000, normalisation = "logNormalize")))
  )

  expect_error(suppressWarnings(run_dataIntegration(scdata, config)))
})
