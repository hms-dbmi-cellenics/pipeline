mock_scdata <- function() {
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )

  scdata <- Seurat::CreateSeuratObject(counts = pbmc_raw)

  # add samples
  scdata$samples <- rep(c("123abc", "123def"), each = 40)
  scdata$cells_id <- 0:79

  # scale and PCA
  scdata <- Seurat::NormalizeData(scdata, normalization.method = "LogNormalize", verbose = FALSE)
  scdata <- Seurat::FindVariableFeatures(scdata, verbose = FALSE)
  scdata <- Seurat::ScaleData(scdata, verbose = FALSE)
  scdata <- Seurat::RunPCA(scdata, verbose = FALSE)

  return(scdata)
}

test_that("Integrate scdata works", {
  scdata <- mock_scdata()
  cells_id <- 0:79
  config <- list(
    dimensionalityReduction = list(numPCs = 2),
    dataIntegration = list(method = "harmony", methodSettings = list(harmony = list(numGenes = 10, normalisation = "logNormalize")))
  )

  scdata <- suppressWarnings(integrate_scdata(scdata, config, "", cells_id, task_name = "dataIntegration"))$data
  expect_s4_class(scdata, "Seurat")
  expect_equal(ncol(scdata), 80)
})

test_that("Integrate scdata filters out cells ids", {
  scdata <- mock_scdata()
  cells_id <- 0:40
  config <- list(
    dimensionalityReduction = list(numPCs = 2),
    dataIntegration = list(method = "harmony", methodSettings = list(harmony = list(numGenes = 10, normalisation = "logNormalize")))
  )

  scdata <- suppressWarnings(integrate_scdata(scdata, config, "", cells_id, task_name = "dataIntegration"))$data
  expect_lt(ncol(scdata), 80)
})

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

test_that("numPCs estimation works", {
  scdata <- mock_scdata()
  config <- list(dimensionalityReduction = list(numPCs = NULL))
  expect_gte(estimate_npcs(scdata, var_threshold = 0.85, max_npcs = 30), 0)
})
