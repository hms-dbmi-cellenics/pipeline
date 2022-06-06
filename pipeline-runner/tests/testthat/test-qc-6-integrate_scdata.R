human_cc_genes <- cc_genes[["human"]]

mock_scdata <- function(rename_genes = c()) {
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )

  if (length(rename_genes) > 0) {
    # rename some genes to match cell cycle genes
    some_genes <- sample(1:nrow(pbmc_raw), length(rename_genes))
    rownames(pbmc_raw)[some_genes] <- rename_genes
  }

  scdata <- Seurat::CreateSeuratObject(counts = pbmc_raw)

  # add samples
  scdata$samples <- rep(c("123abc", "123def"), each = 40)
  scdata$cells_id <- 0:79
  scdata@misc$gene_annotations$input <- rownames(scdata)

  # scale and PCA
  scdata <- Seurat::NormalizeData(scdata, normalization.method = "LogNormalize", verbose = FALSE)
  scdata <- Seurat::FindVariableFeatures(scdata, verbose = FALSE)
  scdata <- Seurat::ScaleData(scdata, verbose = FALSE)
  scdata <- Seurat::RunPCA(scdata, verbose = FALSE)
  scdata@misc[["active.reduction"]] <- "pca"

  return(scdata)
}

test_that("Integrate scdata works",{
  scdata <- mock_scdata()
  cells_id <- 0:79
  config <- list(
    dimensionalityReduction = list(numPCs = 2),
    dataIntegration = list(method = "harmony", methodSettings = list(harmony = list(numGenes = 10, normalisation = "logNormalize")))
  )

  scdata <- suppressWarnings(integrate_scdata(scdata, config, "", cells_id, task_name = "dataIntegration"))$data
  expect_s4_class(scdata, "Seurat")
  expect_equal(ncol(scdata),80)
})

test_that("Integrate scdata filters out cells ids",{
  scdata <- mock_scdata()
  cells_id <- 0:40
  config <- list(
    dimensionalityReduction = list(numPCs = 2),
    dataIntegration = list(method = "harmony", methodSettings = list(harmony = list(numGenes = 10, normalisation = "logNormalize")))
  )

  scdata <- suppressWarnings(integrate_scdata(scdata, config, "", cells_id, task_name = "dataIntegration"))$data
  expect_lt(ncol(scdata),80)
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
  scdata <- suppressWarnings(mock_scdata())
  npcs <- get_npcs(scdata)
  expect_lte(npcs, 30)
})


test_that("build_cc_gene_list correctly makes the list of cc genes when there are matches", {

  n_rename <- 10
  some_cc_genes <- sample(human_cc_genes$symbol, n_rename)
  scdata <- suppressWarnings(mock_scdata(rename_genes = some_cc_genes))
  all_genes <- scdata@misc$gene_annotations$input

  expected_res <- match(some_cc_genes, all_genes)
  res <- build_cc_gene_list(all_genes)

  expect_setequal(res, expected_res)

})


test_that("build_cc_gene_list returns empty int vector when there aren't matches", {

  scdata <- suppressWarnings(mock_scdata())
  all_genes <- scdata@misc$gene_annotations$input

  # empty integer vector
  expected_res <- integer()
  res <- build_cc_gene_list(all_genes)

  expect_equal(res, expected_res)
})


test_that("list_exclude_genes adds custom genes to exclusion", {

  n_rename <- 10
  some_cc_genes <- sample(human_cc_genes$symbol, n_rename)
  scdata <- suppressWarnings(mock_scdata(rename_genes = some_cc_genes))
  all_genes <- scdata@misc$gene_annotations$input

  exclude_custom <- sample(setdiff(all_genes, some_cc_genes), 7)
  exclude_custom_indices <- match(exclude_custom, all_genes)

  expected_res <- c(build_cc_gene_list(all_genes), exclude_custom_indices)

  res <- list_exclude_genes(all_genes, list("cellCycle"), exclude_custom)

  expect_setequal(res, expected_res)

})

test_that("remove_genes removes the correct genes when there are genes to remove", {

  n_rename <- 10
  some_cc_genes <- sample(human_cc_genes$symbol, n_rename)
  scdata <- suppressWarnings(mock_scdata(rename_genes = some_cc_genes))
  all_genes <- scdata@misc$gene_annotations$input

  res <- remove_genes(scdata, exclude_groups = list("cellCycle"))

  # only cc genes
  expect_equal(nrow(res), nrow(scdata) - 10)
  expect_false(any(some_cc_genes %in% rownames(res)))

  exclude_custom <- sample(setdiff(all_genes, some_cc_genes), 7)

  # cc genes and custom
  res <- remove_genes(scdata, exclude_groups = list("cellCycle"), exclude_custom)

  expect_equal(nrow(res), nrow(scdata) - 17)
  expect_false(any(c(some_cc_genes, exclude_custom) %in% rownames(res)))

})


test_that("remove_genes doesn't modify the object when there are no matches", {

  scdata <- suppressWarnings(mock_scdata())

  # empty integer vector
  res <- remove_genes(scdata, exclude_groups = "cellCycle")

  expect_equal(res, scdata)

})
