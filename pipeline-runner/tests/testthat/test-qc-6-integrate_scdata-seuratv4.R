library(Seurat)

mock_prev_out <- function(samples = "sample_a", counts = NULL) {
  if (is.null(counts)) {
    counts <- DropletUtils:::simCounts()
    colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  }

  eout <- DropletUtils::emptyDrops(counts)

  counts_list <- list()
  edrops <- list()
  doublet_scores <- list()

  for (sample in samples) {
    counts_list[[sample]] <- counts
    edrops[[sample]] <- eout
    doublet_scores[[sample]] <- mock_doublet_scores(counts)
  }

  annot <- data.frame(name = row.names(counts), input = row.names(counts))

  # as passed to create_seurat
  prev_out <- list(
    counts_list = counts_list,
    edrops = edrops,
    doublet_scores = doublet_scores,
    annot = annot,
    config = list(name = "project name")
  )

  # call create_seurat to get prev_out to pass to prepare_experiment
  prev_out <- create_seurat(NULL, NULL, prev_out)$output


  prev_out$scdata_list <- add_metadata_to_samples(prev_out$scdata_list, annot = annot, experiment_id = "expid")

  return(prev_out)
}

mock_scdata <- function(rename_genes = c(), n_rep = 1) {
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )
  # replicate matrix columns n times to create a bigger mock dataset
  pbmc_raw <- do.call("cbind", replicate(n_rep, pbmc_raw, simplify = FALSE))

  if (length(rename_genes) > 0) {
    # rename some genes to match cell cycle genes
    some_genes <- sample(1:nrow(pbmc_raw), length(rename_genes))
    rownames(pbmc_raw)[some_genes] <- rename_genes
  }

  sample_1_id <- "123abc"
  sample_2_id <- "123def"

  scdata <- Seurat::CreateSeuratObject(counts = pbmc_raw)

  # add samples
  scdata$samples <- rep(c(sample_1_id, sample_2_id), each = ncol(scdata) / 2)
  scdata$cells_id <- 0:(ncol(scdata) - 1)
  scdata@misc$gene_annotations <- data.frame(input = rownames(scdata), name = paste("SYMBOL -", rownames(scdata)))
  rownames(scdata@misc$gene_annotations) <- rownames(scdata)

  scdata_sample1 <- subset(scdata, samples == sample_1_id)
  scdata_sample2 <- subset(scdata, samples == sample_2_id)


  scdata_list <- list(scdata_sample1, scdata_sample2)
  names(scdata_list) <- c(sample_1_id, sample_2_id)

  return(list(scdata_list, sample_1_id, sample_2_id))
}

mock_ids <- function() {
  return(list("123abc" = 0:39, "123def" = 40:79))
}


scdata_preprocessing <- function(scdata) {

  # scale and PCA
  scdata <- Seurat::NormalizeData(scdata, normalization.method = "LogNormalize", verbose = FALSE)
  scdata <- Seurat::FindVariableFeatures(scdata, verbose = FALSE)
  scdata <- Seurat::ScaleData(scdata, verbose = FALSE)
  scdata <- Seurat::RunPCA(scdata, verbose = FALSE)
  scdata@misc[["active.reduction"]] <- "pca"

  return(scdata)
}

mock_doublet_scores <- function(counts) {
  doublet_scores <- runif(ncol(counts))
  doublet_class <- ifelse(doublet_scores < 0.8, "singlet", "doublet")

  data.frame(
    row.names = colnames(counts),
    barcodes = colnames(counts),
    doublet_class = doublet_class,
    doublet_scores = doublet_scores
  )
}


#' Fix time stamps in seurat object
#'
#' Seurat adds logs to certain command runs, with time stamps that make snapshot
#' tests fail. This function replaces them with a fixed datetime object.
#'
#' @param scdata seurat object
#'
#' @return seurat object with fixed time stamps
#'
clean_timestamp <- function(scdata) {
  fixed_datetime <- as.POSIXct("1991-12-19 05:23:00", tz = "UTC")

  for (slot in names(scdata@commands)) {
    scdata@commands[[slot]]@time.stamp <- fixed_datetime
  }

  if ("ingestionDate" %in% names(scdata@misc)) {
    scdata@misc$ingestionDate <- fixed_datetime
  }

  return(scdata)
}


test_that("SeuratV4 integration doesnt error out with small dataset", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  cells_id <- mock_ids()

  config <- list(
    dimensionalityReduction = list(numPCs = 2, method = "rpca"),
    dataIntegration = list(method = "seuratv4", methodSettings = list(seuratv4 = list(numGenes = 1000, normalisation = "logNormalize")))
  )

  integrated_scdata <- suppressWarnings(integrate_scdata(scdata_list, config, "", cells_id, task_name = "dataIntegration"))$data
  expect_s4_class(integrated_scdata, "Seurat")
})


test_that("SeuratV4 integration works", {
  # mock a bigger dataset to run Seurat v4 integration without skipping it
  c(scdata_list, sample_1_id, sample_2_id) %<-% suppressWarnings(mock_scdata(n_rep = 3))
  cells_id <- list("123abc" = scdata_list$`123abc`$cells_id, "123def" = scdata_list$`123def`$cells_id)

  config <- list(
    dimensionalityReduction = list(method = "rpca"),
    dataIntegration = list(method = "seuratv4", methodSettings = list(seuratv4 = list(numGenes = 1000, normalisation = "logNormalize")))
  )

  integrated_scdata <- suppressWarnings(integrate_scdata(scdata_list, config, "", cells_id, task_name = "dataIntegration"))$data
  expect_s4_class(integrated_scdata, "Seurat")
})


test_that("PCA is computed when RPCA method is selected within SeuratV4 integration", {
  # mock a bigger dataset to run Seurat v4 integration without skipping it
  c(scdata_list, sample_1_id, sample_2_id) %<-% suppressWarnings(mock_scdata(n_rep = 3))
  cells_id <- list("123abc" = scdata_list$`123abc`$cells_id, "123def" = scdata_list$`123def`$cells_id)

  config <- list(
    dimensionalityReduction = list(method = "rpca"),
    dataIntegration = list(method = "seuratv4", methodSettings = list(seuratv4 = list(numGenes = 1000, normalisation = "logNormalize")))
  )

  expect_message(integrate_scdata(scdata_list, config, "", cells_id, task_name = "dataIntegration"), "Running PCA")
})

test_that("PCA is not computed when CCA method is selected within SeuratV4 integration", {
  # mock a bigger dataset to run Seurat v4 integration without skipping it
  c(scdata_list, sample_1_id, sample_2_id) %<-% suppressWarnings(mock_scdata(n_rep = 3))
  cells_id <- list("123abc" = scdata_list$`123abc`$cells_id, "123def" = scdata_list$`123def`$cells_id)

  config <- list(
    dimensionalityReduction = list(method = "cca"),
    dataIntegration = list(method = "seuratv4", methodSettings = list(seuratv4 = list(numGenes = 1000, normalisation = "logNormalize")))
  )

  expect_message(integrate_scdata(scdata_list, config, "", cells_id, task_name = "dataIntegration"), "PCA is not running .*")
})


test_that("SeuratV4 integration finds integration anchors using RPCA method, if method in config is RPCA", {
  # mock a bigger dataset to run Seurat v4 integration without skipping it
  c(scdata_list, sample_1_id, sample_2_id) %<-% suppressWarnings(mock_scdata(n_rep = 3))
  cells_id <- list("123abc" = scdata_list$`123abc`$cells_id, "123def" = scdata_list$`123def`$cells_id)

  config <- list(
    dimensionalityReduction = list(method = "rpca"),
    dataIntegration = list(method = "seuratv4", methodSettings = list(seuratv4 = list(numGenes = 1000, normalisation = "logNormalize")))
  )

  integrated_scdata <- suppressWarnings(integrate_scdata(scdata_list, config, "", cells_id, task_name = "dataIntegration"))$data
  expect_equal(integrated_scdata@commands$FindIntegrationAnchors$reduction, "pca")
})


test_that("SeuratV4 integration finds integration anchors using CCA method, if method in config is CCA", {
  # mock a bigger dataset to run Seurat v4 integration without skipping it
  c(scdata_list, sample_1_id, sample_2_id) %<-% suppressWarnings(mock_scdata(n_rep = 3))
  cells_id <- list("123abc" = scdata_list$`123abc`$cells_id, "123def" = scdata_list$`123def`$cells_id)

  config <- list(
    dimensionalityReduction = list(method = "cca"),
    dataIntegration = list(method = "seuratv4", methodSettings = list(seuratv4 = list(numGenes = 1000, normalisation = "logNormalize")))
  )

  integrated_scdata <- suppressWarnings(integrate_scdata(scdata_list, config, "", cells_id, task_name = "dataIntegration"))$data
  expect_equal(integrated_scdata@commands$FindIntegrationAnchors$reduction, "cca")
})

test_that("SCTransform integration works", {
  # mock a bigger dataset to run Seurat v4 integration without skipping it
  c(scdata_list, sample_1_id, sample_2_id) %<-% suppressWarnings(mock_scdata(n_rep = 3))
  cells_id <- list("123abc" = scdata_list$`123abc`$cells_id, "123def" = scdata_list$`123def`$cells_id)

  config <- list(
    dimensionalityReduction = list(method = "rpca"),
    dataIntegration = list(method = "seuratv4", methodSettings = list(seuratv4 = list(numGenes = 1000, normalisation = "SCT")))
  )

  integrated_scdata <- suppressWarnings(integrate_scdata(scdata_list, config, "", cells_id, task_name = "dataIntegration"))$data
  expect_s4_class(integrated_scdata, "Seurat")
  expect_s4_class(integrated_scdata[["SCT"]], "SCTAssay")
  expect_equal(Seurat::DefaultAssay(integrated_scdata), "integrated")
})

test_that("create_scdata merge only the counts when merge_data is FALSE", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% suppressWarnings(mock_scdata(n_rep = 3))
  cells_id <- list("123abc" = scdata_list$`123abc`$cells_id, "123def" = scdata_list$`123def`$cells_id)

  for (i in 1:length(scdata_list)) {
    scdata_list[[i]] <- Seurat::NormalizeData(scdata_list[[i]], assay = "RNA", normalization.method = "LogNormalize", verbose = FALSE)
    scdata_list[[i]] <- Seurat::FindVariableFeatures(scdata_list[[i]], assay = "RNA", nfeatures = 2000, verbose = FALSE)
    scdata_list[[i]] <- Seurat::ScaleData(scdata_list[[i]], verbose = FALSE)
  }

  merged_data <- create_scdata(scdata_list, cells_id, merge_data = FALSE)
  expect_equal(merged_data@assays$RNA@counts, merged_data@assays$RNA@data)
})

test_that("create_scdata merge only the counts when merge_data is not passed", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% suppressWarnings(mock_scdata(n_rep = 3))
  cells_id <- list("123abc" = scdata_list$`123abc`$cells_id, "123def" = scdata_list$`123def`$cells_id)

  for (i in 1:length(scdata_list)) {
    scdata_list[[i]] <- Seurat::NormalizeData(scdata_list[[i]], assay = "RNA", normalization.method = "LogNormalize", verbose = FALSE)
    scdata_list[[i]] <- Seurat::FindVariableFeatures(scdata_list[[i]], assay = "RNA", nfeatures = 2000, verbose = FALSE)
    scdata_list[[i]] <- Seurat::ScaleData(scdata_list[[i]], verbose = FALSE)
  }

  merged_data <- create_scdata(scdata_list, cells_id)
  expect_equal(merged_data@assays$RNA@counts, merged_data@assays$RNA@data)
})

test_that("create_scdata merge the data slots when merge_data is TRUE", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% suppressWarnings(mock_scdata(n_rep = 3))
  cells_id <- list("123abc" = scdata_list$`123abc`$cells_id, "123def" = scdata_list$`123def`$cells_id)

  for (i in 1:length(scdata_list)) {
    scdata_list[[i]] <- Seurat::NormalizeData(scdata_list[[i]], assay = "RNA", normalization.method = "LogNormalize", verbose = FALSE)
    scdata_list[[i]] <- Seurat::FindVariableFeatures(scdata_list[[i]], assay = "RNA", nfeatures = 2000, verbose = FALSE)
    scdata_list[[i]] <- Seurat::ScaleData(scdata_list[[i]], verbose = FALSE)
  }

  merged_data <- create_scdata(scdata_list, cells_id, merge_data = TRUE)
  expect_false(identical(merged_data@assays$RNA@counts, merged_data@assays$RNA@data))
})


test_that("misc slot is complete after Seurat V4 integration", {

  # mock a bigger dataset to run Seurat v4 integration without skipping it
  c(scdata_list, sample_1_id, sample_2_id) %<-% suppressWarnings(mock_scdata(n_rep = 3))
  cells_id <- list("123abc" = scdata_list$`123abc`$cells_id, "123def" = scdata_list$`123def`$cells_id)

  config <- list(
    dimensionalityReduction = list(method = "rpca"),
    dataIntegration = list(method = "seuratv4", methodSettings = list(seuratv4 = list(numGenes = 1000, normalisation = "logNormalize")))
  )

  integrated_scdata <- suppressWarnings(integrate_scdata(scdata_list, config, "", cells_id, task_name = "dataIntegration"))$data
  integrated_scdata <- clean_timestamp(integrated_scdata)

  expect_s4_class(integrated_scdata, "Seurat")
  expect_snapshot(str(integrated_scdata))
  expect_snapshot(str(integrated_scdata@misc))

})


test_that("misc slot is complete after Seurat V4 integration with geosketch", {

  # mock a bigger dataset to run Seurat v4 integration without skipping it
  c(scdata_list, sample_1_id, sample_2_id) %<-% suppressWarnings(mock_scdata(n_rep = 3))
  cells_id <- list("123abc" = scdata_list$`123abc`$cells_id, "123def" = scdata_list$`123def`$cells_id)

  config <- list(
    dimensionalityReduction = list(numPCs = 10, method = "rpca"),
    dataIntegration = list(method = "seuratv4", methodSettings = list(seuratv4 = list(numGenes = 10, normalisation = "logNormalize"))),
    downsampling = list(method = "geosketch", methodSettings = list(geosketch = list(percentageToKeep = 50))))


  integrated_scdata <- suppressWarnings(integrate_scdata(scdata_list, config, "", cells_id, task_name = "dataIntegration")$data)
  expect_s4_class(integrated_scdata, "Seurat")
  integrated_scdata@misc$ingestionDate <- "fixed_date"

  expected_misc_names <- c("gene_annotations", "color_pool", "ingestionDate", "active.reduction", "numPCs", "geosketch")

  expect_setequal(names(integrated_scdata@misc), expected_misc_names)
})
