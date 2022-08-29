human_cc_genes <- cc_genes[["human"]]

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

test_that("Integrate scdata works", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  cells_id <- mock_ids()
  config <- list(
    dimensionalityReduction = list(numPCs = 2),
    dataIntegration = list(method = "harmony", methodSettings = list(harmony = list(numGenes = 10, normalisation = "logNormalize")))
  )

  integrated_scdata <- suppressWarnings(integrate_scdata(scdata_list, config, "", cells_id, task_name = "dataIntegration"))$data
  expect_s4_class(integrated_scdata, "Seurat")
  expect_equal(ncol(integrated_scdata), 80)
})

test_that("Integrate scdata filters out cells ids", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  cells_id <- mock_ids()
  cells_id[[sample_1_id]] <- cells_id[[sample_1_id]][-c(23:31)]
  cells_id[[sample_2_id]] <- cells_id[[sample_2_id]][-c(40:47)]

  config <- list(
    dimensionalityReduction = list(numPCs = 2),
    dataIntegration = list(method = "harmony", methodSettings = list(harmony = list(numGenes = 10, normalisation = "logNormalize")))
  )

  scdata <- suppressWarnings(integrate_scdata(scdata_list, config, "", cells_id, task_name = "dataIntegration"))$data
  expect_lt(ncol(scdata), 80)
})

test_that("harmony integration works", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  merged_scdata <- create_scdata(scdata_list)
  config <- list(
    dimensionalityReduction = list(numPCs = 2),
    dataIntegration = list(method = "harmony", methodSettings = list(harmony = list(numGenes = 10, normalisation = "logNormalize")))
  )

  integrated_scdata <- suppressWarnings(run_dataIntegration(merged_scdata, config))
  expect_s4_class(integrated_scdata, "Seurat")
})

test_that("SeuratV4 integration doesnt error out with small dataset", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  merged_scdata <- create_scdata(scdata_list)
  config <- list(
    dimensionalityReduction = list(numPCs = 2, method = "rpca"),
    dataIntegration = list(method = "seuratv4", methodSettings = list(seuratv4 = list(numGenes = 1000, normalisation = "logNormalize")))
  )

  integrated_scdata <- suppressWarnings(run_dataIntegration(merged_scdata, config))
  expect_s4_class(integrated_scdata, "Seurat")
})

test_that("Unisample integration works", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  merged_scdata <- create_scdata(scdata_list)
  config <- list(
    dimensionalityReduction = list(numPCs = 2),
    dataIntegration = list(method = "unisample", methodSettings = list(unisample = list(numGenes = 1000, normalisation = "logNormalize")))
  )

  integrated_scdata <- suppressWarnings(run_dataIntegration(merged_scdata, config))
  expect_s4_class(integrated_scdata, "Seurat")
})

test_that("FastMNN is not working", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  merged_scdata <- create_scdata(scdata_list)
  config <- list(
    dimensionalityReduction = list(numPCs = 2),
    dataIntegration = list(method = "fastmnn", methodSettings = list(fastmnn = list(numGenes = 1000, normalisation = "logNormalize")))
  )

  expect_error(suppressWarnings(run_dataIntegration(merged_scdata, config)))
})

test_that("numPCs estimation works", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  merged_scdata <- create_scdata(scdata_list)
  merged_scdata <- suppressWarnings(scdata_preprocessing(merged_scdata))
  npcs <- get_npcs(merged_scdata)
  expect_lte(npcs, 30)
})


test_that("build_cc_gene_list correctly makes the list of cc genes when there are matches", {
  n_rename <- 10
  some_cc_genes <- sample(human_cc_genes$symbol, n_rename)
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata(rename_genes = some_cc_genes)
  merged_scdata <- create_scdata(scdata_list)
  all_genes <- merged_scdata@misc$gene_annotations$input

  expected_res <- match(some_cc_genes, all_genes)
  res <- build_cc_gene_list(all_genes)

  expect_setequal(res, expected_res)
})


test_that("build_cc_gene_list returns empty int vector when there aren't matches", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  merged_scdata <- create_scdata(scdata_list)
  all_genes <- merged_scdata@misc$gene_annotations$input

  # empty integer vector
  expected_res <- integer()
  res <- build_cc_gene_list(all_genes)

  expect_equal(res, expected_res)
})


test_that("list_exclude_genes adds custom genes to exclusion", {
  n_rename <- 10
  some_cc_genes <- sample(human_cc_genes$symbol, n_rename)
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata(rename_genes = some_cc_genes)
  merged_scdata <- create_scdata(scdata_list)
  all_genes <- merged_scdata@misc$gene_annotations$input

  exclude_custom <- sample(setdiff(all_genes, some_cc_genes), 7)
  exclude_custom_indices <- match(exclude_custom, all_genes)

  expected_res <- c(build_cc_gene_list(all_genes), exclude_custom_indices)

  res <- list_exclude_genes(all_genes, list("cellCycle"), exclude_custom)

  expect_setequal(res, expected_res)
})

test_that("remove_genes removes the correct genes when there are genes to remove", {
  n_rename <- 10
  some_cc_genes <- sample(human_cc_genes$symbol, n_rename)
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata(rename_genes = some_cc_genes)
  merged_scdata <- create_scdata(scdata_list)
  all_genes <- merged_scdata@misc$gene_annotations$input

  res <- remove_genes(merged_scdata, exclude_groups = list("cellCycle"))

  # only cc genes
  expect_equal(nrow(res), nrow(merged_scdata) - 10)
  expect_false(any(some_cc_genes %in% rownames(res)))

  exclude_custom <- sample(setdiff(all_genes, some_cc_genes), 7)

  # cc genes and custom
  res <- remove_genes(merged_scdata, exclude_groups = list("cellCycle"), exclude_custom)

  expect_equal(nrow(res), nrow(merged_scdata) - 17)
  expect_false(any(c(some_cc_genes, exclude_custom) %in% rownames(res)))
})


test_that("remove_genes doesn't modify the object when there are no matches", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  merged_scdata <- create_scdata(scdata_list)
  # empty integer vector
  res <- remove_genes(merged_scdata, exclude_groups = "cellCycle")

  expect_equal(res, merged_scdata)
})


test_that("SeuratV4 integration works", {
  # mock a bigger dataset to run Seurat v4 integration without skipping it
  c(scdata_list, sample_1_id, sample_2_id) %<-% suppressWarnings(mock_scdata(n_rep = 3))
  merged_scdata <- create_scdata(scdata_list)

  merged_scdata <- suppressWarnings(scdata_preprocessing(merged_scdata))
  npcs <- get_npcs(merged_scdata)

  config <- list(
    dimensionalityReduction = list(numPCs = npcs, method = "rpca"),
    dataIntegration = list(method = "seuratv4", methodSettings = list(seuratv4 = list(numGenes = 1000, normalisation = "logNormalize")))
  )

  merged_scdata <- suppressWarnings(run_dataIntegration(merged_scdata, config))
  expect_s4_class(merged_scdata, "Seurat")
})


test_that("PCA is computed when RPCA method is selected within SeuratV4 integration", {
  # mock a bigger dataset to run Seurat v4 integration without skipping it
  c(scdata_list, sample_1_id, sample_2_id) %<-% suppressWarnings(mock_scdata(n_rep = 3))
  merged_scdata <- create_scdata(scdata_list)

  merged_scdata <- suppressWarnings(scdata_preprocessing(merged_scdata))
  npcs <- get_npcs(merged_scdata)

  config <- list(
    dimensionalityReduction = list(numPCs = npcs, method = "rpca"),
    dataIntegration = list(method = "seuratv4", methodSettings = list(seuratv4 = list(numGenes = 1000, normalisation = "logNormalize")))
  )

  expect_message(run_dataIntegration(merged_scdata, config), "Running PCA")
})


test_that("PCA is not computed when CCA method is selected within SeuratV4 integration", {
  # mock a bigger dataset to run Seurat v4 integration without skipping it
  c(scdata_list, sample_1_id, sample_2_id) %<-% suppressWarnings(mock_scdata(n_rep = 3))
  merged_scdata <- create_scdata(scdata_list)

  merged_scdata <- suppressWarnings(scdata_preprocessing(merged_scdata))
  npcs <- get_npcs(merged_scdata)

  config <- list(
    dimensionalityReduction = list(numPCs = npcs, method = "cca"),
    dataIntegration = list(method = "seuratv4", methodSettings = list(seuratv4 = list(numGenes = 1000, normalisation = "logNormalize")))
  )

  expect_message(run_dataIntegration(merged_scdata, config), "PCA is not running .*")
})


test_that("SeuratV4 integration finds integration anchors using RPCA method, if method in config is RPCA", {
  # mock a bigger dataset to run Seurat v4 integration without skipping it
  c(scdata_list, sample_1_id, sample_2_id) %<-% suppressWarnings(mock_scdata(n_rep = 3))
  merged_scdata <- create_scdata(scdata_list)

  merged_scdata <- suppressWarnings(scdata_preprocessing(merged_scdata))
  npcs <- get_npcs(merged_scdata)

  config <- list(
    dimensionalityReduction = list(numPCs = npcs, method = "rpca"),
    dataIntegration = list(method = "seuratv4", methodSettings = list(seuratv4 = list(numGenes = 1000, normalisation = "logNormalize")))
  )

  expect_message(merged_scdata <- suppressWarnings(run_dataIntegration(merged_scdata, config)), "Finding integration anchors using RPCA reduction")
  expect_equal(merged_scdata@commands$FindIntegrationAnchors$reduction, "pca")
})


test_that("SeuratV4 integration finds integration anchors using CCA method, if method in config is CCA", {
  # mock a bigger dataset to run Seurat v4 integration without skipping it
  c(scdata_list, sample_1_id, sample_2_id) %<-% suppressWarnings(mock_scdata(n_rep = 3))
  merged_scdata <- create_scdata(scdata_list)

  merged_scdata <- suppressWarnings(scdata_preprocessing(merged_scdata))
  npcs <- get_npcs(merged_scdata)

  config <- list(
    dimensionalityReduction = list(numPCs = npcs, method = "cca"),
    dataIntegration = list(method = "seuratv4", methodSettings = list(seuratv4 = list(numGenes = 1000, normalisation = "logNormalize")))
  )

  expect_message(merged_scdata <- suppressWarnings(run_dataIntegration(merged_scdata, config)), "Finding integration anchors using CCA reduction")
  expect_equal(merged_scdata@commands$FindIntegrationAnchors$reduction, "cca")
})


test_that("merge_scdata_list correctly merges seurat objects", {
  prev_out <- mock_prev_out(samples = c("a", "b", "c"))
  scdata_list <- prev_out$scdata_list


  scdata <- suppressWarnings(merge_scdata_list(scdata_list))

  expect_equal(sum(unlist(lapply(scdata_list, ncol))), ncol(scdata))
  expect_true(all(scdata$samples[1:ncol(scdata_list[[1]])] == "a"))

})
