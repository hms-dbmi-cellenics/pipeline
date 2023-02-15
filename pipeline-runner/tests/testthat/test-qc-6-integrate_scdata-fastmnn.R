library(Seurat)
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


test_that("FastMNN works", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  cells_id <- mock_ids()
  merged_scdata <- create_scdata(scdata_list, cells_id)
  config <- list(
    dimensionalityReduction = list(numPCs = 2),
    dataIntegration = list(method = "fastmnn", methodSettings = list(fastmnn = list(numGenes = 1000, normalisation = "logNormalize")))
  )

  integrated_scdata <- suppressWarnings(run_dataIntegration(merged_scdata, scdata_sketch = NA, config = config))
  expect_s4_class(integrated_scdata, "Seurat")
})


test_that("normalize_data doesn't scale data if integration method is FastMNN", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  cells_id <- mock_ids()
  merged_scdata <- create_scdata(scdata_list, cells_id)
  config <- list(
    dimensionalityReduction = list(numPCs = 5),
    dataIntegration = list(method = "fastmnn", methodSettings = list(fastmnn = list(numGenes = 1000, normalisation = "logNormalize")))
  )

  merged_scdata <- log_normalize(merged_scdata, "LogNormalize", config$dataIntegration$method, config$dataIntegration$methodSettings$fastmnn$numGenes)
  expect_equal(dim(merged_scdata@assays$RNA@scale.data), c(0, 0))
})
