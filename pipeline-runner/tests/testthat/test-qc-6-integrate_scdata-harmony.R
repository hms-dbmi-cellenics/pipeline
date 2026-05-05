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

test_that("harmony integration works", {
  scdata_list <- mock_scdata_list()
  cells_id <- mock_ids(scdata_list)
  config <- list(
    dimensionalityReduction = list(numPCs = 2),
    dataIntegration = list(
      method = "harmony",
      methodSettings = list(
        harmony = list(
          numGenes = 10,
          normalisation = "logNormalize"
        )
      )
    )
  )

  integrated_scdata <- suppressWarnings(
    run_harmony(scdata_list, config, cells_id)
  )
  integrated_scdata <- clean_timestamp(integrated_scdata)
  integrated_scdata <- remove_commands_functions(integrated_scdata)

  expect_s4_class(integrated_scdata, "Seurat")
  skip_if(is_bpcells())
  expect_snapshot(str(integrated_scdata))
})
