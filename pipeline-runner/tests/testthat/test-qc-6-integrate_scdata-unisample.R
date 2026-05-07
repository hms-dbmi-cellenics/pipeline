library(Seurat)
human_cc_genes <- cc_genes[["human"]]

mock_prev_out <- function(samples = "sample_a", counts = NULL, use_bpcells = FALSE) {
  if (is.null(counts)) {
    counts <- DropletUtils:::simCounts()
    colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  }

  eout <- DropletUtils::emptyDrops(counts)

  counts_list <- list()
  edrops <- list()
  doublet_scores <- list()

  for (sample in samples) {

    counts_list[[sample]] <- maybe_bpcells(
      counts,
      withr::local_tempfile(.local_envir = parent.frame()),
      use_bpcells
    )
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

test_that("Unisample integration works", {
  scdata_list <- mock_scdata_list()
  cells_id <- mock_ids(scdata_list)
  config <- list(
    dimensionalityReduction = list(numPCs = 2),
    dataIntegration = list(method = UNISAMPLE, methodSettings = list(unisample = list(numGenes = 1000, normalisation = "logNormalize")))
  )

  integrated_scdata <- suppressWarnings(run_unisample(scdata_list, config, cells_id))
  integrated_scdata <- clean_timestamp(integrated_scdata)
  integrated_scdata <- remove_commands_functions(integrated_scdata)

  expect_s4_class(integrated_scdata, "Seurat")

  skip_on_ci()
  expect_snapshot(str(integrated_scdata))
})
