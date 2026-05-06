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




scdata_preprocessing <- function(scdata) {

  # scale and PCA
  scdata <- Seurat::NormalizeData(
    scdata, normalization.method = "LogNormalize", verbose = FALSE
  )
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

test_that("SeuratV4 integration doesnt error out with small dataset", {
  scdata_list <- mock_scdata_list()
  cells_id <- mock_ids(scdata_list)

  config <- list(
    dimensionalityReduction = list(numPCs = 2, method = "rpca"),
    dataIntegration = list(
      method = "seuratv4",
      methodSettings = list(
        seuratv4 = list(
          numGenes = 1000,
          normalisation = "logNormalize"
        )
      )
    )
  )

  integrated_scdata <- suppressWarnings(
    integrate_scdata(
      scdata_list, config, "", cells_id, task_name = "dataIntegration"
    )$data
  )
  expect_s4_class(integrated_scdata, "Seurat")
})


test_that("SeuratV4 integration works with RPCA method", {
  # mock a bigger dataset to run Seurat v4 integration without skipping it
  scdata_list <- mock_scdata_list(n_rep=3)
  cells_id <- list(
    "123abc" = scdata_list$`123abc`$cells_id,
    "123def" = scdata_list$`123def`$cells_id
  )

  config <- list(
    dimensionalityReduction = list(method = "rpca"),
    dataIntegration = list(
      method = "seuratv4",
      methodSettings = list(
        seuratv4 = list(numGenes = 1000, normalisation = "logNormalize")
      )
    )
  )

  integrated_scdata <- suppressWarnings(
    integrate_scdata(
      scdata_list, config, "", cells_id, task_name = "dataIntegration"
    )$data
  )
  expect_s4_class(integrated_scdata, "Seurat")
})



test_that("SeuratV4 integration works with CCA method", {
  # mock a bigger dataset to run Seurat v4 integration without skipping it
  scdata_list <- mock_scdata_list(n_rep=3)
  cells_id <- list(
    "123abc" = scdata_list$`123abc`$cells_id,
    "123def" = scdata_list$`123def`$cells_id
  )

  config <- list(
    dimensionalityReduction = list(method = "cca"),
    dataIntegration = list(
      method = "seuratv4",
      methodSettings = list(
        seuratv4 = list(numGenes = 1000, normalisation = "logNormalize")
      )
    )
  )

  integrated_scdata <- suppressWarnings(
    integrate_scdata(
      scdata_list, config, "", cells_id, task_name = "dataIntegration"
    )$data
  )
  expect_s4_class(integrated_scdata, "Seurat")
})


test_that("SeuratV4 integration finds integration anchors using RPCA method, if method in config is RPCA", {
  # mock a bigger dataset to run Seurat v4 integration without skipping it
  scdata_list <- mock_scdata_list(n_rep=3)
  cells_id <- list(
    "123abc" = scdata_list$`123abc`$cells_id,
    "123def" = scdata_list$`123def`$cells_id
  )

  config <- list(
    dimensionalityReduction = list(method = "rpca"),
    dataIntegration = list(
      method = "seuratv4",
      methodSettings = list(
        seuratv4 = list(numGenes = 1000, normalisation = "logNormalize")
      )
    )
  )

  integrated_scdata <- suppressWarnings(
    integrate_scdata(
      scdata_list, config, "", cells_id, task_name = "dataIntegration"
    )$data
  )
  expect_contains(
    Seurat::Reductions(integrated_scdata),
    "integrated.rpca"
  )

  expect_equal(
    integrated_scdata@misc$active.reduction,
    "integrated.rpca"
  )
})


test_that("SeuratV4 integration finds integration anchors using CCA method, if method in config is CCA", {
  # mock a bigger dataset to run Seurat v4 integration without skipping it
  scdata_list <- mock_scdata_list(n_rep=3)
  cells_id <- list(
    "123abc" = scdata_list$`123abc`$cells_id,
    "123def" = scdata_list$`123def`$cells_id
  )

  config <- list(
    dimensionalityReduction = list(method = "cca"),
    dataIntegration = list(
      method = "seuratv4",
      methodSettings = list(
        seuratv4 = list(numGenes = 1000, normalisation = "logNormalize")
      )
    )
  )

  integrated_scdata <- suppressWarnings(
    integrate_scdata(
      scdata_list, config, "", cells_id, task_name = "dataIntegration"
    )$data
  )
  expect_contains(
    Seurat::Reductions(integrated_scdata),
    "integrated.cca"
  )

  expect_equal(
    integrated_scdata@misc$active.reduction,
    "integrated.cca"
  )
})

test_that("SCTransform integration works", {
  # mock a bigger dataset to run Seurat v4 integration without skipping it
  scdata_list <- mock_scdata_list(n_rep=3)
  cells_id <- list(
    "123abc" = scdata_list$`123abc`$cells_id,
    "123def" = scdata_list$`123def`$cells_id
  )

  config <- list(
    dimensionalityReduction = list(method = "rpca"),
    dataIntegration = list(
      method = "seuratv4",
      methodSettings = list(
        seuratv4 = list(numGenes = 1000, normalisation = "SCT")
      )
    )
  )

  integrated_scdata <- suppressWarnings(
    integrate_scdata(
      scdata_list, config, "", cells_id, task_name = "dataIntegration"
    )$data
  )
  expect_s4_class(integrated_scdata, "Seurat")
  expect_s4_class(integrated_scdata[["SCT"]], "SCTAssay")
})


test_that("misc slot is complete after Seurat V4 integration", {
  # mock a bigger dataset to run Seurat v4 integration without skipping it
  scdata_list <- mock_scdata_list(n_rep=3)
  cells_id <- list(
    "123abc" = scdata_list$`123abc`$cells_id,
    "123def" = scdata_list$`123def`$cells_id
  )

  config <- list(
    dimensionalityReduction = list(method = "rpca"),
    dataIntegration = list(
      method = "seuratv4",
      methodSettings = list(
        seuratv4 = list(numGenes = 1000, normalisation = "logNormalize")
      )
    )
  )

  integrated_scdata <- suppressWarnings(
    integrate_scdata(
      scdata_list, config, "", cells_id, task_name = "dataIntegration"
    )$data
  )
  integrated_scdata <- clean_timestamp(integrated_scdata)
  integrated_scdata <- remove_commands_functions(integrated_scdata)

  expect_s4_class(integrated_scdata, "Seurat")

  skip_if_bpcells()
  skip_on_ci()
  expect_snapshot(str(integrated_scdata@misc))
  expect_snapshot(str(integrated_scdata))
})


test_that("misc slot is complete after Seurat V4 integration with geosketch", {

  # mock a bigger dataset to run Seurat v4 integration without skipping it
  scdata_list <- mock_scdata_list(n_rep=3)
  cells_id <- list(
    "123abc" = scdata_list$`123abc`$cells_id,
    "123def" = scdata_list$`123def`$cells_id
  )

  config <- list(
    dimensionalityReduction = list(numPCs = 10, method = "rpca"),
    dataIntegration = list(
      method = "seuratv4",
      methodSettings = list(
        seuratv4 = list(numGenes = 1000, normalisation = "logNormalize")
      )
    ),
    downsampling = list(
      method = "geosketch",
      methodSettings = list(geosketch = list(percentageToKeep = 50))
    )
  )


  integrated_scdata <- suppressWarnings(
    integrate_scdata(
      scdata_list,config, "", cells_id, task_name = "dataIntegration"
    )$data
  )
  expect_s4_class(integrated_scdata, "Seurat")
  integrated_scdata@misc$ingestionDate <- "fixed_date"

  expected_misc_names <- c(
    "gene_annotations",
    "color_pool",
    "ingestionDate",
    "active.reduction",
    "numPCs",
    "gene_dispersion"
  )

  expect_setequal(names(integrated_scdata@misc), expected_misc_names)
})


test_that("default assay in the integrated object matches normalisation method after Seurat V4 integration with geosketch", {
  # mock a bigger dataset to run Seurat v4 integration without skipping it
  scdata_list <- mock_scdata_list(n_rep=3)
  cells_id <- list(
    "123abc" = scdata_list$`123abc`$cells_id,
    "123def" = scdata_list$`123def`$cells_id
  )

  normalisation_methods <- c("logNormalize", "SCT")

  for (normalisation_method in normalisation_methods) {
    config <- list(
      dimensionalityReduction = list(numPCs = 10, method = "rpca"),
      dataIntegration = list(
        method = "seuratv4",
        methodSettings = list(
          seuratv4 = list(numGenes = 1000, normalisation = normalisation_method)
        )
      ),
      downsampling = list(
        method = "geosketch",
        methodSettings = list(geosketch = list(percentageToKeep = 50))
      )
    )

    integrated_scdata <- suppressWarnings(
      integrate_scdata(
        scdata_list, config, "", cells_id, task_name = "dataIntegration"
      )$data
    )
    expect_s4_class(integrated_scdata, "Seurat")
  }
})