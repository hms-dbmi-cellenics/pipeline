run_unisample <- function(scdata_list, config, cells_id) {
  settings <- config$dataIntegration$methodSettings[[UNISAMPLE]]
  nfeatures <- settings$numGenes
  normalization <- settings$normalisation

  # grep in case misspelled
  if (grepl("lognorm", normalization, ignore.case = TRUE)) normalization <- "LogNormalize"

  # calculate as many PCs for the PCA as possible, ideally 50, unless few cells
  npcs_for_pca <- min(vapply(scdata_list, ncol, numeric(1)) - 1, 50)
  npcs <- config$dimensionalityReduction$numPCs

  # use the min nPCs of what the user wants and what can be calculated
  if (!is.null(npcs)) {
    npcs <- min(config$dimensionalityReduction$numPCs, npcs_for_pca)
  }

  scdata <- prepare_scdata_for_unisample(scdata_list, config, cells_id)

  # required pre-processing
  # run PCA with npcs_for_pca for the elbow plot and the % of variance explained
  use_geosketch <- is_geosketch(config)
  percent_keep <- config$downsampling$methodSettings$geosketch$percentageToKeep

  pca_reduction <- ifelse(use_geosketch, "pca.sketch", "pca")

  # in unisample we only need to run preprocessing (no integration)
  scdata <- scdata |>
    run_preprocessing(
      normalization = normalization,
      nfeatures = nfeatures,
      use_geosketch = use_geosketch,
      percent_keep = percent_keep
    ) |>
    Seurat::ScaleData(verbose = FALSE) |>
    run_pca(npcs = npcs_for_pca, reduction.name = pca_reduction) |>
    add_dispersions(normalization)


  # estimate number of PCs to be used downstream, for integration and clustering
  if (is.null(npcs)) {
    scdata@misc[["active.reduction"]] <- pca_reduction
    npcs <- get_npcs(scdata)
    message("Estimated number of PCs: ", npcs)
  }

  if (use_geosketch) {
    message("Number of PCs for geosketch: ", npcs)
    scdata <- project_geosketch_integration(
      scdata,
      npcs,
      sketched.reduction = "pca.sketch",
      full.reduction = "pca",
      is.unisample = TRUE
    )
  }

  # NOTE: if sketch, downstream steps use sketch reduction
  scdata@misc[["active.reduction"]] <- pca_reduction
  scdata@misc[["numPCs"]] <- npcs

  return(scdata)
}


prepare_scdata_for_unisample <- function(scdata_list, config, cells_id) {
  scdata <- create_scdata(scdata_list, cells_id)

  exclude_groups <- config$dimensionalityReduction$excludeGeneCategories
  # remove genes groups if required
  if (length(exclude_groups) > 0) {
    scdata <- remove_genes(scdata, exclude_groups)
  }

  return(scdata)
}
