run_unisample <- function(scdata_list, config, cells_id) {
  settings <- config$dataIntegration$methodSettings[[UNISAMPLE]]
  nfeatures <- settings$numGenes
  normalization <- settings$normalisation

  # grep in case misspelled
  if (grepl("lognorm", normalization, ignore.case = TRUE)) normalization <- "LogNormalize"

  # calculate as many PCs for the PCA as possible, ideally 50, unless few cells
  npcs_for_pca <- min(vapply(scdata_list, ncol, integer(1)) - 1, 50)
  npcs <- config$dimensionalityReduction$numPCs

  # use the min nPCs of what the user wants and what can be calculated
  if (!is.null(npcs)) {
    npcs <- min(config$dimensionalityReduction$numPCs, npcs_for_pca)
  }

  scdata <- prepare_scdata_for_unisample(scdata_list, config, cells_id)

  # in unisample we only need to normalize
  scdata <- Seurat::NormalizeData(scdata, normalization.method = normalization, verbose = FALSE) |>
    Seurat::FindVariableFeatures(assay = "RNA", nfeatures = nfeatures, verbose = FALSE) |>
    Seurat::ScaleData(verbose = FALSE)

  scdata <- add_dispersions(scdata, normalization)

  # run PCA with 50 PCs (or less if there are less cells)
  scdata <-
    Seurat::RunPCA(
      scdata,
      npcs = npcs_for_pca,
      features = Seurat::VariableFeatures(object = scdata),
      verbose = FALSE
    )

  # estimate number of PCs to be used downstream, for clustering
  if (is.null(npcs)) {
    scdata@misc[["active.reduction"]] <- "pca"
    npcs <- get_npcs(scdata)
    message("Estimated number of PCs: ", npcs)
  }

  scdata@misc[["active.reduction"]] <- "pca"
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
