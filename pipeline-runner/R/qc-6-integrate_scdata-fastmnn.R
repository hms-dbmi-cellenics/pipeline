# run_fastmnn <- function(scdata, config, npcs) {
#   settings <- config$dataIntegration$methodSettings[["fastmnn"]]
#
#   nfeatures <- settings$numGenes
#   normalization <- settings$normalisation
#
#   # grep in case misspelled
#   if (grepl("lognorm", normalization, ignore.case = TRUE)) normalization <- "LogNormalize"
#
#
#   scdata <- log_normalize(scdata, normalization, "fastmnn", nfeatures)
#   scdata <- add_dispersions(scdata, normalization)
#
#   # @misc slots not preserved so transfer
#   misc <- scdata@misc
#   scdata <- SeuratWrappers::RunFastMNN(object.list = Seurat::SplitObject(scdata, split.by = "samples"), d = 50, get.variance = TRUE)
#   scdata@misc <- misc
#   scdata@misc[["active.reduction"]] <- "mnn"
#
#   return(scdata)
# }


run_fastmnn(scdata_list, config, cells_id) {
  settings <- config$dataIntegration$methodSettings[["fastmnn"]]
  nfeatures <- settings$numGenes
  normalization <- settings$normalisation
  if (grepl("lognorm", normalization, ignore.case = TRUE)) {
    normalization <- "LogNormalize"
  }

  # calculate as many PCs for the PCA as possible, ideally 50, unless few cells
  npcs_for_pca <- min(vapply(scdata_list, ncol, integer(1)) - 1, 50)
  npcs <- config$dimensionalityReduction$numPCs

  # use the min of what the user wants and what can be calculated
  if (!is.null(npcs)) {
    npcs <- min(config$dimensionalityReduction$numPCs, npcs_for_pca)
  }

  scdata <- prepare_scdata_for_fastmnn(scdata_list, config, cells_id)

  # we need RNA assay to compute the integrated matrix
  Seurat::DefaultAssay(scdata) <- "RNA"

  scdata <- scdata |>
    Seurat::NormalizeData(assay = "RNA", normalization.method = normalization_method, verbose = FALSE) |>
    Seurat::FindVariableFeatures(assay = "RNA", nfeatures = nfeatures, verbose = FALSE) |>
    Seurat::ScaleData(verbose = FALSE) |>
    Seurat::RunPCA(npcs = npcs_for_pca, verbose = FALSE)

  # estimate number of PCs to be used downstream, for integration and clustering
  if (is.null(npcs)) {
    scdata@misc[["active.reduction"]] <- "pca"
    npcs <- get_npcs(scdata)
    message("Estimated number of PCs: ", npcs)
  }

  scdata_split <- Seurat::SplitObject(scdata, split.by = "samples")

  scdata <- SeuratWrappers::RunFastMNN(scdata_split, features = nfeatures, d = npcs, get.variance = TRUE)

  scdata <- add_metadata(scdata, scdata_list)
  scdata <- add_dispersions(scdata, normalization)
  scdata@misc[["numPCs"]] <- npcs
  scdata@misc[["active.reduction"]] <- "mnn"

  return(scdata)
}


prepare_scdata_for_fastmnn <- function(scdata_list, config, cells_id) {
  # pre-process
  scdata_list <- order_by_size(scdata_list)
  scdata <- create_scdata(scdata_list, cells_id)

  exclude_groups <- config$dimensionalityReduction$excludeGeneCategories
  # remove genes groups if required
  if (length(exclude_groups) > 0) {
    scdata <- remove_genes(scdata, exclude_groups)
  }

  return(scdata)
}
