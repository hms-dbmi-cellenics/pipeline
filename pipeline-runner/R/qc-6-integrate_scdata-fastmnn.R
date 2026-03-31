run_fastmnn <- function(scdata_list, config, cells_id) {
  settings <- config$dataIntegration$methodSettings[["fastmnn"]]
  nfeatures <- settings$numGenes
  normalization <- settings$normalisation
  if (grepl("lognorm", normalization, ignore.case = TRUE)) {
    normalization <- "LogNormalize"
  }

  # calculate as many PCs for the PCA as possible, ideally 50, unless few cells
  npcs_for_pca <- min(vapply(scdata_list, ncol, numeric(1)) - 1, 50)
  npcs <- config$dimensionalityReduction$numPCs

  # use the min of what the user wants and what can be calculated
  if (!is.null(npcs)) {
    npcs <- min(config$dimensionalityReduction$numPCs, npcs_for_pca)
  }

  scdata <- prepare_scdata_for_fastmnn(scdata_list, config, cells_id)

  # required pre-processing
  # run PCA with npcs_for_pca for the elbow plot and the % of variance explained
  use_geosketch <- is_geosketch(config)
  percent_keep <- config$downsampling$methodSettings$geosketch$percentageToKeep

  pca_reduction <- ifelse(use_geosketch, "pca.sketch", "pca")
  new_reduction <- ifelse(use_geosketch, "mnn.sketch", "mnn")

  # we skip scaling and pca as it's computed by RunFastMNN
  scdata <- scdata |>
    run_preprocessing(
      normalization = normalization,
      nfeatures = nfeatures,
      use_geosketch = use_geosketch,
      percent_keep = percent_keep
    ) |>
    add_dispersions(normalization)


  # NOTE that IntegrateLayers(method = FastMNNIntegration)
  # doesn't work with the current version of SeuratWrappers
  # once it's fixed, we can use it instead
  scdata <- fastmnn_integration(scdata, nfeatures, npcs, new_reduction)

  if (use_geosketch) {
    scdata <- project_geosketch_integration(
      scdata,
      npcs,
      sketched_reduction = "mnn.sketch",
      full_reduction = "mnn"
    )
  }

  scdata@misc[["numPCs"]] <- npcs
  scdata@misc[["active.reduction"]] <- new_reduction

  return(scdata)
}

fastmnn_integration <- function(scdata, nfeatures, npcs, new_reduction) {
  # NOTE: RunFastMNN doesn't preserve misc
  # so we preserve original and add reduction at end
  scdata_orig <- scdata

  # NOTE: RunFastMNN fails if default assay is "sketch"
  # need to subset to sketch cells and remove RNA assay
  if (Seurat::DefaultAssay(scdata) == "sketch") {
    scdata <- scdata[, Seurat::Cells(scdata)]
    scdata[["RNA"]] <- NULL
    message("Running FastMNN on sketch: ", ncol(scdata), " cells")
  }

  scdata_split <- Seurat::SplitObject(scdata, split.by = "samples")

  scdata <-
    SeuratWrappers::RunFastMNN(
      scdata_split,
      features = nfeatures,
      d = npcs,
      get.variance = TRUE
    )

   # add MNN reduction
   scdata_orig[[new_reduction]] <- scdata[["mnn"]]
   
   # add tool with pca info (used by get_explained_variance)
   slot <- "SeuratWrappers::RunFastMNN"
   scdata_orig@tools[[slot]] <- scdata@tools[[slot]]

  return(scdata_orig)
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
