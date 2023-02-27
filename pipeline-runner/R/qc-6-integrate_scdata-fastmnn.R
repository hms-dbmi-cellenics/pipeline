run_fastmnn <- function(scdata_list, config, cells_id) {
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

  use_geosketch <- "downsampling" %in% names(config) && config$downsampling$method == "geosketch"

  scdata <- prepare_scdata_for_fastmnn(scdata_list, config, cells_id)

  # we need RNA assay to compute the integrated matrix
  Seurat::DefaultAssay(scdata) <- "RNA"

  scdata <- scdata |>
    Seurat::NormalizeData(normalization.method = normalization, verbose = FALSE) |>
    Seurat::FindVariableFeatures(nfeatures = nfeatures, verbose = FALSE) |>
    Seurat::ScaleData(verbose = FALSE) |>
    Seurat::RunPCA(npcs = npcs_for_pca, verbose = FALSE)

  # estimate number of PCs to be used downstream, for integration and clustering
  if (is.null(npcs)) {
    scdata@misc[["active.reduction"]] <- "pca"
    npcs <- get_npcs(scdata)
    message("Estimated number of PCs: ", npcs)
  }

  scdata <- add_dispersions(scdata, normalization)
  misc <- scdata@misc
  # remove scale.data slot to avoid errors
  # https://github.com/satijalab/seurat-wrappers/issues/126
  scdata@assays$RNA@scale.data <- as.matrix(0)

  if (use_geosketch) {
    scdata <-
      RunGeosketchFastMNN(
        scdata,
        split.by = "samples",
        features = nfeatures,
        dims = npcs,
        config = config
      )
    misc[["geosketch"]] <- TRUE
  } else {
    scdata_split <- Seurat::SplitObject(scdata, split.by = "samples")
    scdata <-
      SeuratWrappers::RunFastMNN(
        scdata_split,
        features = nfeatures,
        d = npcs,
        get.variance = TRUE
      )
  }

  scdata@misc <- misc
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


RunGeosketchFastMNN <- function(scdata, split.by, features, dims, config) {
  set.seed(RANDOM_SEED)
  perc_num_cells <- config$downsampling$methodSettings$geosketch$percentageToKeep
  geosketch_list <- run_geosketch(
    scdata,
    dims = dims,
    perc_num_cells = perc_num_cells
    )

  scdata_split <- Seurat::SplitObject(geosketch_list$sketch, split.by = split.by)
  scdata_sketch_integrated <-
    SeuratWrappers::RunFastMNN(
      scdata_split,
      features = nfeatures,
      d = dims,
      get.variance = TRUE
    )

  scdata <- learn_from_sketches(
    geosketch_list$scdata,
    geosketch_list$sketch,
    scdata_sketch_integrated,
    dims = dims
  )

  return(scdata)

}
