run_harmony <- function(scdata_list, config, cells_id) {
  settings <- config$dataIntegration$methodSettings[["harmony"]]
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

  scdata <- prepare_scdata_for_harmony(scdata_list, config, cells_id)

  # we need RNA assay to compute the integrated matrix
  Seurat::DefaultAssay(scdata) <- "RNA"

  # required pre-processing
  # run PCA with npcs_for_pca for the elbow plot and the % of variance explained
  scdata <- scdata |>
    RunPreprocessing(
      normalization = normalization,
      nfeatures = nfeatures,
      npcs = npcs_for_pca,
      config = config
    )

  use_geosketch <- is_geosketch(config)
  pca.reduction <- ifelse(use_geosketch, "pca.sketch", "pca")
  new.reduction <- ifelse(use_geosketch, "harmony.sketch", "harmony")

  # estimate number of PCs to be used downstream, for integration and clustering
  if (is.null(npcs)) {
    scdata@misc[["active.reduction"]] <- pca.reduction
    npcs <- get_npcs(scdata)
    message("Estimated number of PCs: ", npcs)
  }

  scdata <- Seurat::IntegrateLayers(
    scdata,
    method = Seurat::HarmonyIntegration,
    orig = pca.reduction,
    new.reduction = new.reduction,
    dims = 1:npcs,
    verbose = FALSE
  )

  if (use_geosketch) {
    scdata <- ProjectGeosketchIntegration(
      scdata,
      npcs,
      sketched.reduction = "harmony.sketch",
      full.reduction = "harmony"
    )
  }

  scdata <- add_dispersions(scdata, normalization)
  scdata@misc[["numPCs"]] <- npcs
  scdata@misc[["active.reduction"]] <- new.reduction

  return(scdata)
}

ProjectGeosketchIntegration <- function(scdata, npcs, sketched.reduction, full.reduction, is.unisample = FALSE) {

  if (!is.unisample) {
    # project full data onto integrated DR from sketch
    # TODO: understand what this does?
    message("Projecting integration...")
    scdata <- Seurat::ProjectIntegration(
      object = scdata,
      assay = "RNA",
      sketched.assay = "sketch",
      reduction = sketched.reduction,
      reduction.name = full.reduction
    )
  }

  # for unisample: projects full PCA
  # for all: called to pre-compute weights and neighbors (stored in tools)
  # for future projections (clustering and UMAP)
  # NOTE: subsequent calls to ProjectData require same npcs as here
  scdata <- Seurat::ProjectData(
    object = scdata,
    assay = "RNA",
    sketched.assay = "sketch",
    full.reduction = full.reduction,
    sketched.reduction = sketched.reduction,
    dims = 1:npcs
  )

  # set to RNA assay prior to upload for worker
  Seurat::DefaultAssay(scdata) <- "RNA"
  return(scdata)
}

is_geosketch <- function(config) {
  "downsampling" %in% names(config) && config$downsampling$method == "geosketch"
}

RunPreprocessing <- function(scdata, normalization, nfeatures, npcs, config) {

  # SketchData takes normalized data and set of variable features
  scdata <-
    Seurat::NormalizeData(
      scdata,
      normalization.method = normalization,
      verbose = FALSE
    )

  # split assay into sample layers
  num_samples <- length(unique(scdata$samples))
  if (num_samples > 1) {
    scdata[["RNA"]] <- split(scdata[["RNA"]], f = scdata$samples)
  }

  scdata <- Seurat::FindVariableFeatures(scdata, nfeatures = nfeatures, verbose = FALSE)

  message("Is this GEOSKETCH? ", is_geosketch(config))

  if (is_geosketch(config)) {
    perc_num_cells <- config$downsampling$methodSettings$geosketch$percentageToKeep
    num_cells <- round(ncol(scdata) * (perc_num_cells / 100) / num_samples)
    message("Sketching data with Geosketch, keeping ", num_cells, " cells per sample.")
    scdata <- Seurat::SketchData(
      object = scdata,
      ncells = num_cells,
      method = "LeverageScore",
      sketched.assay = "sketch"
    )

    # vignettes run FindVariableFeatures before and after SketchData
    Seurat::DefaultAssay(scdata) <- "sketch"
    scdata <- Seurat::FindVariableFeatures(scdata, nfeatures = nfeatures, verbose = FALSE)
  }

  reduction.name <- ifelse(is_geosketch(config), "pca.sketch", "pca")

  scdata <- scdata |>
    Seurat::ScaleData(verbose = FALSE) |>
    RunPCA(npcs = npcs, reduction.name = reduction.name)

  return(scdata)
}

RunPCA <- function(scdata, npcs, reduction.name = "pca") {

  # if BPCells write out operations for faster PCA
  is.bpcells <- is(scdata[['RNA']]$scale.data, "IterableMatrix")

  if (is.bpcells) {
    data_dir <- tempfile()
    scale.data <- scdata[['RNA']]$scale.data
    disk.scale.data <- BPCells::write_matrix_dir(scale.data, data_dir)
    scdata[['RNA']]$scale.data <- disk.scale.data
  }

  # run PCA
  message('Running PCA ...')
  scdata <- Seurat::RunPCA(
    scdata,
    npcs = npcs,
    reduction.name = reduction.name,
    verbose = FALSE
  )
  message('PCA complete.')

  if (is.bpcells) {
    # cleanup restore scale.data to original object (removed disk backing)
    unlink(data_dir, recursive = TRUE)
    scdata[['RNA']]$scale.data <- scale.data
  }

  return(scdata)
}


prepare_scdata_for_harmony <- function(scdata_list, config, cells_id) {
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


RunGeosketchHarmony <- function(scdata,
                                group.by.vars,
                                reduction,
                                npcs,
                                config) {
  set.seed(RANDOM_SEED)
  perc_num_cells <-
    config$downsampling$methodSettings$geosketch$percentageToKeep
  geosketch_list <- run_geosketch(
    scdata,
    dims = npcs,
    perc_num_cells = perc_num_cells,
    reduction = reduction
  )

  scdata_sketch_integrated <-
    harmony::RunHarmony(
      geosketch_list$sketch,
      group.by.vars = "samples",
      reduction.use = reduction,
      dims.use = 1:npcs,
    )
  scdata_sketch_integrated@misc[["active.reduction"]] <- "harmony"

  scdata <- learn_from_sketches(
    geosketch_list$scdata,
    geosketch_list$sketch,
    scdata_sketch_integrated,
    dims = npcs,
    reduction = reduction
  )

  return(scdata)
}
