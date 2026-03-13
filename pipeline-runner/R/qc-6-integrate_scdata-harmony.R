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

  use_geosketch <- "downsampling" %in% names(config) && config$downsampling$method == "geosketch"

  scdata <- prepare_scdata_for_harmony(scdata_list, config, cells_id)

  # we need RNA assay to compute the integrated matrix
  Seurat::DefaultAssay(scdata) <- "RNA"

  # required pre-processing
  # run PCA with npcs_for_pca for the elbow plot and the % of variance explained
  scdata <- scdata |>
    Seurat::NormalizeData(normalization.method = normalization, verbose = FALSE) |>
    Seurat::FindVariableFeatures(nfeatures = nfeatures, verbose = FALSE) |>
    Seurat::ScaleData(verbose = FALSE) |>
    run_pca(npcs = npcs_for_pca, verbose = FALSE)


  # estimate number of PCs to be used downstream, for integration and clustering
  if (is.null(npcs)) {
    scdata@misc[["active.reduction"]] <- "pca"
    npcs <- get_npcs(scdata)
    message("Estimated number of PCs: ", npcs)
  }

  # downsample
  # main function
  if (use_geosketch) {
    scdata <-
      RunGeosketchHarmony(
        scdata,
        group.by.vars = "samples",
        npcs = npcs,
        config
      )
  } else {
    scdata <-
      harmony::RunHarmony(
        scdata,
        group.by.vars = "samples",
        dims.use = 1:npcs
      )
  }

  scdata <- add_dispersions(scdata, normalization)
  scdata@misc[["numPCs"]] <- npcs
  scdata@misc[["active.reduction"]] <- "harmony"

  return(scdata)
}

run_pca <- function(scdata, npcs) {
  scale.data <- scdata[['RNA']]$scale.data
  scale.data_dir <- file.path(tempdir(), 'scale.data_bpcells')

  # if BPCells write out operations for faster PCA
  if (is(scale.data, "IterableMatrix")) {

    BPCells::write_matrix_dir(scale.data, scale.data_dir)
    scdata[['RNA']]$scale.data <- BPCells::open_matrix_dir(scale.data_dir)
  }

  # run PCA
  scdata <- Seurat::RunPCA(
    scdata,
    npcs = npcs_for_pca,
    verbose = FALSE
  )

  # cleanup disk scale.data
  unlink(scale.data_dir, recursive = TRUE)

  # restore scale.data to BPCells operations
  scdata[['RNA']]$scale.data <- scale.data
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
