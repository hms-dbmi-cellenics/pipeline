run_harmony <- function(scdata, config, use_geosketch) {
  settings <- config$dataIntegration$methodSettings[["harmony"]]
  nfeatures <- settings$numGenes
  normalization <- settings$normalisation
  npcs <- config$dimensionalityReduction$numPCs

  # required pre-processing
  scdata <- scdata |>
    Seurat::NormalizeData(normalization.method = normalization, verbose = FALSE) |>
    Seurat::FindVariableFeatures(nfeatures = nfeatures, verbose = FALSE) |>
    Seurat::ScaleData(verbose = FALSE) |>
    Seurat::RunPCA(scdata, npcs = 50, verbose = FALSE)

  # estimate number of PCs to be used downstream after running PCA
  if (is.null(npcs)) {
    npcs <- get_npcs(scdata)
    message("Estimated number of PCs: ", npcs)
  }

  # harmony ignores dims.use. It uses all dims available
  # https://github.com/immunogenomics/harmony/issues/151
  scdata <- Seurat::RunPCA(scdata, npcs = npcs, verbose = FALSE, reduction.name = "pca_for_harmony")

  # downsample
  # main function
  set.seed(RANDOM_SEED)
  if (use_geosketch) {
    geosketch_list <- run_geosketch(
      scdata,
      dims = npcs,
      perc_num_cells = perc_num_cells
    )

    scdata_sketch_integrated <-
      harmony::RunHarmony(
        geosketch_list$scdata,
        group.by.vars = "samples",
        reduction = "pca_for_harmony" ,
        dims.use = npcs
      )

    scdata <- learn_from_sketches(
      geosketch_list$scdata,
      geosketch_list$scdata_sketch,
      scdata_sketch_integrated,
      method,
      npcs
    )
  } else {
    scdata <-
      harmony::RunHarmony(
        scdata,
        group.by.vars = "samples",
        reduction = "pca_for_harmony",
        dims.use = npcs
      )
  }

  scdata <- add_dispersions(scdata, normalization)
  scdata@misc[["active.reduction"]] <- "harmony"

  return(scdata)
}


run_geosketch <- function(scdata, dims, perc_num_cells) {

  reduction <- "pca"
  num_cells <- round(ncol(scdata) * perc_num_cells / 100)

  message("Geosketching to ", num_cells, " cells")

  if (!exists("geosketch")) {
    geosketch <- reticulate::import("geosketch")
  }
  stopifnot(
    "The requested reduction is not present in the Seurat object." = reduction %in% names(scdata@reductions),
    "The number of cells is lower that the number of dimensions." = ncol(scdata@reductions[[reduction]]) >= dims
  )

  embeddings <- scdata@reductions[[reduction]]@cell.embeddings[, 1:dims]
  index <- unlist(geosketch$gs(embeddings, as.integer(num_cells), one_indexed = TRUE))
  sketch <- scdata[, index]
  Seurat::DefaultAssay(sketch) <- "RNA"
  sketch@misc[["active.reduction"]] <- reduction

  return(list(scdata, sketch))
}
