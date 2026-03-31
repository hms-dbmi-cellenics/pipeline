#' Run Seurat v4 integration
#'
#' Integrates two or more samples using the Seurat v4 workflow with IntegrateLayers.
#' It takes into account two different normalization methods: LogNormalize and SCTransform.
#' Supports both CCA and RPCA reduction methods for integration.
#' This function can also be used in combination with Geosketch for downsampling.
#'
#' @param scdata_list list of SeuratObjects
#' @param config list of configuration parameters
#' @param cells_id list of cells ids to keep
#'
#' @return normalized and integrated Seurat object
#' @export
#'
run_seuratv4 <- function(scdata_list, config, cells_id) {
  settings <- config$dataIntegration$methodSettings[["seuratv4"]]
  nfeatures <- settings$numGenes
  normalization <- settings$normalisation

  k_filter <- min(ceiling(sapply(scdata_list, ncol) / 2) - 1, 200)
  message("k_filter: ", k_filter)

  if (grepl("lognorm", normalization, ignore.case = TRUE)) {
    normalization <- "LogNormalize"
  }

  reduction <- config$dimensionalityReduction$method

  # calculate as many PCs for the PCA as possible, ideally 50, unless few cells
  npcs_for_pca <- min(vapply(scdata_list, ncol, numeric(1)) - 1, 50)
  npcs <- config$dimensionalityReduction$numPCs

  # use the min of what the user wants and what can be calculated
  if (!is.null(npcs)) {
    npcs <- min(config$dimensionalityReduction$numPCs, npcs_for_pca)
  }

  scdata <- prepare_scdata_for_seuratv4(scdata_list, config, cells_id)

  is_sct <- normalization == "SCT"
  use_geosketch <- is_geosketch(config)
  percent_keep <- config$downsampling$methodSettings$geosketch$percentageToKeep

  # required pre-processing
  sketch_reduction <- paste0("integrated.", reduction, ".sketch")
  full_reduction <- paste0("integrated.", reduction)

  pca_reduction <- ifelse(use_geosketch, "pca.sketch", "pca")
  new_reduction <- ifelse(use_geosketch, sketch_reduction, full_reduction)

  # Determine the integration method based on reduction type
  integration_method <- if (reduction == "rpca") {
    Seurat::RPCAIntegration
  } else if (reduction == "cca") {
    Seurat::CCAIntegration
  } else {
    stop("Unknown reduction method: ", reduction)
  }

  # Preprocessing handles both SCT and LogNormalize normalization methods
  scdata <- scdata |>
    run_preprocessing(
      normalization = normalization,
      nfeatures = nfeatures,
      use_geosketch = use_geosketch,
      percent_keep = percent_keep
    )

  # Scale data (for LogNormalize; SCT is already scaled)
  if (normalization == "LogNormalize") {
    scdata <- Seurat::ScaleData(scdata, verbose = FALSE)
  }

  scdata <- scdata |>
    run_pca(npcs = npcs_for_pca, reduction_name = pca_reduction) |>
    # NOTE: SCTransform calculates variances per sample so HVFInfo fails
    # see https://github.com/satijalab/seurat/issues/6412
    # instead we add dispersions from log-normalized data
    # dispersions are used for gene list in Data Exploration
    add_dispersions(method = "LogNormalize")

  scdata <- Seurat::IntegrateLayers(
    scdata,
    method = integration_method,
    orig = pca_reduction,
    new.reduction = new_reduction,
    dims = 1:npcs,
    normalization.method = normalization,
    verbose = FALSE,
    k.filter = k_filter
  )

  if (use_geosketch) {
    scdata <- project_geosketch_integration(
      scdata,
      npcs,
      sketched_reduction = sketch_reduction,
      full_reduction = full_reduction,
      assay = ifelse(is_sct, "SCT", "RNA")
    )
  }

  scdata@misc[["numPCs"]] <- npcs
  scdata@misc[["active.reduction"]] <- new_reduction

  return(scdata)
}

#' Prepare scdata for Seurat v4 integration
#'
#' Preprocess the scdata list before integration by merging samples,
#' removing excluded genes, and preparing for integration.
#'
#' @inheritParams run_seuratv4
#' @return merged Seurat object ready for integration
#' @export
#'
prepare_scdata_for_seuratv4 <- function(scdata_list, config, cells_id) {
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
