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

  if (grepl("lognorm", normalization, ignore.case = TRUE)) {
    normalization <- "LogNormalize"
  }

  reduction <- config$dimensionalityReduction$method

  # calculate as many PCs for the PCA as possible, ideally 50, unless few cells
  npcs_for_pca <- min(vapply(scdata_list, ncol, numeric(1)) - 1, 50)
  npcs <- config$dimensionalityReduction$numPCs

  # use the min of what the user wants and what can be calculated
  npcs <- min(config$dimensionalityReduction$numPCs, npcs_for_pca)

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
    rpca_integration
  } else if (reduction == "cca") {
    cca_integration
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

  # adaptive k_filter to avoid errors with small sample sizes
  k_filter <- calculate_k_filter(scdata)

  scdata <- Seurat::IntegrateLayers(
    scdata,
    method = integration_method,
    orig = pca_reduction,
    new.reduction = new_reduction,
    dims = 1:npcs,
    normalization.method = normalization,
    k.filter = k_filter,
    verbose = FALSE
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

calculate_k_filter <- function(scdata) {
  # if sketch, will be kept cells only
  active_cells <- Seurat::Cells(scdata)

  ncells_per_sample <- table(scdata$samples[active_cells])
  k_filter <- min(ceiling(ncells_per_sample / 2) - 1, 200)

  return(k_filter)
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

# reimplementation of Seurat::RPCAIntegration with
# adaptive k.weight calculation to avoid
# errors with small sample sizes and low anchor counts
rpca_integration <- function(
  object = NULL,
  assay = NULL,
  layers = NULL,
  orig = NULL,
  new.reduction = "integrated.dr",
  reference = NULL,
  features = NULL,
  normalization.method = c("LogNormalize", "SCT"),
  dims = 1:30,
  k.filter = NA,
  scale.layer = "scale.data",
  dims.to.integrate = NULL,
  k.weight = 100,
  weight.reduction = NULL,
  sd.weight = 1,
  sample.tree = NULL,
  preserve.order = FALSE,
  verbose = TRUE,
  ...) {

  op <- options(Seurat.object.assay.version = "v3", Seurat.object.assay.calcn = FALSE)
  on.exit(expr = options(op), add = TRUE)
  normalization.method <- match.arg(arg = normalization.method)
  features <- features %||% Seurat:::SelectIntegrationFeatures5(object = object)
  assay <- assay %||% "RNA"
  layers <- layers %||% Seurat::Layers(object = object, search = "data")

  ncells <- sapply(X = layers, FUN = function(x) {
    ncell <- dim(object[x])[2]
    return(ncell)
  })
  if (min(ncells) < max(dims)) {
    abort(message = "At least one layer has fewer cells than dimensions specified, please lower 'dims' accordingly.")
  }
  if (normalization.method == "SCT") {
    groups <- Seurat:::CreateIntegrationGroups(object, layers = layers,
      scale.layer = scale.layer)
    object.sct <- SeuratObject::CreateSeuratObject(counts = object, assay = "SCT")
    object.sct$split <- groups[, 1]
    object.list <- Seurat::SplitObject(object = object.sct, split.by = "split")
    object.list <- Seurat::PrepSCTIntegration(object.list = object.list,
      anchor.features = features)
    object.list <- lapply(X = object.list, FUN = function(x) {
      x <- Seurat::RunPCA(object = x, features = features, verbose = FALSE,
        npcs = max(dims))
      return(x)
    })
  }
  else {
    object.list <- list()
    for (i in seq_along(along.with = layers)) {
      object.list[[i]] <- suppressMessages(suppressWarnings(SeuratObject::CreateSeuratObject(counts = NULL,
        data = object[layers[i]][features, ])))
      Seurat::VariableFeatures(object = object.list[[i]]) <- features
      object.list[[i]] <- suppressWarnings(Seurat::ScaleData(object = object.list[[i]],
        verbose = FALSE))
      object.list[[i]] <- Seurat::RunPCA(object = object.list[[i]],
        verbose = FALSE, npcs = max(dims))
      suppressWarnings(object.list[[i]][["RNA"]]$counts <- NULL)
    }
  }
  anchor <- Seurat::FindIntegrationAnchors(object.list = object.list,
    anchor.features = features, scale = FALSE, reduction = "rpca",
    normalization.method = normalization.method, dims = dims,
    k.filter = k.filter, reference = reference, verbose = verbose,
    ...)

  # Calculate adaptive k.weight based on anchor count and sample count
  # https://github.com/satijalab/seurat/issues/4427#issuecomment-834685413
  nsamples <- length(object.list)
  data_anchors_count <- nrow(anchor@anchors)
  k.weight <- min(floor(data_anchors_count / (3 * nsamples)), k.filter, 100)
  message("k.weight: ", k.weight)

  slot(object = anchor, name = "object.list") <- lapply(X = slot(object = anchor,
    name = "object.list"), FUN = function(x) {
    suppressWarnings(expr = x <- Seurat::DietSeurat(x, features = features[1:2]))
    return(x)
  })
  object_merged <- Seurat::IntegrateEmbeddings(anchorset = anchor,
    reductions = orig, new.reduction.name = new.reduction,
    dims.to.integrate = dims.to.integrate, k.weight = k.weight,
    weight.reduction = weight.reduction, sd.weight = sd.weight,
    sample.tree = sample.tree, preserve.order = preserve.order,
    verbose = verbose)
  output.list <- list(object_merged[[new.reduction]])
  names(output.list) <- c(new.reduction)
  return(output.list)
}

# reimplementation of Seurat::CCAIntegration with
# adaptive k.filter and k.weight calculation to avoid
# errors with small sample sizes and low anchor counts
cca_integration <- function(
  object = NULL,
  assay = NULL,
  layers = NULL,
  orig = NULL,
  new.reduction = "integrated.dr",
  reference = NULL,
  features = NULL,
  normalization.method = c("LogNormalize", "SCT"),
  dims = 1:30,
  k.filter = NA,
  scale.layer = "scale.data",
  dims.to.integrate = NULL,
  k.weight = 100,
  weight.reduction = NULL,
  sd.weight = 1,
  sample.tree = NULL,
  preserve.order = FALSE,
  verbose = TRUE,
  ...) {

  op <- options(Seurat.object.assay.version = "v3", Seurat.object.assay.calcn = FALSE)
  on.exit(expr = options(op), add = TRUE)
  normalization.method <- match.arg(arg = normalization.method)
  features <- features %||% Seurat:::SelectIntegrationFeatures5(object = object)
  assay <- assay %||% "RNA"
  layers <- layers %||% Seurat::Layers(object, search = "data")

  # Calculate adaptive k.filter based on sample sizes
  if (is.na(k.filter)) {
    ncells <- sapply(X = layers, FUN = function(x) {
      ncell <- dim(object[x])[2]
      return(ncell)
    })
    k.filter <- min(ceiling(ncells / 2) - 1, 200)
  }

  if (normalization.method == "SCT") {
    groups <- Seurat:::CreateIntegrationGroups(object, layers = layers,
      scale.layer = scale.layer)
    object.sct <- SeuratObject::CreateSeuratObject(counts = object, assay = "SCT")
    object.sct$split <- groups[, 1]
    object.list <- Seurat::SplitObject(object = object.sct, split.by = "split")
    object.list <- Seurat::PrepSCTIntegration(object.list, anchor.features = features)
  }
  else {
    object.list <- list()
    for (i in seq_along(along.with = layers)) {
      if (inherits(x = object[layers[i]], what = "IterableMatrix")) {
        warning("Converting BPCells matrix to dgCMatrix for integration ",
          "as on-disk CCA Integration is not currently supported",
          call. = FALSE, immediate. = TRUE)
        counts <- as(object = object[layers[i]][features, ],
          Class = "dgCMatrix")
      }
      else {
        counts <- object[layers[i]][features, ]
      }
      object.list[[i]] <- SeuratObject::CreateSeuratObject(counts = counts)
      if (inherits(x = object[scale.layer], what = "IterableMatrix")) {
        scale.data.layer <- as.matrix(object[scale.layer][features,
          Seurat::Cells(object.list[[i]])])
        object.list[[i]][["RNA"]]$scale.data <- scale.data.layer
      }
      else {
        object.list[[i]][["RNA"]]$scale.data <- object[scale.layer][features,
          Seurat::Cells(object.list[[i]])]
      }
      object.list[[i]][["RNA"]]$counts <- NULL
    }
  }
  anchor <- Seurat::FindIntegrationAnchors(object.list = object.list,
    anchor.features = features, scale = FALSE, reduction = "cca",
    normalization.method = normalization.method, dims = dims,
    k.filter = k.filter, reference = reference, verbose = verbose,
    ...)

  # Calculate adaptive k.weight based on anchor count and sample count
  # https://github.com/satijalab/seurat/issues/4427#issuecomment-834685413
  nsamples <- length(object.list)
  data_anchors_count <- nrow(anchor@anchors)
  k.weight <- min(floor(data_anchors_count / (3 * nsamples)), k.filter, 100)
  message("k.weight: ", k.weight)

  suppressWarnings({
    anchor@object.list <- lapply(anchor@object.list, function(x) {
      x <- Seurat::DietSeurat(x, features = features[1:2])
      return(x)
    })
  }, classes = "dimWarning")
  object_merged <- Seurat::IntegrateEmbeddings(anchorset = anchor,
    reductions = orig, new.reduction.name = new.reduction,
    dims.to.integrate = dims.to.integrate, k.weight = k.weight,
    weight.reduction = weight.reduction, sd.weight = sd.weight,
    sample.tree = sample.tree, preserve.order = preserve.order,
    verbose = verbose)
  output.list <- list(object_merged[[new.reduction]])
  names(output.list) <- c(new.reduction)
  return(output.list)
}
