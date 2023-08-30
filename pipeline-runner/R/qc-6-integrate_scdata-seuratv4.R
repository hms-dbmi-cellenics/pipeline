#' Run Seurat v4 integration
#'
#' Integrates two or more samples using the Seurat v4 workflow.
#' It takes into account two different normaization methods: LogNormalize and SCTransform.
#' This function can also be used in combination with Geosketch. In this case,
#' the samples are first merged, then the whole dataset is downsampled using geometric sketching,
#' and the sketches are integrated. The integrated sketches are then used to learn
#' the integration transformation and apply it to the whole dataset.
#'
#' @param scdata_list list of SeuratObjects
#' @param config list of configuration parameters
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
  exclude_groups <- config$dimensionalityReduction$excludeGeneCategories

  use_geosketch <- "downsampling" %in% names(config) && config$downsampling$method == "geosketch"

  # calculate as many PCs for the PCA as possible, ideally 50, unless few cells
  npcs_for_pca <- min(vapply(scdata_list, ncol, integer(1)) - 1, 50)
  # use the min of what the user wants and what can be calculated
  npcs <- min(config$dimensionalityReduction$numPCs, npcs_for_pca)

  scdata_list <- order_by_size(scdata_list)

  # normalize single samples
  for (i in 1:length(scdata_list)) {
    # we need RNA assay to compute the integrated matrix
    Seurat::DefaultAssay(scdata_list[[i]]) <- "RNA"

    # remove cell cycle genes if needed
    if (length(exclude_groups) > 0) {
      scdata_list[[i]] <- remove_genes(scdata_list[[i]], exclude_groups)
    }

    if (normalization == "LogNormalize") {
      scdata_list[[i]] <- scdata_list[[i]] |>
        Seurat::NormalizeData(assay = "RNA", normalization.method = normalization, verbose = FALSE) |>
        Seurat::FindVariableFeatures(assay = "RNA", nfeatures = nfeatures, verbose = FALSE) |>
        Seurat::ScaleData(verbose = FALSE)
    } else if (normalization == "SCT") {
      message("Started normalization using SCTransform")
      # conserve.memory parameter reduces the memory footprint but can significantly increase runtime
      scdata_list[[i]] <- Seurat::SCTransform(
        scdata_list[[i]],
        vst.flavor = "v2",
        conserve.memory = FALSE,
        variable.features.n = nfeatures
      )
    } else {
      stop("No normalization method provided")
    }

    # PCA needs to be run also here
    # otherwise when running FindIntegrationAnchors() with reduction="rpca" it will fail because no "pca" is present
    if (reduction == "rpca") {
      message("Running PCA")
      scdata_list[[i]] <- Seurat::RunPCA(scdata_list[[i]], verbose = FALSE, npcs = npcs)
    } else if (reduction == "cca") {
      message("PCA is not running before integration as CCA method is selected")
    } else {
      stop("No reduction method provided")
    }
  }

  if (!use_geosketch) {
    # if not using geosketch, just integrate
    scdata <- seuratv4_find_and_integrate_anchors(scdata_list, cells_id, reduction, normalization, npcs, npcs_for_pca, nfeatures)
  } else if (use_geosketch) {
    perc_num_cells <- config$downsampling$methodSettings$geosketch$percentageToKeep
    scdata <- seuratv4_geosketch_find_and_integrate_anchors(scdata_list, cells_id, reduction, perc_num_cells, normalization, npcs, npcs_for_pca, nfeatures, use_geosketch)
  }

  scdata@misc[["numPCs"]] <- npcs
  return(scdata)
}



#' Find and integrate anchors
#'
#' This function find and integrate anchors according to the Seurat v4 workflow.
#'
#' @param cells_id list of cells ids to keep
#' @param reduction reduction method
#' @param normalization normalization method
#' @param npcs numer of principal components
#' @param nfeatures number of features
#' @param scdata merged Seurat object
#' @param use_geosketch boolean indicating if geosketch has to be run
#' @param scdata_list list of Seurat objects
#'
#' @return integrated Seurat object
#' @export
#'
seuratv4_find_and_integrate_anchors <-
  function(scdata_list,
           cells_id,
           reduction,
           normalization,
           npcs,
           npcs_for_pca,
           nfeatures,
           scdata = NA,
           use_geosketch = FALSE) {
    k.filter <- min(ceiling(sapply(scdata_list, ncol) / 2) - 1, 200)
    message("k.filter: ", k.filter)
    k.weight <- NA
    tryCatch(
      {
        if (normalization == "SCT") {
          data_anchors <-
            prepare_sct_integration(scdata_list, reduction, normalization, k.filter, npcs)
        } else if (normalization == "LogNormalize") {
          data_anchors <- Seurat::FindIntegrationAnchors(
            object.list = scdata_list,
            dims = 1:npcs,
            k.filter = k.filter,
            normalization.method = normalization,
            verbose = TRUE,
            reduction = reduction
          )
        }
        # correct k.weight for extremely small samples
        # https://github.com/satijalab/seurat/issues/4427#issuecomment-834685413
        nsamples <- length(scdata_list)
        k.weight <- min(floor(nrow(data_anchors@anchors) / (3 * nsamples)), k.filter, 100)
        message("k.weight: ", k.weight)

        scdata <-
          Seurat::IntegrateData(
            anchorset = data_anchors,
            dims = 1:npcs,
            normalization.method = normalization,
            k.weight = k.weight
          )
      },
      error = function(e) {
        # Specifying error message
        # ideally this should be passed to the UI as a error message:
        print(e)
        print(paste("current k.filter:", k.filter))
        print(paste("current k.weight:", k.weight))
        # Should we still continue if data is not integrated? No, right now..
        print("Current number of cells per sample: ")
        print(sapply(scdata_list, ncol))
        warning(
          "Error thrown in IntegrateData: Probably one/many of the samples contain too few cells.\nRule of thumb is that this can happen at around < 100 cells."
        )
        # An ideal solution would be to launch an error to the UI, however, for now, we will skip the integration method.
        print("Skipping integration step")
      }
    )

    if (use_geosketch == FALSE && class(scdata) != "Seurat") {
      message(
        "Merging data because integration was skipped due to one/many samples containing too few cells"
      )
      scdata <- create_scdata(scdata_list, cells_id, merge_data = TRUE)
    }

    # running LogNormalization on SCTransformed data for downstream analyses
    if (normalization == "SCT") {
      Seurat::DefaultAssay(scdata) <- "RNA"
      scdata <- Seurat::NormalizeData(scdata, normalization.method = "LogNormalize", verbose = FALSE)
    }

    if ("integrated" %in% names(scdata@assays)) {
      Seurat::DefaultAssay(scdata) <- "integrated"
    }

    scdata <- Seurat::FindVariableFeatures(scdata, assay = "RNA", nfeatures = nfeatures, verbose = FALSE)
    scdata <- Seurat::ScaleData(scdata, verbose = FALSE)

    # run PCA with 50 PCs (or less if there are less cells)
    scdata <-
      Seurat::RunPCA(
        scdata,
        npcs = npcs_for_pca,
        features = Seurat::VariableFeatures(object = scdata),
        verbose = FALSE
      )
    scdata <- add_metadata(scdata, scdata_list)
    scdata <- add_dispersions(scdata, normalization)
    scdata@misc[["active.reduction"]] <- "pca"

    return(scdata)
  }


#' Downsample and integrate sketches
#'
#' This function uses Geosketch to downsample the dataset, integrates the resulting
#' sketches using the Seurat v4 workflow, then learns
#' the integration transformation and applies it to the whole dataset.
#'
#' @param scdata_list list of Seurat objects
#' @param cells_id list of cells ids to keep
#' @param reduction reduction method
#' @param perc_num_cells percentage of cells to keep when using geosketch
#' @param normalization normalization method
#' @param npcs number of principal components
#' @param nfeatures number of features
#' @param use_geosketch boolean indicating if geosketch has to be run
#' @param npcs_for_pca max number of PCs
#'
#' @return integrated Seurat object
#' @export
#'
seuratv4_geosketch_find_and_integrate_anchors <-
  function(scdata_list,
           cells_id,
           reduction,
           perc_num_cells,
           normalization,
           npcs,
           npcs_for_pca,
           nfeatures,
           use_geosketch) {
    message("Percentage of cells to keep: ", perc_num_cells)
    # merge
    scdata <- create_scdata(scdata_list, cells_id, merge_data = TRUE)
    # geosketch needs PCA to be run
    scdata <- scdata |>
      Seurat::FindVariableFeatures(assay = "RNA", nfeatures = 2000, verbose = FALSE) |>
      Seurat::ScaleData(verbose = FALSE) |>
      Seurat::RunPCA(npcs = npcs, verbose = FALSE)

    scdata@misc[["active.reduction"]] <- "pca"
    # geoesketch
    set.seed(RANDOM_SEED)
    geosketch_list <- run_geosketch(
      scdata = scdata,
      dims = 50,
      perc_num_cells = perc_num_cells
    )
    # split and integrate sketches
    scdata_sketch_split <- Seurat::SplitObject(geosketch_list$sketch, split.by = "samples")
    scdata_sketch_integrated <- seuratv4_find_and_integrate_anchors(
      scdata_sketch_split, cells_id,
      reduction, normalization,
      npcs, npcs_for_pca, nfeatures, scdata, use_geosketch
    )
    # learn from sketches
    message("Learning from sketches")
    scdata <- learn_from_sketches(
      geosketch_list$scdata,
      geosketch_list$sketch,
      scdata_sketch_integrated,
      npcs
    )

    # runPCA with 50 PCs for the elbow plot
    scdata <-
      Seurat::RunPCA(
        scdata,
        npcs = npcs_for_pca,
        features = Seurat::VariableFeatures(object = scdata),
        verbose = FALSE
      )

    return(scdata)
  }


#' Prepare for integration after SCTransform
#'
#' This function runs the steps required to prepare the list of Seurat object normalized with
#' SCTransform for integration, and finds the integration anchors.
#' For further details see the documentation for
#' \code{\link[Seurat:SelectIntegrationFeatures]{Seurat::SelectIntegrationFeatures()}},
#' \code{\link[Seurat:PrepSCTIntegration]{Seurat::PrepSCTIntegration()}},
#' and [sctransform_v2 vignette](https://satijalab.org/seurat/articles/sctransform_v2_vignette.html#perform-integration-using-pearson-residuals-1).
#'
#' @param data.split list of Seurat objects
#' @param reduction reduction method
#' @param normalization normalization method
#' @param k.filter number of neighbors (k) to use when filtering anchors
#' @param npcs number of PCs
#'
#' @return data.anchors to use for integration
#' @export
#'
prepare_sct_integration <- function(data.split, reduction, normalization, k.filter, npcs) {
  features <- Seurat::SelectIntegrationFeatures(object.list = data.split, nfeatures = 3000)
  data.split <- Seurat::PrepSCTIntegration(
    object.list = data.split,
    assay = "SCT",
    anchor.features = features
  )
  data.anchors <- Seurat::FindIntegrationAnchors(
    object.list = data.split,
    dims = 1:npcs,
    k.filter = k.filter,
    verbose = TRUE,
    reduction = reduction,
    normalization.method = normalization,
    anchor.features = features
  )
  return(data.anchors)
}
