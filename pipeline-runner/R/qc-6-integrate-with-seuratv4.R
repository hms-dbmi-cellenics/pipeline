run_seuratv4 <- function(scdata_list, cells_id, exclude_groups, use_geosketch, npcs,
                         nfeatures, normalization, reduction, perc_num_cells) {
  # @misc slots not preserved so transfer. Mics slots are the same for each sample
  misc <- scdata_list[[1]]@misc

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
      scdata_list[[i]] <- Seurat::NormalizeData(scdata_list[[i]], assay = "RNA", normalization.method = normalization, verbose = FALSE)
      scdata_list[[i]] <- Seurat::FindVariableFeatures(scdata_list[[i]], assay = "RNA", nfeatures = nfeatures, verbose = FALSE)
      scdata_list[[i]] <- Seurat::ScaleData(scdata_list[[i]], verbose = FALSE)
    }

    if (normalization == "SCT") {
      message("Started normalization using SCTransform")
      # conserve.memory parameter reduces the memory footprint but can significantly increase runtime
      scdata_list[[i]] <- Seurat::SCTransform(scdata_list[[i]], vst.flavor = "v2", conserve.memory = FALSE)
      message("Finished normalization using SCTransform")
    }

    # PCA needs to be run also here
    # otherwise when running FindIntegrationAnchors() with reduction="rpca" it will fail because no "pca" is present
    if (reduction == "rpca") {
      message("Running PCA")
      scdata_list[[i]] <- Seurat::RunPCA(scdata_list[[i]], verbose = FALSE, npcs = npcs)
    } else {
      message("PCA is not running before integration as CCA method is selected")
    }
  }

  if (!use_geosketch) {
    # if not using geosketch, just integrate
    scdata <- seuratv4_find_and_integrate_anchors(scdata_list, cells_id, reduction, normalization, npcs, misc, nfeatures)
  } else {
    scdata <- integrate_using_geosketch(scdata_list, cells_id, reduction, perc_num_cells, normalization, npcs, misc, nfeatures, use_geosketch)
  }

  return(scdata)
}



seuratv4_find_and_integrate_anchors <-
  function(data_split, cells_id,
           reduction,
           normalization,
           npcs,
           misc,
           nfeatures,
           scdata = NA,
           use_geosketch = FALSE) {
    k.filter <- min(ceiling(sapply(data_split, ncol) / 2), 200)
    tryCatch(
      {
        if (normalization == "SCT") {
          data_anchors <-
            prepare_sct_integration(data_split, reduction, normalization, k.filter, npcs)
        }
        if (normalization == "LogNormalize") {
          data_anchors <- Seurat::FindIntegrationAnchors(
            object.list = data_split,
            dims = 1:npcs,
            k.filter = k.filter,
            normalization.method = normalization,
            verbose = TRUE,
            reduction = reduction
          )
        }
        scdata <-
          Seurat::IntegrateData(
            anchorset = data_anchors,
            dims = 1:npcs,
            normalization.method = normalization
          )
      },
      error = function(e) {
        # Specifying error message
        # ideally this should be passed to the UI as a error message:
        print(e)
        print(paste("current k.filter:", k.filter))
        # Should we still continue if data is not integrated? No, right now..
        print("Current number of cells per sample: ")
        print(sapply(data_split, ncol))
        warning(
          "Error thrown in IntegrateData: Probably one/many of the samples contain too few cells.\nRule of thumb is that this can happen at around < 100 cells."
        )
        # An ideal solution would be to launch an error to the UI, however, for now, we will skip the integration method.
        print("Skipping integration step")
      }
    )

    if (use_geosketch == FALSE && class(scdata) != "Seurat") {
      message(
        "Merging data because integration was skipped due tue one/many samples containing too few cells"
      )
      scdata <- create_scdata(data_split, cells_id)
    }

    # running LogNormalization on SCTransformed data for downstream analyses
    if (normalization == "SCT") {
      Seurat::DefaultAssay(scdata) <- "RNA"
      scdata <- Seurat::NormalizeData(scdata, normalization.method = "LogNormalize", verbose = FALSE)
    }

    if ("integrated" %in% names(scdata@assays)) {
      Seurat::DefaultAssay(scdata) <- "integrated"
    }

    scdata@misc <- misc
    scdata <- Seurat::FindVariableFeatures(scdata, assay = "RNA", nfeatures = nfeatures, verbose = FALSE)
    scdata <- add_dispersions(scdata, normalization)
    scdata <- Seurat::ScaleData(scdata, verbose = FALSE)

    # run PCA
    scdata <-
      Seurat::RunPCA(
        scdata,
        npcs = npcs,
        features = Seurat::VariableFeatures(object = scdata),
        verbose = FALSE
      )

    scdata@misc[["active.reduction"]] <- "pca"

    return(scdata)
  }


integrate_using_geosketch <-
  function(scdata_list,
           cells_id,
           reduction,
           perc_num_cells,
           normalization,
           npcs,
           misc,
           nfeatures,
           use_geosketch) {
    message("Percentage of cells to keep: ", perc_num_cells)
    # merge
    scdata <- create_scdata(scdata_list, cells_id)
    # geosketch needs PCA to be run
    scdata <- run_pca(scdata)
    # geoesketch
    set.seed(RANDOM_SEED)
    c(scdata, scdata_sketch) %<-% run_geosketch(
      scdata = scdata,
      dims = 50,
      perc_num_cells = perc_num_cells
    )
    # split and integrate sketches
    scdata_sketch_split <- Seurat::SplitObject(scdata_sketch, split.by = "samples")
    scdata_sketch_integrated <- seuratv4_find_and_integrate_anchors(
      scdata_sketch_split, cells_id,
      reduction, normalization,
      npcs, misc, nfeatures, scdata, use_geosketch
    )
    # learn from sketches
    message("Learning from sketches")
    scdata <- learn_from_sketches(
      scdata,
      scdata_sketch,
      scdata_sketch_integrated,
      npcs
    )
    message("Finished learning from sketches")

    return(scdata)
  }
