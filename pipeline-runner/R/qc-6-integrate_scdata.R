# STEP 6. DATA INTEGRATION

# Data integration step where batch effect is corrected through data integration methods such as "Seurat V3".
# The data integration step include the normalization and the PCA analysis.

#   "dataIntegration": {
#       "dataIntegration": {
#           "method": "seuratv4",
#           "methodSettings": {
#               "seuratv4": {
#                   "numGenes": 2000,
#                   "normalization": "logNormalise"
#               }
#           }
#       },
#       "dimensionalityReduction": {
#           "method": "rpca",  --> NOT USE!
#           "numPCs": 30,
#           "excludeGeneCategories": []
#       }
#   },

integrate_scdata <- function(scdata_list, config, sample_id, cells_id, task_name = "dataIntegration") {
  # merge the data before integration
  scdata <- merge_scdatas(scdata_list)
  # the sample already contain the metadata, check if the merged object already has the annotations
  # otherwise try to add it from one of the samples
  #ascdata <- add_metadata(scdata, scdata_list$annot, input$experimentId)


  flat_cells_id <- unname(unlist(cells_id))
  scdata <- subset_ids(scdata, flat_cells_id)
  # main function
  set.seed(42)
  scdata.integrated <- run_dataIntegration(scdata, config)

  # get  npcs from the UMAP call in integration functions
  npcs <- length(scdata.integrated@commands$RunUMAP@params$dims)
  message("\nSet config numPCs to npcs used in last UMAP call: ", npcs, "\n")
  config$dimensionalityReduction$numPCs <- npcs

  var_explained <- get_explained_variance(scdata.integrated)

  # This same numPCs will be used throughout the platform.
  scdata.integrated@misc[["numPCs"]] <- config$dimensionalityReduction$numPCs

  scdata.integrated <- colorObject(scdata.integrated)
  cells_order <- rownames(scdata.integrated@meta.data)
  plot1_data <- unname(purrr::map2(scdata.integrated@reductions$umap@cell.embeddings[, 1], scdata.integrated@reductions$umap@cell.embeddings[, 2], function(x, y) {
    c("x" = x, "y" = y)
  }))

  # Adding color and sample id
  plot1_data <- purrr::map2(
    plot1_data,
    unname(scdata.integrated@meta.data[cells_order, "samples"]),
    function(x, y) {
      append(x, list("sample" = y))
    }
  )

  plot1_data <- purrr::map2(
    plot1_data,
    unname(scdata.integrated@meta.data[cells_order, "color_samples"]),
    function(x, y) {
      append(x, list("col" = y))
    }
  )

  plot2_data <- unname(purrr::map2(1:min(50,length(var_explained)), var_explained, function(x, y) {
    c("PC" = x, "percentVariance" = y)
  }))

  plots <- list()
  plots[generate_gui_uuid("", task_name, 0)] <- list(plot1_data)
  plots[generate_gui_uuid("", task_name, 1)] <- list(plot2_data)


  # the result object will have to conform to this format: {data, config, plotData : {plot1, plot2}}
  result <- list(
    data = scdata.integrated,
    new_ids = cells_id,
    config = config,
    plotData = plots
  )

  return(result)
}

merge_scdatas <- function(scdata_list) {
  saveRDS(scdata_list, '/debug/scdata_list.b4-merge.rds')
  # remove the samples slot before merging or it triggers error because it ain't expected
  scdata_list$samples <- NULL
  if (length(scdata_list) == 1) {
    scdata <- scdata_list[[1]]
  } else {
    scdata <- merge(scdata_list[[1]], y = scdata_list[-1])
  }

  scdata@misc <- scdata_list[[1]]@misc

  return(scdata)
}

# This function covers
#   - Integrate the data using the variable "type" (in case of data integration method is selected) and normalize using LogNormalize method.
#   - Compute PCA analysis
#   - To visualize the results of the batch effect, an UMAP with default setting has been made.
run_dataIntegration <- function(scdata, config) {

  # get method and settings
  method <- config$dataIntegration$method
  npcs <- config$dimensionalityReduction$numPCs

  nsamples <- length(unique(scdata$samples))
  if (nsamples == 1) {
    method <- "unisample"
    message("Only one sample detected or method is non integrate.")
  }

  # we need RNA assay to compute the integrated matrix
  Seurat::DefaultAssay(scdata) <- "RNA"

  integration_function <- get(paste0("run_", method))
  scdata <- integration_function(scdata, config)

  if (is.null(npcs)) {
    npcs <- get_npcs(scdata)
    message("Number of PCs: ", npcs)
  }

  # Compute embedding with default setting to get an overview of the performance of the batch correction
  scdata <- Seurat::RunUMAP(scdata, reduction = scdata@misc[["active.reduction"]], dims = 1:npcs, verbose = FALSE)

  return(scdata)
}

run_harmony <- function(scdata, config) {
  settings <- config$dataIntegration$methodSettings[["harmony"]]

  nfeatures <- settings$numGenes
  normalization <- settings$normalisation
  npcs <- config$dimensionalityReduction$numPCs

  # grep in case misspelled
  if (grepl("lognorm", normalization, ignore.case = TRUE)) normalization <- "LogNormalize"

  scdata <- Seurat::NormalizeData(scdata, normalization.method = normalization, verbose = FALSE)
  scdata <- Seurat::FindVariableFeatures(scdata, nfeatures = nfeatures, verbose = FALSE)
  scdata <- Seurat::ScaleData(scdata, verbose = FALSE)
  scdata <- Seurat::RunPCA(scdata, verbose = FALSE)
  scdata <- harmony::RunHarmony(scdata, group.by.vars = "samples")
  scdata <- add_dispersions(scdata)
  scdata@misc[["active.reduction"]] <- "harmony"

  return(scdata)
}

run_seuratv4 <- function(scdata, config) {
  settings <- config$dataIntegration$methodSettings[["seuratv4"]]

  nfeatures <- settings$numGenes
  normalization <- settings$normalisation

  # for data integration
  npcs <- config$dimensionalityReduction$numPCs

  # grep in case misspelled
  if (grepl("lognorm", normalization, ignore.case = TRUE)) normalization <- "LogNormalize"

  # @misc slots not preserved so transfer
  misc <- scdata@misc

  # Currently, we only support Seurat V4 pipeline for the multisample integration
  data.split <- Seurat::SplitObject(scdata, split.by = "samples")
  for (i in 1:length(data.split)) {
    data.split[[i]] <- Seurat::NormalizeData(data.split[[i]], normalization.method = normalization, verbose = FALSE)
    data.split[[i]] <- Seurat::FindVariableFeatures(data.split[[i]], nfeatures = nfeatures, verbose = FALSE)
  }

  # If Number of anchor cells is less than k.filter/2, there is likely to be an error:
  # Note that this is a heuristic and was found to still fail for small data-sets

  # Try to integrate data (catch error most likely caused by too few cells)
  k.filter <- min(ceiling(sapply(data.split, ncol) / 2), 200)
  tryCatch(
    {
      data.anchors <- Seurat::FindIntegrationAnchors(object.list = data.split, dims = 1:npcs, k.filter = k.filter, verbose = TRUE)
      scdata <- Seurat::IntegrateData(anchorset = data.anchors, dims = 1:npcs)
      Seurat::DefaultAssay(scdata) <- "integrated"
    },
    error = function(e) { # Specifying error message
      # ideally this should be passed to the UI as a error message:
      print(table(scdata$samples))
      print(e)
      print(paste("current k.filter:", k.filter))
      # Should we still continue if data is not integrated? No, right now..
      print("Current number of cells per sample: ")
      print(table(scdata$samples))
      warning("Error thrown in IntegrateData: Probably one/many of the samples contain to few cells.\nRule of thumb is that this can happen at around < 100 cells.")
      # An ideal solution would be to launch an error to the UI, however, for now, we will skip the integration method.
      print("Skipping integration step")
      scdata <- Seurat::NormalizeData(scdata, normalization.method = normalization, verbose = FALSE)
    }
  )

  scdata@misc <- misc
  scdata <- Seurat::FindVariableFeatures(scdata, assay = "RNA", nfeatures = nfeatures, verbose = FALSE)
  scdata <- add_dispersions(scdata)
  scdata@misc[["active.reduction"]] <- "pca"

  # scale and run PCA
  scdata <- Seurat::ScaleData(scdata, verbose = FALSE)
  scdata <- Seurat::RunPCA(scdata, npcs = 50, features = Seurat::VariableFeatures(object = scdata), verbose = FALSE)

  return(scdata)
}

run_fastmnn <- function(scdata, config) {
  settings <- config$dataIntegration$methodSettings[["fastmnn"]]

  nfeatures <- settings$numGenes
  normalization <- settings$normalisation
  npcs <- config$dimensionalityReduction$numPCs

  # grep in case misspelled
  if (grepl("lognorm", normalization, ignore.case = TRUE)) normalization <- "LogNormalize"


  scdata <- Seurat::NormalizeData(scdata, normalization.method = normalization, verbose = FALSE)
  scdata <- Seurat::FindVariableFeatures(scdata, nfeatures = nfeatures, verbose = FALSE)
  scdata <- add_dispersions(scdata)

  # @misc slots not preserved so transfer
  misc <- scdata@misc
  scdata <- SeuratWrappers::RunFastMNN(object.list = Seurat::SplitObject(scdata, split.by = "samples"), d = 50, get.variance = TRUE)
  scdata@misc <- misc
  scdata@misc[["active.reduction"]] <- "mnn"

  return(scdata)
}

run_unisample <- function(scdata, config) {
  settings <- config$dataIntegration$methodSettings[["unisample"]]

  nfeatures <- settings$numGenes
  normalization <- settings$normalisation
  npcs <- config$dimensionalityReduction$numPCs

  # grep in case misspelled
  if (grepl("lognorm", normalization, ignore.case = TRUE)) normalization <- "LogNormalize"

  # in unisample we only need to normalize
  scdata <- Seurat::NormalizeData(scdata, normalization.method = normalization, verbose = FALSE)
  scdata <- Seurat::FindVariableFeatures(scdata, assay = "RNA", nfeatures = nfeatures, verbose = FALSE)
  scdata <- add_dispersions(scdata)
  scdata@misc[["active.reduction"]] <- "pca"

  # scale and run PCA
  scdata <- Seurat::ScaleData(scdata, verbose = FALSE)
  scdata <- Seurat::RunPCA(scdata, npcs = 50, features = Seurat::VariableFeatures(object = scdata), verbose = FALSE)

  return(scdata)
}

add_dispersions <- function(scdata) {
  vars <- Seurat::HVFInfo(object = scdata, assay = "RNA", selection.method = "vst")
  annotations <- scdata@misc[["gene_annotations"]]
  vars$SYMBOL <- annotations$name[match(rownames(vars), annotations$input)]
  vars$ENSEMBL <- rownames(vars)
  scdata@misc[["gene_dispersion"]] <- vars
  return(scdata)
}


colorObject <- function(data) {
  if ("color_pool" %in% names(data@misc)) {
    color_pool <- data@misc[["color_pool"]]
  } else { # THIS SHOULD BE REMOVE ONCE THE EXPERIMENT HAS BEEN UPDATED WITH THE NEW VERSION OF THE DATA-INGEST
    color_pool <- get_color_pool()
  }
  data$color_active_ident <- color_pool[as.numeric(data@active.ident)]

  ##########################
  # Coloring samples
  ###########################
  if ("samples" %in% colnames(data@meta.data)) { # In that case we are in multisample experiment
    data@meta.data[, "color_samples"] <- color_pool[as.numeric(as.factor(data$samples))]
  } else {
    data@meta.data[, "color_samples"] <- color_pool[1]
  }
  return(data)
}

get_explained_variance <- function(scdata) {
  # Compute explained variance for plotting and numPCs estimation.
  # It can be computed from pca or other reductions such as mnn
  if (scdata@misc[["active.reduction"]] == "mnn") {
    var_explained <- scdata@tools$`SeuratWrappers::RunFastMNN`$pca.info$var.explained
  } else {
    eig_values <- (scdata@reductions$pca@stdev)^2
    var_explained <- eig_values / sum(eig_values)
  }
  var_explained
}

get_npcs <- function(scdata, var_threshold = 0.85, max_npcs = 30) {
  # estimates the number of PCs to use in data integration and embeddings,
  # using accumulated explained variance
  var_explained <- get_explained_variance(scdata)
  npcs <- min(which(cumsum(var_explained) >= var_threshold))
  min(npcs, max_npcs, na.rm = TRUE)
}


add_metadata <- function(scdata, annot, experiment_id) {

  # Ensure index by rownames in scdata
  annot <- annot[match(rownames(scdata), annot$input), ]
  scdata@misc[["gene_annotations"]] <- annot

  message("Storing cells id...")
  # Keeping old version of ids starting from 0
  scdata$cells_id <- 0:(ncol(scdata) - 1)

  message("Storing color pool...")
  # We store the color pool in a slot in order to be able to access it during configureEmbedding
  scdata@misc[["color_pool"]] <- get_color_pool()
  scdata@misc[["experimentId"]] <- experiment_id
  scdata@misc[["ingestionDate"]] <- Sys.time()

  return(scdata)
}

