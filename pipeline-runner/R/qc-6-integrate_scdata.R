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
  # the following operations give different results depending on sample order
  # make sure they are ordered according to their matrices size
  scdata_list <- order_by_size(scdata_list)
  message("Started create_scdata for sample ", sample_id, "\n")
  scdata <- create_scdata(scdata_list, cells_id)
  scdata <- Seurat::FindVariableFeatures(scdata, assay = "RNA", nfeatures = 2000, verbose = FALSE)
  scdata <- Seurat::ScaleData(scdata, verbose = FALSE)
  scdata <- Seurat::RunPCA(scdata, verbose = FALSE)
  message("Finished create_scdata for sample ", sample_id, "\n")

  # main function
  set.seed(RANDOM_SEED)
  scdata_sketches <- perform_geomsketch(scdata)
  message("Started data integration")
  scdata_integrated <- run_dataIntegration(scdata, scdata_sketches, config)
  message("Finished data integration")

  # get  npcs from the UMAP call in integration functions
  npcs <- length(scdata_integrated@commands$RunUMAP@params$dims)
  message("\nSet config numPCs to npcs used in last UMAP call: ", npcs, "\n")
  config$dimensionalityReduction$numPCs <- npcs

  var_explained <- get_explained_variance(scdata_integrated)

  # This same numPCs will be used throughout the platform.
  scdata_integrated@misc[["numPCs"]] <- config$dimensionalityReduction$numPCs

  scdata_integrated <- colorObject(scdata_integrated)

  plots <- generate_elbow_plot_data(scdata_integrated, config, task_name, var_explained)

  # the result object will have to conform to this format: {data, config, plotData : {plot1, plot2}}
  result <- list(
    data = scdata_integrated,
    new_ids = cells_id,
    config = config,
    plotData = plots
  )

  return(result)
}


#' Create the merged Seurat object
#'
#' This function takes care of merging the sample seurat objects, shuffling
#' and adding the metadata to the complete Seurat object. It does so by calling
#' the corresponding functions.
#'
#' @param scdata_list list of SeuratObjects
#'
#' @return SeuratObject
#' @export
#'
create_scdata <- function(scdata_list, cells_id) {
  scdata_list <- remove_filtered_cells(scdata_list, cells_id)
  merged_scdatas <- merge_scdata_list(scdata_list)
  merged_scdatas <- add_metadata(merged_scdatas, scdata_list)

  return(merged_scdatas)
}

#' For each sample, remove filtered cells from the Seurat object
#'
#' @param scdata_list list of SeuratObjects
#' @param cells_id list of cells ids to keep
#'
#' @return list of filtered SeuratObjects
#' @export
#'
remove_filtered_cells <- function(scdata_list, cells_id) {
  for (sample in names(scdata_list)) {
    flat_cell_ids <- unname(unlist(cells_id[[sample]]))
    scdata_list[[sample]] <- subset_ids(scdata_list[[sample]], flat_cell_ids)
  }

  message("Total cells: ", sum(sapply(scdata_list, ncol)))
  return(scdata_list)
}


#' Merge the list of sample Seurat objects
#'
#' If the list contains more than one seurat object, it merges them all. Else,
#' returns the only element.
#'
#' @param scdata_list list of SeuratObjects
#'
#' @return SeuratObject
#' @export
#'
merge_scdata_list <- function(scdata_list) {
  if (length(scdata_list) == 1) {
    scdata <- scdata_list[[1]]
  } else {
    scdata <- merge(scdata_list[[1]], y = scdata_list[-1], merge.data = FALSE)
  }

  return(scdata)
}


learn_from_sketches <- function (scdata, scdata_sketches, scdata_int) {
  # get embeddings from splitted Seurat object
  embeddings_orig <- list(scdata@reductions[["pca"]]@cell.embeddings[, 1:50])
  embeddings_sketch <- list(scdata_sketches@reductions[["pca"]]@cell.embeddings[, 1:50])
  embeddings_sketch_int <- list(scdata_int@reductions[["harmony"]]@cell.embeddings[, 1:50])

  # use python script to learn integration from sketches and apply to whole dataset
  reticulate::source_python("R/learn-apply-transformation.py")
  learned_int <- apply_transf(embeddings_orig, embeddings_sketch, embeddings_sketch_int)
  rownames(learned_int[[1]]) <- colnames(scdata)

  scdata[["harmony"]] <- Seurat::CreateDimReducObject(embeddings = learned_int[[1]], key = "harmony_", assay = Seurat::DefaultAssay(scdata))
  # scdata <- Seurat::RunUMAP(scdata, reduction = "harmony", dims = 1:50, verbose = FALSE)

  return(scdata)
}


# This function covers
#   - Integrate the data using the variable "type" (in case of data integration method is selected) and normalize using LogNormalize method.
#   - Compute PCA analysis
#   - To visualize the results of the batch effect, an UMAP with default setting has been made.
run_dataIntegration <- function(scdata, scdata_sketches, config) {

  # get method and settings
  # method <- config$dataIntegration$method
  method <- "harmony"
  npcs <- config$dimensionalityReduction$numPCs
  exclude_groups <- config$dimensionalityReduction$excludeGeneCategories


  nsamples <- length(unique(scdata$samples))
  if (nsamples == 1) {
    method <- "unisample"
    message("Only one sample detected or method is non integrate.")
  }

  # we need RNA assay to compute the integrated matrix
  Seurat::DefaultAssay(scdata) <- "RNA"
  Seurat::DefaultAssay(scdata_sketches) <- "RNA"

  # remove cell cycle genes if needed
  if (length(exclude_groups) > 0) {
    message("\n------\n")
    scdata <- remove_genes(scdata, exclude_groups)
    message("\n------\n")
  }

  integration_function <- get(paste0("run_", method))
  scdata_int <- integration_function(scdata_sketches, config)

  message("Started learning from sketched")
  scdata <- learn_from_sketches(scdata, scdata_sketches, scdata_int)
  scdata@misc[["active.reduction"]] <- "harmony"
  message("Finished learning from sketched")

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

  scdata <- normalize_data(scdata, normalization, "harmony", nfeatures)

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

  # get reduction method to find integration anchors
  reduction <- config$dimensionalityReduction$method

  # grep in case misspelled
  if (grepl("lognorm", normalization, ignore.case = TRUE)) normalization <- "LogNormalize"

  # @misc slots not preserved so transfer
  misc <- scdata@misc

  data.split <- Seurat::SplitObject(scdata, split.by = "samples")
  for (i in 1:length(data.split)) {
    data.split[[i]] <- normalize_data(data.split[[i]], normalization, "seuratv4", nfeatures)
    # PCA needs to be run also here
    # otherwise when running FindIntegrationAnchors() with reduction="rpca" it will fail because no "pca" is present
    if (reduction == "rpca") {
      message("Running PCA")
      data.split[[i]] <- Seurat::RunPCA(data.split[[i]], verbose = FALSE, npcs = npcs)
    } else {
      message("PCA is not running before integration as CCA method is selected")
    }
  }

  # If Number of anchor cells is less than k.filter/2, there is likely to be an error:
  # Note that this is a heuristic and was found to still fail for small data-sets

  # Try to integrate data (catch error most likely caused by too few cells)
  k.filter <- min(ceiling(sapply(data.split, ncol) / 2), 200)
  tryCatch(
    {
      if (reduction == "rpca") message("Finding integration anchors using RPCA reduction")
      if (reduction == "cca") message("Finding integration anchors using CCA reduction")
      data.anchors <- Seurat::FindIntegrationAnchors(object.list = data.split, dims = 1:npcs, k.filter = k.filter, verbose = TRUE, reduction = reduction)
      scdata <- Seurat::IntegrateData(anchorset = data.anchors, dims = 1:npcs, normalization.method = normalization)
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
    }
  )

  # seurat v4 requires to call the normalize_data function before (on single objects)
  # and after integration (on the integrated object)
  scdata <- normalize_data(scdata, normalization, "seuratv4", nfeatures)
  scdata@misc <- misc
  scdata <- add_dispersions(scdata)
  scdata@misc[["active.reduction"]] <- "pca"

  # run PCA
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


  scdata <- normalize_data(scdata, normalization, "fastmnn", nfeatures)
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
  scdata <- normalize_data(scdata, normalization, "unisample", nfeatures)
  scdata <- add_dispersions(scdata)
  scdata@misc[["active.reduction"]] <- "pca"

  # run PCA
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
  return(var_explained)
}

get_npcs <- function(scdata, var_threshold = 0.85, max_npcs = 30) {
  # estimates the number of PCs to use in data integration and embeddings,
  # using accumulated explained variance
  var_explained <- get_explained_variance(scdata)
  npcs <- min(which(cumsum(var_explained) >= var_threshold))
  return(min(npcs, max_npcs, na.rm = TRUE))
}


#' Remove genes from downstream analysis
#'
#' This function subsets the seurat object, removing genes to be excluded for
#' integration and downstream analysis.
#'
#' It calls list_exclude_genes to build the list of genes to remove.
#'
#' @param scdata Seurat object
#' @param exclude_groups list of groups to exclude
#' @param exclude_custom list of custom (user provided) genes to exclude
#'
#' @return Seurat object without excluded genes
#' @export
#'
remove_genes <- function(scdata, exclude_groups, exclude_custom = list()) {
  message("Excluding genes...")
  message("Number of genes before excluding: ", nrow(scdata))

  all_genes <- scdata@misc$gene_annotations$input

  # build list of genes to exclude
  exclude_genes <- list_exclude_genes(all_genes, exclude_groups, exclude_custom)

  # only subset if there are genes to remove.
  # Seurat removes reductions when subsetting
  if (length(exclude_genes) > 0) {
    message("Total number of genes to exlude: ", length(exclude_genes))

    # subset using input (either ensID, ensID + sym or sym, depending on dataset)
    # subset.Seurat requires genes to keep.
    keep_genes <- scdata@misc$gene_annotations$input[-exclude_genes]
    scdata <- subset(scdata, features = keep_genes)
  }

  message("Number of genes after excluding: ", nrow(scdata))

  return(scdata)
}


#' Build list of genes to exclude
#'
#' This function builds the union of gene indices to exclude, joining all groups.
#' To do so, it calls all the required list_* functions, getting the exclude gene
#' ids and joins them.
#'
#' @param all_genes character vector, gene_annotations$input
#' @param exclude_groups list of groups to exclude
#' @param exclude_custom list of custom (user provided) genes to exclude
#'
#' @return integer vector of gene indices to exclude
#' @export
#'
list_exclude_genes <- function(all_genes, exclude_groups, exclude_custom) {
  gene_lists <- list(
    "cellCycle" = build_cc_gene_list,
    "ribosomal" = NULL,
    "mitochondrial" = NULL
  )

  exclude_gene_indices <- c()

  for (group in exclude_groups) {
    list_fun <- gene_lists[[group]]
    exclude_gene_indices <- c(exclude_gene_indices, list_fun(all_genes))
  }

  # in case there's a custom list of genes to exclude
  if (length(exclude_custom > 0)) {
    exclude_custom_indices <- na.omit(match(unlist(exclude_custom), all_genes))
    exclude_gene_indices <- c(exclude_gene_indices, exclude_custom_indices)
  }

  # remove duplicates
  return(unique(exclude_gene_indices))
}


#' Make list of cell cycle genes
#'
#' Uses cell cycle gene list from Tirosh et. al. 2016, bundled with Seurat.
#'
#' For now it's only useful for Homo sapiens. But could be easily extended to mice
#' converting the human gene names to their orthologs in mice, as suggested by
#' satijalab in [this github issue](https://github.com/satijalab/seurat/issues/2493#issuecomment-575702274)
#'
#' @param all_genes character vector, gene_annotations$input
#'
#' @return integer vector of cell cycle gene indices
#' @export
#'
build_cc_gene_list <- function(all_genes) {
  message("Excluding Cell Cycle genes...")

  # TODO: change when adding species input
  human_cc_genes <- cc_genes[["human"]]
  mouse_cc_genes <- cc_genes[["mouse"]]

  # match human cc genes
  human_cc_ens_indices <- na.omit(match(na.omit(human_cc_genes[["ensembl_id"]]), all_genes))
  human_cc_sym_indices <- na.omit(match(human_cc_genes[["symbol"]], all_genes))

  # match mouse cc genes
  mouse_cc_ens_indices <- na.omit(match(mouse_cc_genes[["ensembl_id"]], all_genes))
  mouse_cc_sym_indices <- na.omit(match(mouse_cc_genes[["symbol"]], all_genes))

  # questionable bit of code. This should work for human, mice, human + mice
  # and ignore other species, since matching is case sensitive.
  cc_gene_indices <- unique(c(
    human_cc_ens_indices,
    human_cc_sym_indices,
    mouse_cc_ens_indices,
    mouse_cc_sym_indices
  ))

  message(
    "Number of Cell Cycle genes to exclude: ",
    length(cc_gene_indices)
  )

  return(cc_gene_indices)
}

#' Add the metadata present in scdata_list into the merged SeuratObject
#'
#' This function adds metadata, some of which is present in the sample SeuratObjects
#' (the experimentID, the gene annotations), and some that is not (the color pool)
#' to the merged SeuratObject.
#'
#' @param scdata SeuratObject
#' @param scdata_list list of sample SeuratObjects
#'
#' @return SeuratObject with added metadata
#' @export
#'
add_metadata <- function(scdata, scdata_list) {
  # misc data is duplicated in each of the samples and it does not
  # need to be merge so pick the data in the first one and add it to the merged dataset
  scdata@misc <- scdata_list[[1]]@misc
  experiment_id <- scdata_list[[1]]@misc[["experimentId"]]

  annot_list <- list()
  for (sample in names(scdata_list)) {
    annot_list[[sample]] <- scdata_list[[sample]]@misc[["gene_annotations"]]
  }

  scdata@misc[["gene_annotations"]] <- format_annot(annot_list)

  # We store the color pool in a slot in order to be able to access it during configureEmbedding
  scdata@misc[["color_pool"]] <- get_color_pool()
  scdata@misc[["experimentId"]] <- experiment_id
  scdata@misc[["ingestionDate"]] <- Sys.time()

  return(scdata)
}


#' Normalize data according to the specific normalization method
#'
#' This function normalize the data taking into account the integration method.
#'
#' @param scdata SeuratObject
#' @param normalization_method normalization method
#' @param integration_method integration method
#' @param nfeatures number of features to pass to Seurat::FindVariableFeatures()
#'
#' @return normalized and scaled SeuratObject
#' @export
#'
normalize_data <- function(scdata, normalization_method, integration_method, nfeatures) {
  if (normalization_method == "LogNormalize") {
    scdata <- log_normalize(scdata, normalization_method, integration_method, nfeatures)
  }
  return(scdata)
}

#' Perform log normalization
#'
#' #' If the integration method is fastMNN, it will skip ScaleData() because
#' fastMNN already performs its own scaling.
#' If the integration method is SeuratV4, the default assay will be set to "integrated",
#' in this case NormalizeData() will not work (see the [integration vignette](https://satijalab.org/seurat/articles/integration_introduction.html)),
#' so here it's skipped.
#'
#' @param scdata SeuratObject
#' @param normalization_method normalization method
#' @param integration_method integration method
#' @param nfeatures number of features to pass to Seurat::FindVariableFeatures()
#'
#' @return normalized and scaled SeuratObject
#' @export
#'
log_normalize <- function(scdata, normalization_method, integration_method, nfeatures) {
  if (Seurat::DefaultAssay(scdata) == "RNA") {
    scdata <- Seurat::NormalizeData(scdata, normalization.method = normalization_method, verbose = FALSE)
  }
  scdata <- Seurat::FindVariableFeatures(scdata, assay = "RNA", nfeatures = nfeatures, verbose = FALSE)
  if (integration_method != "fastmnn") {
    scdata <- Seurat::ScaleData(scdata, verbose = FALSE)
  }
  return(scdata)
}


#' generate elbow plot data
#'
#' Reshapes table to an UI compatible format for elbow/scree plot.
#'
#' @param scdata_integrated integrated seurat object
#' @param config list
#' @param task_name character
#' @param var_explained numeric
#'
#' @return list of plot data
#' @export
#'
generate_elbow_plot_data <- function(scdata_integrated, config, task_name, var_explained) {

  cells_order <- rownames(scdata_integrated@meta.data)

  # plot1_data is an empty list because it is not used anymore by the UI
  plot1_data <- list()

  plot2_data <- unname(purrr::map2(1:min(50, length(var_explained)), var_explained, function(x, y) {
    c("PC" = x, "percentVariance" = y)
  }))

  plots <- list()
  plots[generate_gui_uuid("", task_name, 0)] <- list(plot1_data)
  plots[generate_gui_uuid("", task_name, 1)] <- list(plot2_data)

  return(plots)
}


perform_geomsketch <- function(scdata) {
  scdata <- Seurat::FindVariableFeatures(scdata, assay = "RNA", nfeatures = 2000, verbose = FALSE)
  scdata <- Seurat::ScaleData(scdata, verbose = FALSE)
  scdata <- Seurat::RunPCA(scdata, verbose = FALSE)
  scdata_sketches <- Geosketch(object = scdata, reduction = "pca", dims = 50, num.cells = round(ncol(scdata)/2))
  return(scdata_sketches)
}

Geosketch <- function(object, reduction, dims, num.cells) {
  if(!exists("geosketch")) {
    geosketch <- reticulate::import("geosketch")
  }
  stopifnot(reduction %in% names(object@reductions))
  stopifnot(ncol(object@reductions[[reduction]]) >= dims)

  embeddings <- object@reductions[[reduction]]@cell.embeddings[, 1:dims]
  index <- unlist(geosketch$gs(embeddings, as.integer(num.cells),  one_indexed = TRUE))
  sketch <- object[, index]
  return(sketch)
}
