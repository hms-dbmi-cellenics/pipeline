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
  # subset each of the samples before merging them
  for (sample in names(scdata_list)) {
    flat_cell_ids <- unname(unlist(cells_id[[sample]]))
    scdata_list[[sample]] <- subset_ids(scdata_list[[sample]], flat_cell_ids)
  }

  message("Total cells: ", sum(sapply(scdata_list, ncol)))
  scdata <- create_scdata(scdata_list)

  # main function
  set.seed(RANDOM_SEED)
  message("running data integration")
  scdata_integrated <- run_dataIntegration(scdata, config)

  # get  npcs from the UMAP call in integration functions
  npcs <- length(scdata_integrated@commands$RunUMAP@params$dims)
  message("\nSet config numPCs to npcs used in last UMAP call: ", npcs, "\n")
  config$dimensionalityReduction$numPCs <- npcs

  var_explained <- get_explained_variance(scdata_integrated)

  # This same numPCs will be used throughout the platform.
  scdata_integrated@misc[["numPCs"]] <- config$dimensionalityReduction$numPCs

  scdata_integrated <- colorObject(scdata_integrated)
  cells_order <- rownames(scdata_integrated@meta.data)
  plot1_data <- unname(purrr::map2(scdata_integrated@reductions$umap@cell.embeddings[, 1], scdata_integrated@reductions$umap@cell.embeddings[, 2], function(x, y) {
    c("x" = x, "y" = y)
  }))

  # Adding color and sample id
  plot1_data <- purrr::map2(
    plot1_data,
    unname(scdata_integrated@meta.data[cells_order, "samples"]),
    function(x, y) {
      append(x, list("sample" = y))
    }
  )

  plot1_data <- purrr::map2(
    plot1_data,
    unname(scdata_integrated@meta.data[cells_order, "color_samples"]),
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
#' @param scdata_list
#'
#' @return SeuratObject
#' @export
#'
create_scdata <- function(scdata_list) {

  merged_scdatas <- merge_scdata_list(scdata_list)
  merged_scdatas <- add_metadata(merged_scdatas, scdata_list)

  return(merged_scdatas)
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
    scdata <- merge(scdata_list[[1]], y = scdata_list[-1])
  }

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
  exclude_groups <- config$dimensionalityReduction$excludeGeneCategories


  nsamples <- length(unique(scdata$samples))
  if (nsamples == 1) {
    method <- "unisample"
    message("Only one sample detected or method is non integrate.")
  }

  # we need RNA assay to compute the integrated matrix
  Seurat::DefaultAssay(scdata) <- "RNA"

  # remove cell cycle genes if needed
  if(length(exclude_groups) > 0) {
    message("\n------\n")
    scdata <- remove_genes(scdata, exclude_groups)
    message("\n------\n")
  }

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

  # get reduction method to find integration anchors
  reduction <- config$dimensionalityReduction$method

  # grep in case misspelled
  if (grepl("lognorm", normalization, ignore.case = TRUE)) normalization <- "LogNormalize"

  # @misc slots not preserved so transfer
  misc <- scdata@misc

  # Currently, we only support Seurat V4 pipeline for the multisample integration
  data.split <- Seurat::SplitObject(scdata, split.by = "samples")
  for (i in 1:length(data.split)) {
    data.split[[i]] <- Seurat::NormalizeData(data.split[[i]], normalization.method = normalization, verbose = FALSE)
    data.split[[i]] <- Seurat::FindVariableFeatures(data.split[[i]], nfeatures = nfeatures, verbose = FALSE)
    # PCA needs to be run also here
    # otherwise when running FindIntegrationAnchors() with reduction="rpca" it will fail because no "pca" is present
    if (reduction == "rpca") {
      message("Running PCA")
      data.split[[i]] <- Seurat::ScaleData(data.split[[i]], verbose = FALSE)
      data.split[[i]] <- Seurat::RunPCA(data.split[[i]], verbose = FALSE, npcs = npcs)
    }
    else {
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

  gene_lists <- list("cellCycle" = build_cc_gene_list,
                     "ribosomal" = NULL,
                     "mitochondrial" = NULL)

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
  cc_gene_indices <- unique(c(human_cc_ens_indices,
                            human_cc_sym_indices,
                            mouse_cc_ens_indices,
                            mouse_cc_sym_indices))

  message("Number of Cell Cycle genes to exclude: ",
                  length(cc_gene_indices))

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

