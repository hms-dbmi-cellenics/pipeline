#' STEP 6. Data Integration
#'
#' Data integration step where batch effect is corrected through data
#' integration methods. The default method is Harmony.
#' The data integration step include ordering the samples according to their
#' matrices size, merging of Seurat objects, normalization and PCA analysis.
#' Optional: performs geometric sketching of merged Seurat object, integrates
#' the sketched data, and learn and apply back the integration transformation
#' to the full data.
#'
#' @param config list containing the following information
#' 		- dataIntegration
#' 			- method: String. Method to be used. "harmony by default.
#' 			- methodSettings: List with the method as key and:
#' 				- numGenes: Numeric. Number of gene to be used.
#' 				- normalisation: String. Normalisation method to be used.
#' 		- dimensionalityReduction
#' 			- method: String. Method to be used. "rpca" by default.
#' 			- numPCs: Numeric. Number of principal components.
#' 			- excludeGeneCategories: List. Categories of genes to be excluded.
#' 		- downsampling
#' 		  - method: String. Method to be used to downsample. "geosketch" by default
#' 		  - methodSettings: List with required parameters for each downsampling method
#' 		    - geosketch
#' 		      - percentageToKeep: percentage of cells to keep
#'
#' @return a list with the integrated seurat object, the cell ids, the config and the plot values.
#' @export
#'
integrate_scdata <- function(scdata_list, config, sample_id, cells_id, task_name = "dataIntegration") {

  message("Started data integration")
  method <- config$dataIntegration$method

  if (length(scdata_list) == 1) {
    method <- UNISAMPLE
    message("Only one sample detected or method is non integrate.")
  }

  # integrate
  integration_function <- get(paste0("run_", method))
  scdata_integrated <- integration_function(scdata_list, config, cells_id)

  message("Finished data integration")

  # Update config numPCs with estimated or user provided nPCs
  config$dimensionalityReduction$numPCs <- scdata_integrated@misc$numPCs

  var_explained <- get_explained_variance(scdata_integrated)

  plots <- generate_elbow_plot_data(scdata_integrated, task_name, var_explained)

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
create_scdata <- function(scdata_list, cells_id, merge_data = FALSE) {
  scdata_list <- remove_filtered_cells(scdata_list, cells_id)
  merged_scdatas <- merge_scdata_list(scdata_list, merge_data)
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
  message("Filtering cells.")
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
merge_scdata_list <- function(scdata_list, merge_data = FALSE) {
  message("Merging Seurat objects.")
  if (length(scdata_list) == 1) {
    scdata <- scdata_list[[1]]
  } else {
    scdata <- merge(scdata_list[[1]], y = scdata_list[-1], merge.data = merge_data)
  }

  return(scdata)
}


add_dispersions <- function(scdata, method) {
  if (method == "SCT" && Seurat::DefaultAssay(scdata) == "integrated") {
    vars <- Seurat::HVFInfo(object = scdata, assay = "integrated", selection.method = "sctransform")
    # change colnames as they are when run with selection.method = "vst", otherwise will break the listGenes worker task
    colnames(vars) <- c("mean", "variance", "variance.standardized")
  } else {
    vars <- Seurat::HVFInfo(object = scdata, assay = "RNA", selection.method = "vst")
  }
  annotations <- scdata@misc[["gene_annotations"]]
  vars$SYMBOL <- annotations$name[match(rownames(vars), annotations$input)]
  vars$ENSEMBL <- rownames(vars)
  scdata@misc[["gene_dispersion"]] <- vars
  return(scdata)
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


get_npcs <- function(scdata, var_threshold = 0.90, max_npcs = 30) {
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

  all_genes <- scdata@misc$gene_annotations[,c("input", "name")]

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
    "ribosomal" = build_ribosomal_gene_list,
    "mitochondrial" = build_mitochondrial_gene_list
  )

  exclude_gene_indices <- c()

  for (group in exclude_groups) {
    list_fun <- gene_lists[[group]]
    exclude_gene_indices <- c(exclude_gene_indices, list_fun(all_genes))
  }

  # in case there's a custom list of genes to exclude
  if (length(exclude_custom > 0)) {
    exclude_custom_indices <- na.omit(match(unlist(exclude_custom), all_genes$name))
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
  all_genes <- all_genes[["input"]]
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


#' Make list of ribosomal genes
#'
#' Matches ribosomal gene symbols using a regular expression that covers most
#' cases across several commonly used species. It also matches 3 extra ribosomal
#' genes for human and mouse, which are not covered by the regex.
#'
#' @inheritParams build_cc_gene_list
#'
#' @return integer vector of ribosomal gene indices
#' @export
#'
build_ribosomal_gene_list <- function(all_genes) {
  all_genes <- all_genes[["name"]]

  ribo_gene_indices <- grep(RIBOSOMAL_REGEX, all_genes, ignore.case = TRUE)

  message(
    "Number of ribosomal genes to exclude: ",
    length(ribo_gene_indices)
  )

  return(ribo_gene_indices)
}


#' Make list of mitochondrial genes
#'
#' Matches mitochondrial genes using the same regular expression as used in QC's
#' mitochondrial filter.
#'
#' @inheritParams build_cc_gene_list
#'
#' @return integer vector of mitochondrial gene indices
#' @export
#'
build_mitochondrial_gene_list <- function(all_genes) {
  all_genes <- all_genes[["name"]]

  mito_gene_indices <- grep(MITOCHONDRIAL_REGEX, all_genes, ignore.case = TRUE)

  message("Number of mitochondrial genes to exclude: ",
          length(mito_gene_indices))

  return(mito_gene_indices)
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
  message("Adding metadata.")
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


#' generate elbow plot data
#'
#' Reshapes table to an UI compatible format for elbow/scree plot.
#'
#' @param scdata_integrated integrated seurat object
#' @param task_name character
#' @param var_explained numeric
#'
#' @return list of plot data
#' @export
#'
generate_elbow_plot_data <- function(scdata_integrated, task_name, var_explained) {
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

