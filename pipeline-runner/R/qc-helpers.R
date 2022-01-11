#' Title
#'
#' @param scdata
#'
#' @return
#' @export
#'
#' @examples
generate_first_step_ids <- function(scdata) {
  cells_id <- list()
  for (sample_id in unique(scdata$samples)) {
    cells_id[[sample_id]] <- scdata$cells_id[scdata$samples == sample_id]
  }
  return(cells_id)
}

remove_cell_ids <- function(pipeline_config, experiment_id) {
  tasks <- list(
      'cellSizeDistribution',
      'mitochondrialContent',
      'numGenesVsNumUmis',
      'doubletScores',
      'dataIntegration',
      'configureEmbedding'
    )
  keys_to_remove <- list()

  s3 <- paws::s3(config = pipeline_config$aws_config)
  for (task_name in tasks) {
    object_list <- s3$list_objects(pipeline_config$cells_id_bucket, Prefix = paste0(experiment_id, "/", task_name, "/"))
    for (object in object_list$Contents) {
      keys_to_remove <- append(keys_to_remove, object$Key)
      s3$delete_object(
        Bucket = pipeline_config$cells_id_bucket,
        Key = object$Key
      )
    }
  }
  message("Cell ids keys deleted: ", keys_to_remove)
}

#
# Returns the positions to keep based on scdata and number of cells to keep
#
get_positions_to_keep <- function(scdata, num_cells_to_downsample) {
  # Downsample plotData
  num_cells_to_downsample <- downsample_plotdata(ncol(scdata), num_cells_to_downsample)
  set.seed(123)
  cells_position_to_keep <- sample(1:ncol(scdata), num_cells_to_downsample, replace = FALSE)
  cells_position_to_keep <- sort(cells_position_to_keep)

  return(cells_position_to_keep)
}

#
# subset_ids subsets a seurat object with the cell ids
#
subset_ids <- function(scdata, cells_id) {
  meta_data_subset <- scdata@meta.data[match(cells_id, scdata@meta.data$cells_id), ]
  current_cells <- rownames(meta_data_subset)
  scdata <- subset_safe(scdata, cells = current_cells)
  return(scdata)
}

#
# Generate GUI plots and data uuid as required by the UI
#
generate_gui_uuid <- function(sample_uuid, task_name, item_idx) {
  if (sample_uuid != "") {
    return(paste(sample_uuid, task_name, item_idx, sep = "-"))
  }

  return(paste(task_name, item_idx, sep = "-"))
}

#
# Subset safe allows us to attempt to subset a seurat object even with an empty list
# this function exists for the case when the user over-filters the object, and we need to return something
# that'd allow the user to realize that they are filtering all the cells, while maintaining certain seurat functionality.
# it's a questionable function and it should be questioned.
#
# IN scdata: object to filter
# IN cells: cell barcodes to subset the object with
#
subset_safe <- function(scdata, cells) {
  if (length(cells) > 0) {
    return(subset(scdata, cells = cells))
  } else {
    return(subset(scdata, cells = colnames(scdata)[1]))
  }
}


#
# down sample plots
#

downsample_plotdata <- function(ncol_sample, max_number_of_cells) {
  nkeep <- min(max_number_of_cells, ncol_sample)
  if (nkeep < ncol_sample) {
    message("sample of size ", ncol_sample, " downsampled to ", nkeep, " cells")
  }

  return(nkeep)
}


handle_debug <- function(scdata, config, task_name, sample_id, debug_config) {
  is_debug <- debug_config$step %in% c(task_name, "all")

  if (is_debug) {
    # variable names used by functions
    scdata <- scdata
    num_cells_to_downsample <- 6000
    sample_str <- ifelse(sample_id == "", "", paste0("_", sample_id))
    fname <- paste0(task_name, sample_str, ".RData")
    fpath_cont <- file.path("/debug", fname)
    fpath_host <- file.path(debug_config$path, fname)
    message(sprintf("⚠️ DEBUG_STEP = %s. Saving arguments.", task_name))
    save(scdata, config, task_name, sample_id, num_cells_to_downsample, file = fpath_cont)
    message(sprintf("⚠️ RUN load('%s') to restore environment.", fpath_host))
  }
}


#' Calculates statistics before/after filter step
#'
#' @param scdata \code{SeuratObject}
#' @param sample_id sample name in \code{scdata$samples} to compute statistics for
#'
#' @return list with \itemize{
#'   \item{"num_cells"}{Number of cells in sample}
#'   \item{"total_genes"}{Number of detected genes in sample}
#'   \item{"median_genes"}{Median number of genes detected per cell}
#'   \item{"median_umis"}{Median number of counts per cell}
#' }
#'
calc_filter_stats <- function(scdata) {
  if (!nrow(scdata)) {
    return(list(
      num_cells = 0,
      total_genes = 0,
      median_genes = 0,
      median_umis = 0
    ))
  }

  # number of counts per gene
  ncount <- Matrix::rowSums(scdata[["RNA"]]@counts)

  list(
    num_cells = ncol(scdata),
    total_genes = sum(ncount > 0),
    median_genes = median(scdata$nFeature_RNA),
    median_umis = median(scdata$nCount_RNA)
  )
}

safeTRUE <- function(x) {
  isTRUE(as.logical(x))
}
