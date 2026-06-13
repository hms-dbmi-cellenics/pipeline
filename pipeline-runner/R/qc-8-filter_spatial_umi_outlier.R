# ── Spatial local-outlier detection ──────────────────────────────────────────
# factored out of SpotSweeper::localOutliers, repeated for each metric.

find_spatial_knn <- function(scdata, workers = 1, n_neighbors = 36) {
  neighborhood_coords <- Seurat::GetTissueCoordinates(scdata)
  neighborhood_coords <- as.matrix(neighborhood_coords[, c("x", "y")])

  BiocNeighbors::findKNN(
    neighborhood_coords,
    k = n_neighbors,
    warn.ties = FALSE,
    BPPARAM = BiocParallel::MulticoreParam(workers = workers)
  )$index
}

# ~10x faster than spatialEco::outliers
outliers_rfast <- function(x, s = 1.4826) {
  med_x <- Rfast::med(x)
  mad_x <- s * Rfast::med(abs(x - med_x))
  0.6745 * (x - med_x) / mad_x
}

#' Compute spatial local-outlier modified z-scores for a metric
#'
#' For each cell, computes a modified z-score describing how much \code{metric}
#' deviates from its spatial neighborhood. Stores a \code{<metric>_z} column (and
#' \code{<metric>_log} when \code{log = TRUE}) in the returned metadata, row-aligned
#' to \code{scdata@meta.data}.
local_outliers <- function(
  scdata, spatial_knn, metric = "percent.mt", log = FALSE
) {

  metadata <- scdata@meta.data
  if (log) {
    metric_log <- paste0(metric, "_log")
    metadata[[metric_log]] <- log1p(metadata[[metric]])
    metric_to_use <- metric_log
  } else {
    metric_to_use <- metric
  }

  mask <- spatial_knn != 0
  row_idx <- row(spatial_knn)[mask]
  col_vals <- spatial_knn[mask]

  values <- metadata[col_vals, metric_to_use]

  neighborhoods <- split(
    values,
    factor(row_idx, levels = seq_len(nrow(spatial_knn)))
  )

  mod_z_matrix <- vapply(neighborhoods, function(x) {
    outliers_rfast(x)[1]
  }, numeric(1))

  mod_z_matrix[!is.finite(mod_z_matrix)] <- 0

  metric_z <- paste0(metric, "_z")
  metadata[[metric_z]] <- mod_z_matrix
  metadata_cols <- c(metric_to_use, metric_z)

  metadata[, metadata_cols]
}

# ── Spatial local-outlier filters (Visium HD) ────────────────────────────────

#' Spatial local-outlier filters (Visium HD)
#'
#' These filters are specific to spatial datasets. Instead of a global threshold,
#' they flag cells whose metric deviates from its spatial neighborhood, using the
#' modified z-scores precomputed by \code{add_spatial_local_outliers} and stored in
#' \code{meta.data} as \code{<metric>_z} (see \code{local_outliers}).
#'
#' A cell is filtered out when its z-score is beyond the chosen \code{cutoff}:
#' \itemize{
#'   \item \code{direction = "lower"}: outlier (removed) if \code{z <= -cutoff}
#'     (used for UMI count and number of genes — low local values are suspicious).
#'   \item \code{direction = "upper"}: outlier (removed) if \code{z >= cutoff}
#'     (used for mitochondrial content — high local values are suspicious).
#' }
#'
#' @param config list with \code{enabled}, \code{auto} and
#'   \code{filterSettings} (\code{cutoff}, \code{direction}). When \code{auto} is
#'   TRUE the default \code{cutoff} of 3 is used.
#' @param metric character name of the \code{meta.data} metric column
#'   (\code{nCount_RNA}, \code{nFeature_RNA} or \code{percent.mt}).
#' @param direction one of \code{"lower"} or \code{"upper"}.
#' @param log_scale logical, whether the metric was log-scaled to compute z-scores;
#'   when TRUE the plot reports the \code{<metric>_log} value.
#'
#' @return a list with the filtered \code{scdata_list}, the updated config and the
#'   plot values (spatial slide data, z-score histogram data, filter statistics).
#'
#' @export
#'
filter_spatial_local_outlier <- function(
  scdata_list, config, sample_id, cells_id,
  task_name, metric, direction, log_scale = FALSE, num_cells_to_downsample = 6000
) {
  sample_cell_ids <- cells_id[[sample_id]]

  if (length(sample_cell_ids) == 0) {
    return(list(
      data = scdata_list[[sample_id]],
      new_ids = cells_id,
      config = config,
      plotData = list()
    ))
  }

  sample_data <- subset_ids(scdata_list[[sample_id]], sample_cell_ids)
  metric_z <- paste0(metric, "_z")

  # spatial z-scores are computed in gem2s (add_spatial_local_outliers)
  if (!metric_z %in% colnames(sample_data@meta.data)) {
    message("Warning! No spatial z-scores (", metric_z, ") computed for this experiment!")
    return(list(data = scdata_list[[sample_id]], config = config, plotData = list()))
  }

  # cutoff: automatic uses the default of 3, manual uses the user-provided value
  cutoff <- config$filterSettings$cutoff
  if (is.null(cutoff) || (exists("auto", where = config) && safeTRUE(config$auto))) {
    cutoff <- 3
  }
  cutoff <- as.numeric(cutoff)
  config$filterSettings$cutoff <- cutoff
  config$filterSettings$direction <- direction

  # Assign updated config to global env so that it can be accessed if there is an error
  config_key <- paste0("config-", task_name, "-", sample_id)
  assign(config_key, config, envir = globalenv())

  z <- sample_data@meta.data[[metric_z]]
  ids <- sample_data@meta.data$cells_id

  if (safeTRUE(config$enabled)) {
    if (direction == "upper") {
      remaining_ids <- ids[z < cutoff]
    } else {
      remaining_ids <- ids[z > -cutoff]
    }
  } else {
    remaining_ids <- sample_cell_ids
  }

  # the slide plot colours by the same (possibly log-scaled) value the z-score
  # was computed on, so the colouring matches the filtering
  value_metric <- if (log_scale) paste0(metric, "_log") else metric
  if (!value_metric %in% colnames(sample_data@meta.data)) value_metric <- metric

  guidata <- list()

  # plot 0: spatial slide data (full per-cell, NOT downsampled — the overlay needs all cells)
  guidata[[generate_gui_uuid(sample_id, task_name, 0)]] <-
    generate_spatial_outlier_plot_data(sample_data, value_metric, metric_z)

  # plot 1: histogram of outlier z-scores (downsampled like other filters)
  nkeep <- downsample_plotdata(ncol(sample_data), num_cells_to_downsample)
  set.seed(RANDOM_SEED)
  cells_position_to_keep <- sort(sample(seq_len(ncol(sample_data)), nkeep, replace = FALSE))
  guidata[[generate_gui_uuid(sample_id, task_name, 1)]] <-
    lapply(unname(z[cells_position_to_keep]), function(x) c("zscore" = x))

  # plot 2: filter statistics
  guidata[[generate_gui_uuid(sample_id, task_name, 2)]] <- list(
    before = calc_filter_stats(sample_data),
    after = calc_filter_stats(subset_ids(sample_data, remaining_ids))
  )

  cells_id[[sample_id]] <- remaining_ids

  list(
    data = scdata_list,
    new_ids = cells_id,
    config = config,
    plotData = guidata
  )
}

# Builds the per-cell array consumed by the spatial slide plot: one record per
# cell with its stable (0-indexed) cells_id, the metric value (for colouring),
# the local-outlier z-score (to recompute the outlier set client-side) and tissue
# coordinates.
generate_spatial_outlier_plot_data <- function(scdata, value_metric, metric_z) {
  md <- scdata@meta.data

  cols <- list(
    cellId = md$cells_id,
    value = unname(md[[value_metric]]),
    zscore = unname(md[[metric_z]])
  )

  coords <- tryCatch(Seurat::GetTissueCoordinates(scdata), error = function(e) NULL)
  if (!is.null(coords) && nrow(coords) == nrow(md)) {
    rownames(coords) <- coords$cell
    coords <- coords[rownames(md), , drop = FALSE]
    cols$x <- unname(coords$x)
    cols$y <- unname(coords$y)
  }

  purrr::transpose(cols)
}

#' STEP 8. Spatial UMI filter (Visium HD)
#'
#' Filters cells whose total UMI count (\code{nCount_RNA}) is a local outlier on
#' the log scale below the chosen z-score cutoff.
#'
#' @inheritParams filter_spatial_local_outlier
#' @export
#'
filter_spatial_umi_outlier <- function(
  scdata_list, config, sample_id, cells_id,
  task_name = "spatialUmiOutlier", num_cells_to_downsample = 6000
) {
  filter_spatial_local_outlier(
    scdata_list, config, sample_id, cells_id, task_name,
    metric = "nCount_RNA", direction = "lower", log_scale = TRUE,
    num_cells_to_downsample = num_cells_to_downsample
  )
}
