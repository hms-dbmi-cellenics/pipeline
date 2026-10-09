# ── Spatial local-outlier detection ──────────────────────────────────────────
# factored out of SpotSweeper::localOutliers, repeated for each metric.

# n_neighbors = 36 is SpotSweeper's default: three rings of a hexagonal grid
# (6 + 12 + 18). Segmented cells are not on a grid, so the square-grid value of
# 48 does not apply either; OSTA's segmented Visium HD workflow also uses 36.
find_spatial_knn <- function(scdata, workers = 1, n_neighbors = 36) {
  neighborhood_coords <- Seurat::GetTissueCoordinates(scdata)
  neighborhood_coords <- as.matrix(neighborhood_coords[, c("x", "y")])

  # BiocNeighbors::findKNN errors if k exceeds the number of other points; cap it
  # so tiny spatial samples degrade gracefully rather than crash the whole run
  k <- min(n_neighbors, nrow(neighborhood_coords) - 1L)

  BiocNeighbors::findKNN(
    neighborhood_coords,
    k = k,
    warn.ties = FALSE,
    BPPARAM = get_bpparam(workers)
  )$index
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

  metric_vals <- metadata[[metric_to_use]]
  n <- nrow(spatial_knn)

  # SpotSweeper::localOutliers computes each spot's modified z-score over the
  # neighbourhood = the focal spot prepended to its spatial kNN, i.e.
  # outliers(metric[c(i, knn_i)])[1]. findKNN excludes the spot itself and
  # returns a 0 index for any missing neighbour (fewer than k available); map
  # those to NA so they drop out of the row-wise median/MAD. The whole thing is
  # vectorised with matrixStats for Visium HD scale.
  idx <- cbind(seq_len(n), spatial_knn) # focal spot is column 1
  idx[idx == 0] <- NA
  neighbourhood <- matrix(metric_vals[idx], nrow = n)

  med <- matrixStats::rowMedians(neighbourhood, na.rm = TRUE)
  mad <- 1.4826 * matrixStats::rowMedians(abs(neighbourhood - med), na.rm = TRUE)
  mod_z <- 0.6745 * (metric_vals - med) / mad # metric_vals == neighbourhood[, 1]

  # zero-MAD (uniform) neighbourhoods produce non-finite scores
  mod_z[!is.finite(mod_z)] <- 0

  metric_z <- paste0(metric, "_z")
  metadata[[metric_z]] <- mod_z
  metadata[, c(metric_to_use, metric_z)]
}

# ── Spatial QC metric filters (Visium HD, Xenium) ────────────────────────────

#' Spatial QC metric filters (Visium HD, Xenium)
#'
#' Filters cells on one per-cell QC metric with one of two methods, chosen by
#' \code{config$filterSettings$method}:
#' \itemize{
#'   \item \code{"localOutlier"}: flags cells whose metric deviates from its
#'     spatial neighborhood, using the modified z-scores precomputed by
#'     \code{add_spatial_local_outliers} and stored in \code{meta.data} as
#'     \code{<metric>_z}. With \code{direction = "lower"} a cell is removed if
#'     \code{z <= -cutoff}; with \code{direction = "upper"} if \code{z >= cutoff}.
#'   \item \code{"absoluteThreshold"}: keeps cells whose raw metric lies within
#'     \code{[min, max]} (inclusive). A NULL bound is not applied. The upper
#'     direction (mitochondrial content, in percent) only has a \code{max}.
#' }
#'
#' In automatic mode the local-outlier cutoff is 3 and the absolute thresholds
#' are 3-MAD global outlier thresholds of this sample (see
#' \code{generate_default_values_spatial_absolute}).
#'
#' @param config list with \code{enabled}, \code{auto} and \code{filterSettings}
#'   (\code{method}, \code{methodSettings}).
#' @param metric character name of the \code{meta.data} metric column
#'   (\code{nCount_RNA}, \code{nFeature_RNA} or \code{percent.mt}).
#' @param direction one of \code{"lower"} or \code{"upper"}.
#' @param log_scale logical, whether the metric was log-scaled to compute z-scores;
#'   when TRUE the plot reports the \code{<metric>_log} value.
#'
#' @return a list with the filtered \code{scdata_list}, the updated config and the
#'   plot values (spatial slide data, histogram data, filter statistics).
#'
#' @export
#'
filter_spatial_qc_metric <- function(
  scdata_list, config, sample_id, cells_id,
  task_name, metric, direction, log_scale = FALSE, num_cells_to_downsample = 6000
) {
  sample_cell_ids <- cells_id[[sample_id]]

  # the step is not configured for this experiment (e.g. a single-cell one)
  if (length(sample_cell_ids) == 0 || is.null(config)) {
    return(list(
      data = scdata_list,
      new_ids = cells_id,
      config = config,
      plotData = list()
    ))
  }

  sample_data <- subset_ids(scdata_list[[sample_id]], sample_cell_ids)
  metadata <- sample_data@meta.data
  metric_z <- paste0(metric, "_z")
  method <- config$filterSettings$method %||% "localOutlier"

  # spatial z-scores are computed in gem2s (add_spatial_local_outliers)
  if (method == "localOutlier" && !metric_z %in% colnames(metadata)) {
    message("Warning! No spatial z-scores (", metric_z, ") computed for this experiment!")
    return(list(
      data = scdata_list,
      new_ids = cells_id,
      config = config,
      plotData = list()
    ))
  }

  auto <- safeTRUE(config$auto)
  settings <- config$filterSettings$methodSettings[[method]]
  raw <- metadata[[metric]]

  if (method == "localOutlier") {
    cutoff <- settings$cutoff
    if (is.null(cutoff) || auto) cutoff <- 3
    settings$cutoff <- as.numeric(cutoff)

    z <- metadata[[metric_z]]
    keep <- if (direction == "upper") z < settings$cutoff else z > -settings$cutoff
  } else if (method == "absoluteThreshold") {
    if (auto) {
      settings[c("min", "max")] <- list(NULL)
      settings <- utils::modifyList(
        settings,
        generate_default_values_spatial_absolute(raw, direction),
        keep.null = TRUE
      )
    }

    min_value <- if (direction == "upper") NULL else settings$min
    keep <- rep(TRUE, length(raw))
    if (!is.null(min_value)) keep <- keep & raw >= min_value
    if (!is.null(settings$max)) keep <- keep & raw <= settings$max
  } else {
    stop("Unknown spatial QC method: ", method)
  }

  config$filterSettings$method <- method
  config$filterSettings$methodSettings[method] <- list(settings)

  # Assign updated config to global env so that it can be accessed if there is an error
  config_key <- paste0("config-", task_name, "-", sample_id)
  assign(config_key, config, envir = globalenv())

  ids <- metadata$cells_id
  remaining_ids <- if (safeTRUE(config$enabled)) ids[keep %in% TRUE] else sample_cell_ids

  # the slide plot colours by the same (possibly log-scaled) value the z-score
  # was computed on, so the colouring matches the local-outlier filtering
  value_metric <- if (log_scale) paste0(metric, "_log") else metric
  if (!value_metric %in% colnames(metadata)) value_metric <- metric
  z_metric <- if (metric_z %in% colnames(metadata)) metric_z else NULL

  guidata <- list()

  # plot 0: spatial slide data (full per-cell, NOT downsampled — the overlay needs all cells)
  guidata[[generate_gui_uuid(sample_id, task_name, 0)]] <-
    generate_spatial_outlier_plot_data(sample_data, value_metric, z_metric, metric)

  # plot 1: histogram data, z-scores and raw values (downsampled like other filters)
  nkeep <- downsample_plotdata(ncol(sample_data), num_cells_to_downsample)
  set.seed(RANDOM_SEED)
  cells_position_to_keep <- sort(sample(seq_len(ncol(sample_data)), nkeep, replace = FALSE))
  histogram_cols <- list(raw = unname(raw[cells_position_to_keep]))
  if (!is.null(z_metric)) {
    histogram_cols$zscore <- unname(metadata[[z_metric]][cells_position_to_keep])
  }
  guidata[[generate_gui_uuid(sample_id, task_name, 1)]] <- purrr::transpose(histogram_cols)

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

#' Default absolute thresholds for a spatial QC metric
#'
#' 3-MAD global outlier thresholds over one sample, as OSTA and the single-cell
#' mitochondrial filter use: the lower threshold of the log metric for
#' \code{direction = "lower"} (UMIs, genes) and the upper threshold of the raw
#' metric for \code{direction = "upper"} (mitochondrial percent). When the MAD is
#' zero the threshold would equal the median, so no bound is set instead.
#'
#' @param values numeric vector of the raw metric
#' @param direction one of \code{"lower"} or \code{"upper"}
#'
#' @return list with \code{min} and \code{max} (lower) or \code{max} (upper);
#'   a NULL bound is not applied
#'
generate_default_values_spatial_absolute <- function(values, direction) {
  values <- values[is.finite(values)]
  median_value <- stats::median(values)
  # the threshold is back-transformed from log space, so compare with tolerance
  at_median <- function(threshold) {
    !is.finite(threshold) || isTRUE(all.equal(threshold, median_value))
  }

  if (direction == "lower") {
    outliers <- scuttle::isOutlier(values, log = TRUE, type = "lower", nmads = 3)
    threshold <- unname(attr(outliers, "thresholds")[["lower"]])
    if (at_median(threshold)) threshold <- NULL
    return(list(min = threshold, max = NULL))
  }

  outliers <- scuttle::isOutlier(values, type = "higher", nmads = 3)
  threshold <- unname(attr(outliers, "thresholds")[["higher"]])
  if (at_median(threshold)) threshold <- NULL
  list(max = threshold)
}

# Builds the per-cell array consumed by the spatial slide plot: one record per
# cell with its stable (0-indexed) cells_id, the metric value (for colouring),
# the raw metric (to recompute absolute-threshold outliers client-side), the
# local-outlier z-score when there is one (to recompute local outliers) and
# tissue coordinates.
generate_spatial_outlier_plot_data <- function(scdata, value_metric, metric_z, raw_metric = value_metric) {
  md <- scdata@meta.data

  cols <- list(
    cellId = md$cells_id,
    value = unname(md[[value_metric]]),
    raw = unname(md[[raw_metric]])
  )
  if (!is.null(metric_z)) cols$zscore <- unname(md[[metric_z]])

  coords <- Seurat::GetTissueCoordinates(scdata)
  rownames(coords) <- coords$cell
  coords <- coords[rownames(md), , drop = FALSE]
  cols$x <- unname(coords$x)
  cols$y <- unname(coords$y)

  purrr::transpose(cols)
}

#' STEP 8. Spatial UMI filter (Visium HD)
#'
#' Filters cells whose total UMI count (\code{nCount_RNA}) is a local outlier on
#' the log scale below the chosen z-score cutoff.
#'
#' @inheritParams filter_spatial_qc_metric
#' @export
#'
filter_spatial_umi_outlier <- function(
  scdata_list, config, sample_id, cells_id,
  task_name = "spatialUmiOutlier", num_cells_to_downsample = 6000
) {
  filter_spatial_qc_metric(
    scdata_list, config, sample_id, cells_id, task_name,
    metric = "nCount_RNA", direction = "lower", log_scale = TRUE,
    num_cells_to_downsample = num_cells_to_downsample
  )
}
