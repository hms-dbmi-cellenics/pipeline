#' STEP 11. Spatial count-density filter (Xenium)
#'
#' The OSTA imaging QC rule for segmented cells. A cell is removed when either:
#' \itemize{
#'   \item its count density, \code{log(nCount_RNA / cell_area)} in counts per
#'     um^2, is below this sample's median by more than \code{nmads} MADs
#'     (\code{scuttle::isOutlier(type = "lower")}), or
#'   \item \code{filterControls} is TRUE and its negative control probe plus
#'     codeword counts are more than \code{maxControlPercent} percent of all its
#'     counts (gene counts plus control counts), as in Vannan et al. 2025 and
#'     Seo et al. 2026.
#' }
#' The density threshold is computed on the cells that reach this step. In
#' automatic mode \code{nmads} is 3 and cells above 5\% control counts are removed.
#'
#' @param config list with \code{enabled}, \code{auto} and \code{filterSettings}
#'   (\code{nmads}, \code{filterControls}, \code{maxControlPercent}, \code{binStep})
#' @inheritParams filter_spatial_qc_metric
#'
#' @return a list with the filtered \code{scdata_list}, the updated config and the
#'   plot values (spatial slide data, histogram data, filter statistics)
#' @export
#'
filter_spatial_count_density <- function(
  scdata_list, config, sample_id, cells_id,
  task_name = "spatialCountDensity", num_cells_to_downsample = 6000
) {
  sample_cell_ids <- cells_id[[sample_id]]
  no_op <- list(
    data = scdata_list,
    new_ids = cells_id,
    config = config,
    plotData = list()
  )

  # the step is not configured for this experiment (e.g. a single-cell one)
  if (length(sample_cell_ids) == 0 || is.null(config)) {
    return(no_op)
  }

  sample_data <- subset_ids(scdata_list[[sample_id]], sample_cell_ids)
  metadata <- sample_data@meta.data

  # cell_area comes from Xenium cells.parquet (read in gem2s)
  if (!"cell_area" %in% colnames(metadata)) {
    message("Warning! No cell areas found for this experiment, skipping count density filter")
    return(no_op)
  }

  settings <- config$filterSettings
  if (safeTRUE(config$auto)) {
    settings$nmads <- 3
    settings$filterControls <- TRUE
    settings$maxControlPercent <- 5
  }
  settings$nmads <- as.numeric(settings$nmads %||% 3)
  settings$filterControls <- isTRUE(settings$filterControls %||% TRUE)
  settings$maxControlPercent <- as.numeric(settings$maxControlPercent %||% 5)

  log_density <- log(metadata$nCount_RNA / metadata$cell_area)
  density_stats <- get_log_density_stats(log_density)
  threshold <- density_stats$median - settings$nmads * density_stats$mad

  # zero-count cells (-Inf) are outliers; cells with no area (NA) are kept, as
  # there is nothing to judge them on
  keep <- !((log_density < threshold) %in% TRUE)

  # empty cells (0 / 0) have no control percent; the density rule removes them
  control_counts <- metadata$control_counts %||% rep(0, nrow(metadata))
  control_percent <- 100 * control_counts / (metadata$nCount_RNA + control_counts)
  if (settings$filterControls) {
    keep <- keep & !((control_percent > settings$maxControlPercent) %in% TRUE)
  }

  config$filterSettings <- settings

  # Assign updated config to global env so that it can be accessed if there is an error
  config_key <- paste0("config-", task_name, "-", sample_id)
  assign(config_key, config, envir = globalenv())

  ids <- metadata$cells_id
  remaining_ids <- if (safeTRUE(config$enabled)) ids[keep] else sample_cell_ids

  guidata <- list()

  # plot 0: spatial slide data (full per-cell, NOT downsampled — the overlay needs all cells)
  coords <- Seurat::GetTissueCoordinates(sample_data)
  rownames(coords) <- coords$cell
  coords <- coords[rownames(metadata), , drop = FALSE]
  guidata[[generate_gui_uuid(sample_id, task_name, 0)]] <- purrr::transpose(list(
    cellId = ids,
    value = unname(finite_or_na(log_density)),
    controlPercent = unname(finite_or_na(control_percent)),
    x = unname(coords$x),
    y = unname(coords$y)
  ))

  # plot 1: log density histogram (downsampled like other filters) together with
  # the median and MAD the UI needs to recompute the threshold for any nmads
  nkeep <- downsample_plotdata(ncol(sample_data), num_cells_to_downsample)
  set.seed(RANDOM_SEED)
  cells_position_to_keep <- sort(sample(seq_len(ncol(sample_data)), nkeep, replace = FALSE))
  sampled_density <- log_density[cells_position_to_keep]
  guidata[[generate_gui_uuid(sample_id, task_name, 1)]] <- list(
    values = lapply(
      unname(sampled_density[is.finite(sampled_density)]),
      function(x) c("value" = x)
    ),
    median = density_stats$median,
    mad = density_stats$mad
  )

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

# Median and MAD (scaled, as scuttle::isOutlier uses) of the finite log
# densities. Cells with zero counts have a log density of -Inf and are left out
# of the statistics; they are not in the data after gem2s drops empty cells.
get_log_density_stats <- function(log_density) {
  finite_density <- log_density[is.finite(log_density)]
  list(
    median = stats::median(finite_density),
    mad = stats::mad(finite_density)
  )
}

finite_or_na <- function(x) {
  x[!is.finite(x)] <- NA
  x
}
