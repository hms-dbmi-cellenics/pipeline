#' STEP 10. Spatial mitochondrial-content outlier filter (Visium HD)
#'
#' Filters cells whose mitochondrial content (\code{percent.mt}) is a local
#' outlier (not log-scaled) above the chosen z-score cutoff.
#'
#' @inheritParams filter_spatial_local_outlier
#' @export
#'
filter_spatial_mito_outlier <- function(
  scdata_list, config, sample_id, cells_id,
  task_name = "spatialMitoOutlier", num_cells_to_downsample = 6000
) {
  filter_spatial_local_outlier(
    scdata_list, config, sample_id, cells_id, task_name,
    metric = "percent.mt", direction = "upper",
    num_cells_to_downsample = num_cells_to_downsample
  )
}
