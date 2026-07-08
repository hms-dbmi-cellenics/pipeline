#' STEP 9. Spatial number-of-genes outlier filter (Visium HD)
#'
#' Filters cells whose number of detected genes (\code{nFeature_RNA}) is a local
#' outlier on the log scale below the chosen z-score cutoff.
#'
#' @inheritParams filter_spatial_local_outlier
#' @export
#'
filter_spatial_num_genes_outlier <- function(
  scdata_list, config, sample_id, cells_id,
  task_name = "spatialNumGenesOutlier", num_cells_to_downsample = 6000
) {
  filter_spatial_local_outlier(
    scdata_list, config, sample_id, cells_id, task_name,
    metric = "nFeature_RNA", direction = "lower", log_scale = TRUE,
    num_cells_to_downsample = num_cells_to_downsample
  )
}
