# factored out of SpotSweeper::localOutliers as repeated for each metric
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

local_outliers <- function(
  scdata, spatial_knn, metric = "percent.mt", direction = "lower",
  log = FALSE, cutoff = 3
) {

  if (!direction %in% c("lower", "higher", "both")) {
    stop("'direction' must be one of 'lower', 'higher', or 'both'.")
  }

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
  metric_outliers <- paste0(metric, "_outliers")
  metadata[[metric_outliers]] <- switch(
    direction,
    higher = sapply(mod_z_matrix, function(x) x > cutoff),
    lower = sapply(mod_z_matrix, function(x) x < -cutoff),
    both = sapply(mod_z_matrix, function(x) {
      x > cutoff | x < -cutoff
    })
  )
  metric_z <- paste0(metric, "_z")
  metadata[[metric_z]] <- mod_z_matrix
  metadata_cols <- paste0(metric, c("_log", "_outliers", "_z"))
  if (!log) metadata_cols <- metadata_cols[-1]

  return(metadata[, metadata_cols])
}

# example usage:
# spatial_knn <- find_spatial_knn(scdata)
# mt_out <- local_outliers(scdata, spatial_knn)