# Helpers for spatial (Visium HD) tests.

# adds the <metric>_z column that add_spatial_local_outliers would have computed
add_zscores <- function(scdata_list, metric, zscores_by_sample) {
  for (sample_id in names(zscores_by_sample)) {
    scdata_list[[sample_id]]@meta.data[[paste0(metric, "_z")]] <-
      zscores_by_sample[[sample_id]]
  }
  scdata_list
}

mock_spatial_config <- function(cutoff = 3, direction = "lower", auto = FALSE) {
  list(
    auto = auto,
    enabled = TRUE,
    filterSettings = list(cutoff = cutoff, direction = direction)
  )
}

# Build a SeuratObject with a tissue-coordinate FOV (centroids), as the spatial
# code paths require (Seurat::GetTissueCoordinates must succeed). Cells are laid
# out on a grid so each cell has well-defined spatial neighbours.
mock_spatial_scdata <- function(ncells = 64, sample_id = "spatial1") {
  set.seed(42)
  ngenes <- 30
  counts <- matrix(
    rpois(ngenes * ncells, lambda = 5),
    nrow = ngenes,
    dimnames = list(
      paste0("gene", seq_len(ngenes)),
      paste0(sample_id, "cell", seq_len(ncells))
    )
  )

  scdata <- Seurat::CreateSeuratObject(counts = as(counts, "dgCMatrix"))
  scdata$cells_id <- seq_len(ncells) - 1
  scdata$samples <- sample_id
  scdata$percent.mt <- runif(ncells, 0, 10)

  # grid coordinates
  side <- ceiling(sqrt(ncells))
  grid <- expand.grid(x = seq_len(side), y = seq_len(side))[seq_len(ncells), ]
  coords <- data.frame(
    x = grid$x,
    y = grid$y,
    cell = colnames(scdata)
  )

  centroids <- SeuratObject::CreateCentroids(coords)
  fov <- SeuratObject::CreateFOV(
    coords = centroids,
    type = "centroids",
    assay = "RNA"
  )
  scdata[[paste0(sample_id, ".fov")]] <- fov

  scdata
}

mock_spatial_scdata_list <- function(ncells = 64, sample_id = "spatial1") {
  setNames(list(mock_spatial_scdata(ncells, sample_id)), sample_id)
}
