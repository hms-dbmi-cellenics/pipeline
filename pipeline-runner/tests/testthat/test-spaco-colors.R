mock_cluster_labels <- function(scdata, n_clusters = 4) {
  ncells <- ncol(scdata)
  set.seed(1)
  labels <- as.character(sample(seq_len(n_clusters) - 1, ncells, replace = TRUE))
  stats::setNames(labels, as.character(scdata$cells_id))
}

is_hex <- function(x) all(grepl("^#[0-9A-Fa-f]{6}$", x))


test_that("resolve_color_pool returns the pool unchanged for non-spatial data", {
  pool <- get_color_pool()

  # NULL scdata, and a spatial mock with its image removed -> non-spatial
  expect_identical(resolve_color_pool(NULL, 0, "0", "0", pool), pool)

  non_spatial <- mock_spatial_scdata(ncells = 16)
  non_spatial[[Seurat::Images(non_spatial)]] <- NULL
  expect_identical(
    resolve_color_pool(
      non_spatial, non_spatial$cells_id, rep("0", ncol(non_spatial)), "0", pool
    ),
    pool
  )
})


test_that("get_spaco_color_map returns NULL for non-spatial data", {
  # a spatial mock with its image removed -> non-spatial
  scdata <- mock_spatial_scdata(ncells = 64)
  scdata[[Seurat::Images(scdata)]] <- NULL
  labels <- stats::setNames(
    rep("0", ncol(scdata)), as.character(scdata$cells_id)
  )

  # no images -> no slices -> NULL, distances are never computed
  mockery::stub(
    get_spaco_color_map, "spatial_distance_r",
    function(...) stop("should not compute distances for non-spatial data")
  )
  expect_null(get_spaco_color_map(scdata, labels))
})


test_that("get_spaco_color_map returns NULL when scdata is NULL", {
  expect_null(get_spaco_color_map(NULL, c("0" = "x")))
})


test_that("get_spaco_color_map auto-generates a color per cluster", {
  scdata <- mock_spatial_scdata(ncells = 100)
  labels <- mock_cluster_labels(scdata, n_clusters = 6)

  color_map <- get_spaco_color_map(scdata, labels)

  # one hex color per cluster, keyed by cluster label
  expect_type(color_map, "character")
  expect_setequal(names(color_map), unique(unname(labels)))
  expect_true(is_hex(color_map))
  expect_equal(length(unique(color_map)), length(color_map))
})


test_that("cluster_color_pool colors a cell_sets frame for spatial data", {
  scdata <- mock_spatial_scdata(ncells = 100)
  labels <- mock_cluster_labels(scdata, n_clusters = 6)
  cell_sets <- data.frame(
    cluster = unname(labels),
    cell_ids = as.integer(names(labels))
  )

  palette <- cluster_color_pool(scdata, cell_sets, get_color_pool())

  # one auto-generated hex color per cluster, in the formatter's sorted order
  categories <- sort_cluster_names(unique(cell_sets$cluster))
  expect_length(palette, length(categories))
  expect_true(is_hex(palette))
})


test_that("cluster_color_pool returns the pool unchanged for non-spatial data", {
  scdata <- mock_spatial_scdata(ncells = 16)
  scdata[[Seurat::Images(scdata)]] <- NULL
  cell_sets <- data.frame(cluster = c("0", "1"), cell_ids = c(0L, 1L))
  pool <- get_color_pool()

  expect_identical(cluster_color_pool(scdata, cell_sets, pool), pool)
})


test_that("spaco_palette aligns colors to the requested category order", {
  scdata <- mock_spatial_scdata(ncells = 100)
  labels <- mock_cluster_labels(scdata, n_clusters = 6)

  categories <- sort_cluster_names(unique(unname(labels)))
  palette <- spaco_palette(scdata, labels, categories)
  color_map <- get_spaco_color_map(scdata, labels)

  expect_length(palette, length(categories))
  # palette[k] must be the auto-generated color for categories[k]
  expect_equal(palette, unname(color_map[as.character(categories)]))
})


test_that("spaco_palette returns NULL when a category is uncolored", {
  scdata <- mock_spatial_scdata(ncells = 100)
  labels <- mock_cluster_labels(scdata, n_clusters = 6)

  # ask for a category that does not exist in the data
  categories <- c(sort_cluster_names(unique(unname(labels))), "not-a-cluster")
  expect_null(spaco_palette(scdata, labels, categories))
})


test_that("spatial_distance_r matches a known checkerboard interlacement", {
  coords <- as.matrix(expand.grid(x = 1:6, y = 1:6))
  labels <- ifelse((coords[, 1] + coords[, 2]) %% 2 == 0, "A", "B")
  m <- spatial_distance_r(coords, labels, radius = 1.5, n_cells = 1L)

  expect_equal(rownames(m), c("A", "B"))
  expect_gt(m["A", "B"], 0)
  expect_equal(m["A", "B"], m["B", "A"])
  expect_equal(unname(diag(m)), c(0, 0))
})


test_that("embed_graph_r returns one distinct hex color per cluster", {
  coords <- as.matrix(expand.grid(x = 1:12, y = 1:12))
  labels <- as.character(rep(1:6, length.out = nrow(coords)))
  m <- spatial_distance_r(coords, labels, radius = 4)

  cols <- embed_graph_r(m)
  expect_named(cols, rownames(m))
  expect_true(is_hex(cols))
  expect_equal(length(unique(cols)), length(cols))
})
