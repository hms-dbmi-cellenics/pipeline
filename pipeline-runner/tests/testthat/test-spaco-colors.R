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


# Xenium and Visium HD call GetTissueCoordinates slightly differently: Xenium is
# scaleless (scale = NULL), Visium HD passes the hires/lowres factor. Neither
# should error. get_spaco_color_map swallows errors and returns NULL, so a
# non-NULL colour map is the proof that the slice was coloured without erroring.

test_that("get_spaco_color_map colours Xenium slices without error (scale = NULL)", {
  scdata <- mock_spatial_scdata(ncells = 100, technology = "xenium")
  labels <- mock_cluster_labels(scdata, n_clusters = 5)

  captured_scale <- "unset"
  mockery::stub(
    get_spaco_color_map, "Seurat::GetTissueCoordinates",
    function(object, image, scale, ...) {
      captured_scale <<- scale
      SeuratObject::GetTissueCoordinates(object, image)
    }
  )

  color_map <- get_spaco_color_map(scdata, labels)

  expect_null(captured_scale) # xenium is scaleless
  expect_false(is.null(color_map))
  expect_true(is_hex(color_map))
})


test_that("get_spaco_color_map colours Visium HD slices without error (scale via get_image_scale)", {
  scdata <- mock_spatial_scdata(ncells = 100, technology = "visium_hd")
  labels <- mock_cluster_labels(scdata, n_clusters = 5)

  # the FOV fixture has no @image, so stub get_image_scale to the factor a real
  # VisiumV2 image would yield; GetTissueCoordinates(scale = "lowres") must work.
  mockery::stub(get_spaco_color_map, "get_image_scale", function(...) "lowres")

  captured_scale <- "unset"
  mockery::stub(
    get_spaco_color_map, "Seurat::GetTissueCoordinates",
    function(object, image, scale, ...) {
      captured_scale <<- scale
      SeuratObject::GetTissueCoordinates(object, image)
    }
  )

  color_map <- get_spaco_color_map(scdata, labels)

  expect_equal(captured_scale, "lowres")
  expect_false(is.null(color_map))
  expect_true(is_hex(color_map))
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


test_that("colouring falls back to the pool when too few clusters to embed", {
  # with 2 clusters n_neighbors collapses to 1 and uwot's spectral init errors;
  # get_spaco_color_map should swallow it and return NULL, and the public entry
  # points should hand back the unchanged colour pool.
  scdata <- mock_spatial_scdata(ncells = 64)
  labels <- mock_cluster_labels(scdata, n_clusters = 2)

  expect_null(get_spaco_color_map(scdata, labels))

  pool <- get_color_pool()
  cell_sets <- data.frame(
    cluster = unname(labels), cell_ids = as.integer(names(labels))
  )
  expect_identical(cluster_color_pool(scdata, cell_sets, pool), pool)
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


# ── Faithfulness to Spaco colorize_mutiple_slices ─────────────────────────────
#
# Our colouring is a native R port of Spaco (BrainStOrmics/Spaco,
# spaco/colorization.py::colorize_mutiple_slices). That function computes a
# per-slice "degree of interlacement" graph (distance.spatial_distance), merges
# them across slices by summing per (cluster_i, cluster_j) pair aligned to the
# union of clusters with missing pairs = 0 (pandas
# stack -> concat -> groupby(["v1","v2"]).sum -> unstack -> fillna(0)), then
# embeds the merged graph into Lab colour space (mapping.embed_graph) to auto-
# pick colours.
#
# spatial_distance and the cross-slice merge are fully deterministic, so we pin
# them exactly with an independent, loop-based re-implementation of the Spaco
# algorithm. The embed_graph stage runs a UMAP, and we use uwot::umap2 rather
# than Spaco's numba UMAP, so the final colours differ numerically by design (we
# only assert colour-level properties for that stage).

# Independent (loop-based) re-implementation of spaco/distance.py::spatial_distance.
# For each cell, take its within-radius kNN (self included, as KDTree/nn2 return
# it first); banish cells with fewer than n_cells same-cluster neighbours by
# keeping only self, so they contribute nothing; otherwise add an inverse-distance
# weight, normalised by the neighbourhood size, to the (cluster_i, cluster_j)
# entry for each non-self neighbour; then symmetrise by max with a zero diagonal.
# (Spaco's spatial_distance also accepts a neighbor_weight arg but never uses it.)
ref_spatial_distance <- function(
  coords, labels, radius = 90, n_neighbors = 16L, n_cells = 3L
) {
  clusters <- sort(unique(labels))
  k <- length(clusters)
  n <- nrow(coords)
  score <- matrix(0, k, k, dimnames = list(clusters, clusters))

  knn <- RANN::nn2(coords, coords, k = min(n_neighbors, n))

  for (i in seq_len(n)) {
    nb <- knn$nn.idx[i, ]
    dd <- knn$nn.dists[i, ]
    keep <- dd <= radius
    nb <- nb[keep]
    dd <- dd[keep]

    if (sum(labels[nb] == labels[i]) < n_cells) next # banished
    size <- length(nb) # neighbourhood size, self included

    ci <- match(labels[i], clusters)
    for (t in seq_along(nb)) {
      if (nb[t] == i) next # skip self
      cj <- match(labels[nb[t]], clusters)
      score[ci, cj] <- score[ci, cj] + (1 / dd[t]) / size
    }
  }

  score <- pmax(score, t(score))
  diag(score) <- 0
  score
}

# colorize_multiple_slices accumulates per-slice graphs aligned to the union of
# clusters across slices.
ref_multi_slice_distance <- function(coords_list, labels_list, clusters) {
  merged <- matrix(0, length(clusters), length(clusters),
    dimnames = list(clusters, clusters)
  )
  for (s in seq_along(coords_list)) {
    m <- ref_spatial_distance(coords_list[[s]], labels_list[[s]])
    rn <- rownames(m)
    merged[rn, rn] <- merged[rn, rn] + m
  }
  merged
}

# random-ish point cloud with cluster labels, spread so the radius actually bites
mock_slice <- function(n, k, seed) {
  set.seed(seed)
  coords <- cbind(x = runif(n, 0, 200), y = runif(n, 0, 200))
  labels <- as.character(sample(seq_len(k) - 1, n, replace = TRUE))
  list(coords = coords, labels = labels)
}


test_that("spatial_distance_r reproduces Spaco distance.spatial_distance", {
  s <- mock_slice(n = 150, k = 4, seed = 11)
  expect_equal(
    spatial_distance_r(s$coords, s$labels),
    ref_spatial_distance(s$coords, s$labels)
  )
})

test_that("spatial_distance_r reproduces Spaco banishing at a tight radius", {
  # a smaller radius and larger n_cells banishes more cells; the vectorised and
  # reference paths must drop exactly the same cells
  s <- mock_slice(n = 200, k = 5, seed = 23)
  expect_equal(
    spatial_distance_r(s$coords, s$labels, radius = 30, n_cells = 4L),
    ref_spatial_distance(s$coords, s$labels, radius = 30, n_cells = 4L)
  )
})

test_that("multi-slice graph matches Spaco colorize_mutiple_slices accumulation", {
  s1 <- mock_slice(n = 120, k = 4, seed = 31)
  s2 <- mock_slice(n = 90, k = 3, seed = 47) # fewer clusters -> exercises alignment
  clusters <- sort(unique(c(s1$labels, s2$labels)))

  ours <- merge_cluster_distances(
    list(
      spatial_distance_r(s1$coords, s1$labels),
      spatial_distance_r(s2$coords, s2$labels)
    ),
    clusters
  )
  expect_equal(
    ours,
    ref_multi_slice_distance(
      list(s1$coords, s2$coords), list(s1$labels, s2$labels), clusters
    )
  )
})

test_that("get_spaco_color_map colours across multiple slices", {
  # two spatial slices with disjoint cells_id, sharing a global cluster set
  s1 <- mock_spatial_scdata(ncells = 80, sample_id = "slice1")
  s2 <- mock_spatial_scdata(ncells = 80, sample_id = "slice2")
  s2$cells_id <- s2$cells_id + ncol(s1) # disjoint ids across slices

  l1 <- mock_cluster_labels(s1, n_clusters = 5)
  l2 <- mock_cluster_labels(s2, n_clusters = 5)
  labels <- c(l1, l2)

  color_map <- get_spaco_color_map(list(s1, s2), labels)

  # one distinct, valid colour per category seen across both slices
  expect_setequal(names(color_map), unique(unname(labels)))
  expect_true(is_hex(color_map))
  expect_equal(length(unique(color_map)), length(color_map))
})
