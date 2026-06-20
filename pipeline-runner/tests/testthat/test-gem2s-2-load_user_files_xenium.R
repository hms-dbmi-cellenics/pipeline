# Tests for the Xenium reader (gem2s-2). The reader is disk-backed: counts come
# from cell_feature_matrix.h5 via the shared 10X H5 reader (Gene Expression only)
# and centroids/boundaries from the raw parquet files. We don't ship real Xenium
# fixtures, so read_10x_h5_sample and arrow::read_parquet are mocked; the contract
# the rest of gem2s depends on is the returned list shape.

# a small disk-backed BPCells matrix, like read_10x_h5_sample produces
mock_bpcells_counts <- function(ngenes = 20, ncells = 30) {
  set.seed(1)
  counts <- matrix(
    rpois(ngenes * ncells, lambda = 3),
    nrow = ngenes,
    dimnames = list(
      paste0("Gene", seq_len(ngenes)),
      paste0("cell", seq_len(ncells))
    )
  )
  matrix_dir <- withr::local_tempdir(.local_envir = parent.frame())
  BPCells::write_matrix_dir(as(counts, "dgCMatrix"), dir = matrix_dir, overwrite = TRUE)
}

# cells.parquet / cell_boundaries.parquet as plain data frames
mock_cells_parquet <- function(ncells = 30, with_method = TRUE) {
  side <- ceiling(sqrt(ncells))
  grid <- expand.grid(x = seq_len(side), y = seq_len(side))[seq_len(ncells), ]
  df <- data.frame(
    cell_id = paste0("cell", seq_len(ncells)),
    x_centroid = grid$x,
    y_centroid = grid$y,
    transcript_counts = rpois(ncells, 50)
  )
  if (with_method) df$segmentation_method <- "cell"
  df
}

mock_boundaries_parquet <- function(ncells = 30, nverts = 4) {
  do.call(rbind, lapply(seq_len(ncells), function(i) {
    data.frame(
      cell_id = rep(paste0("cell", i), nverts),
      vertex_x = i + cos(seq_len(nverts)),
      vertex_y = i + sin(seq_len(nverts))
    )
  }))
}


test_that("read_xenium_sample returns the gem2s list contract", {
  counts <- mock_bpcells_counts()

  mockery::stub(read_xenium_sample, "read_10x_h5_sample", list(
    counts = counts,
    annotations = data.frame(
      input = rownames(counts),
      symbol = rownames(counts)
    ),
    matrix_dir = "/tmp/sample_a_matrix_dir"
  ))
  mockery::stub(read_xenium_sample, "read_xenium_segmentations", list(
    centroids = mock_cells_parquet()[, c("x_centroid", "y_centroid", "cell_id")],
    segmentations = mock_boundaries_parquet(),
    segmentation_method = NULL
  ))

  out <- read_xenium_sample("sample_a", tempdir())

  expect_setequal(
    names(out),
    c("counts", "annotations", "matrix_dir", "segmentations", "edrops", "doublet_scores")
  )

  # counts are disk-backed BPCells (not an in-memory dgCMatrix)
  expect_s4_class(out$counts, "IterableMatrix")

  # 2-col annotations
  expect_equal(colnames(out$annotations), c("input", "symbol"))

  # segmentations is the raw-frame list, FOV assembly deferred to create_seurat
  expect_setequal(
    names(out$segmentations),
    c("centroids", "segmentations", "segmentation_method")
  )

  # no EmptyDrops / Doublets for spatial
  expect_length(out$edrops, 0)
  expect_length(out$doublet_scores, 0)
})


test_that("read_xenium_sample keeps only the Gene Expression feature type", {
  # capture the feature_type passed to the shared 10X H5 reader
  reader <- mockery::mock(list(
    counts = mock_bpcells_counts(ngenes = 15),
    annotations = data.frame(input = character(), symbol = character()),
    matrix_dir = tempdir()
  ))
  mockery::stub(read_xenium_sample, "read_10x_h5_sample", reader)
  mockery::stub(read_xenium_sample, "read_xenium_segmentations", list())

  read_xenium_sample("sample_a", tempdir())

  mockery::expect_called(reader, 1)
  expect_equal(
    mockery::mock_args(reader)[[1]]$feature_type,
    "Gene Expression"
  )
})


test_that("read_xenium_segmentations parses the parquet frames", {
  cells <- mock_cells_parquet(ncells = 12)
  boundaries <- mock_boundaries_parquet(ncells = 12)

  # arrow::read_parquet returns cells then boundaries (call order in the reader)
  reader <- mockery::mock(cells, boundaries)
  mockery::stub(read_xenium_segmentations, "arrow::read_parquet", reader)

  out <- read_xenium_segmentations("/some/sample_dir")

  # centroids: x_centroid/y_centroid/cell_id -> x/y/cell
  expect_equal(colnames(out$centroids), c("x", "y", "cell"))
  expect_equal(out$centroids$cell, cells$cell_id)
  expect_equal(out$centroids$x, cells$x_centroid)

  # boundaries: first three columns -> cell/x/y
  expect_equal(colnames(out$segmentations), c("cell", "x", "y"))
  expect_equal(nrow(out$segmentations), nrow(boundaries))

  # segmentation_method keyed by cell id
  expect_equal(rownames(out$segmentation_method), cells$cell_id)
  expect_true(all(out$segmentation_method$segmentation_method == "cell"))
})


test_that("read_xenium_segmentations tolerates a missing segmentation_method column", {
  cells <- mock_cells_parquet(ncells = 8, with_method = FALSE)
  boundaries <- mock_boundaries_parquet(ncells = 8)

  reader <- mockery::mock(cells, boundaries)
  mockery::stub(read_xenium_segmentations, "arrow::read_parquet", reader)

  out <- read_xenium_segmentations("/some/sample_dir")
  expect_null(out$segmentation_method)
})


test_that("read_xenium_files assembles the per-sample results into the gem2s contract", {
  samples <- c("s1", "s2")
  config <- list(samples = samples)

  per_sample <- function(sample, input_dir) {
    list(
      counts = paste0(sample, "_counts"),
      annotations = data.frame(
        input = c("GeneA", "GeneB"),
        symbol = c("GeneA", "GeneB")
      ),
      matrix_dir = paste0(sample, "_matrix_dir"),
      segmentations = paste0(sample, "_segs"),
      edrops = list(),
      doublet_scores = list()
    )
  }

  # deterministic, fork-free fan-out
  mockery::stub(read_xenium_files, "BiocParallel::bplapply", function(X, FUN, ...) {
    lapply(X, function(s) per_sample(s, input_dir = NULL))
  })

  out <- read_xenium_files(config, "./input")

  expect_setequal(
    names(out),
    c("counts_list", "annot", "matrix_dir_list", "segmentations_list", "edrops", "doublet_scores")
  )
  expect_equal(names(out$counts_list), samples)
  expect_equal(out$segmentations_list[["s2"]], "s2_segs")

  # no EmptyDrops / Doublets for spatial
  expect_length(out$edrops, 0)
  expect_length(out$doublet_scores, 0)
})
