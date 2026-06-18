# Tests for the Xenium reader (gem2s-2). LoadXenium is mocked: we don't ship real
# Xenium files in the test fixtures, and the contract the rest of gem2s depends
# on is the returned list shape, not LoadXenium's internals.

# A Seurat object shaped like LoadXenium's output: a "Xenium" assay (Gene
# Expression counts) plus a centroids FOV image. No tissue image (FOV, microns).
mock_loadxenium_obj <- function(ngenes = 20, ncells = 30) {
  set.seed(1)
  counts <- matrix(
    rpois(ngenes * ncells, lambda = 3),
    nrow = ngenes,
    dimnames = list(
      paste0("Gene", seq_len(ngenes)),
      paste0("cell", seq_len(ncells))
    )
  )

  obj <- SeuratObject::CreateSeuratObject(
    counts = as(counts, "dgCMatrix"),
    assay = "Xenium"
  )

  side <- ceiling(sqrt(ncells))
  grid <- expand.grid(x = seq_len(side), y = seq_len(side))[seq_len(ncells), ]
  coords <- data.frame(x = grid$x, y = grid$y, cell = colnames(obj))
  fov <- SeuratObject::CreateFOV(
    SeuratObject::CreateCentroids(coords),
    type = "centroids",
    assay = "Xenium"
  )
  obj[["fov"]] <- fov

  obj
}


test_that("read_xenium_sample returns the gem2s list contract", {
  obj <- mock_loadxenium_obj()
  mockery::stub(read_xenium_sample, "Seurat::LoadXenium", function(...) obj)

  out <- read_xenium_sample("sample_a", tempdir())

  expect_setequal(
    names(out),
    c("counts", "annotations", "matrix_dir", "segmentations")
  )

  # counts are disk-backed BPCells (not an in-memory dgCMatrix)
  expect_s4_class(out$counts, "IterableMatrix")
  expect_true(dir.exists(out$matrix_dir))

  # 2-col annotations; Xenium is keyed by symbol so input == symbol
  expect_equal(colnames(out$annotations), c("input", "symbol"))
  expect_equal(out$annotations$input, out$annotations$symbol)
  expect_equal(out$annotations$symbol, rownames(out$counts))

  # segmentations is the FOV; GetTissueCoordinates works on it (no image needed)
  expect_no_error(Seurat::GetTissueCoordinates(out$segmentations))
})


test_that("read_xenium_sample keeps only the Gene Expression (Xenium) assay", {
  obj <- mock_loadxenium_obj(ngenes = 15)
  mockery::stub(read_xenium_sample, "Seurat::LoadXenium", function(...) obj)

  out <- read_xenium_sample("sample_a", tempdir())

  # control/codeword assays are dropped: only the 15 Xenium genes survive
  expect_equal(nrow(out$counts), 15)
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
      segmentations = paste0(sample, "_fov")
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
  expect_equal(out$segmentations_list[["s2"]], "s2_fov")

  # no EmptyDrops / Doublets for spatial
  expect_length(out$edrops, 0)
  expect_length(out$doublet_scores, 0)
})
