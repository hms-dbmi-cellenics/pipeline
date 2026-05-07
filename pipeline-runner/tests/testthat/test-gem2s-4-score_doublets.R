mock_counts_dbl <- function(...) {
  set.seed(RANDOM_SEED)
  sce <- scDblFinder::mockDoubletSCE(...)
  counts <- sce@assays@data$counts
  colnames(counts) <- paste0("barcode-", 1:ncol(counts))
  rownames(counts) <- paste0("gene-", 1:nrow(counts))
  return(Matrix::Matrix(counts, sparse = TRUE))
}

mock_input <- function(technology) {
  input <- list()
  input$input$type <- technology
  return(input)
}

test_that("score_doublets works with bpcells", {
  counts <- mock_counts_dbl()
  counts <- counts_to_bpcells(
    counts,
    withr::local_tempfile(.local_envir = parent.frame())
  )
  input <- mock_input("10x")

  prev_out <- list(
    counts_list = list(sample1 = counts),
    config = list(),
    annot = list()
  )
  out <- score_doublets(input, NULL, prev_out)$output

  expect_named(out, c("counts_list", "config", "annot", "doublet_scores"))
  expect_named(out$doublet_scores, c("sample1"))
  expect_true(is.data.frame(out$doublet_scores$sample1))
})

test_that("score_doublets returns expected columns", {
  counts <- mock_counts_dbl()
  input <- mock_input("10x")

  prev_out <- list(
    counts_list = list(sample1 = counts),
    config = list(),
    annot = list()
  )
  out <- score_doublets(input, NULL, prev_out)$output

  expect_setequal(
    colnames(out$doublet_scores$sample1),
    c("barcodes", "doublet_class", "doublet_scores")
  )
})


test_that("score_doublets filters cells to avoid warning of extremely low read counts", {
  counts <- mock_counts_dbl(ncells = c(200, 300, 400, 200, 500, 300))
  counts <- round(counts / 2)

  input <- mock_input("10x")

  expect_warning(scDblFinder:::.checkSCE(counts), "extremely low read counts")

  prev_out <- list(
    counts_list = list(sample1 = counts),
    config = list(),
    annot = list()
  )
  expect_warning(score_doublets(input, NULL, prev_out), NA)
})


test_that("score_doublets fails if prev_out is missing 'config', 'counts_list', or 'annot'", {
  counts <- mock_counts_dbl()
  prev_out <- list(
    counts_list = list(sample1 = counts),
    config = list(),
    annot = list()
  )
  input <- mock_input("10x")

  prev_out$config <- NULL
  expect_error(score_doublets(input, NULL, prev_out), "config is missing")

  prev_out$config <- list()
  prev_out$annot <- NULL
  expect_error(score_doublets(input, NULL, prev_out), "annot is missing")

  prev_out$annot <- data.frame()
  prev_out$counts_list <- NULL
  expect_error(score_doublets(input, NULL, prev_out), "counts_list is missing")
})


test_that("compute_sample_doublet_scores handles technology 'parse' correctly", {
  counts <- mock_counts_dbl()
  input <- mock_input("parse")

  prev_out <- list(
    counts_list = list(sample1 = counts),
    config = list(),
    annot = list()
  )
  out <- suppressWarnings(score_doublets(input, NULL, prev_out)$output)

  expect_setequal(
    colnames(out$doublet_scores$sample1),
    c("barcodes", "doublet_class", "doublet_scores")
  )
})


test_that("compute_sample_doublet_scores uses correct dbr for different Parse kits", {
  counts <- mock_counts_dbl()

  dbr_mini <- DOUBLET_RATE_PARSE[["mini"]]
  set.seed(RANDOM_SEED)
  expected_sce_mini <- suppressWarnings(scDblFinder::scDblFinder(counts, dbr = dbr_mini))

  dbr_wt <- DOUBLET_RATE_PARSE[["WT"]]
  set.seed(RANDOM_SEED)
  expected_sce_wt <- suppressWarnings(scDblFinder::scDblFinder(counts, dbr = dbr_wt))

  dbr_mega <- DOUBLET_RATE_PARSE[["mega"]]
  set.seed(RANDOM_SEED)
  expected_sce_mega <- suppressWarnings(scDblFinder::scDblFinder(counts, dbr = dbr_mega))

  observed_sce_mini <- suppressWarnings(compute_sample_doublet_scores(counts, technology = "parse", parse_kit = "mini"))
  observed_sce_wt <- suppressWarnings(compute_sample_doublet_scores(counts, technology = "parse", parse_kit = "WT"))
  observed_sce_mega <- suppressWarnings(compute_sample_doublet_scores(counts, technology = "parse", parse_kit = "mega"))

  expect_identical(expected_sce_mini$scDblFinder.score, observed_sce_mini$doublet_scores)
  expect_identical(expected_sce_wt$scDblFinder.score, observed_sce_wt$doublet_scores)
  expect_identical(expected_sce_mega$scDblFinder.score, observed_sce_mega$doublet_scores)
})


test_that("compute_sample_doublet_scores stops with an error for invalid parse_kit values", {
  counts <- mock_counts_dbl()

  expect_error(
    compute_sample_doublet_scores(counts, technology = "parse", parse_kit = "invalid_kit"),
    "Invalid parse kit value"
  )
})

# Tests for scdblfinder_bpcells function
test_that("scdblfinder_bpcells returns SingleCellExperiment object", {
  counts <- mock_counts_dbl()

  result <- suppressWarnings(scdblfinder_bpcells(counts, dbr = NULL))

  expect_s4_class(result, "SingleCellExperiment")
  expect_true(
    "scDblFinder.class" %in% colnames(SummarizedExperiment::colData(result))
  )
  expect_true(
    "scDblFinder.score" %in% colnames(SummarizedExperiment::colData(result))
  )
})

test_that("scdblfinder_bpcells uses provided dbr parameter", {
  counts <- mock_counts_dbl()
  dbr <- 0.01

  result <- suppressWarnings(scdblfinder_bpcells(counts, dbr = dbr))

  expect_s4_class(result, "SingleCellExperiment")
  expect_equal(ncol(result), ncol(counts))
})

test_that("scdblfinder_bpcells handles counts matrix correctly", {
  counts <- mock_counts_dbl()

  result <- suppressWarnings(scdblfinder_bpcells(counts, dbr = NULL))

  expect_equal(
    colnames(result),
    colnames(counts)
  )
  expect_equal(
    nrow(result),
    nrow(counts)
  )
})

test_that("scdblfinder_bpcells with disk-backed BPCells matrix", {
  # Create a count matrix and write it to disk using BPCells
  counts <- mock_counts_dbl()

  # Ensure matrix is integer (required by BPCells)
  counts <- as(counts, "dgCMatrix")
  if (mode(counts@x) != "integer") {
    counts@x <- as.integer(counts@x)
  }

  # Write the matrix to disk
  matrix_dir <- withr::local_tempfile()
  bpcells_matrix <- BPCells::write_matrix_dir(
    counts, dir = matrix_dir
  )

  # Run scdblfinder_bpcells with the disk-backed matrix
  result <- suppressWarnings(
    scdblfinder_bpcells(bpcells_matrix, dbr = NULL)
  )

  expect_s4_class(result, "SingleCellExperiment")
  expect_true(
    "scDblFinder.class" %in%
      colnames(SummarizedExperiment::colData(result))
  )
  expect_equal(ncol(result), ncol(counts))

})

test_that("scdblfinder_bpcells produces identical result to scDblFinder", {
  counts <- mock_counts_dbl()

  # Run scdblfinder_bpcells
  set.seed(RANDOM_SEED)
  result_bpcells <- suppressWarnings(scdblfinder_bpcells(counts, dbr = NULL))

  # Run standard scDblFinder directly with default arguments
  set.seed(RANDOM_SEED)
  result_standard <- suppressWarnings(scDblFinder::scDblFinder(counts))

  # Compare doublet scores and classes
  expect_identical(
    result_bpcells$scDblFinder.score,
    result_standard$scDblFinder.score
  )
  expect_identical(
    result_bpcells$scDblFinder.class,
    result_standard$scDblFinder.class
  )
})
