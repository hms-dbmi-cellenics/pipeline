mock_counts <- function(...) {
  set.seed(RANDOM_SEED)
  sce <- scDblFinder::mockDoubletSCE(...)
  counts <- sce@assays@data$counts
  colnames(counts) <- paste0("barcode-", 1:ncol(counts))
  rownames(counts) <- paste0("gene-", 1:nrow(counts))
  return(Matrix::Matrix(counts, sparse = T))
}

mock_input <- function(type) {
  input <- list()
  input$input$type <- type
  return(input)
}

test_that("score_doublets returns expected columns", {
  counts <- mock_counts()
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
  counts <- mock_counts(ncells = c(200, 300, 400, 200, 500, 300))
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
  counts <- mock_counts()
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


test_that("compute_sample_doublet_scores handles type 'parse' correctly", {
  counts <- mock_counts()
  input <- mock_input("parse")

  prev_out <- list(
    counts_list = list(sample1 = counts),
    config = list(),
    annot = list()
  )
  out <-suppressWarnings(score_doublets(input, NULL, prev_out)$output)

  expect_setequal(
    colnames(out$doublet_scores$sample1),
    c("barcodes", "doublet_class", "doublet_scores")
  )
})


test_that("compute_sample_doublet_scores uses correct dbr for different Parse kits", {
  counts <- mock_counts()

  dbr_mini <- DOUBLET_RATE_MINI
  set.seed(RANDOM_SEED)
  expected_sce_mini <- suppressWarnings(scDblFinder::scDblFinder(counts, dbr = dbr_mini))

  dbr_wt <- DOUBLET_RATE_WT
  set.seed(RANDOM_SEED)
  expected_sce_wt <- suppressWarnings(scDblFinder::scDblFinder(counts, dbr = dbr_wt))

  dbr_mega <- DOUBLET_RATE_MEGA
  set.seed(RANDOM_SEED)
  expected_sce_mega <- suppressWarnings(scDblFinder::scDblFinder(counts, dbr = dbr_mega))

  observed_sce_mini <- suppressWarnings(compute_sample_doublet_scores(counts, type = "parse", parse_kit = "mini"))
  observed_sce_wt <- suppressWarnings(compute_sample_doublet_scores(counts, type = "parse", parse_kit = "WT"))
  observed_sce_mega <- suppressWarnings(compute_sample_doublet_scores(counts, type = "parse", parse_kit = "mega"))

  expect_identical(expected_sce_mini$scDblFinder.score, observed_sce_mini$doublet_scores)
  expect_identical(expected_sce_wt$scDblFinder.score, observed_sce_wt$doublet_scores)
  expect_identical(expected_sce_mega$scDblFinder.score, observed_sce_mega$doublet_scores)
})


test_that("compute_sample_doublet_scores stops with an error for invalid parse_kit values", {
  counts <- mock_counts()

  expect_error(
    compute_sample_doublet_scores(counts, type = "parse", parse_kit = "invalid_kit"),
    "Invalid parse kit value"
  )
})
