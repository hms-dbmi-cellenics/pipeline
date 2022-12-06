mock_counts <- function(...) {
  set.seed(RANDOM_SEED)
  sce <- scDblFinder::mockDoubletSCE(...)
  sce@assays@data$counts
}


test_that("score_doublets returns expected columns", {
  counts <- mock_counts()

  prev_out <- list(
    counts_list = list(sample1 = counts),
    config = list(),
    annot = list()
  )
  out <- score_doublets(NULL, NULL, prev_out)$output

  expect_setequal(
    colnames(out$doublet_scores$sample1),
    c("barcodes", "doublet_class", "doublet_scores")
  )
})


test_that("score_doublets filters cells to avoid warning of extremely low read counts", {
  counts <- mock_counts(ncells = c(200, 300, 400, 200, 500, 300))
  counts <- round(counts / 2)

  expect_warning(scDblFinder:::.checkSCE(counts), "extremely low read counts")

  prev_out <- list(
    counts_list = list(sample1 = counts),
    config = list(),
    annot = list()
  )
  expect_warning(score_doublets(NULL, NULL, prev_out), NA)
})


test_that("score_doublets fails if prev_out is missing 'config', 'counts_list', or 'annot'", {
  counts <- mock_counts()
  prev_out <- list(
    counts_list = list(sample1 = counts),
    config = list(),
    annot = list()
  )

  prev_out$config <- NULL
  expect_error(score_doublets(NULL, NULL, prev_out), "config is missing")

  prev_out$config <- list()
  prev_out$annot <- NULL
  expect_error(score_doublets(NULL, NULL, prev_out), "annot is missing")

  prev_out$annot <- data.frame()
  prev_out$counts_list <- NULL
  expect_error(score_doublets(NULL, NULL, prev_out), "counts_list is missing")
})
