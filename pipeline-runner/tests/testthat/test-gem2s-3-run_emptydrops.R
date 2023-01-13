mock_counts <- function() {
  set.seed(RANDOM_SEED)
  DropletUtils:::simCounts()
}


test_that("run_emptydrops is skipped when dataset is pre-filtered", {
  counts <- mock_counts()
  ntot <- Matrix::colSums(counts)
  counts <- counts[, ntot > 100]

  prev_out <- list(
    counts_list = list(sample1 = counts),
    config = list(),
    annot = data.frame()
  )
  expect_message(
    out <- run_emptydrops(NULL, NULL, prev_out)$output,
    "filtered"
  )

  expect_named(out, c("counts_list", "config", "annot", "edrops"))
  expect_equal(out$edrops, list())
})

test_that("run_emptydrops runs when not pre-filtered", {
  counts <- mock_counts()
  prev_out <- list(
    counts_list = list(sample1 = counts, sample2 = counts),
    config = list(),
    annot = data.frame()
  )
  out <- run_emptydrops(NULL, NULL, prev_out)$output

  # names of stored things are right
  expect_named(out, c("counts_list", "config", "annot", "edrops"))
  expect_named(out$edrops, c("sample1", "sample2"))

  # emptyDrops stores something
  expect_s4_class(out$edrops$sample1, "DFrame")

  # preserves counts list
  expect_equal(prev_out$counts_list, out$counts_list)
})


test_that("run_emptydrops fails if prev_out is missing 'config', 'counts_list', or 'annot'", {
  counts <- mock_counts()
  prev_out <- list(
    counts_list = list(sample1 = counts),
    config = list(),
    annot = data.frame()
  )

  prev_out$config <- NULL
  expect_error(run_emptydrops(NULL, NULL, prev_out), "config is missing")

  prev_out$config <- list()
  prev_out$annot <- NULL
  expect_error(run_emptydrops(NULL, NULL, prev_out), "annot is missing")

  prev_out$annot <- data.frame()
  prev_out$counts_list <- NULL
  expect_error(run_emptydrops(NULL, NULL, prev_out), "counts_list is missing")
})
