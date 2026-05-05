mock_config <- function(
  thr = 0.1,
  auto = FALSE,
  enabled = TRUE,
  recomputeDoubletScore = FALSE,
  sampleTechnology = "10x"
) {
  config <- list(
    auto = auto,
    enabled = enabled,
    filterSettings = list(
      probabilityThreshold = thr
    ),
    recomputeDoubletScore = recomputeDoubletScore,
    sampleTechnology = sampleTechnology
  )
  return(config)
}

test_that("filter_doublets filters based on threshold", {
  scdata_list <- mock_scdata_list()
  sample1_id <- names(scdata_list)[1]

  # set first 10 cells to have high doublet scores and doublet class
  scdata_list[[sample1_id]]$doublet_scores[1:10] <- 0.9
  scdata_list[[sample1_id]]$doublet_class[1:10] <- "doublet"

  cells_id <- mock_ids(scdata_list)
  # should filter first 10 cells
  config <- mock_config(0.5)

  out <- filter_doublets(scdata_list, config, sample1_id, cells_id)

  expect_equal(ncol(out$data[[sample1_id]]), 40)
  expect_equal(out$new_ids[[sample1_id]], 10:39)
})

test_that("filter_doublets is sample aware", {
  scdata_list <- mock_scdata_list()
  sample1_id <- names(scdata_list)[1]
  sample2_id <- names(scdata_list)[2]

  # set first 10 cells of sample1 to have high doublet scores and doublet class
  scdata_list[[sample1_id]]$doublet_scores[1:10] <- 0.9
  scdata_list[[sample1_id]]$doublet_class[1:10] <- "doublet"
  # set last 10 cells of sample2 to have high doublet scores and doublet class
  scdata_list[[sample2_id]]$doublet_scores[31:40] <- 0.9
  scdata_list[[sample2_id]]$doublet_class[31:40] <- "doublet"
  cells_id <- mock_ids(scdata_list)
  config <- mock_config(0.5)

  out <- filter_doublets(scdata_list, config, sample1_id, cells_id)
  expect_equal(ncol(out$data[[sample1_id]]), 40)
  expect_equal(out$new_ids[[sample1_id]], 10:39)
  expect_equal(out$new_ids[[sample2_id]], 40:79)

  out <- filter_doublets(out$data, config, sample2_id, out$new_ids)
  expect_equal(ncol(out$data[[sample2_id]]), 40)
  expect_equal(out$new_ids[[sample1_id]], 10:39)
  expect_equal(out$new_ids[[sample2_id]], 40:69)
})

test_that("filter_doublets filters works with auto", {
  scdata_list <- mock_scdata_list()
  sample1_id <- names(scdata_list)[1]

  # set first 10 cells to have high doublet scores and doublet class
  scdata_list[[sample1_id]]$doublet_scores[1:10] <- 0.9
  scdata_list[[sample1_id]]$doublet_class[1:10] <- "doublet"

  cells_id <- mock_ids(scdata_list)
  # should filter first 10 cells
  config <- mock_config(0.001, auto = TRUE)
  out <- filter_doublets(scdata_list, config, sample1_id, cells_id)
  expect_equal(out$new_ids[[sample1_id]], 10:39)

  # reset scdata_list for second test
  scdata_list <- mock_scdata_list()
  scdata_list[[sample1_id]]$doublet_scores[1:10] <- 0.9
  scdata_list[[sample1_id]]$doublet_class[1:10] <- "doublet"
  cells_id <- mock_ids(scdata_list)

  config <- mock_config(0.001, auto = FALSE)
  out <- filter_doublets(scdata_list, config, sample1_id, cells_id)
  expect_equal(length(out$new_ids[[sample1_id]]), 0)
})

test_that("filter_doublets can be disabled", {
  scdata_list <- mock_scdata_list()
  sample1_id <- names(scdata_list)[1]
  sample2_id <- names(scdata_list)[2]

  # set first 10 cells of sample1 to have high doublet scores and doublet class
  scdata_list[[sample1_id]]$doublet_scores[1:10] <- 0.9
  scdata_list[[sample1_id]]$doublet_class[1:10] <- "doublet"

  cells_id <- mock_ids(scdata_list)
  config <- mock_config(0.5, enabled = FALSE)
  out <- filter_doublets(scdata_list, config, sample1_id, cells_id)
  expect_equal(out$new_ids[[sample1_id]], 0:39)
  expect_equal(out$new_ids[[sample2_id]], 40:79)

  config <- mock_config(0.5, enabled = TRUE)
  out <- filter_doublets(scdata_list, config, sample1_id, cells_id)
  expect_equal(out$new_ids[[sample1_id]], 10:39)
  expect_equal(out$new_ids[[sample2_id]], 40:79)
})


test_that("generate_default_values_doubletScores sets threshold to 0 when there are no singlets", {
  scdata_list <- mock_scdata_list()
  sample1_id <- names(scdata_list)[1]

  scdata_list[[1]]$doublet_class <- "doublet"

  expect_equal(generate_default_values_doubletScores(scdata_list[[1]]), 0)
})


test_that("doublet scores are re-computed if the API says so", {

  scdata_list <- mock_scdata_list()
  sample1_id <- names(scdata_list)[1]

  cells_id <- mock_ids(scdata_list)
  config <- mock_config(recomputeDoubletScore = TRUE)

  out <- suppressWarnings(
    filter_doublets(scdata_list, config, sample1_id, cells_id)
  )

  expect_false(
    identical(
      out$data[[sample1_id]]$doublet_scores,
      scdata_list[[sample1_id]]$doublet_scores
    )
  )
})
