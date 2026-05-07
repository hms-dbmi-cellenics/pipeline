mock_config <- function(max_fraction = 0.1) {
  config <- list(
    auto = TRUE,
    enabled = TRUE,
    filterSettings = list(
      method = "absoluteThreshold",
      methodSettings = list(
        absoluteThreshold = list(maxFraction = max_fraction)
      )
    )
  )

  return(config)
}

get_threshold <- function(config) {
  config$filterSettings$methodSettings$absoluteThreshold$maxFraction
}

test_that("filter_high_mito filters based on threshold and works with bpcells", {
  for (use_bpcells in c(FALSE, TRUE)) {
    scdata_list <- mock_scdata_list(use_bpcells = use_bpcells)
    sample1_id <- names(scdata_list)[1]

    # set first 10 cells to have high percent.mt (>10) to be filtered
    scdata_list[[sample1_id]]$percent.mt[1:10] <- 11

    cells_id <- mock_ids(scdata_list)

    # should filter first 10 cells
    config <- mock_config(0.1)
    config$auto <- FALSE
    out <- filter_high_mito(scdata_list, config, sample1_id, cells_id)
    expect_equal(ncol(out$data[[sample1_id]]), 40)
    expect_equal(length(out$new_ids[[sample1_id]]), 30)

    # should filter all cells in first sample
    config <- mock_config(0.01)
    config$auto <- FALSE

    # set all cells to have high percent.mt to filter all out
    scdata_list[[sample1_id]]$percent.mt <- 50

    out <- filter_high_mito(scdata_list, config, sample1_id, cells_id)
    expect_equal(ncol(out$data[[sample1_id]]), 40)
    expect_equal(length(out$new_ids[[sample1_id]]), 0)
  }
})

test_that("filter_high_mito can be disabled", {
  scdata_list <- mock_scdata_list()
  sample1_id <- names(scdata_list)[1]
  sample2_id <- names(scdata_list)[2]

  cells_id <- mock_ids(scdata_list)
  nstart <- ncol(scdata_list[[sample1_id]])

  # should filter first 10 cells
  config <- mock_config(0.1)
  config$enabled <- FALSE
  out <- filter_high_mito(scdata_list, config, sample1_id, cells_id)
  expect_equal(ncol(out$data[[sample1_id]]), nstart)
  expect_equal(out$new_ids[[sample1_id]], 0:39)
  expect_equal(out$new_ids[[sample2_id]], 40:79)
})

test_that("filter_high_mito is sample specific", {
  scdata_list <- mock_scdata_list()
  sample1_id <- names(scdata_list)[1]
  sample2_id <- names(scdata_list)[2]

  # set first 10 cells to have high percent.mt to be filtered
  scdata_list[[sample1_id]]$percent.mt[1:10] <- 11

  cells_id <- mock_ids(scdata_list)

  # should filter first 10 cells of 123abc only
  config <- mock_config(0.1)
  config$auto <- FALSE

  out1 <- filter_high_mito(scdata_list, config, sample1_id, cells_id)
  out2 <- filter_high_mito(scdata_list, config, sample2_id, cells_id)

  expect_equal(ncol(out1$data[[sample1_id]]), 40)
  expect_equal(length(out1$new_ids[[sample1_id]]), 30)
  expect_equal(length(out1$new_ids[[sample2_id]]), 40)
  expect_equal(ncol(out2$data[[sample2_id]]), 40)
  expect_equal(length(out2$new_ids[[sample1_id]]), 40)
  expect_equal(length(out2$new_ids[[sample2_id]]), 40)
})

test_that("filter_high_mito can be set to auto", {
  scdata_list <- mock_scdata_list()
  sample1_id <- names(scdata_list)[1]
  sample2_id <- names(scdata_list)[2]

  cells_id <- mock_ids(scdata_list)

  # would filter all 40 cells in first sample unless auto
  config <- mock_config(0.01)
  out <- filter_high_mito(scdata_list, config, sample1_id, cells_id)

  # check that threshold was updated
  expect_true(get_threshold(out$config) > get_threshold(config))

  # check that updated threshold was used
  nkeep <- length(unlist(out$new_ids[[sample1_id]]))
  expected <- sum(
    scdata_list[[sample1_id]]$percent.mt <= get_threshold(out$config) * 100
  )
  expect_equal(nkeep, expected)

  # didn't subset original data
  expect_equal(ncol(out$data[[sample1_id]]), 40)
})


test_that("data without percent.mt outliers uses the max percentage as the threshold", {
  scdata_list <- mock_scdata_list()
  sample1_id <- names(scdata_list)[1]
  sample2_id <- names(scdata_list)[2]

  set.seed(123)

  scdata_sample_1 <- scdata_list[[sample1_id]]

  scdata_list[[sample1_id]]$percent.mt <- rnorm(ncol(scdata_sample_1), 6)
  cells_id <- mock_ids(scdata_list)

  config <- mock_config()
  out <- filter_high_mito(scdata_list, config, sample1_id, cells_id)

  # nothing filtered
  expect_length(unlist(out$new_ids[[sample1_id]]), ncol(scdata_sample_1))

  # threshold equals max for sample
  expect_equal(
    get_threshold(out$config),
    max(scdata_sample_1$percent.mt / 100)
  )
})
