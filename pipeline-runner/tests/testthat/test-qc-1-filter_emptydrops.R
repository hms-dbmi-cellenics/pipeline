mock_config <- function() {
  config <- list(
    auto = FALSE,
    enabled = TRUE,
    filterSettings = list(
      FDR = 0.01
    )
  )

  return(config)
}


test_that("Initial cells_id are correct", {
  scdata_list <- mock_scdata_list()
  cells_id <- generate_first_step_ids(scdata_list)
  expect_equal(unique(cells_id[["123abc"]]), 0:39)
  expect_equal(unique(cells_id[["123def"]]), 40:79)
})

test_that("filter_emptydrops removes NA with threshold < 1 and works with bpcells", {
  for (use_bpcells in c(FALSE, TRUE)) {
    scdata_list <- mock_scdata_list(use_bpcells = use_bpcells)
    
    # Set first 10 cells to have high FDR (>0.01) to be filtered
    scdata_list[["123abc"]]$emptyDrops_FDR[1:10] <- NA
    scdata_list[["123def"]]$emptyDrops_FDR[1:10] <- NA
    
    config <- mock_config()
    cells_id <- generate_first_step_ids(scdata_list)

    for (sample_id in c("123abc", "123def")){
      out <- filter_emptydrops(scdata_list, config, sample_id, cells_id)
      expect_equal(ncol(out$data[[sample_id]]), 40)
      expect_equal(length(out$new_ids[[sample_id]]), 30)
    }
  }
})

test_that("filter_emptydrops is sample aware", {
  scdata_list <- mock_scdata_list()
  config <- mock_config()
  cells_id <- generate_first_step_ids(scdata_list)

  # NA (empty) drops are in 123def only
  scdata_list[["123abc"]]$emptyDrops_FDR <- 0.009
  scdata_list[["123def"]]$emptyDrops_FDR[1:10] <- NA
  
  out <- filter_emptydrops(scdata_list, config, "123def", cells_id)
  expect_equal(ncol(out$data[["123abc"]]), 40)
  expect_equal(length(out$new_ids[["123abc"]]), 40)
  expect_equal(length(out$new_ids[["123def"]]), 30)
})


test_that("if FDR=1 filter_emptydrops keeps everything", {
  scdata_list <- mock_scdata_list()
  config <- mock_config()
  config$filterSettings$FDR <- 1
  cells_id <- generate_first_step_ids(scdata_list)

  out <- filter_emptydrops(scdata_list, config, "123def", cells_id)
  expect_equal(ncol(out$data[["123def"]]), 40)
  expect_equal(length(out$new_ids[["123abc"]]), 40)
  expect_equal(length(out$new_ids[["123def"]]), 40)
})

test_that("filter_emptydrops can be disabled", {
  scdata_list <- mock_scdata_list()
  config <- mock_config()
  config$enabled <- FALSE
  cells_id <- generate_first_step_ids(scdata_list)

  out <- filter_emptydrops(scdata_list, config, "123def", cells_id)
  expect_equal(ncol(out$data[["123def"]]), 40)
  expect_equal(length(out$new_ids[["123abc"]]), 40)
  expect_equal(length(out$new_ids[["123def"]]), 40)
})


test_that("filter_emptydrops handles missing emptyDrops_FDR", {
  scdata_list <- mock_scdata_list()
  config <- mock_config()
  cells_id <- generate_first_step_ids(scdata_list)

  # remove emptyDrops_FDR field
  scdata_list[["123abc"]]$emptyDrops_FDR <- NULL
  scdata_list[["123def"]]$emptyDrops_FDR <- NULL

  out <- filter_emptydrops(scdata_list, config, "123def", cells_id)
  expect_equal(ncol(out$data[["123def"]]), 40)
  expect_equal(length(out$new_ids[["123abc"]]), 40)
  expect_equal(length(out$new_ids[["123def"]]), 40)
})
