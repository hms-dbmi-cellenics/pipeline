source("qc_mock.R")

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


mock_scdata_list <- function() {
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )

  sample_ids <- c("sample_1", "sample_2")
  scdata_list <- {}
  for (i in seq_along(sample_ids)) {
    message(i)
    scdata <- Seurat::CreateSeuratObject(counts = pbmc_raw)
    n_cells <- ncol(scdata)
    start_idx <- (i-1)*n_cells
    end_idx <- ((i)*n_cells) - 1

    scdata$cells_id <- start_idx:end_idx

    # add empty drops stuff
    scdata$emptyDrops_FDR <- NA
    scdata$emptyDrops_FDR[0:70] <- 0.009

    # add samples
    scdata$samples <- rep(sample_ids[i], each = n_cells)
    # scdata_list <- Seurat::RenameCells(scdata, paste(scdata$samples, colnames(scdata), sep = ""))
    scdata_list[[sample_ids[i]]] <- scdata
  }
  return(scdata_list)
}

test_that("Initial cells_id are correct", {
  scdata_list <- mock_scdata_list()
  cells_id <- generate_first_step_ids(scdata_list)
  expect_equal(unique(cells_id[["sample_1"]]), 0:79)
  expect_equal(unique(cells_id[["sample_2"]]), 80:159)
})

test_that("filter_emptydrops removes NA with threshold < 1", {
  scdata_list <- mock_scdata_list()
  config <- mock_config()
  cells_id <- generate_first_step_ids(scdata_list)


  for (sample_id in c("sample_1", "sample_2")){
    out <- filter_emptydrops(scdata_list, config, sample_id, cells_id)
    expect_equal(ncol(out$data[[sample_id]]), 80)
    expect_equal(length(out$new_ids[[sample_id]]), 70)
  }
})

test_that("filter_emptydrops is sample aware", {
  scdata_list <- mock_scdata_list()
  config <- mock_config()
  cells_id <- generate_first_step_ids(scdata_list)

  # NA (empty) drops are in sample_2 only
  scdata_list[["sample_1"]]$emptyDrops_FDR <- 0.009
  out <- filter_emptydrops(scdata_list, config, "sample_2", cells_id)
  expect_equal(ncol(out$data[["sample_1"]]), 80)
  expect_equal(length(out$new_ids[["sample_1"]]), 80)
  expect_equal(length(out$new_ids[["sample_2"]]), 70)
})


test_that("if FDR=1 filter_emptydrops keeps everything", {
  scdata_list <- mock_scdata_list()
  config <- mock_config()
  config$filterSettings$FDR <- 1
  cells_id <- generate_first_step_ids(scdata_list)

  out <- filter_emptydrops(scdata_list, config, "sample_2", cells_id)
  expect_equal(ncol(out$data[["sample_2"]]), 80)
  expect_equal(length(out$new_ids[["sample_1"]]), 80)
  expect_equal(length(out$new_ids[["sample_2"]]), 80)
})

test_that("filter_emptydrops can be disabled", {
  scdata_list <- mock_scdata_list()
  config <- mock_config()
  config$enabled <- FALSE
  cells_id <- generate_first_step_ids(scdata_list)

  out <- filter_emptydrops(scdata_list, config, "sample_2", cells_id)
  expect_equal(ncol(out$data[["sample_2"]]), 80)
  expect_equal(length(out$new_ids[["sample_1"]]), 80)
  expect_equal(length(out$new_ids[["sample_2"]]), 80)
})


test_that("filter_emptydrops handles missing emptyDrops_FDR", {
  scdata_list <- mock_scdata_list()
  config <- mock_config()
  cells_id <- generate_first_step_ids(scdata_list)

  # remove emptyDrops_FDR field
  scdata_list[["sample_1"]]$emptyDrops_FDR <- NULL
  scdata_list[["sample_2"]]$emptyDrops_FDR <- NULL

  out <- filter_emptydrops(scdata_list, config, "sample_2", cells_id)
  expect_equal(ncol(out$data[["sample_2"]]), 80)
  expect_equal(length(out$new_ids[["sample_1"]]), 80)
  expect_equal(length(out$new_ids[["sample_2"]]), 80)
})
