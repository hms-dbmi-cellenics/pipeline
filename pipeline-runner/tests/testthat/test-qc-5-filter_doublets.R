mock_ids <- function() {
  # TODO: parametrize sample ids
  return(list("sample_1" = 0:39, "sample_2" = 40:79))
}


mock_config <- function(thr = 0.1, auto = FALSE, enabled = TRUE, recomputeDoubletScore = FALSE) {
  config <- list(
    auto = auto,
    enabled = enabled,
    filterSettings = list(
      probabilityThreshold = thr
    ),
    recomputeDoubletScore = recomputeDoubletScore
  )
  return(config)
}

mock_scdata <- function() {
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )

  sample_1_id <- "sample_1"
  sample_2_id <- "sample_2"

  scdata <- Seurat::CreateSeuratObject(counts = pbmc_raw)
  scdata$cells_id <- 0:(ncol(scdata) - 1)

  # add doublet scores and class
  scdata@meta.data$doublet_scores <- 0.01
  scdata@meta.data$doublet_scores[1:10] <- 0.9

  scdata@meta.data$doublet_class <- rep("singlet", 80)
  scdata@meta.data$doublet_class[1:10] <- rep("doublet", 10)

  # add samples
  scdata$samples <- rep(c(sample_1_id, sample_2_id), each = 40)
  scdata <- Seurat::RenameCells(scdata, paste(scdata$samples, colnames(scdata), sep = ""))

  scdata_sample1 <- subset(scdata, samples == sample_1_id)
  scdata_sample2 <- subset(scdata, samples == sample_2_id)

  scdata_list <- list(scdata_sample1, scdata_sample2)
  names(scdata_list) <- c(sample_1_id, sample_2_id)

  return(scdata_list)
}

test_that("filter_doublets filters based on threshold", {
  scdata_list <- mock_scdata()
  sample_1_id <- names(scdata_list)[1]
  sample_2_id <- names(scdata_list)[2]
  cells_id <- mock_ids()
  # should filter first 10 cells
  config <- mock_config(0.5)

  out <- filter_doublets(scdata_list, config, sample_1_id, cells_id)

  expect_equal(ncol(out$data[[sample_1_id]]), 40)
  expect_equal(out$new_ids[[sample_1_id]], 10:39)
})

test_that("filter_doublets is sample aware", {
  scdata_list <- mock_scdata()
  sample_1_id <- names(scdata_list)[1]
  sample_2_id <- names(scdata_list)[2]
  cells_id <- mock_ids()
  scdata_list[[sample_2_id]]@meta.data$doublet_scores[31:40] <- 0.9
  config <- mock_config(0.5)

  out <- filter_doublets(scdata_list, config, sample_1_id, cells_id)
  expect_equal(ncol(out$data[[sample_1_id]]), 40)
  expect_equal(out$new_ids[[sample_1_id]], 10:39)
  expect_equal(out$new_ids[[sample_2_id]], 40:79)

  out <- filter_doublets(out$data, config, sample_2_id, out$new_ids)
  expect_equal(ncol(out$data[[sample_2_id]]), 40)
  expect_equal(out$new_ids[[sample_1_id]], 10:39)
  expect_equal(out$new_ids[[sample_2_id]], 40:69)
})

test_that("filter_doublets filters works with auto", {
  scdata_list <- mock_scdata()
  sample_1_id <- names(scdata_list)[1]
  sample_2_id <- names(scdata_list)[2]
  cells_id <- mock_ids()
  # should filter first 10 cells
  config <- mock_config(0.001, auto = TRUE)
  out <- filter_doublets(scdata_list, config, sample_1_id, cells_id)
  expect_equal(out$new_ids[[sample_1_id]], 10:39)

  config <- mock_config(0.001, auto = FALSE)
  out <- filter_doublets(scdata_list, config, sample_1_id, cells_id)
  expect_equal(length(out$new_ids[[sample_1_id]]), 0)
})

test_that("filter_doublets can be disabled", {
  scdata_list <- mock_scdata()
  sample_1_id <- names(scdata_list)[1]
  sample_2_id <- names(scdata_list)[2]
  cells_id <- mock_ids()
  config <- mock_config(0.5, enabled = FALSE)
  out <- filter_doublets(scdata_list, config, sample_1_id, cells_id)
  expect_equal(out$new_ids[[sample_1_id]], 0:39)
  expect_equal(out$new_ids[[sample_2_id]], 40:79)

  config <- mock_config(0.5, enabled = TRUE)
  out <- filter_doublets(scdata_list, config, sample_1_id, cells_id)
  expect_equal(out$new_ids[[sample_1_id]], 10:39)
  expect_equal(out$new_ids[[sample_2_id]], 40:79)
})


test_that("generate_default_values_doubletScores sets threshold to 0 when there are no singlets", {
  scdata_list <- mock_scdata()
  sample_1_id <- names(scdata_list)[1]
  sample_2_id <- names(scdata_list)[2]
  scdata_list[[1]]$doublet_class <- "doublet"

  expect_equal(generate_default_values_doubletScores(scdata_list[[1]]), 0)
})


test_that("doublet scores are re-computed if the API says so", {

  scdata_list <- mock_scdata()
  sample_1_id <- names(scdata_list)[1]
  cells_id <- mock_ids()
  config <- mock_config(recomputeDoubletScore = TRUE)

  out <- suppressWarnings(filter_doublets(scdata_list, config, sample_1_id, cells_id))

  expect_false(
    identical(
      out$data$sample_1$doublet_scores,
      scdata_list$sample_1$doublet_scores
    )
  )
})
