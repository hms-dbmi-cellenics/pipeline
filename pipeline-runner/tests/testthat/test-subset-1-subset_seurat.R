# required to correctly source SeuratObject dumped R files
library(Seurat)

mock_scdata_list <- function(samples = rep("mock_sample_1_id", 80)) {
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )

  scdata <- Seurat::CreateSeuratObject(counts = pbmc_raw)
  # add samples
  scdata$samples <- samples
  scdata <- Seurat::RenameCells(scdata, paste(scdata$samples, colnames(scdata), sep = ""))

  # add doublet scores
  scdata$doublet_scores <- rep(c(0.01, 0.9), each = 40)
  scdata$doublet_class <- rep(c("singlet", "doublet"), each = 40)

  # add mitochondrial percent
  scdata$percent.mt <- rnorm(ncol(scdata), mean = 6)

  # create an scdata_list with duplicated samples
  scdata_list <- list()
  for (sample_id in scdata$samples) {
    scdata_list[[sample_id]] <- scdata
  }
  return(scdata_list)
}

mock_input <- function(parent_experiment_id, cellset_keys, samples = rep("mock_sample_1_id", 80), sample_ids = c("mock_sample_1_id")) {
  parentProcessingConfig <- construct_qc_config(mock_scdata_list(samples), unfiltered_samples = sample_ids)

  list(
    parentExperimentId = parent_experiment_id,
    experimentId = "mock_subset_experiment_id",
    cellSetKeys = cellset_keys,
    parentProcessingConfig = parentProcessingConfig
  )
}


mock_parent_experiment_data <- function(input, pipeline_config = NULL) {
  parent_experiment_id <- input$parentExperimentId
  paths <- setup_test_paths()
  qc_snaps_path <- file.path(paths$snaps, "qc")

  parent_scdata_name <- paste0(parent_experiment_id, sep = "-", "integrated_scdata", ".R")
  parent_cluster_cellsets_path <- file.path(qc_snaps_path, "cluster_cell_sets.json")

  # load the final test-qc output. it was named `snap_list`
  source(file.path(qc_snaps_path, parent_scdata_name))

  if (parent_experiment_id != snap_list$data@misc$experimentId) {
    stop("experiment id does not match object sourced from gem2s test output")
  }

  parent_cluster_cellsets <-
    list(cellSets = jsonlite::fromJSON(parent_cluster_cellsets_path, flatten = T))

  # set correct structure, since this is an incomplete cellsets file
  parent_cluster_cellsets$cellSets$children <- list(parent_cluster_cellsets$cellSets$children)

  parent_cluster_cellsets <- parse_cellsets(parent_cluster_cellsets)

   return(list(scdata = snap_list$data,
              cellsets = parent_cluster_cellsets))
}


stub_UUID_generate <- function(n) {
  paste0("this-is-not-a-uuid-", 1:n)
}


stubbed_subset_seurat <-
  function(input, pipeline_config, prev_out = NULL) {

    mockery::stub(subset_seurat,
                  "UUIDgenerate",
                  stub_UUID_generate,
                  depth = 2)

    mockery::stub(subset_seurat,
                  "load_parent_experiment_data",
                  mock_parent_experiment_data)

    res <- subset_seurat(input, pipeline_config, prev_out)

    return(res)
  }


test_that("subset_experiment correctly subsets cellsets", {
  parent_experiment_id <- "mock_experiment_id"
  cellset_keys = c("louvain-0", "louvain-1")
  input <- mock_input(parent_experiment_id, cellset_keys)
  parent_data <- mock_parent_experiment_data(input)

  res <- subset_experiment(input, parent_data)

  expect_equal(ncol(res), nrow(parent_data$cellsets[key %in% cellset_keys]))
  expect_setequal(res$cells_id, parent_data$cellsets[key %in% cellset_keys, cell_id])

})


test_that("create_sample_id_map generates a list of the correct size and format", {
  parent_experiment_id <- "mock_experiment_id"
  cellset_keys = c("louvain-0", "louvain-1")
  input <- mock_input(parent_experiment_id, cellset_keys)
  parent_data <- mock_parent_experiment_data(input)

  subset_scdata <- subset_experiment(input, parent_data)

  res <- create_sample_id_map(unique(subset_scdata$samples))

  expect_equal(length(res), length(unique(subset_scdata$samples)))
  expect_setequal(names(res), unique(subset_scdata$samples))

})


test_that("add_new_sample_ids correctly assigns new sample_ids to each cell in a subset experiment", {
  parent_experiment_id <- "mock_experiment_id"
  cellset_keys = c("louvain-0", "louvain-1")
  input <- mock_input(parent_experiment_id, cellset_keys)
  parent_data <- mock_parent_experiment_data(input)

  subset_scdata <- subset_experiment(input, parent_data)
  sample_id_map <- create_sample_id_map(unique(subset_scdata$samples))
  subset_scdata$parent_samples <- subset_scdata$samples

  res <- add_new_sample_ids(subset_scdata, sample_id_map)

  # convert to dt and join to get sample_id vector
  sample_id_map_dt <- data.table::data.table(parent_sample_ids = names(sample_id_map),
                                             sample_ids = unname(unlist(sample_id_map)))

  expected_sample_ids <- sample_id_map_dt[subset_scdata$samples, sample_ids, on = "parent_sample_ids"]

  expect_equal(unname(res$samples), expected_sample_ids)

})


test_that("add_subset_metadata adds required metadata to a subset seurat object", {
  parent_experiment_id <- "mock_experiment_id"
  cellset_keys <-  c("louvain-0", "louvain-1")
  input <- mock_input(parent_experiment_id, cellset_keys)
  parent_data <- mock_parent_experiment_data(input)

  subset_scdata <- subset_experiment(input, parent_data)
  sample_id_map <- create_sample_id_map(unique(subset_scdata$samples))

  res <- add_subset_metadata(input, subset_scdata, sample_id_map)

  expect_true(all(c("parent_samples", "samples") %in% names(res@meta.data)))
  expect_false(setequal(res@meta.data$parent_samples, res@meta.data$samples))
  expect_equal(res@misc$experimentId, input$experimentId)

})


test_that("subset_seurat matches snapshot", {
  parent_experiment_id <- "mock_experiment_id"
  cellset_keys <-  c("louvain-0", "louvain-1")
  input <- mock_input(parent_experiment_id, cellset_keys)
  pipeline_config <- list()

  res <- stubbed_subset_seurat(input, pipeline_config)

  expect_snapshot(res)
})


test_that("filter_low_cell_samples removes samples with cells below the threshold", {
  parent_experiment_id <- "mock_experiment_id"
  cellset_keys = c("louvain-0", "louvain-1")
  input <- mock_input(parent_experiment_id, cellset_keys)
  parent_data <- mock_parent_experiment_data(input)

  min_cells <- 300

  res <- filter_low_cell_samples(parent_data$scdata, min_cells = min_cells)
  expect_false(all(table(parent_data$scdata$samples) > min_cells))
  expect_true(all(table(res$samples) > min_cells))
})

test_that("generate_subset_config works correctly", {
  parent_sample_ids <- c("sample-id-1", "sample-id-2", "sample-id-3", "sample-id-4")
  scdata_list <- mock_scdata_list(samples = rep(parent_sample_ids, 20))

  parent_processing_config <- construct_qc_config(scdata_list, unfiltered_samples = parent_sample_ids)

  # Make some of the configs unique to each sample
  # so we can check that the translation preserves the configs
  # correctly assigned to each sample
  parent_processing_config$cellSizeDistribution$`sample-id-2`$filterSettings$binStep <- 300
  parent_processing_config$cellSizeDistribution$`sample-id-3`$filterSettings$binStep <- 400
  parent_processing_config$mitochondrialContent$`sample-id-1`$filterSettings$absoluteThreshold$maxFraction <- 0.2
  parent_processing_config$mitochondrialContent$`sample-id-3`$filterSettings$absoluteThreshold$maxFraction <- 0.5

  subset_sample_ids <- c("sample-subset-id-1", "sample-subset-id-2", "sample-subset-id-3", "sample-subset-id-4")

  sample_ids_map <- as.list(subset_sample_ids)
  names(sample_ids_map) <- parent_sample_ids

  subset_processing_config <- generate_subset_config(parent_processing_config, sample_ids_map)

  expect_snapshot(subset_processing_config)
})
