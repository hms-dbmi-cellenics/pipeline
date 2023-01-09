library(Seurat)

mock_parent_experiment_data <- function(parent_experiment_id) {
  paths <- setup_test_paths()
  qc_snaps_path <- file.path(paths$snaps, "qc")

  parent_scdata_name <-
    paste0(paste(parent_experiment_id, "integrated_scdata", sep = "-"),
           ".R")
  parent_cluster_cellsets_path <-
    file.path(qc_snaps_path, "cluster_cell_sets.json")

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


mock_input <- function(parent_experiment_id, cellset_keys) {
  list(
    parentExperimentId = parent_experiment_id,
    experimentId = "mock_subset_experiment_id",
    cellSetKeys = cellset_keys
  )
}


stub_load_parent_experiment_data <- function(input, pipeline_config) {
  mock_parent_experiment_data(input$parentExperimentId)
}


stubbed_subset_seurat <-
  function(input, pipeline_config, prev_out = NULL) {
    mockery::stub(subset_seurat,
                  "load_parent_experiment_data",
                  stub_load_parent_experiment_data)

    res <- subset_seurat(input, pipeline_config, prev_out)

    return(res)
  }


test_that("subset_experiment correctly subsets cellsets", {
  parent_experiment_id <- "mock_experiment_id"
  parent_data <- mock_parent_experiment_data(parent_experiment_id)
  cellset_keys = c("louvain-0", "louvain-1")
  input <- mock_input(parent_experiment_id, cellset_keys)

  res <- subset_experiment(input, parent_data)

  expect_equal(ncol(res), nrow(parent_data$cellsets[key %in% cellset_keys]))
  expect_setequal(res$cells_id, parent_data$cellsets[key %in% cellset_keys, cell_id])

})


test_that("create_sample_id_map generates a list of the correct size and format", {
  parent_experiment_id <- "mock_experiment_id"
  parent_data <- mock_parent_experiment_data(parent_experiment_id)
  cellset_keys = c("louvain-0", "louvain-1")
  input <- mock_input(parent_experiment_id, cellset_keys)

  subset_scdata <- subset_experiment(input, parent_data)

  res <- create_sample_id_map(unique(subset_scdata$samples))

  expect_equal(length(res), length(unique(subset_scdata$samples)))
  expect_setequal(names(res), unique(subset_scdata$samples))

})


test_that("add_new_sample_ids correctly assigns new sample_ids to each cell in a subset experiment", {
  parent_experiment_id <- "mock_experiment_id"
  parent_data <- mock_parent_experiment_data(parent_experiment_id)
  cellset_keys = c("louvain-0", "louvain-1")
  input <- mock_input(parent_experiment_id, cellset_keys)

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
  parent_data <- mock_parent_experiment_data(parent_experiment_id)
  cellset_keys <-  c("louvain-0", "louvain-1")
  input <- mock_input(parent_experiment_id, cellset_keys)

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
  expect_snapshot(stubbed_subset_seurat(input, pipeline_config))

})
