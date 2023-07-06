mock_source_data <- function(sample_ids, from_experiment_id) {
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )

  scdata_list <- list()
  for (i in seq_along(sample_ids)) {
    message(i)
    scdata <- Seurat::CreateSeuratObject(counts = pbmc_raw)
    n_cells <- ncol(scdata)
    start_idx <- (i-1)*n_cells
    end_idx <- ((i)*n_cells) - 1

    scdata$cells_id <- start_idx:end_idx


    scdata@meta.data$samples <- rep(sample_ids[i], each = n_cells)
    scdata@misc$experimentId <- from_experiment_id

    # add samples
    scdata$samples <- rep(sample_ids[i], each = n_cells)
    scdata_list[i] <- scdata
  }
  return(scdata_list)
}

get_mock_processed_scdata <- function(sample_ids) {
  data("pbmc_small", package = "SeuratObject", envir = environment())
  pbmc_small$cells_id <- 0:(ncol(pbmc_small) - 1)
  pbmc_small@misc$gene_annotations <- data.frame(
    input = paste0("ENSG", seq_len(nrow(pbmc_small))),
    name = row.names(pbmc_small),
    row.names = paste0("ENSG", seq_len(nrow(pbmc_small)))
  )

  pbmc_small$samples <- rep(sample_ids, each = ncol(pbmc_small) / 2)
  return(pbmc_small)
}

get_mock_params <- function() {
  from_sample_ids <- list("sample-id-1", "sample-id-2")
  to_sample_ids <- list("cloned-sample-id-1", "cloned-sample-id-2")

  from_experiment_id <- "mock-from-experiment-id"
  to_experiment_id <- "mock-to-experiment-id"

  sample_ids_map <- to_sample_ids
  names(sample_ids_map) <- from_sample_ids


  input <- list(
    fromExperimentId = from_experiment_id,
    toExperimentId = to_experiment_id,
    sampleIdsMap = sample_ids_map
  )

  pipeline_config <- list(aws_config = "")

  return(list(input = input, pipeline_config = pipeline_config))
}

get_mock_cell_sets <- function(flatten = TRUE) {
  cell_sets_path <- file.path(setup_test_paths()$mock_data, "cell_sets")
  path <- file.path(cell_sets_path, "cell_sets_2_samples.json")
  cell_sets <- jsonlite::fromJSON(path, flatten = flatten)

  return(cell_sets)
}

test_that("copy_s3_objects calls the correct functions", {
  params <- get_mock_params()
  input <- params$input
  pipeline_config <- params$pipeline_config

  mock_copy_processed_rds <- mockery::mock()
  mock_copy_cell_sets <- mockery::mock()
  mock_copy_filtered_cells <- mockery::mock()
  mock_copy_source <- mockery::mock()

  mockery::stub(copy_s3_objects, "copy_processed_rds", mock_copy_processed_rds)
  mockery::stub(copy_s3_objects, "copy_cell_sets", mock_copy_cell_sets)
  mockery::stub(copy_s3_objects, "copy_filtered_cells", mock_copy_filtered_cells)
  mockery::stub(copy_s3_objects, "copy_source", mock_copy_source)

  result <- copy_s3_objects(input, pipeline_config)

  expect_called(mock_copy_processed_rds, 1)
  expect_called(mock_copy_cell_sets, 1)
  expect_called(mock_copy_filtered_cells, 1)
  expect_called(mock_copy_source, 1)
})

test_that("copy_processed_rds works correctly", {
  params <- get_mock_params()
  pipeline_config <- params$pipelineConfig
  from_experiment_id <- params$input$fromExperimentId
  to_experiment_id <- params$input$toExperimentId
  sample_ids_map <- params$input$sampleIdsMap

  processed_scdata <- get_mock_processed_scdata(names(sample_ids_map))

  mock_load_processed_scdata <- mock(processed_scdata)
  mock_upload_matrix_to_s3 <- mock()
  mockery::stub(copy_processed_rds, "load_processed_scdata", mock_load_processed_scdata)
  mockery::stub(copy_processed_rds, "upload_matrix_to_s3", mock_upload_matrix_to_s3)

  copy_processed_rds(from_experiment_id, to_experiment_id, sample_ids_map, pipeline_config)

  expect_called(mock_load_processed_scdata, 1)
  expect_called(mock_upload_matrix_to_s3, 1)
})

test_that("copy_cell_sets works correctly", {
  params <- get_mock_params()
  pipeline_config <- params$pipelineConfig
  from_experiment_id <- params$input$fromExperimentId
  to_experiment_id <- params$input$toExperimentId
  sample_ids_map <- params$input$sampleIdsMap

  cell_sets <- get_mock_cell_sets()

  mock_load_cellsets <- mock(cell_sets)
  mock_put_object_in_s3 <- mock()
  mockery::stub(copy_cell_sets, "load_cellsets", mock_load_cellsets)
  mockery::stub(copy_cell_sets, "put_object_in_s3", mock_put_object_in_s3)

  copy_cell_sets(from_experiment_id, to_experiment_id, sample_ids_map, pipeline_config)

  expect_called(mock_load_cellsets, 1)
  expect_called(mock_put_object_in_s3, 1)
})

test_that("copy_filtered_cells works correctly", {
  params <- get_mock_params()
  pipeline_config <- params$pipelineConfig
  from_experiment_id <- params$input$fromExperimentId
  to_experiment_id <- params$input$toExperimentId
  sample_ids_map <- params$input$sampleIdsMap


  filtered_cells_rds <- list("sample-id-1" = c(0, 1), "sample-id-2" = c(60, 61))

  mock_get_s3_rds <- mock(filtered_cells_rds, cycle = TRUE)
  mock_put_s3_rds <- mock(cycle = TRUE)
  mockery::stub(copy_filtered_cells, "get_s3_rds", mock_get_s3_rds)
  mockery::stub(copy_filtered_cells, "put_s3_rds", mock_put_s3_rds)

  copy_filtered_cells(from_experiment_id, to_experiment_id, sample_ids_map, pipeline_config)

  # expect 10 calls: #filter_steps(5) * #samples (2) = 10
  expect_called(mock_get_s3_rds, 10)
  # expect 10 calls: #filter_steps(5) * #samples (2) = 10
  expect_called(mock_put_s3_rds, 10)
})

test_that("copy_source works correctly", {
  params <- get_mock_params()
  pipeline_config <- params$pipelineConfig
  from_experiment_id <- params$input$fromExperimentId
  to_experiment_id <- params$input$toExperimentId
  sample_ids_map <- params$input$sampleIdsMap

  source_rds_list <- mock_source_data(names(sample_ids_map), from_experiment_id)

  # mockery needs us to set each of the responses in different params and in order
  mock_get_s3_rds <- mock(source_rds_list[[1]], source_rds_list[[2]], cycle = FALSE)
  mock_put_s3_rds <- mock(cycle = TRUE)
  mockery::stub(copy_source, "get_s3_rds", mock_get_s3_rds)
  mockery::stub(copy_source, "put_s3_rds", mock_put_s3_rds)

  copy_source(from_experiment_id, to_experiment_id, sample_ids_map, pipeline_config)

  # expect 2 calls, 1 for each sample
  expect_called(mock_get_s3_rds, 2)
  # expect 2 calls, 1 for each sample
  expect_called(mock_put_s3_rds, 2)
})
