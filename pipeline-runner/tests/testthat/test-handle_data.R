mock_publish <- NULL
mock_sns <- NULL

before_each <- function() {
  mock_publish <<- mockery::mock(
    list(MessageId = "ok"),
    cycle = TRUE
  )

  mock_sns <<- function(config) {
    return(list(publish = mock_publish))
  }
}

get_gem2s_mock_cellsets <- function() {
  # get a snapshot cellsets json
  paths <- setup_test_paths()
  jsonlite::fromJSON(
    file.path(
      paths$snaps,
      "gem2s",
      "gem2s-7-mock_experiment_id-cellsets.json"
    ),
    flatten = TRUE
  )
}


get_mock_cell_sets <- function(flatten = TRUE) {
  cell_sets_path <- file.path(setup_test_paths()$mock_data, "cell_sets")
  path <- file.path(cell_sets_path, "cell_sets_2_samples.json")
  cell_sets <- jsonlite::fromJSON(path, flatten = flatten)

  return(cell_sets)
}

mock_counts <- function(ncells = 200, ngenes = 200) {
  sce <- scuttle::mockSCE(ncells, ngenes)
  counts <- SummarizedExperiment::assay(sce, "counts")
  rownames(counts) <- gsub("_", "-", rownames(counts))
  counts <- as(counts, "dgCMatrix")
  maybe_bpcells(counts)
}


test_that("send_pipeline_update_to_api completes successfully", {
  before_each()

  pipeline_config <- list(
    sns_topic = "ExampleTopic",
    aws_config = NULL
  )

  mockery::stub(send_pipeline_update_to_api, "paws::sns", mock_sns)

  response <- send_pipeline_update_to_api(pipeline_config,
    experiment_id = "dfgdfg",
    task_name = "dsfdsdf",
    data = 1:5,
    input = list(auth_JWT = "ayylmao"),
    string_value = "GEM2SResponse"
  )


  expect_true(response == "ok")
})

stub_put_object_in_s3_multipart <- function(
  pipeline_config, bucket, object, key
) {
  NULL
}

test_that("upload_matrix_to_s3 completes successfully", {
  before_each()

  # mock things
  data <- matrix()
  pipeline_config <- list(processed_bucket = "processed-bucket")
  experiment_id <- "1234"

  mockery::stub(
    upload_matrix_to_s3,
    "put_object_in_s3_multipart",
    stub_put_object_in_s3_multipart
  )

  # generates correct S3 key
  key <- upload_matrix_to_s3(pipeline_config, experiment_id, data)
  expect_equal(key, "1234/r.qs")
})

test_that("send_output_to_api completes successfully", {
  before_each()

  c(pipeline_config, input, plot_data_keys, output) %<-%
    send_output_to_api_mock_data

  mockery::stub(send_output_to_api, "put_object_in_s3", NULL)
  mockery::stub(send_output_to_api, "paws::sns", mock_sns)

  response <- send_output_to_api(pipeline_config, input, plot_data_keys, output)

  expect_true(response == "ok")
})


test_that("safe_cbind returns empty data.table when binding an empty data.table with a vector", {
  before_each()

  dt_empty <- data.table::data.table()
  col <- c(a_col = "a_value")

  res <- safe_cbind(dt_empty, col)

  expect_identical(res, dt_empty)
})


test_that("safe_cbind adds a column to a non-empty data.table", {
  before_each()

  dt <- data.table::data.table(col1 = 1:10, col2 = 11:20)
  values <- seq(1, 20, 2)
  res <- safe_cbind(dt, bound_col = values)

  expect_identical(res[, bound_col], values)
  expect_equal(ncol(res), ncol(dt) + 1)
})


test_that("safe_cbind names bound column as expected", {
  before_each()

  dt <- data.table::data.table(col1 = 1:10, col2 = 11:20)
  values <- seq(1, 20, 2)
  res <- safe_cbind(dt, my_expected_column_name = values)

  expect_true("my_expected_column_name" %in% names(res))
  expect_identical(res[, my_expected_column_name], values)
})


test_that("safe_cbind binds more than one column and names accordingly", {
  before_each()

  dt <- data.table::data.table(col1 = 1:10, col2 = 11:20)
  values_1 <- seq(1, 20, 2)
  values_2 <- values_1 + 2

  res <- safe_cbind(
    dt,
    an_interesting_variable = values_1,
    an_interesting_variable_plus_2 = values_2
  )

  expect_true("an_interesting_variable" %in% names(res))
  expect_identical(res[, an_interesting_variable], values_1)

  expect_true("an_interesting_variable_plus_2" %in% names(res))
  expect_identical(res[, an_interesting_variable_plus_2], values_2)
})


test_that("cbind_cellset_type names the bound column correctly", {
  before_each()

  dt <- data.table::data.table(col1 = 1:10, col2 = 11:20)
  values <- seq(1, 20, 2)

  res <- cbind_cellset_type(dt, values)

  expect_true("cellset_type" %in% names(res))
  expect_identical(res[, cellset_type], values)
})


test_that("parse_cellsets parses a cellset object", {
  before_each()

  cellsets <- get_gem2s_mock_cellsets()

  res <- parse_cellsets(cellsets)

  expect_s3_class(res, "data.table")
  expect_identical(names(res), c("key", "name", "type", "cell_id"))
})

stub_s3_put_object <- function(Bucket, Key, Body, Tagging) {
  response <- list(
    Expiration = character(),
    ETag = "this_is_not_an_etag",
    ServerSideEncryption = character(),
    VersionId = character(),
    SSECustomerAlgorithm = character(),
    SSECustomerKeyMD5 = character(),
    SSEKMSKeyId = character(),
    SSEKMSEncryptionContext = character(),
    BucketKeyEnabled = logical(),
    RequestCharged = character()
  )

  return(response)
}

test_that("put_object_in_s3 works", {
  before_each()

  mockery::stub(put_object_in_s3, "s3$put_object", stub_s3_put_object)

  pipeline_config <- mock_pipeline_config()
  bucket <- "mock_bucket"
  key <- "mock_key"
  object <- "something"
  key <- "a_key"

  expect_message(put_object_in_s3(pipeline_config, bucket, object, key),
    regexp = "a_key"
  )
})


test_that("put_object_in_s3 retries if s3$put_object throws an error", {
  before_each()
  mockery::stub(
    put_object_in_s3,
    "s3$put_object",
    mockery::mock(stop("an error"), stub_s3_put_object)
  )

  pipeline_config <- mock_pipeline_config()
  bucket <- "mock_bucket"
  key <- "mock_key"
  object <- "something"
  key <- "a_key"

  expect_message(put_object_in_s3(pipeline_config, bucket, object, key),
    regexp = ".*Retrying \\(1/2\\).*"
  )
})

test_that("is_uuid detects uuids correctly", {
  before_each()

  expect_true(is_uuid(uuid::UUIDgenerate()))
  expect_false(is_uuid("not-a-uuid"))
})


test_that("get_cellset_types correctly gets cellset types", {
  before_each()

  key <-
    c(
      "louvain",
      "scratchpad",
      "sample",
      "metadata_1",
      "metadata_2",
      "faa74528-c4d7-11ed-9fda-0242ac120003"
    )

  type <-
    c(
      "cellSets",
      "cellSets",
      "metadataCategorical",
      "metadataCategorical",
      "metadataCategorical",
      "cellSets"
    )

  expected_cellset_types <- c(
    "cluster",
    "scratchpad",
    "sample",
    "metadata",
    "metadata",
    "sctype"
  )

  expect_identical(
    purrr::map2_chr(key, type, get_cellset_type),
    expected_cellset_types
  )
})

test_that("send_pipeline_fail_update handles a gem2s call successefully", {
  before_each()
  mockery::stub(send_pipeline_fail_update, "paws::sns", mock_sns)

  pipeline_config <- list(
    aws_config = list(),
    results_bucket = "test_bucket",
    api_url = "test_url",
    sns_topic = "test_topic"
  )
  input <- list(
    processName = "gem2s",
    experimentId = "test_experiment",
    taskName = "test_task"
  )
  error_message <- "test_error"

  result <- send_pipeline_fail_update(pipeline_config, input, error_message)

  # completes successfully
  expect_equal(result, list(MessageId = "ok"))

  # Check that sns$publish was called with the correct parameters
  expect_snapshot(str(mockery::mock_args(mock_publish)))
})

test_that("send_pipeline_fail_update handles a qc call successfully with no global_env config", {
  before_each()
  mockery::stub(send_pipeline_fail_update, "paws::sns", mock_sns)
  mockery::stub(send_pipeline_fail_update, "ids::uuid", "mock-uuid")

  pipeline_config <- list(
    aws_config = list(),
    results_bucket = "test_bucket",
    api_url = "test_url",
    sns_topic = "test_topic"
  )
  input <- list(
    processName = "qc"
  )
  error_message <- "test_error"

  result <- send_pipeline_fail_update(pipeline_config, input, error_message)

  # completes successfully
  expect_equal(result, list(MessageId = "ok"))

  # Check that sns$publish was called with the correct parameters
  expect_snapshot(mockery::mock_args(mock_publish))
})

test_that("send_pipeline_fail_update handles a qc call successfully with global_env config", {
  before_each()
  mock_put_object_in_s3 <- mockery::mock(NULL, cycle = TRUE)

  mockery::stub(send_pipeline_fail_update, "paws::sns", mock_sns)
  mockery::stub(send_pipeline_fail_update, "ids::uuid", "mock-uuid")
  mockery::stub(
    send_pipeline_fail_update, "put_object_in_s3", mock_put_object_in_s3
  )

  c(pipeline_config, input, plot_data_keys, output) %<-%
    send_output_to_api_mock_data
  c(config, plot_data = plotData) %<-% output

  pipeline_config <- list(
    aws_config = list(),
    results_bucket = "test_bucket",
    api_url = "test_url",
    sns_topic = "test_topic"
  )
  input <- list(
    processName = "qc",
    sampleUuid = "00000000-0000-0000-000000000111",
    taskName = "mitochondrialContent"
  )
  error_message <- "test_error"

  # Set up a global env config
  config_key <- paste0("config-", input$taskName, "-", input$sampleUuid)
  assign(config_key, config, envir = globalenv())

  result <- send_pipeline_fail_update(pipeline_config, input, error_message)

  # completes successfully
  expect_equal(result, list(MessageId = "ok"))

  # Check that sns$publish was called with the correct parameters
  expect_snapshot(mockery::mock_args(mock_publish))

  # Check that put_object_in_s3 was called with the correct parameters
  expect_snapshot(mockery::mock_args(mock_put_object_in_s3))
})

test_that("unflatten_cell_sets works", {
  # Get flattened cell sets
  flattened_cell_sets <- get_mock_cell_sets()

  # Unflatten them
  cell_sets <- unflatten_cell_sets(flattened_cell_sets$cellSets)

  # They are unflattened
  expect_snapshot(cell_sets)
})

# Tests for handle_data helper functions

test_that("get_nnzero returns sum of nFeature_RNA", {
  # Create a mock Seurat object
  set.seed(RANDOM_SEED)
  counts <- mock_counts()

  # Create a Seurat object
  scdata <- Seurat::CreateSeuratObject(counts = counts)

  # Call get_nnzero
  result <- get_nnzero(scdata)

  expect_true(is.numeric(result))
  expect_true(result > 0)
  expect_equal(as.integer(result), sum(scdata$nFeature_RNA))
})

test_that("get_nnzero with Seurat object backed by BPCells matrices", {
  # Create a Seurat object using BPCells matrix storage
  set.seed(RANDOM_SEED)
  counts <- mock_counts()

  # Ensure counts are integer
  counts <- as(counts, "dgCMatrix")
  mode(counts@x) <- "integer"

  # Write counts to disk and load as BPCells matrix
  matrix_dir <- file.path(withr::local_tempdir(),
    paste0("nnzero_test_", sample(100000, 1))
  )

  bpcells_counts <- BPCells::write_matrix_dir(counts, dir = matrix_dir)

  # Create Seurat object with BPCells-backed counts
  scdata <- Seurat::CreateSeuratObject(counts = bpcells_counts)

  # Call get_nnzero
  result <- get_nnzero(scdata)

  expect_type(result, "double")
  expect_true(result > 0)
})

test_that("order_by_size orders scdata_list by size", {
  # Create multiple mock Seurat objects of different sizes
  set.seed(RANDOM_SEED)
  counts1 <- mock_counts(ncells = 100)
  counts2 <- mock_counts(ncells = 200)
  counts3 <- mock_counts(ncells = 150)

  scdata1 <- Seurat::CreateSeuratObject(counts = counts1)
  scdata2 <- Seurat::CreateSeuratObject(counts = counts2)
  scdata3 <- Seurat::CreateSeuratObject(counts = counts3)

  scdata_list <- list(
    sample1 = scdata2,
    sample2 = scdata1,
    sample3 = scdata3
  )

  # Call order_by_size
  ordered_list <- order_by_size(scdata_list)

  # Check that they are ordered by size
  sizes <- sapply(ordered_list, get_nnzero)
  expect_true(sizes[1] <= sizes[2])
  expect_true(sizes[2] <= sizes[3])
})

test_that("group_files_by_sample groups files correctly", {
  # Create mock experiment files
  experiment_files <- list(
    list(Key = "exp1/sample1/file1.rds"),
    list(Key = "exp1/sample1/file2.txt"),
    list(Key = "exp1/sample2/file1.rds"),
    list(Key = "exp1/sample2/file2.txt")
  )

  # Call group_files_by_sample
  result <- group_files_by_sample(experiment_files)

  # Check the structure
  expect_true("sample1" %in% names(result))
  expect_true("sample2" %in% names(result))
  expect_equal(length(result$sample1), 2)
  expect_equal(length(result$sample2), 2)
  expect_equal(result$sample1[["file1.rds"]], "exp1/sample1/file1.rds")
  expect_equal(result$sample2[["file2.txt"]], "exp1/sample2/file2.txt")
})

test_that("replace_matrix_dir_paths updates MatrixDir path", {
  # Create a mock BPCells matrix directory
  set.seed(42)
  counts <- Matrix::Matrix(
    sample(0:10, 100, replace = TRUE),
    nrow = 10,
    ncol = 10,
    sparse = TRUE
  )
  mode(counts@x) <- "integer"
  counts <- as(counts, "dgCMatrix")

  # Create initial matrix dir
  tmp_dir <- withr::local_tempdir()
  original_dir <- file.path(
    tmp_dir,
    paste0("test_matrix_", sample(100000, 1))
  )
  bpcells_matrix <- BPCells::write_matrix_dir(counts, dir = original_dir)

  # Create Seurat object with BPCells matrix
  seurat_obj <- SeuratObject::CreateSeuratObject(
    counts = bpcells_matrix
  )

  # Save RDS file
  rds_file <- withr::local_tempfile(fileext = ".rds")
  saveRDS(seurat_obj, rds_file)

  # Read RDS and verify matrix access works
  loaded_obj <- readRDS(rds_file)
  expect_equal(nrow(loaded_obj), 10)
  expect_equal(ncol(loaded_obj), 10)

  # Compute reference data before moving directory
  subset_before_data <- as.matrix(loaded_obj[["RNA"]]$counts[1:3, 1:3])
  col_sums_before <- Matrix::colSums(loaded_obj[["RNA"]]$counts[1:3, 1:3])

  # Move matrix directory to new location
  new_dir <- file.path(
    tmp_dir,
    paste0("test_matrix_new_", sample(100000, 1))
  )
  dir.create(new_dir)

  # Copy matrix dir contents to new location
  file.copy(
    list.files(original_dir, full.names = TRUE),
    new_dir,
    recursive = TRUE
  )

  # Delete the original matrix directory
  unlink(original_dir, recursive = TRUE)

  # Get the matrix object from loaded Seurat object
  loaded_matrix <- loaded_obj[["RNA"]]$counts

  # Verify that operations fail when matrix dir is in wrong location
  # (original location no longer exists)
  expect_error(
    Matrix::colSums(loaded_matrix[1:3, 1:3]),
    regexp = "Missing directory|cannot open"
  )

  # Update matrix dir paths using replace_matrix_dir_paths
  updated_matrix <- replace_matrix_dir_paths(loaded_matrix, new_dir)

  # Verify matrix access still works after path update
  expect_equal(nrow(updated_matrix), 10)
  expect_equal(ncol(updated_matrix), 10)
  subset_after <- updated_matrix[1:3, 1:3]
  expect_equal(nrow(subset_after), 3)
  expect_equal(ncol(subset_after), 3)

  # Verify matrix dimensions are preserved
  col_sums_after <- Matrix::colSums(subset_after)
  expect_equal(col_sums_before, col_sums_after)

  # Verify we can extract features from updated matrix
  feature_names <- rownames(updated_matrix)[1:3]
  expect_equal(length(feature_names), 3)
})

test_that("load_source_rds retrieves and deserializes RDS from S3", {
  # Create a test RDS object
  test_obj <- list(data = data.frame(x = 1:5, y = 6:10))
  rds_raw <- serialize(test_obj, NULL)

  # Mock S3 get_object response
  mock_s3_get <- mockery::mock(
    list(
      body = rds_raw,
      Metadata = NULL
    )
  )

  # Call function with mocked S3
  result <- load_source_rds(
    list(get_object = mock_s3_get),
    "test_key",
    "test_bucket"
  )

  # Verify result
  expect_type(result, "list")
  expect_equal(names(result), "data")
  expect_equal(nrow(result$data), 5)
  expect_equal(result$data$x, 1:5)

  # Verify S3 was called correctly
  expect_called(mock_s3_get, 1)
})

test_that("load_source_matrix_dir extracts matrix dir from S3", {
  # Create a temporary matrix directory
  set.seed(42)
  counts <- Matrix::Matrix(
    sample(0:10, 100, replace = TRUE),
    nrow = 10,
    ncol = 10,
    sparse = TRUE
  )
  mode(counts@x) <- "integer"
  counts <- as(counts, "dgCMatrix")

  # Create and tar the matrix directory
  temp_base <- withr::local_tempdir()
  sample_id <- "test_sample"
  original_dir <- file.path(temp_base, paste0(sample_id, "_matrix_dir"))
  BPCells::write_matrix_dir(counts, dir = original_dir)

  # Create tar.zst file
  tarfile <- file.path(temp_base, paste0(sample_id, "_matrix_dir.tar.zst"))

  # Create tar file
  current_dir <- getwd()
  setwd(temp_base)
  system(paste("tar --zstd -cf", tarfile, basename(original_dir)))
  setwd(current_dir)

  # Create mock S3 that writes tarfile when download_file is called
  mock_s3_func <- function(Bucket, Key, Filename) {
    file.copy(tarfile, Filename)
  }
  mock_s3 <- list(download_file = mock_s3_func)

  # Call function
  result <- load_source_matrix_dir(
    mock_s3,
    "test_key",
    "test_bucket",
    sample_id
  )

  # Verify result
  expect_true(dir.exists(result))
  expect_match(basename(result), paste0(sample_id, "_matrix_dir"))

})

test_that(
  "upload_matrix_dir_to_s3 uploads tarred matrix directory",
  {
    # Create a mock Seurat object with BPCells matrix
    set.seed(42)
    counts <- Matrix::Matrix(
      sample(0:10, 100, replace = TRUE),
      nrow = 10,
      ncol = 10,
      sparse = TRUE
    )
    mode(counts@x) <- "integer"
    counts <- as(counts, "dgCMatrix")

    # Create matrix directory
    matrix_dir <- file.path(
      withr::local_tempdir(),
      paste0("upload_test_", sample(100000, 1))
    )
    bpcells_matrix <- BPCells::write_matrix_dir(counts, dir = matrix_dir)

    # Create Seurat object (single object, not list)
    scdata <- SeuratObject::CreateSeuratObject(
      counts = bpcells_matrix
    )

    # Create mock function to capture upload info
    uploaded_info <- list()
    mock_put_func <- function(pipeline_config, bucket,
                              file_path, key) {
      uploaded_info$bucket <<- bucket
      uploaded_info$key <<- key
      uploaded_info$file <<- file_path
    }

    # Create pipeline config
    pipeline_config <- list(
      processed_bucket = "test_bucket",
      aws_config = list()
    )

    # Patch the function and call
    local_mocked_bindings(
      put_object_in_s3_multipart = mock_put_func,
      .package = "pipeline"
    )

    # Call function directly with a Seurat object
    # (Note: this tests with a single Seurat object)
    upload_matrix_dir_to_s3(
      pipeline_config,
      "test_exp_id",
      scdata
    )

    # Verify the upload was called with correct parameters
    expect_equal(uploaded_info$bucket, "test_bucket")
    expect_equal(
      uploaded_info$key,
      "test_exp_id/matrix_dir.tar.zst"
    )
  }
)
