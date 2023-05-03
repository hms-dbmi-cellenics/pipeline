mock_publish <- NULL
mock_sns <- NULL

before_each <- function() {
  mock_publish <<- mockery::mock(
    list(MessageId = 'ok'),
    cycle = TRUE
  )

  mock_sns <<- function(config) {
    return(list(publish = mock_publish))
  }
}

mock_cellsets <- function(){
  # get a snapshot cellsets json
  paths <- setup_test_paths()
  jsonlite::fromJSON(file.path(paths$snaps, "gem2s", "gem2s-7-mock_experiment_id-cellsets.json"), flatten = TRUE)

}

test_that("send_gem2s_update_to_api completes successfully", {
    before_each()

    pipeline_config <- list(
        sns_topic = 'ExampleTopic',
        aws_config = NULL
    )

    mockery::stub(send_gem2s_update_to_api, 'paws::sns', mock_sns)

    response <- send_gem2s_update_to_api(pipeline_config,
                            experiment_id = 'dfgdfg',
                            task_name = 'dsfdsdf',
                            data = 1:5,
                            input = list(auth_JWT='ayylmao'))


    expect_true(response == 'ok')
})

stub_put_object_in_s3_multipart <- function(pipeline_config, bucket, object, key) {
    return(NULL)
}

stub_put_object_in_s3 <- function(Bucket, Key, Body) {
  return(NULL)
}

test_that("upload_debug_folder_to_s3 completes successfully", {
    before_each()

    # create fake logs and dump file
    debug_path <- tempdir()
    debug_timestamp <- 'now'
    experiment_id <- '1234'
    pipeline_config <- list(debug_bucket = 'pipeline-debug-development')

    debug_prefix <- file.path(experiment_id, debug_timestamp)

    debug_dir <- file.path(debug_path, debug_prefix)
    dir.create(debug_dir, recursive = TRUE)
    file.create(file.path(debug_dir, c('logs.txt', 'dump.rda')))

    mockery::stub(upload_debug_folder_to_s3, 'put_object_in_s3_multipart', stub_put_object_in_s3_multipart)

    # generates correct prefix
    expect_message(
        upload_debug_folder_to_s3(debug_prefix, pipeline_config),
        regexp = "with prefix 1234/now"
    )

    # puts in right bucket
    expect_message(
        upload_debug_folder_to_s3(debug_prefix, pipeline_config),
        regexp = pipeline_config$debug_bucket
    )

    # cleanup
    unlink(list.files(debug_path), recursive = TRUE)
})


test_that("upload_matrix_to_s3 completes successfully", {
    before_each()

    # mock things
    data <- matrix()
    pipeline_config <- list(processed_bucket = 'processed-bucket')
    experiment_id <- '1234'

    mockery::stub(upload_matrix_to_s3, 'put_object_in_s3_multipart', stub_put_object_in_s3_multipart)

    # generates correct S3 key
    key <- upload_matrix_to_s3(pipeline_config, experiment_id, data)
    expect_equal(key, '1234/r.rds')
})

test_that("send_output_to_api completes successfully", {
    before_each()

    c(pipeline_config, input, plot_data_keys, output) %<-% send_output_to_api_mock_data

    mockery::stub(send_output_to_api, 'put_object_in_s3', NULL)
    mockery::stub(send_output_to_api, 'paws::sns', mock_sns)

    response <- send_output_to_api(pipeline_config, input, plot_data_keys, output)

    expect_true(response == 'ok')
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

  expect_identical(res[,bound_col], values)
  expect_equal(ncol(res), ncol(dt) + 1)
})


test_that("safe_cbind names bound column as expected", {
  before_each()

  dt <- data.table::data.table(col1 = 1:10, col2 = 11:20)
  values <- seq(1, 20, 2)
  res <- safe_cbind(dt, my_expected_column_name = values)

  expect_true("my_expected_column_name" %in% names(res))
  expect_identical(res[,my_expected_column_name], values)


})


test_that("safe_cbind binds more than one column and names accordingly", {
  before_each()

  dt <- data.table::data.table(col1 = 1:10, col2 = 11:20)
  values_1 <- seq(1, 20, 2)
  values_2 <- values_1 + 2

  res <- safe_cbind(dt, an_interesting_variable = values_1, an_interesting_variable_plus_2 = values_2)

  expect_true("an_interesting_variable" %in% names(res))
  expect_identical(res[,an_interesting_variable], values_1)

  expect_true("an_interesting_variable_plus_2" %in% names(res))
  expect_identical(res[,an_interesting_variable_plus_2], values_2)

})


test_that("cbind_cellset_type names the bound column correctly", {
  before_each()

  dt <- data.table::data.table(col1 = 1:10, col2 = 11:20)
  values <- seq(1, 20, 2)

  res <- cbind_cellset_type(dt, values)

  expect_true("cellset_type" %in% names(res))
  expect_identical(res[,cellset_type], values)

})


test_that("parse_cellsets parses a cellset object", {
  before_each()

  cellsets <- mock_cellsets()

  res <- parse_cellsets(cellsets)

  expect_s3_class(res, "data.table")
  expect_identical(names(res), c("key", "name", "type", "cell_id"))

})

stub_s3_put_object <- function(Bucket, Key, Body, Tagging) {
  response <- list(Expiration = character(),
                   ETag = "this_is_not_an_etag",
                   ServerSideEncryption = character(),
                   VersionId = character(),
                   SSECustomerAlgorithm = character(),
                   SSECustomerKeyMD5 = character(),
                   SSEKMSKeyId = character(),
                   SSEKMSEncryptionContext = character(),
                   BucketKeyEnabled = logical(),
                   RequestCharged = character())

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
                 regexp = "Putting a_key in mock_bucket")
})


test_that("put_object_in_s3 retries if s3$put_object throws an error", {
  before_each()
  mockery::stub(put_object_in_s3,
                "s3$put_object",
                mockery::mock(stop("an error"), stub_s3_put_object))

  pipeline_config <- mock_pipeline_config()
  bucket <- "mock_bucket"
  key <- "mock_key"
  object <- "something"
  key <- "a_key"

  expect_message(put_object_in_s3(pipeline_config, bucket, object, key),
                 regexp = ".*Retrying \\(1/2\\).*")
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

  expected_cellset_types <- c("cluster",
                              "scratchpad",
                              "sample",
                              "metadata",
                              "metadata",
                              "sctype")

  expect_identical(purrr::map2_chr(key, type, get_cellset_type),
                   expected_cellset_types)
})

test_that("send_pipeline_fail_update handles a gem2s call successefully", {
  before_each()
  mockery::stub(send_pipeline_fail_update, 'paws::sns', mock_sns)

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
  mockery::stub(send_pipeline_fail_update, 'paws::sns', mock_sns)
  mockery::stub(send_pipeline_fail_update, 'ids::uuid', "mock-uuid")

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

  mockery::stub(send_pipeline_fail_update, 'paws::sns', mock_sns)
  mockery::stub(send_pipeline_fail_update, 'ids::uuid', "mock-uuid")
  mockery::stub(send_pipeline_fail_update, 'put_object_in_s3', mock_put_object_in_s3)

  c(pipeline_config, input, plot_data_keys, output) %<-% send_output_to_api_mock_data
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
