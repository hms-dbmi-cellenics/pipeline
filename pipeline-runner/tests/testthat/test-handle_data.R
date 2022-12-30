mock_sns <- function(config) {
    return(
        list(publish = function (Message, TopicArn, MessageAttributes) {
            return (list(MessageId = 'ok'))
        }
    ))
}

mock_cellsets <- function(){
  # get a snapshot cellsets json
  paths <- setup_test_paths()
  jsonlite::fromJSON(file.path(paths$snaps, "gem2s", "gem2s-7-mock_experiment_id-cellsets.json"), flatten = TRUE)

}

test_that("send_gem2s_update_to_api completes successfully", {
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
    c(pipeline_config, input, plot_data_keys, output) %<-% send_output_to_api_mock_data

    mockery::stub(send_output_to_api, 's3$put_object', stub_put_object_in_s3)
    mockery::stub(send_output_to_api, 'paws::sns', mock_sns)

    response <- send_output_to_api(pipeline_config, input, plot_data_keys, output)

    expect_true(response == 'ok')
})


test_that("safe_cbind returns empty data.table when binding an empty data.table with a vector", {

  dt_empty <- data.table::data.table()
  col <- c(a_col = "a_value")

  res <- safe_cbind(dt_empty, col)

  expect_identical(res, dt_empty)

})


test_that("safe_cbind adds a column to a non-empty data.table", {
  dt <- data.table::data.table(col1 = 1:10, col2 = 11:20)
  values <- seq(1, 20, 2)
  res <- safe_cbind(dt, bound_col = values)

  expect_identical(res[,bound_col], values)
  expect_equal(ncol(res), ncol(dt) + 1)
})


test_that("safe_cbind names bound column as expected", {
  dt <- data.table::data.table(col1 = 1:10, col2 = 11:20)
  values <- seq(1, 20, 2)
  res <- safe_cbind(dt, my_expected_column_name = values)

  expect_true("my_expected_column_name" %in% names(res))
  expect_identical(res[,my_expected_column_name], values)


})


test_that("safe_cbind binds more than one column and names accordingly", {

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

  dt <- data.table::data.table(col1 = 1:10, col2 = 11:20)
  values <- seq(1, 20, 2)

  res <- cbind_cellset_type(dt, values)

  expect_true("cellset_type" %in% names(res))
  expect_identical(res[,cellset_type], values)

})


test_that("parse_cellsets parses a cellset object", {

  cellsets <- mock_cellsets()

  res <- parse_cellsets(cellsets)

  expect_s3_class(res, "data.table")
  expect_identical(names(res), c("key", "name", "type", "cell_id"))

})
