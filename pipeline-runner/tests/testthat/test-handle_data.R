mock_sns <- function(config) {
    return(
        list(publish = function (Message, TopicArn, MessageAttributes) {
            return (list(MessageId = 'ok'))
        }
    ))
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
