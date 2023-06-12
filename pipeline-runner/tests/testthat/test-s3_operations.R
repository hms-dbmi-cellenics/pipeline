# mock params
get_test_params <- function() {
  source_prefix <- "src_prefix"
  source_keys <- c("key-1", "key-2")

  params <- list(
    source_bucket = "src_bucket",
    source_prefix = source_prefix,
    destination_bucket = "dst_bucket",
    destination_prefix = "dst_prefix",
    aws_config = list(fake = TRUE),
    source_keys = source_keys
  )

  mock_list_objects_response <- list(
    Contents = list(
      list(Key = paste0(source_prefix, source_keys[[1]])),
      list(Key = paste0(source_prefix, source_keys[[2]]))
    )
  )

  # s3 mocked functions
  mock_list_objects <- mockery::mock(mock_list_objects_response)
  mock_copy_object <- mockery::mock(cycle = TRUE)
  mock_paws_s3 <- c(
    list_objects = mock_list_objects,
    copy_object = mock_copy_object
  )

  params$mock_paws_s3 <- mock_paws_s3
  params$mock_list_objects <- mock_list_objects
  params$mock_copy_object <- mock_copy_object

  return(params)
}

test_that("s3_copy_by_prefix works using default transform (identity))", {
  params <- get_test_params()

  mockery::stub(s3_copy_by_prefix, "paws::s3", params$mock_paws_s3)

  s3_copy_by_prefix(
    params$source_bucket,
    params$source_prefix,
    params$destination_bucket,
    params$destination_prefix,
    params$aws_config
  )

  expect_args(params$mock_list_objects, 1, params$source_bucket, params$source_prefix)
  expect_args(
    params$mock_copy_object,
    1,
    paste0(params$source_bucket, "/", params$source_prefix, params$source_keys[[1]]),
    params$destination_bucket,
    paste0(params$destination_prefix, params$source_keys[[1]])
  )

  expect_args(
    params$mock_copy_object,
    2,
    paste0(params$source_bucket, "/", params$source_prefix, params$source_keys[[2]]),
    params$destination_bucket,
    paste0(params$destination_prefix, params$source_keys[[2]])
  )
})

test_that("s3_copy_by_prefix works passing a custom transform", {
  params <- get_test_params()
  mockery::stub(s3_copy_by_prefix, "paws::s3", params$mock_paws_s3)

  mock_transform <- mockery::mock("dest1", "dest2")

  s3_copy_by_prefix(
    params$source_bucket,
    params$source_prefix,
    params$destination_bucket,
    params$destination_prefix,
    params$aws_config,
    mock_transform
  )

  expect_args(params$mock_list_objects, 1, params$source_bucket, params$source_prefix)
  expect_args(
    params$mock_copy_object,
    1,
    paste0(params$source_bucket, "/", params$source_prefix, params$source_keys[[1]]),
    params$destination_bucket,
    "dest1"
  )

  expect_args(
    params$mock_copy_object,
    2,
    paste0(params$source_bucket, "/", params$source_prefix, params$source_keys[[2]]),
    params$destination_bucket,
    "dest2"
  )

  expect_args(mock_transform, 1, paste0(params$destination_prefix, params$source_keys[[1]]))
  expect_args(mock_transform, 2, paste0(params$destination_prefix, params$source_keys[[2]]))
})
