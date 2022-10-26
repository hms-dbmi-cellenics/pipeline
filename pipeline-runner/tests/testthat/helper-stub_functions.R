stub_s3_list_objects <- function(Bucket, Prefix) {
  # this workaround is the lesser evil ("bucket/./project/sample")
  Prefix <- gsub("^./", "", Prefix)

  # returns list with structure like s3$list_objects, but with mocked paths
  files <-
    list.files(file.path(Bucket, Prefix),
               full.names = TRUE,
               recursive = TRUE)
  l <- as.list(files)
  names(l) <- NULL
  l2 <- lapply(l, list)

  for (i in seq_along(l2)) {
    names(l2[[i]]) <- "Key"
  }
  list(Contents = l2)
}


stub_s3_get_object <- function(Bucket, Key) {
  # returns a list of the raw file read from the mocked s3 bucket.
  out <- list(body = readBin(file.path(Bucket, Key),
                             what = "raw", n = 120000000L),
              rest = list())
  return(out)
}


stub_file.path <- function(...) {
  file.path(".", ...)
}


stubbed_download_user_files <-
  function(input, pipeline_config, prev_out = list()) {
    # helper to simplify calls to the stubbed function

    mockedS3 <- list(list_objects = stub_s3_list_objects,
                     get_object = stub_s3_get_object)

    # where makes sure where we are stubbing the what calls.
    mockery::stub(where = download_user_files,
                  what = "paws::s3",
                  how = mockedS3)
    mockery::stub(get_gem2s_file, "s3$list_objects", mockedS3$list_objects)

    mockery::stub(download_user_files, "file.path", stub_file.path)
    mockery::stub(get_gem2s_file, "file.path", stub_file.path)

    mockery::stub(download_and_store, "s3$get_object", mockedS3$get_object)

    res <-
      download_user_files(input,
                          pipeline_config,
                          input_dir = "./input",
                          prev_out = prev_out)
    # download_user_files creates a "/input" folder in the pod. defer deleting
    # it during tests.
    withr::defer(unlink("./input", recursive = TRUE), envir = parent.frame())

    res
  }


stub_put_object_in_s3 <-
  function(pipeline_config, bucket, object, key) {
    if (!dir.exists(bucket))
      dir.create(bucket)
    writeBin(object, file.path(bucket, key))
  }


stub_remove_bucket_folder <-
  function(pipeline_config, bucket, folder) {
    # we cleanup inside the corresponding tests
  }


stub_tempdir <- function() {
  # stub to write always to the same place and be able to capture written files
  # consistently
  base_path <- ifelse(basename(getwd()) == "pipeline-runner",
                      "./tests/testthat",
                      ".")

  mock_path <- file.path(base_path,
                         "mock_data")

  temp_path <- file.path(mock_path, "temp")

  if (!dir.exists(temp_path))
    dir.create(temp_path, recursive = TRUE)

  return(temp_path)
}


stub_put_object_in_s3_multipart <-
  function(pipeline_config, bucket, object, key) {
    # if we do not test the raw RDSs uploaded by upload_to_aws, we can remove
    # this stub function
    dir_path <- file.path(bucket, dirname(key))
    if (!dir.exists(dir_path))
      dir.create(dir_path, recursive = TRUE)
    file.copy(object, file.path(bucket, key))
  }


stubbed_upload_to_aws <-
  function(input, pipeline_config, prev_out) {
    mockery::stub(upload_to_aws, "put_object_in_s3", stub_put_object_in_s3)
    mockery::stub(upload_to_aws,
                  "remove_bucket_folder",
                  stub_remove_bucket_folder)
    mockery::stub(upload_to_aws, "tempdir", stub_tempdir)
    mockery::stub(upload_to_aws,
                  "put_object_in_s3_multipart",
                  stub_put_object_in_s3_multipart)

    upload_to_aws(input, pipeline_config, prev_out)

  }


path_setup <- function() {
  base_path <- ifelse(basename(getwd()) == "pipeline-runner",
                      "./tests/testthat",
                      ".")

  mock_path <- file.path(base_path,
                         "mock_data")

  snaps_path <- file.path(base_path, "_snaps")

  return(list(base = base_path, mock_data = mock_path, snaps = snaps_path))

}


mock_pipeline_config <-
  function(development_aws_server = "mock_aws_server") {
    local_envvar(list("AWS_ACCOUNT_ID" = "000000000000",
                      "ACTIVITY_ARN" = "mock_arn"))

    pipeline_config <- load_config(development_aws_server)

    paths <- path_setup()

    # replace buckets with the local path
    for (bucket in grep("bucket$", names(pipeline_config), value = TRUE)) {
      pipeline_config[[bucket]] <- file.path(paths$mock_data, bucket)
    }

    return(pipeline_config)

  }
