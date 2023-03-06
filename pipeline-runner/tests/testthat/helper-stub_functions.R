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


#' (stubbed) download user files
#'
#' Stubs interactions with S3 to be able to test the original function. Cleans
#' up after itself.
#'
#' @inheritParams download_user_files
#'
#' @return list
#' @export
#'
stubbed_download_user_files <- function(input, pipeline_config, prev_out = list()) {
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

    res <- download_user_files(input,
                               pipeline_config,
                               input_dir = "./input",
                               prev_out = prev_out)
    # download_user_files creates a "/input" folder in the pod. defer deleting
    # it during tests.
    withr::defer(unlink("./input", recursive = TRUE), envir = parent.frame(2))

    return(res)
}

#' (stubbed) Read user input files
#'
#' @inheritSection load_user_files description params return
#'
#' @export
#'
stubbed_load_user_files <- function(input, pipeline_config, prev_out, input_dir) {
  INPUT_DIR <- "./input"
  return(load_user_files(input, pipeline_config, prev_out, INPUT_DIR))
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


stub_put_object_in_s3_multipart <-
  function(pipeline_config, bucket, object, key) {
    # if we do not test the raw RDSs uploaded by upload_to_aws, we can remove
    # this stub function
    dir_path <- file.path(bucket, dirname(key))
    if (!dir.exists(dir_path))
      dir.create(dir_path, recursive = TRUE)
    file.copy(object, file.path(bucket, key))
  }


#' setup useful paths for tests
#'
#' This returns the correct path to the tests/testhat folder, independent of
#' running the tests interactively or in a container.
#'
#' @return list of paths
#'
setup_test_paths <- function() {
  base_path <- ifelse(basename(getwd()) == "pipeline-runner",
                      "./tests/testthat",
                      ".")

  mock_data_path <- file.path(base_path,
                              "mock_data")

  snaps_path <- file.path(base_path, "_snaps")

  return(list(
    base = base_path,
    mock_data = mock_data_path,
    snaps = snaps_path
  ))
}


#' create a mock pipeline config
#'
#' locally sets environment variables. Replaces the bucket paths with the paths
#' to potentially mocked buckets in the `tests/testthat/mock_data` folder
#'
#' @param development_aws_server character
#'
#' @return list - pipeline config
#'
mock_pipeline_config <- function(development_aws_server = "mock_aws_server") {
  withr::local_envvar(list("AWS_ACCOUNT_ID" = "000000000000",
                           "ACTIVITY_ARN" = "mock_arn"))

  pipeline_config <- load_config(development_aws_server)

  paths <- setup_test_paths()

  # replace buckets with the local path
  for (bucket in grep("bucket$", names(pipeline_config), value = TRUE)) {
    pipeline_config[[bucket]] <- file.path(paths$mock_data, bucket)
  }

  return(pipeline_config)
}


#' stub tempdir
#'
#' creates the temp dir to the same place to be able to capture written files
#' consistently
#'
#' @return character - temp dir path
#'
stub_tempdir <- function() {

  paths <- setup_test_paths()
  temp_path <- file.path(paths$mock_data, "temp")

  if (!dir.exists(temp_path))
    dir.create(temp_path, recursive = TRUE)

  return(temp_path)
}


#' (stubbed) Upload Seurat and cellsets objects to aws
#'
#' Stubs interactions with S3 to be able to test the original function.
#'
#' @inheritParams upload_to_aws
#'
#' @return list with experiment parameters
#' @export
#'
stubbed_upload_to_aws <- function(input, pipeline_config, prev_out) {
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


stub_update_sets_through_api <- function(cell_sets_object,
                                         api_url,
                                         experiment_id,
                                         cell_set_key,
                                         auth_JWT,
                                         ignore_ssl_cert = FALSE) {
  cellsets_bucket <- "./mock_data/cell_sets_bucket"
  if (!dir.exists(cellsets_bucket)) {
    dir.create(cellsets_bucket)
  }
  cellset_json <- RJSONIO::toJSON(cell_sets_object)
  write(cellset_json, file.path(cellsets_bucket, "cluster_cellsets.json"))
}


#' (Stubbed) Run clustering and embedding
#'
#' Stubs interactions with the API to be able to test the original function.
#'
#' @inheritParams embed_and_cluster
#'
#' @return list with scdata, cell_ids and config
#' @export
#'
stubbed_embed_and_cluster <- function(scdata, config, sample_id, cells_id, task_name, ignore_ssl_cert) {

  mockery::stub(embed_and_cluster,
                "update_sets_through_api",
                stub_update_sets_through_api)

  embed_and_cluster(scdata, config, sample_id, cells_id, task_name, ignore_ssl_cert)
}
