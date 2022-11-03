#' stub_s3_list_objects <- function(Bucket, Prefix) {
#'   # this workaround is the lesser evil ("bucket/./project/sample")
#'   Prefix <- gsub("^./", "", Prefix)
#'
#'   # returns list with structure like s3$list_objects, but with mocked paths
#'   files <-
#'     list.files(file.path(Bucket, Prefix),
#'                full.names = TRUE,
#'                recursive = TRUE)
#'   l <- as.list(files)
#'   names(l) <- NULL
#'   l2 <- lapply(l, list)
#'
#'   for (i in seq_along(l2)) {
#'     names(l2[[i]]) <- "Key"
#'   }
#'   list(Contents = l2)
#' }
#'
#'
#' stub_s3_get_object <- function(Bucket, Key) {
#'   # returns a list of the raw file read from the mocked s3 bucket.
#'   out <- list(body = readBin(file.path(Bucket, Key),
#'                              what = "raw", n = 120000000L),
#'               rest = list())
#'   return(out)
#' }
#'
#'
#' stub_file.path <- function(...) {
#'   file.path(".", ...)
#' }
#'
#'
#' stubbed_download_user_files <- function(input, pipeline_config, prev_out = list()) {
#'   # helper to simplify calls to the stubbed function
#'
#'   mockedS3 <- list(list_objects = stub_s3_list_objects,
#'                    get_object = stub_s3_get_object)
#'
#'   # where makes sure where we are stubbing the what calls.
#'   mockery::stub(where = download_user_files,
#'                 what = "paws::s3",
#'                 how = mockedS3)
#'   mockery::stub(get_gem2s_file, "s3$list_objects", mockedS3$list_objects)
#'
#'   mockery::stub(download_user_files, "file.path", stub_file.path)
#'   mockery::stub(get_gem2s_file, "file.path", stub_file.path)
#'
#'   mockery::stub(download_and_store, "s3$get_object", mockedS3$get_object)
#'
#'   res <- download_user_files(input,
#'                              pipeline_config,
#'                              input_dir = "./input",
#'                              prev_out = prev_out)
#'   # download_user_files creates a "/input" folder in the pod. defer deleting
#'   # it during tests.
#'   withr::defer(unlink("./input", recursive = TRUE), envir = parent.frame(2))
#'
#'   return(res)
#' }
#'
#' stubbed_load_user_files <- function(input, pipeline_config, prev_out, input_dir) {
#'   INPUT_DIR <- "./input"
#'   return(load_user_files(input, pipeline_config, prev_out, INPUT_DIR))
#' }
#'
#'
#' stub_put_object_in_s3 <-
#'   function(pipeline_config, bucket, object, key) {
#'     if (!dir.exists(bucket))
#'       dir.create(bucket)
#'     writeBin(object, file.path(bucket, key))
#'   }
#'
#'
#' stub_remove_bucket_folder <-
#'   function(pipeline_config, bucket, folder) {
#'     # we cleanup inside the corresponding tests
#'   }
#'
#'
#' stub_put_object_in_s3_multipart <-
#'   function(pipeline_config, bucket, object, key) {
#'     # if we do not test the raw RDSs uploaded by upload_to_aws, we can remove
#'     # this stub function
#'     dir_path <- file.path(bucket, dirname(key))
#'     if (!dir.exists(dir_path))
#'       dir.create(dir_path, recursive = TRUE)
#'     file.copy(object, file.path(bucket, key))
#'   }
#'
#'
#' #' setup useful paths for tests
#' #'
#' #' This returns the correct path to the tests/testhat folder, independent of
#' #' running the tests interactively or in a container.
#' #'
#' #' @return list of paths
#' #'
#' setup_test_paths <- function() {
#'   base_path <- ifelse(basename(getwd()) == "pipeline-runner",
#'                       "./tests/testthat",
#'                       ".")
#'
#'   mock_data_path <- file.path(base_path,
#'                               "mock_data")
#'
#'   snaps_path <- file.path(base_path, "_snaps")
#'
#'   return(list(
#'     base = base_path,
#'     mock_data = mock_data_path,
#'     snaps = snaps_path
#'   ))
#' }
#'
#'
#' #' create a mock pipeline config
#' #'
#' #' locally sets environment variables. Replaces the bucket paths with the paths
#' #' to potentially mocked buckets in the `tests/testthat/mock_data` folder
#' #'
#' #' @param development_aws_server character
#' #'
#' #' @return list - pipeline config
#' #'
#' mock_pipeline_config <- function(development_aws_server = "mock_aws_server") {
#'   withr::local_envvar(list("AWS_ACCOUNT_ID" = "000000000000",
#'                            "ACTIVITY_ARN" = "mock_arn"))
#'
#'   pipeline_config <- load_config(development_aws_server)
#'
#'   paths <- setup_test_paths()
#'
#'   # replace buckets with the local path
#'   for (bucket in grep("bucket$", names(pipeline_config), value = TRUE)) {
#'     pipeline_config[[bucket]] <- file.path(paths$mock_data, bucket)
#'   }
#'
#'   return(pipeline_config)
#' }
#'
#'
#' #' stub tempdir
#' #'
#' #' creates the temp dir to the same place to be able to capture written files
#' #' consistently
#' #'
#' #' @return character - temp dir path
#' #'
#' stub_tempdir <- function() {
#'
#'   paths <- setup_test_paths()
#'   temp_path <- file.path(paths$mock_data, "temp")
#'
#'   if (!dir.exists(temp_path))
#'     dir.create(temp_path, recursive = TRUE)
#'
#'   return(temp_path)
#' }
#'
#'
#' stubbed_upload_to_aws <- function(input, pipeline_config, prev_out) {
#'   mockery::stub(upload_to_aws, "put_object_in_s3", stub_put_object_in_s3)
#'   mockery::stub(upload_to_aws,
#'                 "remove_bucket_folder",
#'                 stub_remove_bucket_folder)
#'   mockery::stub(upload_to_aws, "tempdir", stub_tempdir)
#'   mockery::stub(upload_to_aws,
#'                 "put_object_in_s3_multipart",
#'                 stub_put_object_in_s3_multipart)
#'
#'   upload_to_aws(input, pipeline_config, prev_out)
#' }




##############



load_experiment_input <- function(mock_data_path, experiment_id) {
  # use RJSONIO because the pipeline parses with RJSONIO
  RJSONIO::fromJSON(file.path(
    mock_data_path,
    "input",
    paste(experiment_id, "input.json", sep = "-")
  ))

}


make_snapshot_name <- function(step_n, experiment_id, output_name) {
  paste("gem2s", step_n, experiment_id, output_name, sep = "-")
}


snapshot_final_output <- function(res, experiment_id) {
  step_n <- 6

  qc_config <- res$output$qc_config
  snap_name <- make_snapshot_name(step_n, experiment_id, "qc_config.R")
  withr::with_tempfile("tf_qc_config", {
    dump("qc_config", tf_qc_config)
    expect_snapshot_file(tf_qc_config, name = snap_name)
  })

  # fully snapshot final gem2s file (step 7 does not return the scdata_list)
  snap_name <- make_snapshot_name(step_n, experiment_id, "out.R")
  withr::with_tempfile("tf", {
    dump("res", tf)
    expect_snapshot_file(tf, name = snap_name)
  })

}


snapshot_cellsets_file <- function(experiment_id, pipeline_config) {
  step_n <- 7
  cellset_snap_name <- make_snapshot_name(step_n, experiment_id, "cellsets.json")
  expect_snapshot_file(file.path(pipeline_config$cell_sets_bucket, experiment_id),
                       name = cellset_snap_name)
}


#' Run gem2s sequentially and snapshot test each step
#'
#' The snapshots saved are only the structure of each steps' output and the hash,
#' except for the final object (the output to prepare experiment) which is
#' completely saved as a text representation (using `dump`). That can be sourced
#' and used as input for similar tests in QC.
#'
#' To run on other experiments two things are required:
#'
#' 1. add sample files to the `mock_data/originals_bucket` folder
#' 2. create an `input.json` file, named like `<experiment_id>-input.json` and put it
#'    in the `mock_data/input` folder. It should contain the correct sample file
#'    ids (the name of the sample files in the `originals_bucket` folder) and
#'    `experiment_id`
#'
#' @param experiment_id character
#'
test_gem2s <- function(experiment_id) {
  test_that("gem2s is reproducible", {
    paths <- setup_test_paths()
    input <- load_experiment_input(paths$mock_data, experiment_id)
    pipeline_config <- mock_pipeline_config()

    tasks <- lapply(GEM2S_TASK_LIST, get)

    tasks$downloadGem <- stubbed_download_user_files
    tasks$preproc <- stubbed_load_user_files
    tasks$uploadToAWS <- stubbed_upload_to_aws

    res <- list()

    for (task_name in names(tasks)) {
      res <- run_gem2s_step(res$output,
                            input,
                            pipeline_config,
                            tasks,
                            task_name)
      expect_snapshot({
        task_name
        rlang::hash(res)
        str(res)
      })

      if (task_name == "prepareExperiment") {
        snapshot_final_output(res, experiment_id)
      }
    }

    snapshot_cellsets_file(experiment_id, pipeline_config)

    # cleanup
    withr::defer(unlink(pipeline_config$cell_sets_bucket, recursive = TRUE))
    withr::defer(unlink(pipeline_config$source_bucket, recursive = TRUE))
    withr::defer(unlink(file.path(paths$mock_data, "temp"), recursive = TRUE))

  })

}


test_gem2s("mock_experiment_id")
