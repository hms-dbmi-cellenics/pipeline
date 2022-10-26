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


mock_pipeline_config <-
  function(development_aws_server = "mock_aws_server") {
    local_envvar(list("AWS_ACCOUNT_ID" = "000000000000",
                      "ACTIVITY_ARN" = "mock_arn"))

    pipeline_config <- load_config(development_aws_server)

    base_path <- ifelse(basename(getwd()) == "pipeline-runner",
                        "./tests/testthat",
                        ".")


    mock_path <- file.path(base_path,
                           "mock_data")

    # replace buckets with the local path
    for (bucket in grep("bucket$", names(pipeline_config), value = TRUE)) {
      pipeline_config[[bucket]] <- file.path(mock_path, bucket)
    }

    return(pipeline_config)

  }


load_experiment_input <- function(mock_data_path, experiment_id) {
  RJSONIO::fromJSON(file.path(
    mock_data_path,
    "input",
    paste(experiment_id, "input.json", sep = "-")
  ))

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


make_snapshot_name <- function(step_n, experiment_id, output_name) {
  paste("gem2s", step_n, experiment_id, output_name, sep = "-")
}


gem2s_setup <- function(experiment_id) {

  paths <- path_setup()
  input <- load_experiment_input(paths$mock_data, experiment_id)
  pipeline_config <- mock_pipeline_config()

  output <-  list(paths = paths,
                  input = input,
                  pipeline_config = pipeline_config)

  return(output)
}


test_gem2s_v2 <- function(experiment_id) {
  test_that("some gem2s steps work", {
    setup <- gem2s_setup(experiment_id)

    #step_n <- 1
    res <- stubbed_download_user_files(setup$input, setup$pipeline_config)

    # snapshot download_user_files output
    #snapshot_name <- make_snapshot_name(step_n, experiment_id, "out")
    expect_snapshot(str(res, vec.len = 20))

    # gem2s 2

    # folder where download_user_files puts things in
    input_path <- file.path(setup$paths$base, "input")
    res <- load_user_files(setup$input, NULL, res$output, input_path)
    expect_snapshot(str(res, vec.len = 20))

    # gem2s 3

    res <- run_emptydrops(setup$input, NULL, res$output)
    expect_snapshot(str(res, vec.len = 20))

    # gem2s 4

    res <- score_doublets(setup$input, NULL, res$output)
    expect_snapshot(str(res, vec.len = 20))

    # gem2s 5

    res <- create_seurat(setup$input, NULL, res$output)
    expect_snapshot(str(res, vec.len = 20))

    # gem2s 6
    step_n <- 6
    res <- prepare_experiment(setup$input, NULL, res$output)
    expect_snapshot(str(res, vec.len = 20))

    # snapshot qc_config independently as it is used frequently in QC tests
    qc_config <- res$output$qc_config
    qc_config_snapshot_name <- make_snapshot_name(step_n, experiment_id, "qc_config.R")
    withr::with_tempfile("tf_qc_config", {
      dump("qc_config", tf_qc_config)
      expect_snapshot_file(tf_qc_config, name = qc_config_snapshot_name)
    })

    # fully snapshot final gem2s file (step 7 does not return the scdata_list)
    snapshot_name <- make_snapshot_name(step_n, experiment_id, "out.R")
    withr::with_tempfile("tf", {
      dump("res", tf)
      expect_snapshot_file(tf, name = snapshot_name)
    })


    #gem2s 7
    step_n <- 7

    res <- stubbed_upload_to_aws(setup$input, setup$pipeline_config, res$output)
    expect_snapshot(str(res, vec.len = 20))

    # cellsets file
    cellset_snapshot_name <- make_snapshot_name(step_n, experiment_id, "cellsets.json")
    expect_snapshot_file(
      file.path(setup$pipeline_config$cell_sets_bucket, experiment_id),
      name = cellset_snapshot_name
    )


    # cleanup
    withr::defer(unlink(setup$pipeline_config$cell_sets_bucket, recursive = TRUE))
    withr::defer(unlink(setup$pipeline_config$source_bucket, recursive = TRUE))
    withr::defer(unlink(file.path(setup$paths$mock_data, "temp"), recursive = TRUE))
  })

}

test_gem2s_v2("mock_experiment_id")
