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
      res <- run_pipeline_step(res$output,
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
