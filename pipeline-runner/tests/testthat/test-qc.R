# it's only required for tests in this file, so it's not worth it loading
# seurat in the setup.R file
library(Seurat)

qc_setup <- function(experiment_id) {
  paths <- setup_test_paths()
  gem2s_snaps_path <- file.path(paths$snaps, "gem2s")

  gem2s_snap_basename <- paste("gem2s-6", experiment_id, sep = "-")
  gem2s_out_name <- paste(gem2s_snap_basename, "out.R", sep = "-")
  gem2s_qc_config_name <-
    paste(gem2s_snap_basename, "qc_config.R", sep = "-")

  # load QC config from test-gem2s output
  source(file.path(gem2s_snaps_path,  gem2s_qc_config_name))
  # load the final test-gem2s output. it was named `res`
  source(file.path(gem2s_snaps_path, gem2s_out_name))

  if (experiment_id != res$output$scdata_list[[1]]@misc$experimentId) {
    stop("experiment id does not match object sourced from gem2s test output")
  }

  return(list(prev_out = res,
              qc_config = qc_config))

}


#' Fix time stamps in seurat object
#'
#' Seurat adds logs to certain command runs, with time stamps that make snapshot
#' tests fail. This function replaces them with a fixed datetime object.
#'
#' @param scdata seurat object
#'
#' @return seurat object with fixed time stamps
#'
clean_timestamp <- function(scdata) {
  fixed_datetime <- as.POSIXct("1991-12-19 05:23:00", tz = "GMT")

  for (slot in names(scdata@commands)) {
    scdata@commands[[slot]]@time.stamp <- fixed_datetime
  }

  if ("ingestionDate" %in% names(scdata@misc)) {
    scdata@misc$ingestionDate <- fixed_datetime
  }

  return(scdata)
}


snapshot_qc_output <- function(scdata_list, cells_id, guidata) {
  # scdata_list is returned for every sample, so we'll have duplicated stuff here.
  if (is.list(scdata_list)) {
    clean_scdata_list <- lapply(scdata_list, clean_timestamp)
  }
  else {
    clean_scdata_list <- clean_timestamp(scdata_list)
  }

  expect_snapshot(rlang::hash(clean_scdata_list))
  expect_snapshot(str(clean_scdata_list, vec.len = 20))

  expect_snapshot(rlang::hash(cells_id))
  expect_snapshot(str(cells_id, vec.len = 20))

  expect_snapshot(rlang::hash(guidata))
  expect_snapshot(str(guidata, vec.len = 20))

}


test_qc <- function(experiment_id) {
  test_that("QC results are reproducible", {

    # initial setup
    setup <- qc_setup(experiment_id)
    prev_out <- setup$prev_out
    pipeline_config <- mock_pipeline_config()

    # load task functions, and use required stubs
    QC_TASK_LIST$configureEmbedding <- "stubbed_embed_and_cluster"
    tasks <- lapply(QC_TASK_LIST, get)

    filtered_cells_id <- generate_first_step_ids(prev_out$output$scdata_list)
    expect_snapshot(filtered_cells_id)

    # QC stores everything as global variables. This replicates it
    global_vars <- list(data = prev_out$output$scdata_list,
                        new_ids = filtered_cells_id,
                        config = list(),
                        plotData = list())


    for (task_name in names(QC_TASK_LIST)) {
      # the config changes at the integration step (no more independent samples)
      if (which(task_name == names(QC_TASK_LIST)) < which(names(QC_TASK_LIST) == "dataIntegration")) {

        for (sample_id in names(prev_out$output$scdata_list)) {

          config <- setup$qc_config[[task_name]][[sample_id]]

          global_vars <- run_qc_step(
            scdata = global_vars$data,
            config = config,
            tasks = tasks,
            task_name = task_name,
            cells_id = global_vars$new_ids,
            sample_id = sample_id,
            debug_config =  pipeline_config$debug_config
          )
        }

      } else {
        config <- setup$qc_config[[task_name]]
        global_vars <- run_qc_step(
          scdata = global_vars$data,
          config = config,
          tasks = tasks,
          task_name = task_name,
          cells_id = global_vars$new_ids,
          sample_id = "",
          debug_config =  pipeline_config$debug_config
        )
      }

      snapshot_qc_output(global_vars$data,
                         global_vars$new_ids,
                         global_vars$plotData)

    }
    expect_snapshot_file(file.path(pipeline_config$cell_sets_bucket, "cluster_cellsets.json"),
                        name = "cluster_cell_sets.json")
    withr::defer(unlink(pipeline_config$cell_sets_bucket, recursive = TRUE))
  })

}


test_qc("mock_experiment_id")
