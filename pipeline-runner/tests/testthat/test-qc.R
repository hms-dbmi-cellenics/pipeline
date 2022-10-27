test_that("qc-1 filters empty drops",{

  task_name <- "classifier"

  base_path <- ifelse(basename(getwd()) == "pipeline-runner",
                      "./tests/testthat",
                      ".")

  mock_path <- file.path(base_path,
                         "mock_data")

  pipeline_config <- mock_pipeline_config()

  prev_out <-
    readRDS(file.path(base_path, "_snaps/pipeline", "gem2s-6-out.rds"))

  prev_out <- prev_out$output

  scdata_list <- prev_out$scdata_list
  config <- prev_out$qc_config[[task_name]]
  experiment_id <- prev_out$scdata_list[[1]]@misc$experimentId


  cells_id <-
    readRDS(file.path(
      base_path,
      "_snaps/pipeline",
      paste("qc-0", experiment_id, "filtered_cells.rds", sep = "-")
    ))

  guidata <- list()
  out_config <- list()
  out_data <- list()
  for (sample_id in names(scdata_list)) {
    res <- filter_emptydrops(scdata_list, config[[sample_id]], sample_id, cells_id)
    out_data[[sample_id]] <- res$data
    cells_id[[sample_id]] <- res$new_ids
    guidata[[sample_id]] <- res$plotData
    out_config[[sample_id]] <- res$config
  }

  # scdata_list is returned for every sample, so we'll have duplicated stuff here.
  withr::with_tempfile("temp_out_data", {
    saveRDS(out_data, temp_out_data)
    expect_snapshot_file(temp_out_data, name = paste("qc-1", experiment_id, "out_data.rds", sep = "-"))
  })

  # filtered cells are saved as a global variable and overwritten in each pipeline step
  withr::with_tempfile("temp_cells_id", {
    saveRDS(cells_id, temp_cells_id)
    expect_snapshot_file(temp_cells_id, name = paste("qc-1", experiment_id, "filtered_cells.rds", sep = "-"))
  })

  withr::with_tempfile("temp_guidata", {
    saveRDS(guidata, temp_guidata)
    expect_snapshot_file(temp_guidata, name = paste("qc-1", experiment_id, "guidata.rds", sep = "-"))
  })

  withr::with_tempfile("temp_out_config", {
    saveRDS(out_config, temp_out_config)
    expect_snapshot_file(temp_out_config, name = paste("qc-1", experiment_id, "out_config.rds", sep = "-"))
  })

})


test_qc <- function(experiment_id) {
  test_that("sequentially test qc", {

    paths <- setup_test_paths()

    pipeline_config <- mock_pipeline_config()
    gem2s_out_name <- paste("gem2s-6", experiment_id, "out.R", sep = "-")

    # the variable was called res
    source(file.path(paths$base,
                     "_snaps/gem2s",
                     gem2s_out_name))
    if (experiment_id != res$output$scdata_list[[1]]@misc$experimentId) {
      stop("experiment id does not match object sourced from gem2s test output")
    }

    filtered_cells_id <- generate_first_step_ids(res$output$scdata_list)

    withr::with_tempfile("tf", {
      dump("filtered_cells_id", tf)
      expect_snapshot_file(tf, name = paste("qc-0", experiment_id, "filtered_cells.R", sep = "-"))
    })

    # qc 1

    task_name <- "classifier"

    scdata_list <- prev_out$scdata_list
    config <- prev_out$qc_config[[task_name]]
    experiment_id <- prev_out$scdata_list[[1]]@misc$experimentId

    guidata <- list()
    out_config <- list()
    out_data <- list()
    for (sample_id in names(scdata_list)) {
      res <- filter_emptydrops(scdata_list, config[[sample_id]], sample_id, cells_id)
      out_data[[sample_id]] <- res$data
      cells_id[[sample_id]] <- res$new_ids
      guidata[[sample_id]] <- res$plotData
      out_config[[sample_id]] <- res$config
    }



  })

}
