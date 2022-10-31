library(Seurat)
fdump <- function(var, file) {
  dump(deparse(substitute(var, env = environment())), file = file)
}

qc_step <- function(task_name, scdata_list, qc_step_config, cells_id) {
  qc_function <- QC_TASK_LIST[[task_name]]
  guidata <- list()
  for (sample_id in names(scdata_list)) {
    res <-
      qc_function(scdata_list, qc_step_config[[sample_id]], sample_id, cells_id)
    # update the cells_id "global variable"
    cells_id <- res$new_ids
    # update the scdata_list "global variable"
    scdata_list <- res$data

    # the plot data needs to be saved by sample_id
    guidata[[sample_id]] <- res$plotData

    # nothing is done with the config after it's assigned to rest_of_results
    # out_config <- res$config
  }

  return(list(
    scdata_list = scdata_list,
    cells_id = cells_id,
    guidata = guidata
  ))

}

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

  expect_snapshot(str(clean_scdata_list, vec.len = 20))
  expect_snapshot(str(cells_id, vec.len = 20))
  expect_snapshot(str(guidata, vec.len = 20))

}

qc_setup <- function(experiment_id) {
  paths <- setup_test_paths()
  gem2s_snaps_path <- file.path(paths$snaps, "gem2s")

  gem2s_snap_basename <- paste("gem2s-6", experiment_id, sep = "-")
  gem2s_out_name <- paste(gem2s_snap_basename, "out.R", sep = "-")
  gem2s_qc_config_name <- paste(gem2s_snap_basename, "qc_config.R", sep = "-")

  # load QC config from test-gem2s output
  source(file.path(gem2s_snaps_path,  gem2s_qc_config_name))
  # load the final test-gem2s output. it was named `res`
  source(file.path(gem2s_snaps_path, gem2s_out_name))

  if (experiment_id != res$output$scdata_list[[1]]@misc$experimentId) {
    stop("experiment id does not match object sourced from gem2s test output")
  }

  return(
    list(
      prev_out = res,
      qc_config = qc_config
    )
  )

}


test_qc <- function(experiment_id) {
  test_that("sequentially test qc", {

    # initial setup
    setup <- qc_setup(experiment_id)
    prev_out <- setup$prev_out
    pipeline_config <- mock_pipeline_config()
    message(pipeline_config$cell_sets_bucket)


    filtered_cells_id <- generate_first_step_ids(prev_out$output$scdata_list)
    expect_snapshot(filtered_cells_id)


    global_vars <- list(data = prev_out$output$scdata_list,
                        new_ids = filtered_cells_id,
                        config = list(),
                        plotData = list())

    QC_TASK_LIST$configureEmbedding <- "stubbed_embed_and_cluster"
    tasks <- lapply(QC_TASK_LIST, get)

    for (task_name in names(QC_TASK_LIST)) {

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

        snapshot_qc_output(global_vars$data,
                           global_vars$new_ids,
                           global_vars$plotData)

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

        snapshot_qc_output(clean_timestamp(global_vars$data),
                           global_vars$new_ids,
                           global_vars$plotData)
      }
    }
    expect_snapshot_file(file.path(pipeline_config$cell_sets_bucket, "cluster_cellsets.json"),
                        name = "cluster_cell_sets.json")
    withr::defer(unlink(pipeline_config$cell_sets_bucket, recursive = T))
  })

}


test_qc("mock_experiment_id")


