#' Build activity ARN
#'
#' @param aws_region character
#' @param aws_account_id character
#' @param activity_id character
#'
#' @return string of activity ARN
build_activity_arn <- function(aws_region, aws_account_id, activity_id) {
  if (is.na(activity_id)) {
    return(NA)
  }

  activity_arn <- sprintf(
    "arn:aws:states:%s:%s:activity:%s",
    aws_region, aws_account_id, activity_id
  )
  return(activity_arn)
}


#' Load the pipeline step config
#'
#' Waits for the activity ARN to be assigned. Once it is, it creates a config
#' list with pod, AWS, api, and debug information.
#'
#' @param development_aws_server
#'
#' @return list with config parameters
load_config <- function(development_aws_server) {

  # running in linux needs the IP of the host to work. If it is set as an
  # environment variable (by makefile) honor it instead of the provided
  # parameter
  overriden_server <- Sys.getenv("HOST_IP", "")
  if (overriden_server != "") {
    development_aws_server <- overriden_server
  }

  label_path <- "/etc/podinfo/labels"
  aws_account_id <- Sys.getenv("AWS_ACCOUNT_ID", unset = "242905224710")
  aws_region <- Sys.getenv("AWS_DEFAULT_REGION", unset = "eu-west-1")
  running_in_batch <- Sys.getenv("BATCH", "false")
  domain_name <- Sys.getenv("DOMAIN_NAME")

  activity_arn <- NA

  repeat {
    # if /etc/podinfo/labels exists we are running in a remote aws pod
    if (file.exists(label_path)) {
      labels <- read.csv(label_path,
        sep = "=",
        row.names = 1,
        header = FALSE
      )
      activity_id <- labels["activityId", ]
      activity_arn <- build_activity_arn(
        aws_region,
        aws_account_id,
        activity_id
      )
    }

    # if we didn't find an activity in podinfo try to get the from the env
    # (means we are running locally)
    if (is.na(activity_arn)) {
      activity_arn <- Sys.getenv("ACTIVITY_ARN", unset = NA)
    }

    if (is.na(activity_arn)) {
      message("No activity ARN label set yet, waiting...")
      Sys.sleep(5)
    } else {
      message("Welcome to the pipeline")

      message(activity_arn)
      break
    }
  }

  sandbox <- Sys.getenv("SANDBOX_ID", "default")
  config <- list(
    cluster_env = Sys.getenv("CLUSTER_ENV", "development"),
    sandbox_id = sandbox,
    aws_account_id = aws_account_id,
    aws_region = aws_region,
    aws_config = list(region = aws_region),
    pod_name = Sys.getenv("K8S_POD_NAME", "local"),
    activity_arn = activity_arn,
    api_url = sprintf("http://%s:3000", development_aws_server),
    api_version = "v2",
    debug_config = list(
      step = Sys.getenv("DEBUG_STEP", ""),
      path = Sys.getenv("DEBUG_PATH", "")
    )
  )

  if (config$cluster_env == "staging") {
    config$api_url <- paste0("https://api-", sandbox, ".", domain_name)
  }
  if (config$cluster_env == "production") {
    config$api_url <- paste0("https://api.", domain_name)
  }

  if (config$cluster_env == "development") {
    # DOCKER_GATEWAY_HOST
    config$aws_config[["endpoint"]] <- sprintf(
      "http://%s:4566",
      development_aws_server
    )
    config$aws_config[["credentials"]] <- list(
      creds = list(
        access_key_id = "mock-access-key",
        secret_access_key = "mock-secure-acces-key"
      )
    )
    config$aws_account_id <- "000000000000"

    # This fixes a bug where paws would try to connect to InfraMock as if it
    # was a proper AWS endpoint, i.e. to
    # http://bucket-name.host.docker.internal
    assignInNamespace("update_endpoint_for_s3_config", function(request) {
      return(request)
    }, ns = "paws.common")
  }

  bucket_list <- lapply(bucket_list,
    paste, config$cluster_env, config$aws_account_id,
    sep = "-"
  )

  config <- append(config, bucket_list)

  config[["sns_topic"]] <- paste(
    paste("arn:aws:sns", config$aws_region, config$aws_account_id, "work-results", sep = ":"),
    config$cluster_env, config$sandbox_id, config$api_version,
    sep = "-"
  )

  return(config)
}


#' Run qc pipeline step
#'
#' Calls the corresponding task_name QC pipeline step function.
#'
#' @param scdata Seurat Object
#' @param config list of configuration parameters
#' @param tasks list of pipeline tasks
#' @param task_name character
#' @param cells_id list of filtered cell ids
#' @param sample_id character
#' @param debug_config list of debug parameters
#'
#' @return list of task results
#'
run_qc_step <- function(scdata, config, tasks, task_name, cells_id, sample_id, ignore_ssl_cert, debug_config) {
  if (!task_name %in% names(tasks)) {
    stop("Invalid task: ", task_name)
  }

  handle_debug(scdata, config, task_name, sample_id, debug_config)

  # print info
  task <- tasks[[task_name]]
  message("Running: ", task_name)
  message("Config:")
  str(config)

  # run task and time it
  tstart <- Sys.time()

  # The httr package used in Configure embedding requires ssl configuration
  # to know whether to verify/not ssl configuration. Other qc functions
  # do not use this param, so it'd be useless to pass this param into the other functions
  if(task_name == "configureEmbedding") {
    out <- task(scdata, config, sample_id, cells_id, task_name, ignore_ssl_cert)
  } else {
    out <- task(scdata, config, sample_id, cells_id, task_name)
  }

  ttask <- format(Sys.time() - tstart, digits = 2)
  message(
    "⏱️ Time to complete ", task_name,
    " for sample ", sample_id, ": ", ttask, "\n"
  )

  return(out)
}


#' Run pipeline step
#'
#' Calls the corresponding `task_name` pipeline  step function.
#'
#' The input list only contains experiment level parameters, such as project ID,
#' and sample names and it's only used for downloading user files.
#'
#' @param task_name character
#' @param input list
#'  - sampleIds character
#'  - sampleNames character
#'  - projectId character
#' @param pipeline_config list as defined by load_config
#' @param prev_out list output from previous step
#'
#' @return list of task results
#'
run_pipeline_step <- function(prev_out, input, pipeline_config, tasks, task_name) {
  if (!task_name %in% names(tasks)) {
    stop("Invalid task name given: ", task_name)
  }
  task <- tasks[[task_name]]

  tstart <- Sys.time()
  out <- task(input, pipeline_config, prev_out)
  ttask <- format(Sys.time() - tstart, digits = 2)
  message("⏱️ Time to complete ", task_name, ": ", ttask, "\n")

  return(out)
}


#' Call gem2s
#'
#' Runs step `task_name` of the GEM2S pipeline, sends output message to the API
#'
#' @param task_name character name of the step
#' @param input list containing
#'   - experimentID
#'   - sample IDs, names and S3 paths
#' @param pipeline_config list as defined by load_config
#'
#' @return character message id
#'
call_gem2s <- function(task_name, input, pipeline_config) {
  experiment_id <- input$experimentId

  if (!exists("prev_out")) {
    remove_cell_ids(pipeline_config, experiment_id)
    assign("prev_out", NULL, pos = ".GlobalEnv")
  }

  check_input(input)
  tasks <- lapply(GEM2S_TASK_LIST, get)

  c(data, task_out) %<-% run_pipeline_step(prev_out, input, pipeline_config, tasks, task_name)
  assign("prev_out", task_out, pos = ".GlobalEnv")

  message_id <- send_pipeline_update_to_api(pipeline_config, experiment_id, task_name, data, input, 'GEM2SResponse')

  return(message_id)
}


#' Call subset seurat
#'
#' Runs step `task_name` of the subset seurat pipeline, sends output message to the API
#'
#' @param task_name character name of the step
#' @param input list containing
#'   - parentExperimentId
#'   - childExperimentId
#'   - sample IDs, and names
#' @param pipeline_config list as defined by load_config
#'
#' @return character message id
#'
call_subset <- function(task_name, input, pipeline_config) {
  experiment_id <- input$experimentId

  if (!exists("prev_out")) {
    remove_cell_ids(pipeline_config, experiment_id)
    assign("prev_out", NULL, pos = ".GlobalEnv")
  }

  check_input(input)
  tasks <- lapply(SUBSET_SEURAT_TASK_LIST, get)

  c(data, task_out) %<-% run_pipeline_step(prev_out, input, pipeline_config, tasks, task_name)
  assign("prev_out", task_out, pos = ".GlobalEnv")

  message_id <- send_pipeline_update_to_api(pipeline_config, experiment_id, task_name, data, input, 'GEM2SResponse')

  return(message_id)
}

call_seurat <- function(task_name, input, pipeline_config) {

  experiment_id <- input$experimentId

  # initial step
  if (!exists("prev_out")) {
    remove_cell_ids(pipeline_config, experiment_id)
    assign("prev_out", NULL, pos = ".GlobalEnv")
  }

  check_input(input)
  tasks <- lapply(SEURAT_TASK_LIST, get)

  c(data, task_out) %<-% run_pipeline_step(prev_out, input, pipeline_config, tasks, task_name)
  assign("prev_out", task_out, pos = ".GlobalEnv")

  message_id <- send_pipeline_update_to_api(pipeline_config, experiment_id, task_name, data, input, 'SeuratResponse')
  return(message_id)
}

#' Call copy
#'
#' Runs step `task_name` of the copy pipeline, sends output message to the API
#'
#' @param task_name character name of the step
#' @param input list containing
#'   - parentExperimentId
#'   - childExperimentId
#'   - sample IDs, and names
#' @param pipeline_config list as defined by load_config
#'
#' @return character message id
#'
call_copy <- function(task_name, input, pipeline_config) {
  experiment_id <- input$experimentId

  if (!exists("prev_out")) {
    remove_cell_ids(pipeline_config, experiment_id)
    assign("prev_out", NULL, pos = ".GlobalEnv")
  }

  check_input(input)
  tasks <- lapply(COPY_TASK_LIST, get)

  c(data, task_out) %<-% run_pipeline_step(prev_out, input, pipeline_config, tasks, task_name)
  assign("prev_out", task_out, pos = ".GlobalEnv")

  message_id <- send_pipeline_update_to_api(pipeline_config, experiment_id, task_name, data, input, 'GEM2SResponse')

  return(message_id)
}


#' Call QC pipeline
#'
#' Runs step `task_name` of the data processing pipeline, sends plot data to s3
#' and the output message to the API
#'
#' @param task_name character name of the step
#' @param input list containing:
#'   - step parameters for all samples
#'   - current sample UUID
#'   - uploadCountMatrix (whether or not to upload matrix after step)
#' @param pipeline_config list as defined by load_config
#'
#' @return character message id
#'
call_qc <- function(task_name, input, pipeline_config) {
  experiment_id <- input$experimentId
  config <- input$config
  upload_count_matrix <- input$uploadCountMatrix
  ignore_ssl_cert <- input$ignoreSslCert
  sample_id <- input$sampleUuid
  debug_config <- pipeline_config$debug_config

  if (sample_id != "") {
    config <- config[[sample_id]]
    input$config <- config
  }

  tasks <- lapply(QC_TASK_LIST, get)

  # need this for embed_and_cluster
  config$api_url <- pipeline_config$api_url
  config$auth_JWT <- input$authJWT

  if (!exists("scdata")) {
    message("No single-cell data has been loaded, reloading from S3...")

    # assign it to the global environment so we can
    # persist it across runs of the wrapper
    assign("scdata",
      reload_data_from_s3(pipeline_config, experiment_id, task_name, tasks),
      pos = ".GlobalEnv"
    )

    message("Single-cell data loaded.")
  }

  if (!exists("cells_id")) {
    message("No filtered cell ids have been loaded, loading from S3...")
    if (task_name == names(tasks)[1]) {
      assign("cells_id", generate_first_step_ids(scdata), pos = ".GlobalEnv")
    } else if (task_name %in% names(tasks)) {
      samples <- names(scdata)
      assign("cells_id",
        load_cells_id_from_s3(pipeline_config, experiment_id, task_name, tasks, samples),
        pos = ".GlobalEnv")

      # won't be cells_id in S3 for uploaded Seurat object that is being subsetted
      if (!length(cells_id)) cells_id <- generate_first_step_ids(scdata)

    } else {
      stop("Invalid task name given: ", task_name)
    }
    message("Cells id loaded.")
  }


  # call function to run and update global variable
  c(
    data, new_ids, ...rest_of_results
  ) %<-% run_qc_step(scdata, config, tasks, task_name, cells_id, sample_id, ignore_ssl_cert, debug_config)


  assign("cells_id", new_ids, pos = ".GlobalEnv")

  task_names <- names(tasks)
  integration_index <- match("dataIntegration", task_names)
  task_index <- match(task_name, task_names)
  if (task_index < integration_index) {
    message(
      "Filtered cell ids from ", length(cells_id[[sample_id]]),
      " to ", length(new_ids[[sample_id]])
    )
    next_task <- names(tasks)[[task_index + 1]]
    object_key <- paste0(experiment_id, "/", next_task, "/", sample_id, ".rds")
    upload_cells_id(pipeline_config, object_key, cells_id)
  }

  # upload plot data result to S3
  tstart <- Sys.time()
  plot_data_keys <- send_plot_data_to_s3(pipeline_config, experiment_id, rest_of_results)

  # Upload count matrix data
  if (upload_count_matrix) {
    assign("scdata", data, pos = ".GlobalEnv")
    object_key <- upload_matrix_to_s3(pipeline_config, experiment_id, scdata)
    message(
      "Count matrix uploaded to ", pipeline_config$processed_bucket,
      " with key ", object_key
    )
  }

  # send result to API
  message_id <- send_output_to_api(pipeline_config, input, plot_data_keys, rest_of_results)
  ttask <- format(Sys.time() - tstart, digits = 2)
  message(
    "⏱️ Time to upload ", task_name,
    " objects for sample ", sample_id,
    ": ", ttask
  )

  return(message_id)
}


#' Run heartbeat in the background
#'
#' Sends a heartbeat to the state machine every 'wait_time' seconds Once the task
#' is completed the heartbeat will fail accordingly with a task timeout and exit
#' the loop and a new heartbeat will be set up by next task. This method is
#' invoked with `callr::r_bg` which creates a new process which does not inherit
#' the current workspace or memory, only the provided parameters; that's why we
#' need to re-import `tryCatchLog` & initialize states again.
#'
#' @param task_token character authorization token
#' @param aws_config list of parameters to access AWS step functions
#'
pipeline_heartbeat <- function(task_token, aws_config) {
  library(tryCatchLog)
  message("Starting heartbeat")
  states <- paws::sfn(config = aws_config)

  keep_running <- TRUE
  # amount of time to wait between heartbeats
  wait_time <- 30
  i <- 0
  while (keep_running) {
    tryCatchLog(
      {
        states$send_task_heartbeat(
          taskToken = task_token
        )
        message("Heartbeat sent: ", i)
      },
      error = function(e) {
        message("Send task heartbeat failed: ", e$message)
        message("Stopping heartbeat after ", i + 1)
        keep_running <- FALSE
      }
    )
    i <- i + 1
    # sleep until next heartbeat
    Sys.sleep(wait_time)
  }
}


#' Start heartbeat as a background process
#'
#' messages inside the background process will ONLY be printed into
#' `/tmp/[out|err]`. To see them:
#'  1. log into the R container
#'  2. `cat /tmp/out` or `tail -f /tmp/out`
#'
#' @inheritParams pipeline_heartbeat
#'
start_heartbeat <- function(task_token, aws_config) {
  message("Starting heartbeat")

  heartbeat_proc <- callr::r_bg(
    func = pipeline_heartbeat, args = list(
      task_token, aws_config
    ),
    stdout = "/tmp/out",
    stderr = "/tmp/err"
  )
  return(heartbeat_proc)
}

handlers <- c(
  qc = call_qc,
  gem2s = call_gem2s,
  seurat = call_seurat,
  subset = call_subset,
  copy = call_copy
)

#' calls the appropriate process, QC or gem2s
#'
#' @param input list with parsed input json. should contain:
#'   - task name
#'   - current sample_id
#'   - extra config parameters
#' @param pipeline_config list as generated by load_config
#'
#' @return character message id
#'
wrapper <- function(input, pipeline_config) {
  task_name <- input$taskName
  message("\n------\nStarting task: ", task_name, "\n")
  message("Input:")
  # remove config from print to avoid huge redundant logs
  str(input[names(input) != "config"])

  # common to gem2s and data processing
  server <- input$server
  input <- input[names(input) != "server"]
  process_name <- input$processName

  if (!process_name %in% names(handlers)) {
    stop("Process name not recognized")
  }

  message_id <- handlers[[process_name]](task_name, input, pipeline_config)

  return(message_id)
}

get_user_error <- function(msg) {
  # check if error is a defined code with a corresponding message in the UI
  if (msg %in% errors) return(msg)

  return("We had an issue while processing your data.")
}

#' run the pipeline
#'
#' Loads configurations and repeats the wrapper call until no more messages are
#' received
#'
init <- function() {
  pipeline_config <- load_config("host.docker.internal")
  message("Loaded pipeline config")
  states <- paws::sfn(config = pipeline_config$aws_config)
  message("Loaded step function")

 print(sessionInfo())

  futile.logger::flog.layout(futile.logger::layout.simple)
  futile.logger::flog.threshold(futile.logger::ERROR)

  message("Waiting for tasks")

  repeat {
    c(task_token, input_json) %<-% states$get_activity_task(
      activityArn = pipeline_config$activity_arn,
      workerName = pipeline_config$pod_name
    )

    if (is.null(task_token) || !length(task_token) || task_token == "") {
      message("No input received during last poll, shutting down...")
      quit("no")
    }

    # parse data from state machine input
    input <- RJSONIO::fromJSON(input_json, simplify = FALSE)

    # save logs to file
    debug_prefix <- file.path(input$experimentId, debug_timestamp)
    dump_folder <- file.path(DEBUG_PATH, debug_prefix)
    futile.logger::flog.appender(
      futile.logger::appender.tee(file.path(dump_folder, "logs.txt"))
    )

    heartbeat_proc <- start_heartbeat(task_token, pipeline_config$aws_config)

    tryCatchLog(
      {
        # Refresh pipeline_config with the new task input
        pipeline_config <- load_config(input$server)

        wrapper(input, pipeline_config)

        message("Send task success\n------\n")
        states$send_task_success(
          taskToken = task_token,
          output = "{}"
        )
      },
      error = function(e) {
        futile.logger::flog.error("🚩 ---------")
        sample_text <- ifelse(is.null(input$sampleUuid),
                              "",
                              paste0(" for sample ", input$sampleUuid)
        )

        error_txt <- paste0(
          "R error at filter step ",
          input$taskName, sample_text, "! : ", e$message
        )

        message(error_txt)
        states$send_task_failure(
          taskToken = task_token,
          error = get_user_error(e$message),
          cause = error_txt
        )

        send_pipeline_fail_update(pipeline_config, input, error_txt)
        message("Sent task failure to state machine task: ", task_token)

        if (pipeline_config$cluster_env != "development") {
          upload_debug_folder_to_s3(debug_prefix, pipeline_config)
        }

        message("recovered from error:", e$message)
      },
      write.error.dump.file = TRUE,
      write.error.dump.folder = dump_folder
    )

    # kill heartbeat process
    heartbeat_proc$kill()
  }
}
