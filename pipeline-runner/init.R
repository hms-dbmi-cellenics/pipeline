# if Seurat not attached can cause errors when accessing metadata
library(Seurat)
library(zeallot)
library(tryCatchLog)
library(futile.logger)
library(magrittr)
library(callr)

# time stamp used for directory to store log/dump files in event of error
debug_timestamp <- format(Sys.time(), format = "%Y-%m-%d_at_%H-%M-%OS3")

# increase maxSize from the default of 500MB to 32GB
options(future.globals.maxSize = 32 * 1024 * 1024^2)

# show line numbers for tryCatchLog
options(keep.source.pkgs = TRUE)

for (f in list.files('R', '.R$', full.names = TRUE)) source(f, keep.source = TRUE)
load('R/sysdata.rda') # constants

buildActivityArn <- function(aws_region, aws_account_id, activity_id) {
    if(is.na(activity_id)) {
        return(NA)
    }

    activity_arn <- sprintf("arn:aws:states:%s:%s:activity:%s", aws_region, aws_account_id, activity_id)
    return(activity_arn)
}

load_config <- function(development_aws_server) {
    label_path <- "/etc/podinfo/labels"
    aws_account_id <- Sys.getenv("AWS_ACCOUNT_ID", unset="242905224710")
    aws_region <- Sys.getenv("AWS_DEFAULT_REGION", unset="eu-west-1")

    activity_arn <- NA

    repeat {
        # if /etc/podinfo/labels exists we are running in a remote aws pod
        if(file.exists(label_path)) {
            labels <- read.csv(label_path, sep="=", row.names=1, header=FALSE)
            activity_id <- labels["activityId", ]
            activity_arn <- buildActivityArn(aws_region, aws_account_id, activity_id)
        }

        # if we didn't find an activity in podinfo try to get the from the env (means we are running locally)
        if(is.na(activity_arn)) {
            activity_arn <- Sys.getenv("ACTIVITY_ARN", unset = NA)
        }

        if(is.na(activity_arn)) {
            message("No activity ARN label set yet, waiting...")
            Sys.sleep(5)
        } else {
            message(paste("Welcome to Biomage R pipeline, activity arn", activity_arn))
            break
        }
    }

    sandbox <- Sys.getenv("SANDBOX_ID", "default")
    config <- list(
        cluster_env = Sys.getenv("CLUSTER_ENV", "development"),
        sandbox_id = sandbox,
        aws_account_id = aws_account_id,
        aws_region = aws_region,
        pod_name = Sys.getenv("K8S_POD_NAME", "local"),
        activity_arn = activity_arn,
        api_url = paste0("http://api-",sandbox,".api-",sandbox,".svc.cluster.local:3000"),
        api_version = "v2",
        debug_config = list(
            step = Sys.getenv("DEBUG_STEP", ""),
            path = Sys.getenv("DEBUG_PATH", "")
        )
    )

    config[["aws_config"]] <- list(region = config$aws_region)

    # running in linux needs the IP of the host to work. If it is set as an environment variable (by makefile) honor it instead of the
    # provided parameter
    overriden_server = Sys.getenv("HOST_IP", "")
    if (overriden_server != "") {
        development_aws_server = overriden_server
    }

    if(config$cluster_env == 'development') {
        config$api_url <- sprintf("http://%s:3000", development_aws_server)
        config$aws_config[['endpoint']] <- sprintf("http://%s:4566", development_aws_server) # DOCKER_GATEWAY_HOST
        config$aws_config[['credentials']] <- list(
            creds = list(
                access_key_id = "mock-access-key",
                secret_access_key = "mock-secure-acces-key"
            )
        )
        config$aws_account_id <- '000000000000'

        # This fixes a bug where paws would try to connect to InfraMock as if it was
        # a proper AWS endpoint, i.e. to http://bucket-name.host.docker.internal
        assignInNamespace("update_endpoint_for_s3_config", function(request) {
            return(request)
        }, ns="paws.common")

    }

    config[["originals_bucket"]] <- paste("biomage-originals", config$cluster_env, config$aws_account_id, sep = "-")
    config[["source_bucket"]] <- paste("biomage-source", config$cluster_env, config$aws_account_id, sep = "-")
    config[["processed_bucket"]] <- paste("processed-matrix", config$cluster_env, config$aws_account_id, sep = "-")
    config[["results_bucket"]] <- paste("worker-results", config$cluster_env, config$aws_account_id, sep = "-")
    config[["cells_id_bucket"]] <- paste("biomage-filtered-cells", config$cluster_env, config$aws_account_id, sep = "-")
    config[["plot_data_bucket"]] <- paste("plots-tables", config$cluster_env, config$aws_account_id, sep = "-")
    config[["cell_sets_bucket"]] <- paste("cell-sets", config$cluster_env, config$aws_account_id, sep = "-")
    config[["debug_bucket"]] <- paste("biomage-pipeline-debug", config$cluster_env, config$aws_account_id, sep = "-")
    config[["sns_topic"]] <- paste(
        paste("arn:aws:sns", config$aws_region, config$aws_account_id, "work-results", sep = ":"),
        config$cluster_env, config$sandbox_id, config$api_version, sep = "-"
    )

    return(config)
}

run_processing_step <- function(scdata, config, tasks,task_name, cells_id,sample_id, debug_config) {
    if (!task_name %in% names(tasks)) stop('Invalid task: ', task_name)


    handle_debug(scdata, config, task_name, sample_id, debug_config)

    # print info
    task <- tasks[[task_name]]
    message("Running: ", task_name)
    message("Config:")
    str(config)

    # run task and time it
    tstart <- Sys.time()

    out <- task(scdata, config, sample_id, cells_id, task_name)
    ttask <- format(Sys.time()-tstart, digits = 2)
    message("⏱️ Time to complete ", task_name, " for sample ", sample_id, ": ", ttask, '\n')

    return(out)
}

#
# input needs:
# sampleIds for all samples
# sampleNames for all samples
# projectId
# pipeline_config as defined by load_config
#

run_gem2s_step <- function(task_name, input, pipeline_config, prev_out) {
    # list of task functions named by task name
    tasks <- list(
        "downloadGem" = download_user_files,
        "preproc" = load_user_files,
        "emptyDrops" = run_emptydrops,
        "doubletScores" = score_doublets,
        "createSeurat" = create_seurat,
        "prepareExperiment" = prepare_experiment,
        "uploadToAWS" = upload_to_aws
    )

    if (!task_name %in% names(tasks)) stop("Invalid task name given: ", task_name)
    task <- tasks[[task_name]]

    tstart <- Sys.time()
    res <- task(input, pipeline_config, prev_out)
    ttask <- format(Sys.time()-tstart, digits = 2)
    message("⏱️ Time to complete ", task_name, ": ", ttask, '\n')

    return(res)
}


call_gem2s <- function(task_name, input, pipeline_config) {
    experiment_id <- input$experimentId

    if (!exists("prev_out")) {
        remove_cell_ids(pipeline_config, experiment_id)
        assign("prev_out", NULL, pos = ".GlobalEnv")
    }

    check_input(input)

    c(data, task_out) %<-% run_gem2s_step(task_name, input, pipeline_config, prev_out)
    assign("prev_out", task_out, pos = ".GlobalEnv")

    message_id <- send_gem2s_update_to_api(pipeline_config, experiment_id, task_name, data, input)

    return(message_id)
}

#
# call_data_processing
# Runs step @task_name of the data processing pipeline, send plot data to s3 and output message to api.
# IN task_name: str, name of the step.
# IN input: json str, config settings for all samples, current sample uuid, uploadCountMatrix (whether or not to upload after step)
# IN pipeline_config: json str, environment config resulting from the load_config function
#
call_data_processing <- function(task_name, input, pipeline_config) {
    experiment_id <- input$experimentId
    config <- input$config
    upload_count_matrix <- input$uploadCountMatrix
    sample_id <- input$sampleUuid
    debug_config <- pipeline_config$debug_config

    if (sample_id != "") {
        config <- config[[sample_id]]
        input$config <- config
    }
    # vector of task functions named by task name
    tasks <- list(
        'classifier' = filter_emptydrops,
        'cellSizeDistribution' = filter_low_cellsize,
        'mitochondrialContent' = filter_high_mito,
        'numGenesVsNumUmis' = filter_gene_umi_outlier,
        'doubletScores' = filter_doublets,
        'dataIntegration' = integrate_scdata,
        'configureEmbedding' = embed_and_cluster
    )

    # need this for embed_and_cluster
    config$api_url <- pipeline_config$api_url
    config$auth_JWT <- input$authJWT

    if (!exists("scdata")) {
        message("No single-cell data has been loaded, reloading from S3...")

        # assign it to the global environment so we can
        # persist it across runs of the wrapper
        assign("scdata", reload_data_from_s3(pipeline_config, experiment_id, task_name, tasks), pos = ".GlobalEnv")

        message("Single-cell data loaded.")
    }

    if (!exists("cells_id")) {
        message("No filtered cell ids have been loaded, loading from S3...")
        if(task_name == names(tasks)[1]){
            assign("cells_id", generate_first_step_ids(scdata), pos = ".GlobalEnv")
        }else if(task_name %in% names(tasks)){
            samples <- names(scdata)
            assign("cells_id", load_cells_id_from_s3(pipeline_config, experiment_id, task_name, tasks, samples), pos = ".GlobalEnv")
        }else{
            stop("Invalid task name given: ", task_name)
        }
        message("Cells id loaded.")
    }

    # call function to run and update global variable
    c(
        data, new_ids,...rest_of_results
    ) %<-% run_processing_step(scdata, config, tasks, task_name, cells_id, sample_id, debug_config)

    assign("cells_id", new_ids, pos = ".GlobalEnv")

    task_names <- names(tasks)
    integration_index <- match("dataIntegration", task_names)
    task_index <- match(task_name, task_names)
    if(task_index < integration_index){
        message("Filtered cell ids from ", length(cells_id[[sample_id]]), " to ", length(new_ids[[sample_id]]))
        next_task <- names(tasks)[[task_index+1]]
        object_key <- paste0(experiment_id,"/",next_task,"/",sample_id,".rds")
        upload_cells_id(pipeline_config,object_key,cells_id)
    }

    # upload plot data result to S3
    tstart <- Sys.time()
    plot_data_keys <- send_plot_data_to_s3(pipeline_config, experiment_id, rest_of_results)

    # Upload count matrix data
    if(upload_count_matrix) {
        assign("scdata", data, pos = ".GlobalEnv")
        object_key <- upload_matrix_to_s3(pipeline_config, experiment_id, scdata)
        message('Count matrix uploaded to ', pipeline_config$processed_bucket, ' with key ',object_key)
    }

    # send result to API
    message_id <- send_output_to_api(pipeline_config, input, plot_data_keys, rest_of_results)
    ttask <- format(Sys.time()-tstart, digits = 2)
    message("⏱️ Time to upload ", task_name, " objects for sample ", sample_id, ": ", ttask)

    return(message_id)
}

#
# start_heartbeat(task_token, aws_config)
# IN task_token, aws_config
#
# Sends a heartbeat to the state machine every 'wait_time' seconds
# Once the task is completed the heartbeat will fail accordingly with a
# task timeout and exit the loop and a new heartbeat will be set up by next task.
# This method is invoked with `callr::r_bg` which creates a new process which does not inherit
# the current workspace or memory, only the provided parameters; that's why we need to
# reimport tryCatchLog & initialize states again.
#
start_heartbeat <- function(task_token, aws_config) {
    library(tryCatchLog)
    message("Starting heartbeat")
    states <- paws::sfn(config=aws_config)

    keep_running <- TRUE
    # amount of time to wait between heartbeats
    wait_time <- 30
    i <- 0
    while (keep_running) {
        tryCatchLog({
            states$send_task_heartbeat(
                taskToken = task_token
            )
            message("Heartbeat sent: ", i)
        },
        error = function(e) {
            message("Send task heartbeat failed: ", e$message)
            message("Stopping heartbeat after ", i+1)
            keep_running <- FALSE
        })
        i <- i + 1
        # sleep until next heartbeat
        Sys.sleep(wait_time)

    }
}

#
# Wrapper(input_json)
# IN input_json: json input from message. Input should have:
# taskname, server, extra config parameters.
#
# Calls the appropiate process: data processing pipeline or gem2s.
#
wrapper <- function(input, pipeline_config) {
    task_name <- input$taskName
    message("------\nStarting task: ", task_name, '\n')
    message("Input:")
    str(input)
    message("")

    # common to gem2s and data processing
    server <- input$server
    input <- input[names(input) != "server"]

    process_name <- input$processName

    if (process_name == 'qc') {
        message_id <- call_data_processing(task_name, input, pipeline_config)
    } else if (process_name == 'gem2s') {
        message_id <- call_gem2s(task_name, input, pipeline_config)
    } else {
        stop("Process name not recognized.")
    }

    return(message_id)
}

#
# Init()
# Loads configuration and repeats the wrapper call until no more messages are received.
#
init <- function() {
    pipeline_config <- load_config('host.docker.internal')
    message("Loaded pipeline config")
    states <- paws::sfn(config=pipeline_config$aws_config)
    message("Loaded step function")

    print(sessionInfo())

    flog.layout(layout.simple)
    flog.threshold(ERROR)

    message("Waiting for tasks")

    repeat {
        c(task_token, input_json) %<-% states$get_activity_task(
            activityArn = pipeline_config$activity_arn,
            workerName = pipeline_config$pod_name
        )

        if(is.null(task_token) || !length(task_token) || task_token == "") {
            message('No input received during last poll, shutting down...')
            quit('no')
        }

        # parse data from state machine input
        input <- RJSONIO::fromJSON(input_json, simplify = FALSE)
        message('Input json: ', input)

        # save logs to file
        debug_prefix <- file.path(input$experimentId, debug_timestamp)
        dump_folder <- file.path(DEBUG_PATH, debug_prefix)
        flog.appender(appender.tee(file.path(dump_folder, "logs.txt")))

        # start heartbeat as a different process in the background
        message("Starting heartbeat")
        # message inside r_bg will ONLY be printed into /tmp/[out|err]
        # to see them
        # 1. log into the R container
        # 2. cat /tmp/out or tail -f /tmp/out
        heartbeat_proc <- r_bg(func=start_heartbeat, args=list(
            task_token, pipeline_config$aws_config),
            stdout = "/tmp/out",
            stderr = "/tmp/err"
        )

        tryCatchLog({
                # Refresh pipeline_config with the new task input
                pipeline_config <- load_config(input$server)

                wrapper(input, pipeline_config)

                message('Send task success\n------\n')
                states$send_task_success(
                    taskToken = task_token,
                    output = "{}"
                )
        },
            error = function(e) {
                flog.error("🚩 ---------")
                sample_text <- ifelse(is.null(input$sampleUuid),
                                      "",
                                      paste0(" for sample ", input$sampleUuid))

                error_txt <- paste0("R error at filter step ",
                                    input$taskName, sample_text, "! : ", e$message)

                message(error_txt)
                states$send_task_failure(
                    taskToken = task_token,
                    error = "We had an issue while processing your data.",
                    cause = error_txt
                )

                send_pipeline_fail_update(pipeline_config, input, error_txt)
                message("Sent task failure to state machine task: ", task_token)

                if (pipeline_config$cluster_env != 'development') {
                    upload_debug_folder_to_s3(debug_prefix, pipeline_config)
                }

                message("recovered from error:", e$message)
            },
        write.error.dump.file = TRUE,
        write.error.dump.folder = dump_folder)

        # kill heartbeat process
        heartbeat_proc$kill()
    }
}

init()
