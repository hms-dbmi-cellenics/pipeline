remove_bucket_folder <- function(pipeline_config, bucket, folder) {
  keys_to_remove <- list()

  s3 <- paws::s3(config = pipeline_config$aws_config)
  object_list <- s3$list_objects(bucket, Prefix = paste0(folder, "/"))
  for (object in object_list$Contents) {
    keys_to_remove <- append(keys_to_remove, object$Key)
    s3$delete_object(
      Bucket = bucket,
      Key = object$Key
    )
  }

  message("Removed files from ", bucket, ": ", paste0(keys_to_remove, sep=' '))
}

remove_cell_ids <- function(pipeline_config, experiment_id) {
  remove_bucket_folder(pipeline_config, pipeline_config$cells_id_bucket, experiment_id)
}

upload_cells_id <- function(pipeline_config, object_key, cells_id) {
  object_data <- tempfile()
  saveRDS(cells_id, file = object_data)
  put_object_in_s3(pipeline_config, pipeline_config$cells_id_bucket, object_data, object_key)
  return(object_key)
}

load_processed_scdata <- function(s3, pipeline_config, experiment_id) {
  bucket <- pipeline_config$processed_bucket
  message("Loading processed scdata")
  message(bucket)
  message(paste(experiment_id, "r.rds", sep = "/"))

  c(body, ...rest) %<-% s3$get_object(
    Bucket = bucket,
    Key = paste(experiment_id, "r.rds", sep = "/")
  )
  conn <- gzcon(rawConnection(body))
  obj <- readRDS(conn)
  close(conn)
  return(obj)
}

# get_nnzero will return how many non-zero counts the count matrix has
# it is used to order samples according to their size
get_nnzero <- function (x) {
  return(length(x@assays[["RNA"]]@counts@i))
}

order_by_size <- function(scdata_list) {
    return(scdata_list[order(sapply(scdata_list, get_nnzero))])
}

load_source_scdata_list <- function (s3, pipeline_config, experiment_id) {
  bucket <- pipeline_config$source_bucket
  objects <- s3$list_objects(
    Bucket = bucket,
    Prefix = experiment_id
  )
  samples <- objects$Contents

  scdata_list <- list()
  for (sample in samples) {
    key <- sample$Key

    c(body, ...rest) %<-% s3$get_object(
      Bucket = bucket,
      Key = paste(key, sep = "/")
    )
    conn <- gzcon(rawConnection(body))
    obj <- readRDS(conn)
    sample_id <- strsplit(key, "/")[[1]][[2]]
    scdata_list[[sample_id]] <- obj
    # close connection explicitly or R will run out of available connections and fail
    close(conn)
  }

  # order samples according to their size to make the merge independent of samples order in the UI
  return(order_by_size(scdata_list))
}

# reload_data_from_s3 will reload:
# * scdata_list for all steps before integration (included)
# * scdata file for all steps after data integration
reload_data_from_s3 <- function(pipeline_config, experiment_id, task_name, tasks) {
  task_names <- names(tasks)
  integration_index <- match("dataIntegration", task_names)
  s3 <- paws::s3(config = pipeline_config$aws_config)

  # TODO: remove if block
  # this never runs, because embed and cluster runs in the worker if modified.
  # If the task is after data integration, we need to get scdata from processed_matrix
  if (match(task_name, task_names) > integration_index) {
    return(load_processed_scdata(s3, pipeline_config, experiment_id))
  }

  # Otherwise, return scdata_list
  return(load_source_scdata_list(s3, pipeline_config, experiment_id))

}

load_cells_id_from_s3 <- function(pipeline_config, experiment_id, task_name, tasks, samples) {
  s3 <- paws::s3(config = pipeline_config$aws_config)
  object_list <- s3$list_objects(
                    Bucket = pipeline_config$cells_id_bucket,
                    Prefix = paste0(experiment_id, "/", task_name, "/")
                  )
  message(pipeline_config$cells_id_bucket)
  message(paste(experiment_id, "r.rds", sep = "/"))
  cells_id <- list()
  message("Total of ", length(object_list$Contents), " samples.")
  task_names <- names(tasks)
  integration_index <- match("dataIntegration", task_names)

  # after data integration the cell ids are no longer used because there is no filtering
  # so they are not loaded
  if (match(task_name, task_names) <= integration_index) {
    for (object in object_list$Contents) {
      key <- object$Key
      sample_id <- tools::file_path_sans_ext(basename(key))

      if (!sample_id %in% samples) {
        message("Unexpected filtered ids object for sample ", sample_id, ". Filtered cell ids' objects for removed samples
        should be removed in the first step of the QC pipeline after gem2s is triggered (due to sample removal).
        Removing it now.")
        s3$delete_object(
          Bucket = pipeline_config$cells_id_bucket,
          Key = key
        )
        next
      }

      c(body, ...rest) %<-% s3$get_object(
        Bucket = pipeline_config$cells_id_bucket,
        Key = key
      )
      id_file <- tempfile()
      writeBin(body, con = id_file)
      sample <- readRDS(id_file)
      cells_id[[sample_id]] <- sample[[sample_id]]
      message("Sample ", sample_id, " with ", length(cells_id[[sample_id]]), " cells")
    }
  }
  return(cells_id)
}

build_qc_response <- function(id, input, error, pipeline_config) {
  msg <- list(
      experimentId = input$experimentId,
      taskName = input$taskName,
      input = input,
      response = list(
        error = error
      ),
      pipelineVersion = pipeline_version,
      apiUrl = pipeline_config$api_url
    )

  if (!is.null(id)) {
    msg$output <- list(
      bucket = pipeline_config$results_bucket,
      key = id
    )
  }

  return(msg)
}

send_output_to_api <- function(pipeline_config, input, plot_data_keys, output) {
  c(config, plot_data = plotData) %<-% output

  config <- config[!names(config) %in% c("auth_JWT", "api_url")]

  # upload output
  id <- ids::uuid()
  output <- RJSONIO::toJSON(
    list(
      config = config,
      plotDataKeys = plot_data_keys,
      plotData = plot_data
    )
  )

  message("Uploading results to S3 bucket", pipeline_config$results_bucket, " at key ", id, "...")
  put_object_in_s3(pipeline_config, pipeline_config$results_bucket, charToRaw(output), id)

  message("Sending to SNS topic ", pipeline_config$sns_topic)
  sns <- paws::sns(config = pipeline_config$aws_config)

  message("Building the message")
  msg <- build_qc_response(id, input, FALSE, pipeline_config)

  message("Publishing the message")
  result <- sns$publish(
    Message = RJSONIO::toJSON(msg),
    TopicArn = pipeline_config$sns_topic,
    MessageAttributes = list(
      type = list(
        DataType = "String",
        StringValue = "PipelineResponse",
        BinaryValue = NULL
      )
    )
  )
  message("Done publishing")

  return(result$MessageId)
}

send_gem2s_update_to_api <- function(pipeline_config, experiment_id, task_name, data, input) {
  message("Sending to SNS topic ", pipeline_config$sns_topic)
  sns <- paws::sns(config = pipeline_config$aws_config)
  job_id <- Sys.getenv("AWS_BATCH_JOB_ID", unset = "")

  # TODO -REMOVE DUPLICATE AUTHJWT IN RESPONSE
  msg <- c(
    data,
    taskName = list(task_name),
    experimentId = list(experiment_id),
    jobId = list(job_id),
    authJWT = list(input$auth_JWT),
    input = list(input),
    apiUrl = pipeline_config$api_url
  )

  result <- sns$publish(
    Message = RJSONIO::toJSON(msg),
    TopicArn = pipeline_config$sns_topic,
    MessageAttributes = list(
      type = list(
        DataType = "String",
        StringValue = "GEM2SResponse",
        BinaryValue = NULL
      )
    )
  )

  return(result$MessageId)
}

send_pipeline_fail_update <- function(pipeline_config, input, error_message) {
  sns <- paws::sns(config = pipeline_config$aws_config)

  string_value <- ""
  response <- ""

  process_name <- input$processName
  if (process_name == "qc") {
    string_value <- "PipelineResponse"

    global_env_config <- NULL
    # We don't have dynamic config set for non-sample steps so we can ignore it in that cases
    if (!is.null(input$sampleUuid)) {
      config_key <- paste0("config-", input$taskName, "-", input$sampleUuid)
      global_env_config <- get0(config_key, envir = globalenv(), ifnotfound = NULL)
    }

    # If the step didn't backup any config, don't upload anything
    if (is.null(global_env_config)) {
      response <- build_qc_response(NULL, input, process_name, pipeline_config)
    } else {
      # upload output
      id <- ids::uuid()
      output <- RJSONIO::toJSON(
        list(
          config = global_env_config
        )
      )

      message("Uploading config to S3 bucket", pipeline_config$results_bucket, " at key ", id, "...")
      put_object_in_s3(pipeline_config, pipeline_config$results_bucket, charToRaw(output), id)

      response <- build_qc_response(id, input, process_name, pipeline_config)
    }


  } else if (process_name == "gem2s") {
    string_value <- "GEM2SResponse"

    # TODO - REMOVE THE DUPLICATE EXPERIMENT ID FROM INPUT RESPONSE
    response <- list(
      experimentId <- input$experimentId,
      taskName <- input$taskName,
      input <- input,
      apiUrl <- pipeline_config$api_url,
      response <- list(
        error <- process_name
      )
    )
  } else {
    message(paste("Invalid process_name given: ", process_name))
    return()
  }

  result <- sns$publish(
    Message = RJSONIO::toJSON(response),
    TopicArn = pipeline_config$sns_topic,
    MessageAttributes = list(
      type = list(
        DataType = "String",
        StringValue = string_value,
        BinaryValue = NULL
      )
    )
  )
}

send_plot_data_to_s3 <- function(pipeline_config, experiment_id, output) {
  plot_data <- output$plotData

  s3 <- paws::s3(config = pipeline_config$aws_config)

  plot_names <- names(plot_data)

  plot_data_keys <- list()

  for (plot_data_name in names(plot_data)) {
    id <- ids::uuid()

    plot_data_keys[[plot_data_name]] <- id

    output <- RJSONIO::toJSON(
      list(
        plotData = plot_data[[plot_data_name]]
      )
    )

    message("Uploading plotData to S3 bucket", pipeline_config$plot_data_bucket, " at key ", id, "...")

    tags <- paste0("experimentId=", experiment_id , "&plotUuid=" , plot_data_name)
    put_object_in_s3(pipeline_config, pipeline_config$plot_data_bucket, charToRaw(output), id, tags)
  }

  return(plot_data_keys)
}

upload_matrix_to_s3 <- function(pipeline_config, experiment_id, data) {
  object_key <- paste0(experiment_id, "/r.rds")

  count_matrix <- tempfile()
  saveRDS(data, file = count_matrix)

  message("Count matrix file size : ", format(object.size(count_matrix)))
  message("Uploading updated count matrix to S3 bucket ", pipeline_config$processed_bucket, " at key ", object_key, "...")

  put_object_in_s3_multipart(pipeline_config, pipeline_config$processed_bucket, count_matrix, object_key)

  return(object_key)
}

upload_debug_folder_to_s3 <- function(debug_prefix, pipeline_config) {
  fnames <- list.files(file.path(DEBUG_PATH, debug_prefix))
  bucket <- pipeline_config$debug_bucket

  message("Uploading logs and dump file to S3 bucket ", bucket, " with prefix ", debug_prefix, "...")
  for (fname in fnames) {
    fpath <- file.path(DEBUG_PATH, debug_prefix, fname)
    key <- file.path(debug_prefix, fname)
    put_object_in_s3_multipart(pipeline_config, bucket, fpath, key)
  }

  return(NULL)
}

put_object_in_s3 <- function(pipeline_config, bucket, object, key, tagging = NULL) {
  message(sprintf("Putting %s in %s", key, bucket))
  s3 <- paws::s3(config = pipeline_config$aws_config)

  retry_count <- 0
  max_retries <- 2

  while (retry_count < max_retries) {
    tryCatch({
      response <- s3$put_object(
        Bucket = bucket,
        Key = key,
        Body = object,
        Tagging = switch(!is.null(tagging), tagging)
      )
      message("Object successfully uploaded to S3.")
      return(response)
    }, error = function(e) {
      retry_count <<- retry_count + 1
      message(sprintf("Upload failed. Retrying (%d/%d)...", retry_count, max_retries))
      # exponential back off, preventing sending a new request too fast
      Sys.sleep(2 ^ retry_count)
    })
  }
  stop("Failed to upload object to S3 after maximum number of retries.")
}
#' Upload a file to S3 using multipart upload
#'
#' @param pipeline_config A Paws S3 config object, e.g. from `paws::s3()`.
#' @param object The path to the file to be uploaded.
#' @param bucket The name of the S3 bucket to be uploaded to, e.g. `my-bucket`.
#' @param key The name to assign to the file in the S3 bucket, e.g. `path/to/file`.
put_object_in_s3_multipart <- function(pipeline_config, bucket, object, key) {
  message(sprintf("Putting %s in %s from object %s", key, bucket, object))

  s3 <- paws::s3(config = pipeline_config$aws_config)
  multipart <- s3$create_multipart_upload(
    Bucket = bucket,
    Key = key
  )
  resp <- NULL
  on.exit({
    if (is.null(resp) || inherits(resp, "try-error")) {
      s3$abort_multipart_upload(
        Bucket = bucket,
        Key = key,
        UploadId = multipart$UploadId
      )
    }
  })
  message("Uploading multiparts")
  resp <- try({
    parts <- upload_multipart_parts(s3, bucket, object, key, multipart$UploadId)
    s3$complete_multipart_upload(
      Bucket = bucket,
      Key = key,
      MultipartUpload = list(Parts = parts),
      UploadId = multipart$UploadId
    )
  })
  return(resp)
}

# The limit for part numbers is 10000, so we have a limit on 50GB objects.
# This can be increased by changing the number of megabytes in the part_size parameter
upload_multipart_parts <- function(s3, bucket, object, key, upload_id) {
  file_size <- file.size(object)
  megabyte <- 2^20
  part_size <- 5 * megabyte
  num_parts <- ceiling(file_size / part_size)

  con <- base::file(object, open = "rb")
  on.exit({
    close(con)
  })
  parts <- list()
  for (i in 1:num_parts) {
    part <- readBin(con, what = "raw", n = part_size)
    part_resp <- s3$upload_part(
      Body = part,
      Bucket = bucket,
      Key = key,
      PartNumber = i,
      UploadId = upload_id
    )
    parts <- c(parts, list(list(ETag = part_resp$ETag, PartNumber = i)))
  }

  return(parts)
}


#' Load cellsets object from s3 and flatten it for easier manipulation in R
#'
#' @param s3 paws::s3 object
#' @param pipeline_config list
#' @param experiment_id character
#'
#' @return cellsets list
#' @export
#'
load_cellsets <- function(s3, pipeline_config, experiment_id) {
  message("loading cellsets file")

  bucket <- pipeline_config$cell_sets_bucket

  c(body, ...rest) %<-% s3$get_object(
    Bucket = bucket,
    Key = experiment_id
  )

  obj <- jsonlite::fromJSON(rawConnection(body), flatten = T)
  return(obj)

}

#' Load cellsets object flattened (as is done by load_cellsets)
#' and change it back to the original shape
#' This is useful to prepare the cell sets to upload back into s3
#'
#' @param s3 paws::s3 object
#' @param pipeline_config list
#' @param experiment_id character
#'
#' @return cellsets list
#' @export
#'
unflatten_cell_sets <- function(cell_sets) {
  formatted <- list()

  for (i in seq_along(cell_sets$key)) {
    cell_class <- list(
      key = cell_sets$key[[i]],
      name = cell_sets$name[[i]],
      rootNode = cell_sets$rootNode[[i]],
      type = cell_sets$type[[i]],
      children = list()
    )

    for (j in seq_along(cell_sets$children[[i]]$key)) {
      children <- cell_sets$children[[i]]

      cell_set <- list(
        key = children$key[[j]],
        name = children$name[[j]],
        rootNode = FALSE,
        color = children$color[[j]],
        # Take the type from the parent, some of the child cell sets don't have "type"
        # Some others do, but all parents do
        type = cell_sets$type[[i]],
        cellIds = ensure_is_list_in_json(children$cellIds[[j]])
      )

      cell_class$children <- append(cell_class$children, list(cell_set))
    }

    formatted <- append(formatted, list(cell_class))
  }

  return(formatted)
}


#' Bind columns not creating rows if there's an empty data.table
#'
#' `cbind` on `data.table` adds a row if binding an empty data.table to a non-empty
#' one. We do not want that behavior when parsing cellsets, because it implies
#' the "creation" of a cell that does not exists (i.e. when binding scratchpad
#' cellsets slots of an experiment without custom cellsets)
#'
#' @param dt data.table
#' @param ... columns to add
#'
#' @return data.table with new columns
#' @export
#'
safe_cbind <- function(dt, ...) {
  if (nrow(dt) > 0) {
    dt <- cbind(dt, ...)
  }
  return(dt)
}


#' add cellset type column to cellsets data.table
#'
#' helper to correctly name the cellset type column. some cellsets already
#' contain a "type" slot, which complicates matters, so we chose `cellset_type`,
#'
#' @param dt data.table
#' @param col string of corresponding cellset type
#'
#' @return data.table with cellset_type column
#' @export
#'
cbind_cellset_type <- function(dt, col) {
  dt <- safe_cbind(dt, cellset_type = col)
}


is_uuid <- function(x) {
  uuid_regex <- "^\\b[0-9a-f]{8}-[0-9a-f]{4}-[1-5][0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}\\b$"
  return(grepl(uuid_regex, x))
}


get_cellset_type <- function(key, type) {
  cellset_type <- switch(key,
                         louvain = "cluster",
                         scratchpad = "scratchpad",
                         sample = "sample",
                         ifelse(is_uuid(key), "sctype", "metadata")
  )

  return(cellset_type)
}


#' Parse cellsets object to data.table
#'
#' Gets the cellsets list and converts it to a tidy data.table
#'
#' @param cellsets list
#'
#' @return data.table of cellset keys, names and corresponding cell_ids
#' @export
#'
parse_cellsets <- function(cellsets) {

  dt_list <- cellsets$cellSets$children
  cellset_types <- purrr::map2_chr(cellsets$cellSets$key, cellsets$cellSets$type, get_cellset_type)

  lapply(dt_list, data.table::setDT)

  dt_list <- purrr::map2(dt_list, cellset_types, cbind_cellset_type)

  # fill columns in case there are empty cellset classes
  dt <- data.table::rbindlist(dt_list, fill = TRUE)

  # unnest, and change column name
  dt <- dt[, setNames(.(unlist(cellIds)), "cell_id"), by = .(key, name, cellset_type)]
  data.table::setnames(dt, "cellset_type", "type")
  return(dt)
}

get_s3_rds <- function(bucket, key, aws_config) {
  s3 <- paws::s3(config = aws_config)

  c(body, ...rest) %<-% s3$get_object(
    Bucket = bucket,
    Key = key
  )

  conn <- gzcon(rawConnection(body))
  object <- readRDS(conn)
  
  close(conn)

  return(object)
}

put_s3_rds <- function(bucket, key, aws_config, rds) {
  s3 <- paws::s3(config = aws_config)

  file <- tempfile()
  saveRDS(rds, file = file)

  s3$put_object(
    Bucket = bucket,
    Key = key,
    Body = file
  )
}
