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

reload_scdata_from_s3 <- function (s3, pipeline_config, experiment_id) {
  bucket <- pipeline_config$processed_bucket
  message(bucket)
  message(paste(experiment_id, "r.rds", sep = "/"))

  c(body, ...rest) %<-% s3$get_object(
    Bucket = bucket,
    Key = paste(experiment_id, "r.rds", sep = "/")
  )
  obj <- readRDS(gzcon(rawConnection(body)))
  return(obj)
}

# get_nnzero will return how many non-zero counts the count matrix has
# it is used to order samples according to their size
get_nnzero <- function (x) {
  return(length(x@assays[["RNA"]]@counts@i))
}

order_by_size <- function(scdata_list) {
    return(scdata_list <- scdata_list[ order( sapply(scdata_list, get_nnzero)) ])
}

reload_scdata_list_from_s3 <- function (s3, pipeline_config, experiment_id) {
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
    obj <- readRDS(gzcon(rawConnection(body)))
    sample_id <- strsplit(key, "/")[[1]][[2]]
    scdata_list[[sample_id]] <- obj
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

  # If the task is after data integration, we need to get scdata from processed_matrix
  if (match(task_name, task_names) > integration_index) {
    return(reload_scdata_from_s3(s3, pipeline_config, experiment_id))
  }

  # Otherwise, return scdata_list
  return(reload_scdata_list_from_s3(s3, pipeline_config, experiment_id))

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

send_output_to_api <- function(pipeline_config, input, plot_data_keys, output) {
  c(config, plot_data = plotData) %<-% output

  config <- config[!names(config) %in% c("auth_JWT", "api_url")]

  # upload output
  s3 <- paws::s3(config = pipeline_config$aws_config)
  id <- ids::uuid()
  output <- RJSONIO::toJSON(
    list(
      config = config,
      plotDataKeys = plot_data_keys,
      plotData = plot_data
    )
  )

  message("Uploading results to S3 bucket", pipeline_config$results_bucket, " at key ", id, "...")
  s3$put_object(
    Bucket = pipeline_config$results_bucket,
    Key = id,
    Body = charToRaw(output)
  )

  message("Sending to SNS topic ", pipeline_config$sns_topic)
  sns <- paws::sns(config = pipeline_config$aws_config)

  message("Building the message")
  msg <- list(
    experimentId = input$experimentId,
    taskName = input$taskName,
    input = input,
    output = list(
      bucket = pipeline_config$results_bucket,
      key = id
    ),
    response = list(
      error = FALSE
    ),
    pipelineVersion = pipeline_version,
    apiUrl = pipeline_config$public_api_url
  )

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
  # TODO -REMOVE DUPLICATE AUTHJWT IN RESPONSE
  msg <- c(
    data, 
    taskName = list(task_name),
    experimentId = list(experiment_id),
    authJWT = list(input$auth_JWT),
    input = list(input),
    apiUrl = pipeline_config$public_api_url
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
  process_name <- input$processName

  error_msg <- list()

  # TODO - REMOVE THE DUPLICATE EXPERIMENT ID FROM INPUT RESPONSE

  error_msg$experimentId <- input$experimentId
  error_msg$taskName <- input$taskName
  error_msg$response$error <- process_name
  error_msg$input <- input
  error_msg$apiUrl <- pipeline_config$public_api_url
  sns <- paws::sns(config = pipeline_config$aws_config)

  string_value <- ""
  if (process_name == "qc") {
    string_value <- "PipelineResponse"
  } else if (process_name == "gem2s") {
    string_value <- "GEM2SResponse"
  } else {
    message(paste("Invalid process_name given: ", process_name))
    return()
  }

  result <- sns$publish(
    Message = RJSONIO::toJSON(error_msg),
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

    s3$put_object(
      Bucket = pipeline_config$plot_data_bucket,
      Key = id,
      Body = charToRaw(output),
      Tagging = tags
    )
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

put_object_in_s3 <- function(pipeline_config, bucket, object, key) {
  message(sprintf("Putting %s in %s", key, bucket))

  s3 <- paws::s3(config = pipeline_config$aws_config)
  s3$put_object(
    Bucket = bucket,
    Key = key,
    Body = object
  )
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
