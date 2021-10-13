reload_scdata_from_s3 <- function(pipeline_config, experiment_id) {
  s3 <- paws::s3(config = pipeline_config$aws_config)
  message(pipeline_config$source_bucket)
  message(paste(experiment_id, "r.rds", sep = "/"))

  c(body, ...rest) %<-% s3$get_object(
    Bucket = pipeline_config$source_bucket,
    Key = paste(experiment_id, "r.rds", sep = "/")
  )
  obj <- readRDS(rawConnection(body))
  return(obj)
}

send_output_to_api <- function(pipeline_config, input, plot_data_keys, output) {
  c(config, plot_data = plotData) %<-% output

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

  msg <- list(
    experimentId = input$experimentId,
    input = input,
    output = list(
      bucket = pipeline_config$results_bucket,
      key = id
    ),
    response = list(
      error = FALSE
    )
  )

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

  return(result$MessageId)
}

send_gem2s_update_to_api <- function(pipeline_config, experiment_id, task_name, data, auth_JWT) {
  message("Sending to SNS topic ", pipeline_config$sns_topic)
  sns <- paws::sns(config = pipeline_config$aws_config)

  msg <- c(data, taskName = list(task_name), experimentId = list(experiment_id), authJWT = list(auth_JWT))

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

  # todo - remove the duplicate experiment id from input response
  error_msg$experimentId <- input$experimentId
  error_msg$response$error <- process_name
  error_msg$input <- input
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
    s3$put_object(
      Bucket = pipeline_config$plot_data_bucket,
      Key = id,
      Body = charToRaw(output)
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
