require("RJSONIO")
require("paws")
require("zeallot")
require("ids")

reload_scdata_from_s3 <- function(pipeline_config, experiment_id) {
    s3 <- paws::s3(config=pipeline_config$aws_config)
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
    s3 <- paws::s3(config=pipeline_config$aws_config)
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
    sns <- paws::sns(config=pipeline_config$aws_config)

    msg = list(
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

send_gem2s_update_to_api <- function(pipeline_config, experiment_id, task_name, data) {
    message("Sending to SNS topic ", pipeline_config$sns_topic)
    sns <- paws::sns(config=pipeline_config$aws_config)

    msg = c(data, taskName = list(task_name))
    msg = c(msg, experimentId = list(experiment_id))

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

send_pipeline_fail_update <- function(pipeline_config, experiment_id, process_name, error_message) {
    error_msg = list()
    error_msg$experimentId = experiment_id
    error_msg$response$error = error_message

    sns <- paws::sns(config=pipeline_config$aws_config)

    string_value = ""
    if (process_name == 'qc') {
        string_value = "PipelineResponse"
    } else if (process_name == 'gem2s') {
        string_value = "GEM2SResponse"
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

    s3 <- paws::s3(config=pipeline_config$aws_config)

    plot_names <- names(plot_data)

    plot_data_keys = list()

    for(plot_data_name in names(plot_data)){
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
    s3 <- paws::s3(config=pipeline_config$aws_config)

    object_key <- paste0(experiment_id, '/r.rds')

    count_matrix <- tempfile()
    saveRDS(data, file=count_matrix)

    message("Uploading updated count matrix to S3 bucket ", pipeline_config$processed_bucket, " at key ", object_key, "...")
    s3$put_object(
        Bucket = pipeline_config$processed_bucket,
        Key = object_key,
        Body = count_matrix
    )

    return(object_key)
}

put_object_in_s3 <- function(pipeline_config, bucket, object, key) {

    print(sprintf("Putting %s in %s", key, bucket))

    s3 <- paws::s3(config=pipeline_config$aws_config)
    s3$put_object(
        Bucket = bucket,
        Key = key,
        Body = object
    )
}

