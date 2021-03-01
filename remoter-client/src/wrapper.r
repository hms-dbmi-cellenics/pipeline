require("RJSONIO")
require("paws")
require("zeallot")

reload_from_s3 <- function(experiment_id) {
    s3 <- paws::s3()
    bucket_name <- paste("biomage-source", "staging", sep = "-")

    c(body, ...rest) %<-% s3$get_object(
        Bucket = bucket_name,
        Key = paste(experiment_id, "r.rds", sep = "/")
    )

    obj <- readRDS(rawConnection(body))
    return(obj)
}


run_step <- function(task_name, scdata, config) {
    switch(task_name,
        test_fn = {
            import::from("test_fn.r", task)
        },
        stop(paste("Invalid task name given:", task_name))
    )

    out <- task(scdata, config)
    return(out)
}

send_output_to_api <- function(input, output) {
    c(config, plot_data = plotData) %<-% output

    cluster_env <- Sys.getenv("CLUSTER_ENV", "staging")
    sandbox_id <- Sys.getenv("SANDBOX_ID", "demo")
    topic_name <- paste("work-results", cluster_env, sandbox_id, sep = "-")

    message("Sending to SNS topic", topic_name)
    sns <- paws::sns()
    result <- sns$publish(
        Message = RJSONIO::toJSON(
            list(
                input = input,
                config = config,
                plot_data = plot_data
            )
        ),
        TopicArn = paste(
            "arn:aws:sns:eu-west-1:242905224710",
            topic_name,
            sep = ":"
        ),
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


wrapper <- function(input_json) {
    # Get data from state machine input.
    input <- RJSONIO::fromJSON(input_json)
    c(
        experiment_id = experimentId,
        task_name = taskName,
        config = config
    ) %<-% input

    if (!exists("scdata")) {
        message("No single-cell data has been loaded, reloading from S3...")

        # assign it to the global environment so we can
        # persist it across runs of the wrapper
        assign("scdata", reload_from_s3(experiment_id), pos = ".GlobalEnv")

        message("Single-cell data loaded.")
    }

    # call function to run and update global variable
    c(
        data, ...rest_of_results
    ) %<-% run_step(task_name, scdata, config)

    assign("scdata", data, pos = ".GlobalEnv")

    # send result to API
    message_id <- send_output_to_api(input, rest_of_results)

    return(message_id)
}

wrapper(
    "{'experimentId': \"5928a56c7cbff9de78974ab50765ed20\", 'taskName': \"test_fn\", 'config': {'auto': true, 'enabled': true, 'name': \"a\", 'limit': 200}}"
)