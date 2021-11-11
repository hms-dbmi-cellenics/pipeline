
mock_sns <- function(config) {
    return(
        list(publish = function (Message, TopicArn, MessageAttributes) {
            return (list(MessageId = 'ok'))
        }
    ))
}

test_that("send_gem2s_update_to_api completes successfully", {
    pipeline_config <- list(
        sns_topic = 'ExampleTopic',
        aws_config = NULL
    )

    mockery::stub(send_gem2s_update_to_api, 'paws::sns', mock_sns)

    response <- send_gem2s_update_to_api(pipeline_config,
                            experiment_id = 'dfgdfg',
                            task_name = 'dsfdsdf',
                            data = 1:5,
                            input = list(auth_JWT='ayylmao'))


    expect_true(response == 'ok')
})
