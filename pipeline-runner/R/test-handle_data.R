sns <- paws::sns()

# Create a topic to which we can send notifications.
topic <- sns$create_topic("ExampleTopic")



pipeline_config <- list(
    sns_topic = 'ExampleTopic',
    aws_config = NULL
)

mock_sns <- function(){
    return(
        list(publish=function(){
            return(list(MessageId = 'ok'))
        }
    ))
}

test_that("send_gem2s_update_to_api completes successfully", {
    stub(send_gem2s_update_to_api, 'paws::sns', mock_sns)

    response <- send_gem2s_update_to_api(pipeline_config,
                            experiment_id = 'dfgdfg',
                            task_name = 'dsfdsdf',
                            data = 1:5,
                            input = list(auth_JWT='ayylmao'))
    expect_true(response$MessageId=='ok')
})
