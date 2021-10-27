sns <- paws::sns()

# Create a topic to which we can send notifications.
topic <- sns$create_topic("ExampleTopic")



pipeline_config <- list(
    sns_topic = 'ExampleTopic',
    aws_config = NULL
)


send_gem2s_update_to_api(pipeline_config,
                         experiment_id = 'dfgdfg',
                         task_name = 'dsfdsdf',
                         data = 1:5,
                         input = 'blah')


test_that("filter_gene_umi_outlier updates gam p.level in config if auto", {
    scdata <- mock_scdata()
    config <- mock_config()
    config$filterSettings$regressionTypeSettings$gam$p.level <- 1
    out <- filter_gene_umi_outlier(scdata, config, '123def')
    new <- out$config$filterSettings$regressionTypeSettings$gam$p.level

    expect_lt(new, 1)
})
