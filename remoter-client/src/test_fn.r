some_other_stuff <- function(a) {
    a + 5
}

task <- function(input_data, input_config) {
    config <- input_config

    config$limit <- 202

    # do some dummy transformation
    result <- list(
        data = input_data,
        config = config,
        plotData = list(
            plot1 = c(some_other_stuff(6), 2, 3),
            plot2 = c(4, 5, 6)
        )
    )

    return(result)
}