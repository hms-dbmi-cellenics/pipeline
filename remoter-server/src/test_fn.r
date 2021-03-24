# This is a dummy function.
# Its only purpose is to demonstrate the input and output format
# that you will expect to your pipeline steps.

source('utils.r')

# some dummy function to check that imports do not import additional
# functions from a file into the namespace
some_other_stuff <- function(a) {
    a + 5
}

# a sample task
task <- function(input_data, input_config, task_name, sample_id) {

    # example where after coming up with sensible defaults the configuration
    # will be changed to a different number, say, 202
    config <- input_config
    config$filterSettings[["minCellSize"]] <- 420
    plots <- list()
    plots[generate_plotuuid(sample_id, task_name, 0)] = c(1, 2, 3)
    plots[generate_plotuuid(sample_id, task_name, 1)] = c(4, 5, 6)

    # the result object will have to conform to this format.
    result <- list(
        data = input_data,
        config = config,
        plotData = plots
    )

    return(result)
}