# This is a dummy function.
# Its only purpose is to demonstrate the input and output format
# that you will expect to your pipeline steps.

# some dummy function to check that imports do not import additional
# functions from a file into the namespace
some_other_stuff <- function(a) {
    a + 5
}

# a sample task
task <- function(input_data, input_config, task_name, sample_uuid) {

    # example where after coming up with sensible defaults the configuration
    # will be changed to a different number, say, 202
    config <- input_config
    config$filterSettings[["minCellSize"]] <- 420

    plotData <- list()
    plotData[[paste(sample_uuid, task_name, 0, collapse="-")]] <- c(1, 2, 3)
    plotData[[paste(sample_uuid, task_name, 1, collapse="-")]] <- c(4, 5, 6)

    # the result object will have to conform to this format.
    result <- list(
        data = input_data,
        config = config,
        plotData = plotData
    )

    return(result)
}