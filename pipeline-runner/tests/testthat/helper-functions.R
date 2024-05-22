#' Fix time stamps in seurat object
#'
#' Seurat adds logs to certain command runs, with time stamps that make snapshot
#' tests fail. This function replaces them with a fixed datetime object.
#'
#' @param scdata seurat object
#'
#' @return seurat object with fixed time stamps
#'
clean_timestamp <- function(scdata) {
  fixed_datetime <- as.POSIXct("1991-12-19 05:23:00", tz = "UTC")

  for (slot in names(scdata@commands)) {
    scdata@commands[[slot]]@time.stamp <- fixed_datetime
  }

  if ("ingestionDate" %in% names(scdata@misc)) {
    scdata@misc$ingestionDate <- fixed_datetime
  }

  return(scdata)
}


# Seurat object commands are sometimes functions which will have a new environment
# for every call which breaks snapshots
remove_commands_functions <- function(data) {
  command_names <-  names(data@commands)
  for (command_name in command_names) {
    param_names <- names(data@commands[[command_name]]@params)

    for (param_name in param_names) {
      param <- data@commands[[command_name]]@params[[param_name]]

      if (methods::is(param, 'function'))
        data@commands[[command_name]]@params[[param_name]] <- NULL
    }
  }
  return(data)
}
