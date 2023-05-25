copy_s3_objects <- function(input, pipeline_config, prev_out = NULL) {
  from_experiment_id <- input$fromExperimentId
  to_experiment_id <- input$toExperimentId

  original_data <- load_from_experiment_data(from_experiment_id, pipeline_config)

  original_scdata <- original_data$scdata
  original_cell_sets <- original_data$cell_sets
  sample_ids_map <- input$sampleIdsMap

  translated_scdata <- translate_sample_ids(original_scdata, sample_ids_map)
  translated_scdata@misc$experimentId <- to_experiment_id

  translated_cell_sets <- list()
  translated_cell_sets$cellSets <- format_cell_sets(
    translate_cell_sets(original_cell_sets$cellSets, sample_ids_map)
  )

  # cell sets file to s3
  put_object_in_s3(pipeline_config,
    bucket = pipeline_config$cell_sets_bucket,
    object = charToRaw(RJSONIO::toJSON(translated_cell_sets)),
    key = to_experiment_id
  )

  upload_matrix_to_s3(pipeline_config, to_experiment_id, translated_scdata)

  return()
}


#' load parent experiment data
#'
#' Loads the processed rds and cellsets file from the parent experiment from s3.
#'
#' @param input list of input parameters
#' @param pipelne_config list of pipeline parameters
#'
#' @return list with scdata and parsed cellsets
#' @export
#'
load_from_experiment_data <- function(experiment_id, pipeline_config) {
  # load parent processed scdata and cellsets
  s3 <- paws::s3(config = pipeline_config$aws_config)
  parent_scdata <- load_processed_scdata(s3, pipeline_config, experiment_id)
  # parent_cellsets <- parse_cellsets(load_cellsets(s3, pipeline_config, experiment_id))
  parent_cellsets <- load_cellsets(s3, pipeline_config, experiment_id)

  return(list(scdata = parent_scdata, cell_sets = parent_cellsets))
}

#' Add new sample ids to the subset Seurat Object
#'
#' @param scdata Seurat Object
#' @param sample_id_map data.table of parent/subset sample id map
#'
#' @return SeuratObject with new sample ids
#' @export
#'
translate_sample_ids <- function(scdata, sample_id_map) {
  sample_map_idx <- match(scdata$samples, names(sample_id_map))
  scdata$samples <- unname(unlist(sample_id_map[sample_map_idx]))

  return(scdata)
}

translate_cell_sets <- function(cell_sets, sample_ids_map) {
  # Loop through each element in the 'children' list
  for (i in seq_along(cell_sets$children)) {
    # Get the 'key' value
    key <- cell_sets$children[[i]]$key

    # Apply translate() to the values in the list that exist in the vector
    cell_sets$children[[i]]$key <- lapply(
      key,
      function(current_key) {
        ifelse(
          current_key %in% names(sample_ids_map),
          sample_ids_map[[current_key]],
          current_key
        )
      }
    )
  }

  return(cell_sets)
}

format_cell_sets <- function(cell_sets) {
  formatted <- list()

  for (i in seq_along(cell_sets$key)) {
    cell_class <- list(
      key = cell_sets$key[[i]],
      name = cell_sets$name[[i]],
      rootNode = cell_sets$rootNode[[i]],
      type = cell_sets$type[[i]],
      children = list()
    )

    for (j in seq_along(cell_sets$children[[i]]$key)) {
      children <- cell_sets$children[[i]]

      cell_set <- list(
        key = children$key[[j]],
        name = children$name[[j]],
        rootNode = children$rootNode[[j]],
        color = children$color[[j]],
        type = children$type[[j]],
        cellIds = ensure_is_list_in_json(children$cellIds[[j]])
      )

      cell_class$children <- append(cell_class$children, list(cell_set))
    }

    formatted <- append(formatted, list(cell_class))
  }

  return(formatted)
}