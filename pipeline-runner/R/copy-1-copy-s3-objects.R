get_sample_ids_replacer <- function(sample_ids_map) {

  replacer <- function(key) {
    new_key <- key
    for (sample_id in names(sample_ids_map)) {
      new_key <- gsub(sample_id, sample_ids_map[[sample_id]], new_key)
    }

    return(new_key)
  }

  return(replacer)
}

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
  translated_cell_sets$cellSets <- unflatten_cell_sets(
    translate_cell_sets(original_cell_sets$cellSets, sample_ids_map)
  )

  # cell sets file to s3
  put_object_in_s3(
    pipeline_config,
    bucket = pipeline_config$cell_sets_bucket,
    object = charToRaw(RJSONIO::toJSON(translated_cell_sets)),
    key = to_experiment_id
  )

  upload_matrix_to_s3(pipeline_config, to_experiment_id, translated_scdata)

  replace_sample_ids <- get_sample_ids_replacer(sample_ids_map)

  filtered_cells_migrate(from_experiment_id, to_experiment_id, sample_ids_map, pipeline_config)

  s3_copy_by_prefix("biomage-source-development-000000000000", from_experiment_id, "biomage-source-development-000000000000", to_experiment_id, pipeline_config$aws_config, key_transform = replace_sample_ids)

  message("\n")
  message("Copy s3 objects to experiment clone step complete.")
  return(list(
    data = list(),
    output = prev_out
  ))
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

  parent_cellsets <- load_cellsets(s3, pipeline_config, experiment_id)

  return(list(scdata = parent_scdata, cell_sets = parent_cellsets))
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

# TODO: Change this name for a good one
filtered_cells_migrate <- function(from_experiment_id, to_experiment_id, sample_ids_map, pipeline_config) {
  s3 <- paws::s3(config = pipeline_config$aws_config)

  steps <- c(
    "cellSizeDistribution",
    "dataIntegration",
    "doubletScores",
    "mitochondrialContent",
    "numGenesVsNumUmis"
  )

  for (step in steps) {
    c(body, ...rest) %<-% s3$get_object(
      Bucket = "biomage-filtered-cells-development-000000000000",
      Key = paste0(from_experiment_id, "/", step, "/", names(sample_ids_map)[[1]], ".rds")
    )

    conn <- gzcon(rawConnection(body))
    object <- readRDS(conn)

    translated <- translate_filtered_cells(object, sample_ids_map)

    rds_translated <- tempfile()
    saveRDS(translated, file = rds_translated)

    for (sample_id in unlist(sample_ids_map)) {
      s3$put_object(
        Bucket = "biomage-filtered-cells-development-000000000000",
        Key = paste0(to_experiment_id, "/", step, "/", sample_id, ".rds"),
        Body = rds_translated
      )
    }
  }

}

translate_filtered_cells <- function(filtered_cells, sample_ids_map) {
  sample_map_idx <- match(names(filtered_cells), names(sample_ids_map))
  names(filtered_cells) <- unname(unlist(sample_ids_map[sample_map_idx]))

  return(filtered_cells)
}

#' Replaces sample ids with their matches values in sample_ids_map
#'
#' @param scdata Seurat Object
#' @param sample_ids_map data.table of parent/subset sample id map
#'
#' @return SeuratObject with new sample ids
#' @export
#'
translate_sample_ids <- function(scdata, sample_ids_map) {
  sample_map_idx <- match(scdata$samples, names(sample_ids_map))
  scdata$samples <- unname(unlist(sample_ids_map[sample_map_idx]))

  return(scdata)
}