copy_s3_objects <- function(input, pipeline_config, prev_out = NULL) {
  from_experiment_id <- input$fromExperimentId
  to_experiment_id <- input$toExperimentId
  sample_ids_map <- input$sampleIdsMap

  message("Copying processed rds.")
  copy_processed_rds(from_experiment_id, to_experiment_id, sample_ids_map, pipeline_config)
  message("Copying cell sets.")
  copy_cell_sets(from_experiment_id, to_experiment_id, sample_ids_map, pipeline_config)
  message("Copying filtered cells.")
  copy_filtered_cells(from_experiment_id, to_experiment_id, sample_ids_map, pipeline_config)
  message("Copying source.")
  copy_source(from_experiment_id, to_experiment_id, sample_ids_map, pipeline_config)

  message("\n")
  message("Copy s3 objects step complete.")
  return(list(data = list(), output = prev_out))
}

copy_processed_rds <- function(
  from_experiment_id, to_experiment_id, sample_ids_map, pipeline_config
) {
  s3 <- paws::s3(config = pipeline_config$aws_config)
  scdata <- load_processed_scdata(s3, pipeline_config, from_experiment_id)

  scdata$samples <- translate_sample_ids(scdata$samples, sample_ids_map)
  scdata@misc$experimentId <- to_experiment_id

  upload_matrix_to_s3(pipeline_config, to_experiment_id, scdata)
}

copy_cell_sets <- function(from_experiment_id, to_experiment_id, sample_ids_map, pipeline_config) {
  s3 <- paws::s3(config = pipeline_config$aws_config)
  original_cell_sets <- load_cellsets(s3, pipeline_config, from_experiment_id)

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
}

copy_source <- function(from_experiment_id, to_experiment_id, sample_ids_map, pipeline_config) {
  bucket <- pipeline_config$source_bucket

  for (original_sample_id in names(sample_ids_map)) {
    copy_sample_id <- sample_ids_map[[original_sample_id]]

    original_key <- paste0(from_experiment_id, "/", original_sample_id, "/", "r.rds")
    copy_key <- paste0(to_experiment_id, "/", copy_sample_id, "/", "r.rds")

    source_data <- get_s3_rds(bucket, original_key, pipeline_config$aws_config)

    # Translate relevant parts
    source_data@meta.data$samples <- translate_sample_ids(
      source_data@meta.data$samples,
      sample_ids_map
    )

    source_data@misc$experimentId <- to_experiment_id
    # Finish translating relevant parts

    message("Uploading ", original_key, "to ", copy_key, " with samples translated")
    put_s3_rds(bucket, copy_key, pipeline_config$aws_config, source_data)
  }
}

# TODO: Once the filtered_cells objects are actually separated by sample (not identical)
#  then we can use s3_copy_by_prefix for this function instead of all this code,
#  so basically replace this function by something similar to copy_source
copy_filtered_cells <- function(
  from_experiment_id, to_experiment_id, sample_ids_map, pipeline_config
) {
  steps <- c(
    "cellSizeDistribution",
    "dataIntegration",
    "doubletScores",
    "mitochondrialContent",
    "numGenesVsNumUmis"
  )

  bucket <- pipeline_config$cells_id_bucket
  aws_config <- pipeline_config$aws_config

  for (step in steps) {
    for (original_sample_id in names(sample_ids_map)) {
      copy_sample_id <- sample_ids_map[[original_sample_id]]

      original_key <- paste0(from_experiment_id, "/", step, "/", original_sample_id, ".rds")
      copy_key <- paste0(to_experiment_id, "/", step, "/", copy_sample_id, ".rds")

      filtered_cells <- get_s3_rds(bucket, original_key, aws_config)

      names(filtered_cells) <- translate_sample_ids(names(filtered_cells), sample_ids_map)

      message("Uploading ", original_key, "to ", copy_key, " with samples translated")
      put_s3_rds(bucket, copy_key, aws_config, filtered_cells)
    }
  }
}

#' Replaces sample cell sets keys with their matched values in sample_ids_map
#'
#' @param cell_sets original cell sets with old sample ids as keys
#' @param sample_ids_map data.table of parent/subset sample id map
#'
#' @return cell sets with new sample ids as keys instead of old ones
#' @export
#'
translate_cell_sets <- function(cell_sets, sample_ids_map) {
  for (i in seq_along(cell_sets$children)) {
    key <- cell_sets$children[[i]]$key

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

#' Replaces sample ids with their matched values in sample_ids_map
#'
#' @param old_sample_ids old sample ids vector
#' @param sample_ids_map data.table of parent/subset sample id map
#'
#' @return vector with new sample ids matching the order of old_sample_ids
#' @export
#'
translate_sample_ids <- function(old_sample_ids, sample_ids_map) {
  sample_map_idx <- match(old_sample_ids, names(sample_ids_map))
  new_sample_ids <- unname(unlist(sample_ids_map[sample_map_idx]))

  return(new_sample_ids)
}
