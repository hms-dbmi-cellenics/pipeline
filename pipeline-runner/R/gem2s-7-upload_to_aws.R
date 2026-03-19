
#' Upload Seurat and cellsets objects to aws
#'
#' @param input The input object from the request
#' @param pipeline_config result of \code{load_config}
#' @param prev_out list with results appended in each gem2s task
#'
#' @return list with experiment parameters
#' @export
#'
upload_to_aws <- function(input, pipeline_config, prev_out) {
  message("Uploading to AWS ...")
  check_names <- c("config", "scdata_list", "qc_config", "disable_qc_filters")
  check_prev_out(prev_out, check_names)

  experiment_id <- input$experimentId
  project_id <- input$projectId

  # de-structure what need from prev_out
  scdata_list <- prev_out$scdata_list
  config <- prev_out$config
  qc_config <- prev_out$qc_config
  matrix_dir_list <- prev_out$matrix_dir_list
  default_qc_config <- prev_out$default_qc_config
  disable_qc_filters <- prev_out$disable_qc_filters

  # TODO: replace with subset_experiment flag when available
  if (disable_qc_filters == FALSE) {
    message("Constructing cell sets ...")
    cell_sets <- get_cell_sets(scdata_list, input)
  } else {
    message("Constructing cell sets for subset experiment ...")
    cell_sets <- get_subset_cell_sets(scdata_list, input, prev_out, disable_qc_filters)
  }

  # cell sets file to s3
  cell_sets_data <- RJSONIO::toJSON(cell_sets)

  put_object_in_s3(pipeline_config,
    bucket = pipeline_config$cell_sets_bucket,
    object = charToRaw(cell_sets_data),
    key = experiment_id
  )

  # remove previous existing data
  remove_bucket_folder(pipeline_config, pipeline_config$source_bucket, experiment_id)

  for (sample in names(scdata_list)) {
    message("Uploading sample ", sample, " object to S3 ...")
    fpath <- file.path(tempdir(), "experiment.rds")
    saveRDS(scdata_list[[sample]], fpath)

    # can only upload up to 50Gb because multipart uploads work by uploading
    # smaller chunks (parts) and the max amount of parts is 10,000.
    put_object_in_s3_multipart(pipeline_config,
      bucket = pipeline_config$source_bucket,
      object = fpath,
      key = file.path(experiment_id, sample, "r.rds")
    )

    # if this is bpcells also upload matrix_dir
    matrix_dir <- matrix_dir_list[[sample]]

    if (!is.null(matrix_dir)) {
      message("Uploading sample ", sample, " bpcells matrix_dir to S3 ...")
      tarfile <- tar_matrix_dir(sample, matrix_dir)
      put_object_in_s3_multipart(pipeline_config,
        bucket = pipeline_config$source_bucket,
        object = tarfile,
        key = file.path(experiment_id, sample, "matrix_dir.tar")
      )
    }
  }

  message("Samples uploaded")

  cluster_env <- pipeline_config$cluster_env

  experiment_data <- list(
    apiVersion = "2.0.0-data-ingest-seurat-rds-automated",
    experimentId = experiment_id,
    experimentName = config$name,
    meta = list(
      organism = config$organism,
      type = config$input$type
    ),
    processingConfig = qc_config,
    defaultProcessingConfig = default_qc_config
  )

  res <- list(
    data = list(
      item = experiment_data,
      table = pipeline_config$experiments_table
    ),
    output = list()
  )

  message("\nUpload to AWS step complete.")
  return(res)
}

tar_matrix_dir <- function(id, matrix_dir, suffix='') {

  # where we write the tar file
  tarfile <- file.path(tempdir(), paste0(id, '_matrix_dir', suffix, '.tar.zst'))

  # move to parent of matrix_dir and tar it
  matrix_dir_parent <- dirname(matrix_dir)

  current_dir <- getwd()
  setwd(matrix_dir_parent)
  tar_zstd(tarfile, files = basename(matrix_dir))
  setwd(current_dir)
  return(tarfile)
}

# TODO: newer R utils::tar supports zstd
tar_zstd <- function(tarfile, files) {
  system(paste("tar --zstd -cf", tarfile, files))
}

untar_zstd <- function(tarfile, exdir) {
  system(paste("tar --zstd -xf", tarfile, "-C", exdir))
}

#' Create initial cell sets object
#'
#' @param scdata list of Seurat objects
#' @param input The input object from the request
#'
#' @return cell set object
#' @export
#'
get_cell_sets <- function(scdata, input) {
  scratchpad <- list(
    key = "scratchpad",
    name = "Custom cell sets",
    rootNode = TRUE,
    children = list(),
    type = "cellSets"
  )

  color_pool <- get_color_pool()
  sample_cellsets <- build_sample_cellsets(input, scdata, color_pool)

  # remove used colors from pool
  color_pool <- remove_used_colors(sample_cellsets, color_pool)

  # Design cell_set metadata
  cell_sets <- c(list(scratchpad), list(sample_cellsets))

  if ("metadata" %in% names(input)) {
    cell_sets <- c(cell_sets, build_metadata_cellsets(input, scdata, color_pool))
  }

  cell_sets <- list(cellSets = cell_sets)
}

#' Create cell set using Seurat samples information
#'
#' @param input The input object from the request
#' @param scdata list of Seurat objects or a single Seurat object
#' @param color_pool list of colors to use
#'
#' @return cell set filled with samples information
#' @export
#'
build_sample_cellsets <- function(input, scdata, color_pool) {
  cell_set <- list(
    key = "sample",
    name = "Samples",
    rootNode = TRUE,
    children = list(),
    type = "metadataCategorical"
  )

  sample_ids <- unlist(input$sampleIds)
  sample_names <- unlist(input$sampleNames)

  for (i in seq_along(sample_ids)) {
    sample_id <- sample_ids[i]
    sample_name <- sample_names[i]

    if (is.list(scdata)) {
      cell_ids <- scdata[[sample_id]]$cells_id
    } else {
      cell_ids <- scdata$cells_id[scdata$samples == sample_id]
    }

    cell_set$children[[i]] <- list(
      key = sample_id,
      name = sample_name,
      color = color_pool[i],
      cellIds = ensure_is_list_in_json(unname(cell_ids))
    )
  }

  return(cell_set)
}


#' Create cellsets from user-supplied metadata
#'
#' This function creates the cellsets for the user-supplied metadata (in data
#' management module). It also takes care of assigning a color from the color pool
#' to each cellset.
#'
#' @param input list
#' @param scdata list of Seurat objects or a single Seurat object
#' @param color_pool list of colors to use
#'
#' @return list of cellsets
#' @export
#'
build_metadata_cellsets <- function(input, scdata, color_pool, disable_qc_filters = FALSE, subset_cellsets = NA) {
  cell_set_list <- c()
  user_metadata <- lapply(input$metadata, unlist)

  # user-supplied metadata track names
  if (disable_qc_filters == TRUE) {
    user_metadata <- extract_subset_user_metadata(subset_cellsets)
    } else {
    user_metadata <- lapply(input$metadata, unlist)
  }

  metadata_names <- names(user_metadata)

  # user supplied metadata names must be made syntactically valid, as seurat does.
  # We use them to access the cell ids stored in the seurat object, to create the
  # cellsets. The same names as used in construct_metadata including internal
  # 'samples' column which is dropped after making the names.
  valid_metadata_names <- make.names(c("samples", metadata_names), unique = TRUE)[-1]

  color_index <- 1
  for (i in seq_along(metadata_names)) {
    user_metadata_name <- metadata_names[i]
    valid_metadata_name <- valid_metadata_names[i]

    cell_set <- list(
      "key" = user_metadata_name,
      "name" = user_metadata_name,
      "rootNode" = TRUE,
      "children" = c(),
      "type" = "metadataCategorical"
    )

    # values of current metadata track
    values <- unique(user_metadata[[i]])

    for (j in seq_along(values)) {
      value <- values[j]

      if (is.list(scdata)) {
        cell_ids <- get_metadata_cell_ids(scdata, valid_metadata_name, value)
      } else {
        # need @meta.data for factor columns
        cells_in_value <- scdata@meta.data[[valid_metadata_name]] == value
        cell_ids <- scdata$cells_id[cells_in_value]
      }

      cell_set$children[[j]] <- list(
        "key" = paste(user_metadata_name, value, sep = "-"),
        "name" = value,
        "color" = color_pool[color_index],
        "cellIds" = unname(cell_ids)
      )

      color_index <- color_index + 1
    }
    cell_set_list <- c(cell_set_list, list(cell_set))
  }
  return(cell_set_list)
}

# get metadata cell ids from a list of Seurat object
get_metadata_cell_ids <- function(scdata_list, valid_metadata_name, value) {
  cell_ids <- list()
  for (scdata in scdata_list) {
    cells_in_value <- scdata[[valid_metadata_name]] == value
    cell_ids <- append(cell_ids,  scdata$cells_id[cells_in_value])
  }

  return(cell_ids)
}


#' create list of user supplied metadata tracks and values
#'
#' @param subset_cellsets data.table of cellsets
#'
#' @return list of metadata tracks and values
#' @export
#'
extract_subset_user_metadata <- function(subset_cellsets) {
  # metadata keys are the <track_name>-<value>, and name are the values alone
  metadata <- unique(subset_cellsets[type == "metadata"], by = "key")

  # Extract track name by removing the last component after hyphen (the value)
  # Keys have format "TrackName-Value"
  track_names <- sapply(metadata$key, function(k) {
    parts <- strsplit(k, "-")[[1]]
    if (length(parts) > 1) {
      paste(parts[-length(parts)], collapse = "-")
    } else {
      parts[1]
    }
  }, USE.NAMES = FALSE)

  metadata[, metadata_track := track_names]

  metadata <- metadata[, c("metadata_track", "name")]
  data.table::setnames(metadata, "name", "metadata_value")

  user_metadata <- list()
  for (track in unique(metadata$metadata_track)) {
    user_metadata[[track]] <- metadata[metadata_track == track, metadata_value]
  }

  return(user_metadata)
}


#' Remove used colors from pool
#'
#' @param cellsets cell set object
#' @param color_pool list of colors
#'
#' @return color pool with used colors removed
#' @export
#'
remove_used_colors <- function(cellsets, color_pool) {
  used <- seq_along(cellsets$children)
  color_pool <- color_pool[-used]
  return(color_pool)
}


#' Create cell set object for subset experiment
#'
#' @param scdata_list list of Seurat objects
#' @param input The input object from the request
#' @param prev_out list with results appended in each gem2s task
#' @param disable_qc_filters bool indicating if the data derives from the
#' subsetting of another experiment
#'
#' @return cell set object
#' @export
#'
get_subset_cell_sets <- function(scdata_list, input, prev_out, disable_qc_filters) {

  sample_id_map <- prev_out$sample_id_map
  parent_cellsets <- prev_out$parent_cellsets
  cell_ids_to_keep <- unlist(lapply(scdata_list, function(x) {
    x$cells_id
  }))

  subset_cellsets <- filter_parent_cellsets(parent_cellsets, cell_ids_to_keep)

  # replace old sample ids with new sample ids in the new cellsets
  for (i in 1:length(sample_id_map)) {
    subset_cellsets[key %in% names(sample_id_map)[i], key := unname(sample_id_map[i])]
  }

  input$sampleIds <- names(scdata_list)
  input$sampleNames <- subset_cellsets[input$sampleIds, name, on = "key", mult = "first"]


  # convert back cellsets to list format
  color_pool <- get_color_pool()
  sample_cellsets <- build_sample_cellsets(input, scdata_list, color_pool)
  color_pool <- remove_used_colors(sample_cellsets, color_pool)

  cell_sets <- c(list(sample_cellsets))

  if ("metadata" %in% unique(subset_cellsets$type)) {
    message("adding metadata cellsets to subset experiment")
    metadata_cellsets <- build_metadata_cellsets(input, scdata_list, color_pool, disable_qc_filters, subset_cellsets)
    color_pool <- remove_used_colors(metadata_cellsets[[1]], color_pool)
    cell_sets <- c(cell_sets, metadata_cellsets)
  }

  if ("scratchpad" %in% unique(subset_cellsets$type)) {
    message("adding custom cellsets to subset experiment")
    scratchpad_cellsets <- build_scratchpad_cellsets(color_pool, subset_cellsets)
  } else {
    scratchpad_cellsets <- list(
      key = "scratchpad",
      name = "Custom cell sets",
      rootNode = TRUE,
      children = list(),
      type = "cellSets"
    )
  }
  cell_sets <- c(cell_sets, list(scratchpad_cellsets))

  cell_sets <- list(cellSets = cell_sets)

  return(cell_sets)
}


#' Filter parent cellsets
#'
#' This function filters the parent cellsets removing clusters, which will be
#' recalculated in the subset experiment. It also filters cell ids not present
#' after the subset.
#'
#' @param parent_cellsets data.table with cellsets from the parent experiment
#' @param cell_ids_to_keep integer vector of cell ids to keep
#'
#' @return filtered cellsets data.table
#' @export
#'
filter_parent_cellsets <- function(parent_cellsets, cell_ids_to_keep) {
  # filter out all clustering cellsets
  subset_cellsets <- parent_cellsets[type != "cluster"]

  # filter out cells from cell_sets_original scratchpad and metadata
  subset_cellsets <- subset_cellsets[which(subset_cellsets$cell_id %in% cell_ids_to_keep)]

  return(subset_cellsets)
}


#' Create cellsets using scratchpad information from parent cellsets
#'
#' This function creates cellsets with scratchpad information including only
#' the cell ids in the subset experiment. If all cell ids from a custom cellset
#' of the parent experiments are filtered out in the subset experiment,
#' then that custom cellset is not included in the subset cell set.
#'
#' @param color_pool list of colors to use
#' @param subset_cellsets cell set resulting from parent cell set filtering
#'
#' @return cell set filled with scratchpad information from the parent cell set
#' @export
#'
build_scratchpad_cellsets <- function(color_pool, subset_cellsets) {
  scratchpad <- list(
    key = "scratchpad",
    name = "Custom cell sets",
    rootNode = TRUE,
    children = list(),
    type = "cellSets"
  )

  scratchpad_cellsets <- unique(subset_cellsets[type == "scratchpad", .(key, name)], by = "key")
  scratchpad_ids <- scratchpad_cellsets[, key]
  scratchpad_names <- scratchpad_cellsets[, name]

  for (i in seq_along(scratchpad_ids)) {

    scratchpad$children[[i]] <- list(
      key = scratchpad_ids[i],
      name = scratchpad_names[i],
      color = color_pool[i],
      cellIds =  ensure_is_list_in_json(subset_cellsets[key == scratchpad_ids[i], cell_id])
    )
  }

  return(scratchpad)
}


#' Synchronize metadata values from cellsets
#'
#' Updates metadata column values in the Seurat object to match the current
#' cellsets. This handles the case where metadata labels were renamed in the UI
#' and the Seurat object's metadata needs to be updated.
#'
#' @param scdata seurat object
#' @param cellsets data.table of cellsets
#'
#' @return scdata with synchronized metadata values
#'
sync_metadata_from_cellsets <- function(scdata, cellsets) {
  # Extract metadata cellsets only
  metadata <- cellsets[type == "metadata"]

  if (nrow(metadata) == 0) {
    # No metadata cellsets to sync
    return(scdata)
  }

  # Get the cell IDs from scdata
  scdata_cell_ids <- scdata$cells_id

  # Build a mapping of cell_id -> values for each metadata track
  track_mappings <- list()

  for (i in seq_len(nrow(metadata))) {
    row <- metadata[i]
    key <- row$key
    cell_id <- row$cell_id
    value <- row$name

    # Extract track name from key format "TrackName-Value"
    # Split by hyphen and use all but the last part as the track name
    parts <- strsplit(key, "-")[[1]]
    if (length(parts) > 1) {
      track_name <- paste(parts[-length(parts)], collapse = "-")
    } else {
      track_name <- parts[1]
    }

    # Make column name syntactically valid
    valid_track_name <- make.names(track_name)

    # Initialize track mapping if needed
    if (!valid_track_name %in% names(track_mappings)) {
      track_mappings[[valid_track_name]] <- list()
    }

    # Store the new value for this cell_id
    track_mappings[[valid_track_name]][[as.character(cell_id)]] <- value
  }

  # Create or update each metadata track in the Seurat object
  for (track_name in names(track_mappings)) {
    # Initialize with NA for all cells
    col_values <- rep(NA_character_, length(scdata_cell_ids))

    # Populate with values from cellsets
    for (i in seq_len(length(scdata_cell_ids))) {
      cell_id <- scdata_cell_ids[i]
      cell_id_str <- as.character(cell_id)

      # Look up the value for this cell in cellsets
      if (cell_id_str %in% names(track_mappings[[track_name]])) {
        col_values[i] <- track_mappings[[track_name]][[cell_id_str]]
      }
    }

    # Add or update the column in the Seurat object
    scdata@meta.data[[track_name]] <- col_values
  }

  return(scdata)
}
