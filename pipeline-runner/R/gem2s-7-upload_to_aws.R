
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

  # destructure what need from prev_out
  scdata_list <- prev_out$scdata_list
  config <- prev_out$config
  qc_config <- prev_out$qc_config
  disable_qc_filters <- prev_out$disable_qc_filters

  if("sample_id_map" %in% names(prev_out)) {
    input$sampleIds <- names(scdata_list)
    input$sampleNames <- names(scdata_list)
  }

  message("Constructing cell sets ...")
  cell_sets <- get_cell_sets(scdata_list, input)

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
    processingConfig = qc_config
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

# creates initial cell sets object
get_cell_sets <- function(scdata_list, input) {
  scratchpad <- list(
    key = "scratchpad",
    name = "Custom cell sets",
    rootNode = TRUE,
    children = list(),
    type = "cellSets"
  )

  color_pool <- get_color_pool()
  sample_cellsets <- build_sample_cellsets(input, scdata_list, color_pool)

  # remove used colors from pool
  used <- seq_along(sample_cellsets$children)
  color_pool <- color_pool[-used]

  # Design cell_set meta_data for DynamoDB
  cell_sets <- c(list(scratchpad), list(sample_cellsets))

  if ("metadata" %in% names(input)) {
    cell_sets <- c(cell_sets, build_metadata_cellsets(input, scdata_list, color_pool))
  }

  cell_sets <- list(cellSets = cell_sets)
}

# cell_sets fn for seurat samples information
build_sample_cellsets <- function(input, scdata_list, color_pool) {
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
    scdata <- scdata_list[[sample_id]]
    sample_name <- sample_names[i]

    cell_ids <- scdata$cells_id

    cell_set$children[[i]] <- list(
      key = sample_id,
      name = sample_name,
      color = color_pool[i],
      cellIds = unname(cell_ids)
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
#' @param scdata_list list of Seurat objects
#' @param color_pool list of colors to use
#'
#' @return list of cellsets
#' @export
#'
build_metadata_cellsets <- function(input, scdata_list, color_pool) {
  cell_set_list <- c()
  user_metadata <- lapply(input$metadata, unlist)

  # user-supplied metadata track names
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

      cell_ids <- list()
      for (scdata in scdata_list) {
        cells_in_value <- scdata[[valid_metadata_name]] == value
        cell_ids <- append(cell_ids,  scdata$cells_id[cells_in_value])
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
