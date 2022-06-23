
upload_to_aws <- function(input, pipeline_config, prev_out) {
  message("Uploading to AWS ...")
  check_names <- c("config", "counts_list", "annot", "doublet_scores", "scdata_list", "qc_config")
  check_prev_out(prev_out, check_names)

  experiment_id <- input$experimentId
  project_id <- input$projectId

  # destructure what need from prev_out
  scdata_list <- prev_out$scdata_list
  config <- prev_out$config
  config_dataProcessing <- prev_out$qc_config

  message("Constructing cell sets ...")
  # scdata <- scdata_list[[sample]]
  cell_sets <- get_cell_sets(scdata_list, input)

  # cell sets file to s3
  cell_sets_data <- RJSONIO::toJSON(cell_sets)

  put_object_in_s3(pipeline_config,
    bucket = pipeline_config$cell_sets_bucket,
    object = charToRaw(cell_sets_data),
    key = experiment_id
  )

  message('ExperimentID: ', experiment_id)
 for (sample in names(scdata_list)) {
    message('Uploading: ', sample)
    message("scdata_list[[sample]]: ")
    message(scdata_list[[sample]])
    # seurat object to s3
    message("Uploading Seurat Object to S3 ...")
    fpath <- file.path(tempdir(), "experiment.rds")
    saveRDS(scdata_list[[sample]], fpath, compress = FALSE)
    message("Created fpath")
    # can only upload up to 50Gb because part numbers can be any number from 1 to 10,000, inclusive.
    put_object_in_s3_multipart(pipeline_config,
      bucket = pipeline_config$source_bucket,
      object = fpath,
      key = file.path(experiment_id, sample, "r.rds")
    )
  }

  quit(status=1)

  cluster_env <- pipeline_config$cluster_env

  experiment_data <- list(
    apiVersion = "2.0.0-data-ingest-seurat-rds-automated",
    experimentId = experiment_id,
    experimentName = config$name,
    meta = list(
      organism = config$organism,
      type = config$input$type
    ),
    processingConfig = config_dataProcessing
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

  # cells_sample <- scdata$samples
  # cells_sample <- names(scdata)
  # saveRDS(input, '/debug/input.alpha.gem7.rds')
  sample_ids <- unlist(input$sampleIds)
  sample_names <- unlist(input$sampleNames)

  # now the scdata contains only one sample so we don't really have to loop anymore
  # I've added this to keep the same cellsets structure but the nestig of $children
  # can probably be removed in the future
  # sample_ids <- scdata@meta.data$samples[[1]]
  # sample_ids <- names(scdata_list)
  # assuming that accessing by index & key (ie sample_id) yields the same order might be dangerous
  # gotta check if it's true in R
  for (i in seq_along(sample_ids)) {
    scdata <- scdata_list[[i]]
    sample_id <- sample_ids[i]
    sample_name <- sample_names[i]
    # saveRDS(scdata, '/debug/scdata_cells_id.rds')
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


# cell_sets fn for seurat metadata information
build_metadata_cellsets <- function(input, scdata_list, color_pool) {
  cell_set_list <- c()
  meta <- lapply(input$metadata, unlist)

  # user-supplied metadata track names
  keys <- names(meta)

  # syntactically valid metadata names as stored in scdata
  # same names as used in construct_metadata including internal 'samples' column (dropped)
  seurat_keys <- make.names(c("samples", keys), unique = TRUE)[-1]

  color_index <- 1
  for (i in seq_along(keys)) {
    key <- keys[i]
    seurat_key <- seurat_keys[i]

    cell_set <- list(
      "key" = key,
      "name" = key,
      "rootNode" = TRUE,
      "children" = c(),
      "type" = "metadataCategorical"
    )

    # values of current metadata track
    values <- unique(meta[[i]])

    for (j in seq_along(values)) {
      value <- values[j]

      cell_ids <- list()
      for (scdata in scdata_list) {
        cell_ids <- append(cell_ids,  scdata$cells_id[scdata[[seurat_key]] == value])
      }

      cell_set$children[[j]] <- list(
        "key" = paste(key, value, sep = "-"),
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
