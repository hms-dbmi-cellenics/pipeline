
upload_to_aws <- function(input, pipeline_config, prev_out) {
  message('Uploading to AWS ...')

  experiment_id <- input$experimentId
  project_id <- input$projectId

  # destructure what need from prev_out
  scdata <- prev_out$scdata
  config <- prev_out$config
  config_dataProcessing <- prev_out$qc_config


  message("Constructing cell sets ...")
  cell_sets <- get_cell_sets(scdata, input)

  # cell sets file to s3
  cell_sets_data <- RJSONIO::toJSON(cell_sets)

  put_object_in_s3(pipeline_config,
                   bucket = pipeline_config$cell_sets_bucket,
                   object = charToRaw(cell_sets_data),
                   key = experiment_id
  )

  # seurat object to s3
  message("Uploading Seurat Object to S3 ...")
  fpath <- file.path(tempdir(), 'experiment.rds')
  saveRDS(scdata, fpath, compress = FALSE)

  put_object_in_s3_multipart(pipeline_config,
                             bucket = pipeline_config$source_bucket,
                             object = fpath,
                             key = file.path(experiment_id, "r.rds")
  )

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
      table = pipeline_config$experiments_table),
    output = list())

  return(res)
}

# creates initial cell sets object
get_cell_sets <- function(scdata, input) {

  scratchpad <- list(
    key = "scratchpad",
    name = "Custom cell sets",
    rootNode = TRUE,
    children = list(),
    type = "cellSets"
  )

  color_pool <- get_color_pool()
  samples_set <- samples_sets(input, scdata, color_pool)

  # remove used colors from pool
  used <- seq_along(samples_set$children)
  color_pool <- color_pool[-used]

  # Design cell_set meta_data for DynamoDB
  cell_sets <- c(list(scratchpad), list(samples_set))

  # TODO: refactor meta_sets
  if ("metadata" %in% names(input)) {
    cell_sets <- c(cell_sets, meta_sets(output_dir, color_pool))
  }

  cell_sets <- list(cellSets = cell_sets)
}

# cell_sets fn for seurat samples information
samples_sets <- function(input, scdata, color_pool) {


  cell_set <- list(
    key = "sample",
    name = "Samples",
    rootNode = TRUE,
    children = list(),
    type = "metadataCategorical"
  )

  cell_sample <- scdata$samples
  samples <- unique(cell_sample)

  for (i in seq_along(samples)) {
    sample <- samples[i]
    cell_ids <- scdata$cells_id[cell_sample == sample]

    cell_set$children[[i]] <- list(
      key = paste0(sample),
      name = toString(input$sampleNames[[match(sample, input$sampleIds)]]),
      color = color_pool[i],
      cellIds = unname(cell_ids)
    )
  }

  return(cell_set)
}


# cell_sets fn for seurat metadata information
meta_sets <- function(output_dir, color_pool) {
  meta_annotations <- read.csv(file.path(output_dir, "metadata-cells.csv"), sep = "\t")

  cell_set_list <- c()

  # The first column is the cells_id, the rest is the metadata information
  for (i in seq(2, ncol(meta_annotations))) {
    key <- name <- toString(colnames(meta_annotations)[i])

    cell_set <- list(
      "key" = key,
      "name" = name,
      "rootNode" = TRUE,
      "children" = c(),
      "type" = "metadataCategorical"
    )

    annot <- meta_annotations[[i]]
    values <- unique(annot)

    for (i in seq_along(values)) {
      value <- toString(values[i])
      view <- meta_annotations[which(annot == value), "cells_id"]

      cell_set$children[[i]] <- list(
        "key" = paste(key, value, sep = "-"),
        "name" = value,
        "color" = color_pool[i],
        "cellIds" = view
      )
    }
    cell_set_list <- c(cell_set_list, list(cell_set))
  }
  return(cell_set_list)
}
