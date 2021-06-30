upload_to_aws <- function(input, pipeline_config) {
  # set to local dirs for interactive dev
  input_dir <- "/input"
  output_dir <- "/output"

  experiment_id <- input$experimentId
  project_id <- input$projectId

  # save experiment_id for record-keeping
  writeLines(experiment_id, file.path(output_dir, "experiment_id.txt"))

  # read experiment config
  config <- RJSONIO::fromJSON(file.path(input_dir, "/meta.json"))

  # read config related to QC pipeline
  config_dataProcessing <- RJSONIO::fromJSON(file.path(output_dir, "config_dataProcessing.json"))

  # Design cell_set scratchpad for DynamoDB
  scratchpad <- list(
    key = "scratchpad",
    name = "Scratchpad",
    rootNode = TRUE,
    children = list(),
    type = "cellSets"
  )

  color_pool <- get_color_pool()
  samples_set <- samples_sets(input, output_dir, color_pool)

  # remove used colors from pool
  used <- seq_along(samples_set$children)
  color_pool <- color_pool[-used]

  # Design cell_set meta_data for DynamoDB
  cell_sets <- c(list(scratchpad), list(samples_set))

  if ("metadata" %in% names(config)) {
    cell_sets <- c(cell_sets, meta_sets(output_dir, color_pool))
  }

  cell_sets <- list(cellSets = cell_sets)

  print(paste("Experiment name is", config$name))

  experiment_data <- list(
    apiVersion = "2.0.0-data-ingest-seurat-rds-automated",
    experimentId = experiment_id,
    experimentName = config$name,
    meta = list(
      organism = config$organism,
      type = config$input[["type"]]
    ),
    processingConfig = config_dataProcessing
  )

  # cell sets file to s3
  cell_sets_data <- RJSONIO::toJSON(cell_sets)

  put_object_in_s3(pipeline_config,
    bucket = pipeline_config$cell_sets_bucket,
    object = charToRaw(cell_sets_data),
    key = experiment_id
  )

  # seurat object to s3
  put_object_in_s3_multipart(pipeline_config,
    bucket = pipeline_config$source_bucket,
    object = file.path(output_dir, "experiment.rds"),
    key = file.path(experiment_id, "r.rds")
  )

  cluster_env <- pipeline_config$cluster_env
  print(sprintf("Experiment ID: %s uploaded to %s.", experiment_id, cluster_env))

  data <- list(
    item = experiment_data,
    table = pipeline_config$experiments_table
  )

  if (cluster_env == "production") {
    print(sprintf("https://scp.biomage.net/experiments/%s/data-exploration", experiment_id))
  }

  return(data)
}

samples_sets <- function(config, output_dir, color_pool) {
  sample_annotations <- read.csv(file.path(output_dir, "samples-cells.csv"),
    sep = "\t",
    col.names = c("Cells_ID", "Value"),
    na.strings = "None"
  )

  cell_set <- list(
    key = "sample",
    name = "Samples",
    rootNode = TRUE,
    children = list(),
    type = "metadataCategorical"
  )

  samples <- unique(sample_annotations$Value)

  for (i in seq_along(samples)) {
    sample <- toString(samples[i])
    view <- sample_annotations[sample_annotations$Value == sample, "Cells_ID"]

    cell_set$children[[i]] <- list(
      key = paste0(sample),
      name = toStrng(config$sampleNames[[match(sample, config$sampleIds)]]),
      color = color_pool[i],
      cellIds = view
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
