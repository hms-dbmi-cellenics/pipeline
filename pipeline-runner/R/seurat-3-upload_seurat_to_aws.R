upload_seurat_to_aws <- function(input, pipeline_config, prev_out) {
  message("Uploading to AWS ...")

  experiment_id <- input$experimentId

  # destructure what need from prev_out
  scdata <- prev_out$scdata
  config <- prev_out$config

  scdata <- format_seurat(scdata, experiment_id)

  # change sample ids/names so that get sample cell sets
  input <- add_samples_to_input(scdata, input)
  input <- add_metadata_to_input(scdata, input)
  scdata <- change_sample_names_to_ids(scdata, input)
  cell_sets <- get_cell_sets(scdata, input)

  # add louvain clusters
  cluster_sets <- data.frame(
    cluster = scdata$seurat_clusters,
    cell_ids = scdata$cells_id
  ) %>%
    format_cell_sets_object('louvain', scdata@misc$color_pool)

  # cell sets file to s3
  cell_sets$cellSets <- c(list(cluster_sets), cell_sets$cellSets)
  cell_sets_data <- RJSONIO::toJSON(cell_sets)

  put_object_in_s3(pipeline_config,
                   bucket = pipeline_config$cell_sets_bucket,
                   object = charToRaw(cell_sets_data),
                   key = experiment_id)

  # replicate qc config for simplicity
  # could also create a 'seurat_config' column in experiment table and change the ui/api around more
  qc_config <- construct_qc_config(list(one = scdata), FALSE, FALSE)
  qc_config$configureEmbedding$embeddingSettings$useSaved <- TRUE
  qc_config$configureEmbedding$embeddingSettings$method <- SeuratObject::DefaultDimReduc(scdata)

  # seurat object to s3
  object_key <- upload_matrix_to_s3(pipeline_config, experiment_id, scdata)
  message('Count matrix uploaded to ', pipeline_config$processed_bucket, ' with key ',object_key)

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

add_samples_to_input <- function(scdata, input) {
  samples <- unique(scdata$samples)
  input$sampleNames <- samples
  input$sampleIds <- ids::uuid(n = length(samples))
  return(input)
}

change_sample_names_to_ids <- function(scdata, input) {
  sample_ids <- input$sampleIds
  names(sample_ids) <- input$sampleNames
  scdata$samples <- sample_ids[scdata$samples]
  return(scdata)
}

add_metadata_to_input <- function(scdata, input) {
  group_cols <- find_group_columns(scdata)

  metadata <- list()
  meta_vals <- scdata@meta.data[!duplicated(scdata$samples), ]

  for (col in group_cols) {
    metadata[[col]] <- meta_vals[, col]
  }

  if (length(group_cols)) {
    input$metadata <- metadata
  }

  return(input)
}

# get column names that are consistent with sample groups
find_group_columns <- function(scdata) {
  meta <- scdata@meta.data

  ndistinct_sample <- meta |>
    dplyr::group_by(samples) |>
    dplyr::summarise_all(dplyr::n_distinct)

  ndistinct <- meta |>
    dplyr::summarise_all(dplyr::n_distinct)

  nsamples <- length(unique(scdata$samples))

  # group columns must:
  # - have fewer than the number of samples
  # - have at least two values
  # - have only one value per sample
  too.many <- ndistinct >= nsamples
  too.few <- ndistinct <= 1
  one.per.sample <- apply(ndistinct_sample, 2, function(x) all(x == 1))
  group.cols <- names(ndistinct)[!too.many & !too.few & one.per.sample]

  return(group.cols)
}



# add 'cells_id'
# 'samples' must be already added
# current input$metadata not yet implemented
format_seurat <- function(scdata, experiment_id) {

  scdata <- add_samples_col(scdata)
  scdata$cells_id <- seq_len(ncol(scdata))-1

  # other
  scdata@misc$experimentId <- experiment_id
  scdata@misc$color_pool <- get_color_pool()
  scdata@misc$ingestionDate <- Sys.time()

  # need to mock processing config
  metadata_cols <- list('percent.mt' = 0, 'doublet_scores' = 0, 'doublet_class' = 'singlet')
  scdata <- mock_metadata(scdata, metadata_cols)

  return(scdata)
}

# use 'samples' or 'sample' if present, otherwise assume one sample
add_samples_col <- function(scdata) {
  in.meta <- c('samples', 'sample') %in% colnames(scdata@meta.data)

  if (!any(in.meta)) {
    scdata$samples <- 'NA'
  } else {
    sample_col <- c('samples', 'sample')[which(in.meta)[1]]
    scdata$samples <- scdata@meta.data[[sample_col]]
  }

  return(scdata)
}

mock_metadata <- function(scdata, metadata_cols) {
  for (col in names(metadata_cols)) {
    if (is.null(scdata@meta.data[[col]]))
      scdata[[col]] <- metadata_cols[[col]]
  }

  return(scdata)
}
