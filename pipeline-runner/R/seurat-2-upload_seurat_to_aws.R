upload_seurat_to_aws <- function(input, pipeline_config, prev_out) {
  message("Uploading to AWS ...")

  experiment_id <- input$experimentId
  project_id <- input$projectId

  # destructure what need from prev_out
  scdata <- prev_out$scdata
  config <- prev_out$config

  scdata <- format_seurat(scdata, experiment_id)

  # TODO: change sample ids/names in so that get appropriate cell sets
  cell_sets <- get_cell_sets(scdata, input)

  # TODO: add louvain clusters

  # cell sets file to s3
  cell_sets_data <- RJSONIO::toJSON(cell_sets)

  put_object_in_s3(pipeline_config,
                   bucket = pipeline_config$cell_sets_bucket,
                   object = charToRaw(cell_sets_data),
                   key = experiment_id
  )


  # embedding settings
  qc_config <- construct_saved_config(scdata)



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



# add 'cells_id'
# 'samples' must be already added
# current input$metadata not yet implemented
format_seurat <- function(scdata, experiment_id) {
  scdata$cells_id <- seq_len(ncol(scdata))-1

  # add gene annotations
  rns <- row.names(scdata)
  scdata@misc$gene_annotations <- data.frame(input = rns, name = rns, original_name = rns)

  # add dispersion
  scdata <- add_dispersions(scdata)

  # other
  scdata@misc$experimentId <- experiment_id
  scdata@misc$color_pool <- get_color_pool()
  scdata@misc$ingestionData <- Sys.time()

  stopifnot('pca' %in% Seurat::Reductions(scdata))
  scdata@misc$active.reduction <- 'pca'

  return(scdata)
}

construct_saved_config <- function(scdata) {

  # embedding
  config <- list(
    configureEmbedding <- list(
      embeddingSettings = list(
        method = "umap",
        useSaved = TRUE,
        methodSettings = list(
          umap = list(
            minimumDistance = 0.3,
            distanceMetric = "cosine"
          ),
          tsne = list(
            perplexity = min(30, ncol(scdata) / 100),
            learningRate = max(200, ncol(scdata) / 12)
          )
        )
      ),
      clusteringSettings = list(
        method = "louvain",
        useSaved = TRUE,
        methodSettings = list(louvain = list(resolution = 0.8))
      )
    )
  )

  return(config)
}
