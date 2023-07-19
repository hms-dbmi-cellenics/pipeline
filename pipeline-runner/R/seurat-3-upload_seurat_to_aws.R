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

  # detect cluster metadata
  cluster_columns <- find_cluster_columns(scdata)

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
  qc_config <- construct_qc_config(list(one = scdata), unfiltered_samples = 'one')
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

find_cluster_columns <- function(scdata) {

  group_cols <- find_group_columns(scdata)
  exclude_cols <- c(group_cols, 'samples')
  metadata <- scdata@meta.data
  check_cols <- setdiff(colnames(metadata), exclude_cols)

  cluster_cols <- c()
  for (check_col in check_cols) {
    check_vals <- metadata[[check_col]]

    # skip if col is same as samples or group column
    is_sample_col <- FALSE
    for (exclude_col in exclude_cols)  {
      exclude_vals <- metadata[[exclude_col]]
      if (test_groups_equal(check_vals, exclude_vals)) {
        is_sample_col <- TRUE
        break
      }
    }
    if (is_sample_col) next()

    # skip if too few or way too many values
    n.vals <- length(table(check_vals))
    if (n.vals < 2) next()
    if (n.vals > 1000) next()

    # skip if values are numeric but non-integer
    if (is.numeric(check_vals) &&
        !all(as.integer(check_vals) == check_vals)) next()

    # skip if most values are repeated fewer than 3 times
    ordered_nreps <- names(sort(table(table(check_vals)), decreasing = TRUE))
    if (any(1:3 %in% head(ordered_nreps, 3))) next()

    # skip if values are boolean?
    cat("\ncheck_col:", check_col, "--> nvals:", n.vals)

    # add remains to cluster_cols
    cluster_cols <- c(cluster_cols, check_col)
  }

  # remove duplicate columns
  cluster_meta <- metadata[, cluster_cols, drop = FALSE]
  cluster_cols <- cluster_cols[!duplicated(as.list(cluster_meta))]
}

test_groups_equal <- function(vals1, vals2) {
  vals1 <- as.numeric(factor(vals1))
  vals2 <- as.numeric(factor(vals2))

  all(vals1 == vals2)
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
    dplyr::summarise_all(dplyr::n_distinct) |>
    dplyr::select(colnames(meta))

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
