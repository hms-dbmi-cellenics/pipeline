upload_obj2s_to_aws <- function(input, pipeline_config, prev_out) {
  message("Uploading to AWS ...")

  experiment_id <- input$experimentId

  # destructure what need from prev_out
  scdata <- prev_out$scdata
  config <- prev_out$config

  scdata <- format_obj2s(scdata, experiment_id)

  # change sample ids/names so that get sample cell sets
  input <- add_samples_to_input(scdata, input)
  input <- add_metadata_to_input(scdata, input)
  scdata <- change_sample_names_to_ids(scdata, input)
  cell_sets <- get_cell_sets(scdata, input)

  # detect cluster metadata
  cluster_columns <- find_cluster_columns(scdata)
  if (!length(cluster_columns))
    stop(errors$ERROR_OBJ2S_CLUSTERS, call. = FALSE)

  # add clusters
  cluster_sets <- list()
  for (i in seq_along(cluster_columns)) {
    col <- cluster_columns[i]

    # first cluster column used as 'louvain' (default) clusters
    method <- ifelse(i == 1, 'louvain', col)

    cluster_sets[[i]] <- data.frame(
      cluster = scdata@meta.data[[col]],
      cell_ids = scdata$cells_id
    ) |>
      format_cluster_cellsets(method, scdata@misc$color_pool, name = col)
  }

  # cell sets file to s3
  cell_sets$cellSets <- c(cluster_sets, cell_sets$cellSets)
  cell_sets_data <- RJSONIO::toJSON(cell_sets)

  put_object_in_s3(pipeline_config,
                   bucket = pipeline_config$cell_sets_bucket,
                   object = charToRaw(cell_sets_data),
                   key = experiment_id)

  # replicate qc config for simplicity
  # could also create a 'obj2s_config' column in experiment table and change the ui/api around more
  qc_config <- construct_qc_config(list(one = scdata), unfiltered_samples = 'one', technology = config$input$type)
  qc_config$configureEmbedding$embeddingSettings$useSaved <- TRUE
  qc_config$configureEmbedding$embeddingSettings$method <- SeuratObject::DefaultDimReduc(scdata)

  # seurat object to s3
  object_key <- upload_matrix_to_s3(pipeline_config, experiment_id, scdata)
  message('Count matrix uploaded to ', pipeline_config$processed_bucket, ' with key ',object_key)

  # images for spatial to s3
  upload_images_to_s3(pipeline_config, input, experiment_id, scdata)

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

# rename active.ident (not familiar except to bioinformaticians)
# to the user supplied column name
# ensures it will still be the default clustering
rename_active_ident <- function(meta) {

  # determine columns that are duplicates of active.ident
  test_vals <- as.character(meta$active.ident)
  dup.cols <- c()
  for (col in setdiff(colnames(meta), 'active.ident')) {
    col_vals <- as.character(meta[[col]])
    if (identical(test_vals, col_vals))
      dup.cols <- c(dup.cols, col)
  }

  # remove active.ident and move originating column first column (makes default)
  if (length(dup.cols)) {
    meta$active.ident <- NULL
    meta <- meta |> dplyr::relocate(dup.cols[1])
  }

  return(meta)
}


find_cluster_columns <- function(scdata) {
  meta <- scdata@meta.data

  # exclude columns with NA values
  meta <- meta |> dplyr::select(where(~ !any(is.na(.))))

  # exclude all group columns, including duplicates
  group_cols <- find_group_columns(meta, remove.dups = FALSE)
  group_cols <- c(group_cols, 'samples')
  scdblfinder_cols <- grep('^scDblFinder', colnames(meta), value = TRUE)

  # order meta to indicate preference for default clusters
  louvain_cols <- c('louvain', 'active.ident', 'cell_type', 'seurat_clusters')
  meta <- meta |> dplyr::relocate(dplyr::any_of(louvain_cols))

  # use user supplied name for active.ident if possible
  # unless active.ident is a group column
  if ('active.ident' %in% setdiff(colnames(meta), group_cols))
    meta <- rename_active_ident(meta)

  check_cols <- setdiff(colnames(meta), c(scdblfinder_cols, group_cols))

  cluster_cols <- c()
  for (check_col in check_cols) {
    check_vals <- meta[[check_col]]

    # skip if too few or way too many values
    value_counts <- table(check_vals)
    n.vals <- length(value_counts)
    if (n.vals < 2) next()
    if (n.vals > 1000) next()
    cat("candidate cluster columns:", check_col, "--> nvals:", n.vals, '\n')

    # skip if values are numeric but non-integer
    if (is.numeric(check_vals) &&
        !all(as.integer(check_vals) == check_vals)) next()

    # skip if more than 1/3 of values are repeated fewer than 4 times
    nreps_lt4 <- sum(value_counts < 4)
    if (nreps_lt4 > n.vals/3) next()

    # skip if col is same as samples or group column
    is_sample_col <- FALSE
    for (group_col in group_cols)  {
      group_vals <- meta[[group_col]]
      if (test_groups_equal(check_vals, group_vals)) {
        is_sample_col <- TRUE
        break
      }
    }
    if (is_sample_col) next()

    # TODO: decide if should skip if values are boolean?
    # add remains to cluster_cols
    cluster_cols <- c(cluster_cols, check_col)
  }

  return(cluster_cols)
}

make_vals_numeric <- function(vals) {
  suppressWarnings({
    vals <- as.character(vals)
    as.numeric(factor(vals, levels = unique(vals)))
  })
}

test_groups_equal <- function(vals1, vals2) {
  vals1 <- make_vals_numeric(vals1)
  vals2 <- make_vals_numeric(vals2)

  all(vals1 == vals2)
}


add_samples_to_input <- function(scdata, input) {
  samples <- unique(scdata$samples)
  input$sampleNames <- samples
  input$obj2sSampleId <- input$sampleIds
  input$sampleIds <- ids::uuid(n = length(samples))
  return(input)
}

change_sample_names_to_ids <- function(scdata, input) {
  sample_ids <- input$sampleIds
  names(sample_ids) <- input$sampleNames
  scdata$samples <- unname(sample_ids[scdata$samples])
  return(scdata)
}

add_metadata_to_input <- function(scdata, input) {
  group_cols <- find_group_columns(scdata@meta.data)

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

get_n_distinct_per_sample <- function(metadata) {
  metadata |>
    dplyr::group_by(samples) |>
    dplyr::summarise_all(dplyr::n_distinct) |>
    dplyr::select(colnames(metadata))
}


# get column names that are consistent with sample groups
find_group_columns <- function(metadata, remove.dups = TRUE) {

  ndistinct_sample <- get_n_distinct_per_sample(metadata)

  ndistinct <- metadata |>
    dplyr::summarise_all(dplyr::n_distinct)

  nsamples <- length(unique(metadata$samples))


  # group columns must:
  # - have fewer than the number of samples
  # - have at least two values
  # - have only one value per sample
  too.many <- ndistinct >= nsamples
  too.few <- ndistinct <= 1
  one.per.sample <- apply(ndistinct_sample, 2, function(x) all(x == 1))
  group.cols <- names(ndistinct)[!too.many & !too.few & one.per.sample]

  # remove duplicated group columns
  if (remove.dups) {
    metadata[group.cols] <- lapply(metadata[group.cols], make_vals_numeric)
    dups <- duplicated(as.list(metadata))
    group.cols <- group.cols[group.cols %in% colnames(metadata)[!dups]]
  }

  return(group.cols)
}

# add 'cells_id'
# 'samples' must be already added
# current input$metadata not yet implemented
format_obj2s <- function(scdata, experiment_id) {

  scdata <- add_samples_col(scdata)
  scdata$cells_id <- seq_len(ncol(scdata))-1

  # other
  scdata@misc$experimentId <- experiment_id
  scdata@misc$color_pool <- get_color_pool()
  scdata@misc$ingestionDate <- Sys.time()

  # need to mock processing config
  metadata_cols <- list('percent.mt' = 0, 'doublet_scores' = 0, 'doublet_class' = 'singlet')
  scdata <- mock_metadata(scdata, metadata_cols)

  # need that logcounts and counts have same nrow
  common.genes <- intersect(row.names(scdata[['RNA']]$counts),
                            row.names(scdata[['RNA']]$data))

  scdata <- scdata[common.genes, ]
  scdata@misc$gene_annotations <- scdata@misc$gene_annotations[common.genes, ]

  return(scdata)
}

# use 'samples' or 'sample' if present, otherwise assume one sample
add_samples_col <- function(scdata) {
  samples_cols <- c('samples', 'sample')
  in.meta <- samples_cols %in% colnames(scdata@meta.data)

  if (!any(in.meta)) {
    scdata$samples <- 'NA'
  } else {
    sample_col <- samples_cols[which(in.meta)[1]]
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
