load_seurat <- function(input, pipeline_config, prev_out, input_dir = "/input") {

  config <- prev_out$config
  dataset_dir <- config$samples[1]
  dataset_fpath <- file.path(input_dir, dataset_dir, 'r.rds')

  scdata <- reconstruct_seurat(dataset_fpath)

  prev_out$scdata <- scdata

  res <- list(
    data = list(),
    output = prev_out
  )
  return(res)
}


reconstruct_seurat <- function(dataset_fpath) {

  # read it
  tryCatch(
    user_scdata <- readRDS(dataset_fpath),
    error = function(e) {
      stop("ERROR_SEURAT_RDS", call. = FALSE)
    })


  # get counts
  tryCatch({
    counts <- user_scdata[['RNA']]@counts
    test_user_sparse_mat(counts)
    rns <- row.names(counts)
    check_type_is_safe(rns)
  },
  error = function(e) {
    stop('ERROR_SEURAT_COUNTS', call. = FALSE)
  })

  gene_annotations <- data.frame(input = rns, name = rns, original_name = rns, row.names = rns)
  user_scdata@misc$gene_annotations <- gene_annotations

  # add dispersions
  tryCatch({
    user_scdata <- add_dispersions(user_scdata)
    dispersions <- user_scdata@misc$gene_dispersion
    test_user_df(dispersions)
  },
  error = function(e) {
    stop('ERROR_SEURAT_HVFINFO', call. = FALSE)
  })

  # get meta data
  tryCatch({
    metadata <- user_scdata@meta.data
    test_user_df(metadata)

    # check for clusters
  },
  error = function(e) {
    stop('ERROR_SEURAT_METADATA', call. = FALSE)
  })

  # check for clusters
  if (!'seurat_clusters' %in% colnames(metadata))
    stop('ERROR_SEURAT_METADATA', call. = FALSE)


  # reconstruct seurat object
  scdata <- SeuratObject::CreateSeuratObject(
    counts,
    meta.data = metadata,
  )

  # add dispersions and gene annotations to new scdata
  scdata@misc$gene_dispersion <- dispersions
  scdata@misc$gene_annotations <- gene_annotations


  # add dimensionality reduction
  tryCatch({
    red_name <- SeuratObject::DefaultDimReduc(user_scdata)
    check_type_is_safe(red_name)
    embedding <- user_scdata@reductions[[red_name]]@cell.embeddings
    test_user_df(embedding)
    red <- SeuratObject::CreateDimReducObject(
      embeddings = embedding,
      assay = 'RNA'
    )
    scdata@reductions[[red_name]] <- red
  },
  error = function(e) {
    stop('ERROR_SEURAT_REDUCTION', call. = FALSE)
  })

  # add logcounts
  tryCatch({
    data <- user_scdata[['RNA']]@data
    test_user_sparse_mat(data)
    scdata[['RNA']]@data <- data
  },
  error = function(e) {
    stop('ERROR_SEURAT_LOGCOUNTS', call. = FALSE)
  })

  return(scdata)
}


test_user_df <- function(df) {
  for (col in colnames(df)) {
    check_type_is_safe(df[, col])
  }
}

test_user_sparse_mat <- function(counts) {
  stopifnot(isS4(counts))

  for (slot in methods::slotNames(counts)) {
    check_type_is_safe(methods::slot(counts, slot))
  }
}

check_type_is_safe <- function(x) {

  # typeof determines the R internal type or storage mode of any object
  # whereas methods::is can be tricked by setting class inappropriately (e.g. disguising a function as a numeric)
  safe.types <- c('character', 'double', 'integer', 'logical')

  # recurse into lists until reach node
  if (typeof(x) == 'list') {
    lapply(x, check_type_is_safe)
  } else if (!typeof(x) %in% safe.types) {
    stop('Unexpected data type in uploaded .rds file.')
  }
}
