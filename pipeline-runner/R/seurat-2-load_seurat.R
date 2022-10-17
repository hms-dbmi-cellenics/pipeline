load_seurat <- function(input, pipeline_config, prev_out, input_dir = "/input") {

  config <- prev_out$config
  dataset_dir <- config$samples[1]
  dataset_fpath <- file.path(input_dir, dataset_dir, 'r.rds')

  scdata <- RAppArmor::eval.secure(
    reconstruct_seurat(dataset_fpath),
    profile = 'read-seurat')

  prev_out$scdata <- scdata

  res <- list(
    data = list(),
    output = prev_out
  )
  return(res)
}

reconstruct_seurat <- function(dataset_fpath) {

  user_scdata <- readRDS(dataset_fpath)
  user_scdata <- add_dispersions(user_scdata)
  dispersions <- user_scdata@misc$gene_dispersion
  test_user_df(dispersions)


  counts <- user_scdata[['RNA']]@counts
  metadata <- user_scdata@meta.data
  test_user_df(metadata)
  test_user_sparse_mat(counts)

  scdata <- SeuratObject::CreateSeuratObject(
    counts,
    meta.data = metadata
  )

  red_name <- SeuratObject::DefaultDimReduc(user_scdata)
  check_type_is_safe(red_name)

  embedding <- user_scdata@reductions[[red_name]]@cell.embeddings
  test_user_df(embedding)

  red <- SeuratObject::CreateDimReducObject(
    embeddings = embedding,
    assay = 'RNA'
  )

  data <- user_scdata[['RNA']]@data
  test_user_sparse_mat(data)
  scdata[['RNA']]@data <- data
  scdata@reductions[[red_name]] <- red

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
  } else {
    stopifnot(typeof(x) %in% safe.types)
  }
}
