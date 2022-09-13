load_seurat <- function(input, pipeline_config, prev_out, input_dir = "/input") {

  config <- prev_out$config
  dataset_dir <- config$samples[1]
  dataset_fpath <- file.path(input_dir, dataset_dir, 'r.rds')

  scdata <- readRDS(dataset_fpath)
  prev_out$scdata <- scdata

  res <- list(
    data = list(),
    output = prev_out
  )
  return(res)
}
