#' Download user files from S3
#'
#' @param input The input object from the request
#' @param pipeline_config result of \code{load_config}
#' @param prev_out \code{list()} that is added to in each with gem2s task.
#'
#' @return list where 'output' slot has config used by \code{load_cellranger}
#' @export
#'
download_user_files <- function(input, pipeline_config, prev_out = list()) {
  project_id <- input$projectId
  sample_names <- input$sampleNames
  sample_uuids <- input$sampleIds

  s3 <- paws::s3(config = pipeline_config$aws_config)

  cellranger_fnames <- c("features.tsv.gz", "barcodes.tsv.gz", "matrix.mtx.gz")
  unlink("/input", recursive = TRUE)

  for (sample in sample_uuids) {
    message("\nSample --> ", sample)

    res <- s3$list_objects(
      Bucket = pipeline_config$originals_bucket,
      Prefix = file.path(project_id, sample)
    )

    keys <- unlist(lapply(res$Contents, `[[`, "Key"))

    for (gem_key in keys) {
      message("GEM key: ", gem_key)
      fname <- basename(gem_key)

      sample_name <- sample_names[[match(sample, sample_uuids)]]
      # Preparing directories
      local_dir <- file.path("/input", sample)
      dir.create("/input")
      dir.create(local_dir)
      local_fpath <- file.path(local_dir, fname)

      # Download the file and store the output in a variable
      c(body, ...rest) %<-% s3$get_object(
        Bucket = pipeline_config$originals_bucket,
        Key = gem_key
      )

      # Write output to file
      writeBin(body, con = local_fpath)
    }
  }

  config <- list(
    name = input$experimentName,
    samples = input$sampleIds,
    organism = input$organism,
    input = list(type = input$input$type)
  )

  if ("metadata" %in% names(input)) {
    config$metadata <- input$metadata
  }

  prev_out$config <- config
  res <- list(
    data = list(),
    ouput = prev_out
  )

  message("\nDownloading of cellranger files step complete.")
  return(res)
}
