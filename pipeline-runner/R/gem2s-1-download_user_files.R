#' Download user files from S3
#'
#' @param input The input object from the request
#' @param pipeline_config result of \code{load_config}
#' @param prev_out \code{list()} that is added to in each with gem2s task.
#'
#' @return list where 'output' slot has config used by \code{load_user_files}
#' @export
#'
download_user_files <- function(input, pipeline_config, prev_out = list()) {
  project_id <- input$projectId
  sample_uuids <- unlist(input$sampleIds)
  sample_uuids <- sample_uuids[order(sample_uuids)]

  s3 <- paws::s3(config = pipeline_config$aws_config)

  input_dir <- "/input"

  unlink(input_dir, recursive = TRUE)

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

      # Preparing directories
      local_fpath <- file.path(input_dir, sample, fname)
      fs::dir_create(dirname(local_fpath))

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
    output = prev_out
  )

  message("\nDownloading of user files step complete.")
  return(res)
}
