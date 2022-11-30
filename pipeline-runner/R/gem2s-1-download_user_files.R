# Download the file and stores the output in a file
download_and_store <- function(bucket, key, file_path, s3) {
  fs::dir_create(dirname(file_path))

  # Download the file and store the output in a variable
  c(body, ...rest) %<-% s3$get_object(
    Bucket = bucket,
    Key = key
  )

  # Write output to file
  writeBin(body, con = file_path)
}

#' Download user files from S3
#'
#' @param input The input params received from the api
#' @param originals_bucket The bucket where the sample files are
#' @param aws_config The input object from the request
#' @param input_dir A string, where the s3 object is going to be stored
#' @param s3 S3 aws client
#'
#' @export
#'
get_gem2s_file <- function(input, originals_bucket, input_dir, s3) {
  project_id <- input$projectId
  sample_ids <- input$sampleIds
  sample_s3_paths <- input$sampleS3Paths
  technology <- input$input$type

  unlink(input_dir, recursive = TRUE)

  for (sample_id in sample_ids) {
    for (file_type in file_types_by_technology[[technology]]) {
      s3_path <- sample_s3_paths[[sample_id]][[file_type]]

      local_fpath <- file.path(input_dir, sample_id, file_names[[file_type]])
      download_and_store(originals_bucket, s3_path, local_fpath, s3)
    }
  }
}

#' Download user files from S3
#'
#' @param input The input object from the request
#' @param pipeline_config result of \code{load_config}
#' @param input_dir A string, where the s3 object is going to be stored
#' @param prev_out \code{list()} that is added to in each with gem2s task.
#'
#' @return list where 'output' slot has config used by \code{load_user_files}
#' @export
#'
download_user_files <- function(input, pipeline_config, prev_out = list(), input_dir = "/input") {
  s3 <- paws::s3(config = pipeline_config$aws_config)

  get_gem2s_file(input, pipeline_config$originals_bucket, input_dir, s3)

  config <- list(
    name = input$experimentName,
    samples = input$sampleIds,
    organism = input$organism,
    input = list(type = input$input$type),
    sampleOptions = input$sampleOptions
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
