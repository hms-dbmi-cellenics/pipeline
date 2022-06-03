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
#' @param input_dir A string, where the s3 object is going to be stored
#' @param s3 S3 aws client
#'
#' @export
#'
get_gem2s_file_v1 <- function(input, originals_bucket, input_dir, s3) {
  project_id <- input$projectId
  sample_ids <- unlist(input$sampleIds)
  sample_ids <- sort(sample_ids)


  unlink(input_dir, recursive = TRUE)

  for (sample_id in sample_ids) {
    message("\nSample --> ", sample_id)

    res <- s3$list_objects(
      Bucket = originals_bucket,
      Prefix = file.path(project_id, sample_id)
    )

    keys <- unlist(lapply(res$Contents, `[[`, "Key"))

    for (gem_key in keys) {
      message("GEM key: ", gem_key)
      fname <- basename(gem_key)

      # Preparing directories
      local_fpath <- file.path(input_dir, sample_id, fname)

      download_and_store(originals_bucket, gem_key, local_fpath, s3)
    }
  }
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
get_gem2s_file_v2 <- function(input, originals_bucket, input_dir, s3) {
  project_id <- input$projectId
  sample_ids <- input$sampleIds
  sample_s3_paths <- input$sampleS3Paths
  technology <- input$input$type

  unlink(input_dir, recursive = TRUE)

  for (sample_id in sample_ids) {
    for (file_type in file_types_by_technology[[technology]]) {
      s3_path <- sample_s3_paths[[sample_id]][[file_type]]

      local_fpath <- file.path(input_dir, sample_id, file_names_v1[[file_type]])
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
download_user_files <- function(input, pipeline_config, prev_out = list(), input_dir= "/input") {
  s3 <- paws::s3(config = pipeline_config$aws_config)

  if (input$apiVersion == "v1") {
    get_gem2s_file_v1(input, pipeline_config$originals_bucket, input_dir, s3)
  } else if (input$apiVersion == "v2") {
    get_gem2s_file_v2(input, pipeline_config$originals_bucket, input_dir, s3)
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
