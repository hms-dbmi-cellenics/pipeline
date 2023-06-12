#' Copy recursively all s3 objects matching prefix
#'
#' @param source_bucket string
#' @param source_prefix string
#' @param destination_bucket string
#' @param destination_prefix string
#' @param aws_config
#' @param key_transform function to transform keys of objects
#' before uploading, identity by default
#'
#' @return
#'
s3_copy_by_prefix <- function(
  source_bucket,
  source_prefix,
  destination_bucket,
  destination_prefix,
  aws_config,
  key_transform = identity
) {
  s3 <- paws::s3(config = aws_config)

  list_objects_response <- s3$list_objects(
    Bucket = source_bucket,
    Prefix = source_prefix
  )

  objects <- list_objects_response$Contents

  for (object in objects) {
    source_key <- object$Key
    destination_key <- gsub(source_prefix, destination_prefix, source_key, fixed = TRUE)

    destination_key <- key_transform(destination_key)

    s3$copy_object(
      CopySource = paste0(source_bucket, "/", source_key),
      Bucket = destination_bucket,
      Key = destination_key
    )

    message("Copied object:", source_key, "to", destination_key, "\n")
  }

  message("Copy by prefix completed successfully\n")
}
