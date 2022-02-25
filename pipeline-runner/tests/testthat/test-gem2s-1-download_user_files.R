mock_cellranger_files <- function(sample_dir) {

  counts <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )

  features <- data.frame(
    ensid = paste0("ENSFAKE", seq_len(nrow(counts))),
    symbol = row.names(counts)
  )

  # save features
  features_path <- file.path(sample_dir, "features.tsv")
  write.table(features, features_path, col.names = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
  R.utils::gzip(features_path)

  # save barcodes
  barcodes <- colnames(counts)
  barcodes_path <- file.path(sample_dir, "barcodes.tsv")
  writeLines(barcodes, barcodes_path)
  R.utils::gzip(barcodes_path)

  # save Matrix
  if (is(counts, "data.frame")) counts <- as.matrix(counts)
  sparse.mat <- Matrix::Matrix(counts, sparse = TRUE)
  matrix_path <- file.path(sample_dir, "matrix.mtx")
  Matrix::writeMM(sparse.mat, matrix_path)
  R.utils::gzip(matrix_path)
}

mock_counts <- function() {

}

mock_features <- function() {


  features
}

create_samples <- function(bucket, project, samples) {

  for (id in names(samples)) {
    f <- fs::path(bucket, project, id, samples[[id]])
    fs::dir_create(dirname(f))


    fs::file_create(f)
  }
}

local_create_samples <- function(project, samples, env = parent.frame()) {
  # creates and deletes sample folders in mocked s3 bucket
  bucket <- fs::path("./a_fake_bucket")
  fs::dir_create(bucket)

  create_samples(bucket, project, samples)
  withr::defer(fs::dir_delete(bucket), envir = env)

  bucket
}


# test_that("mocking bucket works", {
#     projectId <- "projectID"
#   files_10x <- c("features.tsv", "barcodes.tsv", "matrix.mtx")
#   samples <- list(sample1 = files_10x, sample2 = files_10x)
#
#   bucket <- local_create_samples(projectId, samples)
#
#   print(dir_ls(fs::path(bucket, projectId), full.names = T))
#   #print(list.files(bucket, recursive = T, full.names = T))
#
#   expect_true(is.character(bucket))
# })


stub_s3_list_objects <- function(Bucket, Prefix) {

  # this workaround is the lesser evil, trust me ("bucket/./project/sample")
  Prefix <- gsub("^./", "", Prefix)

  # returns list with the same structure as s3$list_objects, but with mocked paths
  files <- fs::dir_ls(fs::path(Bucket, Prefix), type = "file", recurse = TRUE)
  l <- as.list(files)
  names(l) <- NULL
  l2 <- lapply(l, list)

  for (i in seq_along(l2)) {
    names(l2[[i]]) <- "Key"
  }
  list(Contents = l2)
}

stub_s3_get_objects <- function(Bucket, Key){
  # returns a list of the raw file read from the mocked s3 bucket.
    list(body = readBin(Key, what="raw"), rest = list())
}

# stub_dir_create <- function(path) {
#   fs::dir_create(fs::path("./", path))
# }

stub_file.path <- function(...) {
  file.path(".", ...)
}

stubbed_up_download_user_files <- function(input, pipeline_config, prev_out = list()) {
  # helper to simplify calls to the stubbed function

  mockery::stub(where = download_user_files, "paws::s3", how = NA)
  mockery::stub(where = download_user_files, "s3$list_objects", how = stub_s3_list_objects)
  mockery::stub(where = download_user_files, "s3$get_object", how = stub_s3_get_objects)
  #mockery::stub(where = download_user_files, "fs::dir_create", how = stub_dir_create)
  mockery::stub(where = download_user_files, "file.path", how = stub_file.path)

  download_user_files(input, pipeline_config, prev_out)
}

test_that("download_user_files downloads user's files.", {

    files_10x <- c("features.tsv", "barcodes.tsv", "matrix.mtx")
    samples <- list(sample1 = files_10x, sample2 = files_10x)

    input <- list(projectId = "projectID",
                  sampleNames = as.list(names(samples)),
                  sampleIds = as.list(names(samples)),
                  experimentName = "test_exp",
                  input = list(type = "techno"))


    bucket <- local_create_samples(input$projectId, samples)
    pipeline_config <- list(originals_bucket = bucket)



    res <- stubbed_up_download_user_files(input, pipeline_config)
    withr::defer(fs::dir_delete("./input"), envir = parent.frame())
})
