mock_cellranger_files <- function(sample_dir, compressed, api_version, sample_id) {
  if (api_version == "v1") {
    features_path <- file.path(sample_dir, "features.tsv")
    barcodes_path <- file.path(sample_dir, "barcodes.tsv")
    matrix_path <- file.path(sample_dir, "matrix.mtx")
  } else if (api_version == "v2") {
    features_path <- file.path(sample_dir, paste("features10x", "path", sample_id, sep = "-"))
    barcodes_path <- file.path(sample_dir, paste("barcodes10x", "path", sample_id, sep = "-"))
    matrix_path <- file.path(sample_dir, paste("matrix10x", "path", sample_id, sep = "-"))
  }

  counts <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )

  features <- data.frame(
    ensid = paste0("ENSFAKE", seq_len(nrow(counts))),
    symbol = row.names(counts)
  )

  # save features
  write.table(features,
    features_path,
    col.names = FALSE,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )

  # save barcodes
  barcodes <- colnames(counts)
  writeLines(barcodes, barcodes_path)

  # save Matrix
  if (is(counts, "data.frame")) counts <- as.matrix(counts)
  sparse.mat <- Matrix::Matrix(counts, sparse = TRUE)
  Matrix::writeMM(sparse.mat, matrix_path)

  files <- c(features_path, barcodes_path, matrix_path)

  if (compressed && api_version == "v1") {
    R.utils::gzip(features_path)
    R.utils::gzip(barcodes_path)
    R.utils::gzip(matrix_path)

    files <- paste0(files, ".gz")
  }

  return(files)
}

create_samples <- function(bucket, project, samples, compressed, env, api_version) {
  # helper to create samples in project in bucket, like S3
  files <- c()

  for (id in samples) {
    if(api_version == "v1") {
      f <- file.path(bucket, project, id)
      dir.create(f, recursive = TRUE)

    } else if (api_version == "v2") {
      f <- paste(bucket, id, sep = "-")
      dir.create(f)
    }

    these_files <- mock_cellranger_files(f, compressed, api_version, id)
    files <- c(files, these_files)
  }

  return(files)
}

local_create_samples <- function(project, samples, compressed = FALSE, env = parent.frame(), api_version = "v1") {
  # calls creates_samples but makes them "local" (in withr speech), deleting
  # created stuff after the test finishes.
  bucket <- "./a_fake_bucket"
  dir.create(bucket)


  files <- create_samples(bucket, project, samples, compressed, env, api_version)

  if (api_version == "v2") {
    # api v2 creates several directories (one per sample), so we have to remove
    # then all, passing  a character vector with all dirnames
    withr::defer(unlink(unique(dirname(files)), recursive = TRUE), envir = env)
  }

  withr::defer(unlink(bucket, recursive = TRUE), envir = env)

  list(bucket = bucket, files = files)
}

stub_s3_list_objects <- function(Bucket, Prefix) {

  # this workaround is the lesser evil ("bucket/./project/sample")
  Prefix <- gsub("^./", "", Prefix)

  # returns list with structure like s3$list_objects, but with mocked paths
  files <- list.files(file.path(Bucket, Prefix), full.names = TRUE, recursive = TRUE)
  l <- as.list(files)
  names(l) <- NULL
  l2 <- lapply(l, list)

  for (i in seq_along(l2)) {
    names(l2[[i]]) <- "Key"
  }
  list(Contents = l2)
}

stub_s3_get_object <- function(Bucket, Key) {
  # returns a list of the raw file read from the mocked s3 bucket.
  list(body = readBin(Key, what = "raw"), rest = list())
}

stub_file.path <- function(...) {
  file.path(".", ...)
}

stubbed_download_user_files <- function(input, pipeline_config, prev_out = list()) {
  # helper to simplify calls to the stubbed function

  mockedS3 <- list(
    list_objects=stub_s3_list_objects,
    get_object=stub_s3_get_object
  )

  # where makes sure where we are stubbing the what calls.
  mockery::stub(where = download_user_files, what = "paws::s3", how=mockedS3)
  mockery::stub(get_gem2s_file_v1, "s3$list_objects", mockedS3$list_objects)

  mockery::stub(download_user_files, "file.path", stub_file.path)
  mockery::stub(get_gem2s_file_v1, "file.path", stub_file.path)
  mockery::stub(get_gem2s_file_v2, "file.path", stub_file.path)

  mockery::stub(download_and_store, "s3$get_object", mockedS3$get_object)

  res <- download_user_files(input, pipeline_config, input_dir = "./input", prev_out = prev_out)
  # download_user_files creates a "/input" folder in the pod. defer deleting
  # it during tests.
  withr::defer(unlink("./input", recursive = TRUE), envir = parent.frame())

  res
}

mock_sample_ids <- function(n_samples = 1) {
  paste0("sample_", seq_len(n_samples))
}

mock_input_v1 <- function(samples) {
  input <- list(
    projectId = "projectID",
    sampleIds = as.list(samples),
    sampleNames = as.list(paste0(samples, "_name")),
    experimentName = "test_exp",
    input = list(type = "techno"),
    apiVersion = "v1"
  )
}

mock_input_v2 <- function(samples) {
  sample_s3_paths <- list()

  for (sample_id in samples) {
    sample_s3_paths[sample_id] <- list(
      "barcodes10x" = paste("barcodes10x", "path", sample_id, sep = "-"),
      "features10x" = paste("features10x", "path", sample_id, sep = "-"),
      "matrix10x" = paste("matrix10x", "path", sample_id, sep = "-")
    )
  }

  input <- list(
    projectId = "projectID",
    sampleIds = as.list(samples),
    sampleNames = as.list(paste0(samples, "_name")),
    sampleS3Paths = sample_s3_paths,
    experimentName = "test_exp",
    input = list(type = "techno"),
    apiVersion = "v2"
  )
}

test_that("download_user_files downloads user's files. one sample", {
  samples <- mock_sample_ids()
  input <- mock_input_v1(samples)
  s3_stuff <- local_create_samples(input$projectId, samples)
  pipeline_config <- list(originals_bucket = s3_stuff$bucket)

  res <- stubbed_download_user_files(input, pipeline_config)

  # download_user_files does not return the paths. So have to build them
  downloaded_file_paths <- gsub(
    file.path(s3_stuff$bucket, input$projectId),
    "./input", s3_stuff$files
  )

  # read the downloaded files as raw files
  expected_files <- lapply(s3_stuff$files, readBin, what = "raw")
  downloaded_files <- lapply(downloaded_file_paths, readBin, what = "raw")

  expect_identical(expected_files, downloaded_files)
})


test_that("download_user_files downloads user's files. 3 samples", {
  samples <- mock_sample_ids(n_samples = 3)
  input <- mock_input_v1(samples)
  s3_stuff <- local_create_samples(input$projectId, samples)
  pipeline_config <- list(originals_bucket = s3_stuff$bucket)

  res <- stubbed_download_user_files(input, pipeline_config)

  # download_user_files does not return the paths. So have to build them
  downloaded_file_paths <- gsub(
    file.path(s3_stuff$bucket, input$projectId),
    "./input", s3_stuff$files
  )

  # read the downloaded files as raw files
  expected_files <- lapply(s3_stuff$files, readBin, what = "raw")
  downloaded_files <- lapply(downloaded_file_paths, readBin, what = "raw")

  expect_identical(expected_files, downloaded_files)
})

test_that("metadata is passed over correctly", {
  samples <- mock_sample_ids(n_samples = 3)
  input <- mock_input_v1(samples)
  s3_stuff <- local_create_samples(input$projectId, samples)
  pipeline_config <- list(originals_bucket = s3_stuff$bucket)

  # metadata shape does not matter here
  expected_metadata <- data.frame(
    a = seq_len(30),
    b = paste0("meta_", seq_len(30))
  )

  input$metadata <- expected_metadata

  res <- stubbed_download_user_files(input, pipeline_config)

  downloaded_metadata <- res$output$config$metadata

  expect_identical(expected_metadata, downloaded_metadata)
})


test_that("download_user_files correctly downloads compressed files", {
  samples <- mock_sample_ids(n_samples = 2)
  input <- mock_input_v1(samples)
  s3_stuff <- local_create_samples(input$projectId, samples, compressed = TRUE)
  pipeline_config <- list(originals_bucket = s3_stuff$bucket)

  res <- stubbed_download_user_files(input, pipeline_config)

  # download_user_files does not return the paths. So have to build them
  downloaded_file_paths <- gsub(
    file.path(s3_stuff$bucket, input$projectId),
    "./input", s3_stuff$files
  )
  # read the downloaded files as raw files
  expected_files <- lapply(s3_stuff$files, readBin, what = "raw")
  downloaded_files <- lapply(downloaded_file_paths, readBin, what = "raw")

  expect_identical(expected_files, downloaded_files)
})

test_that("download_user_files downloads user's files. one sample in v2", {
  samples <- mock_sample_ids()
  input <- mock_input_v2(samples)
  s3_stuff <- local_create_samples(input$projectId, samples, api_version = "v2")
  pipeline_config <- list(originals_bucket = s3_stuff$bucket)

  res <- stubbed_download_user_files(input, pipeline_config)

  # download_user_files does not return the paths. So have to build them
  downloaded_file_paths <- gsub(
    file.path(s3_stuff$bucket, input$projectId),
    "./input", s3_stuff$files
  )

  # read the downloaded files as raw files
  expected_files <- lapply(s3_stuff$files, readBin, what = "raw")
  downloaded_files <- lapply(downloaded_file_paths, readBin, what = "raw")

  expect_identical(expected_files, downloaded_files)
})


test_that("download_user_files downloads user's files. 3 samples in v2", {
  samples <- mock_sample_ids(n_samples = 3)
  input <- mock_input_v2(samples)
  s3_stuff <- local_create_samples(input$projectId, samples, api_version = "v2")
  pipeline_config <- list(originals_bucket = s3_stuff$bucket)

  res <- stubbed_download_user_files(input, pipeline_config)

  # download_user_files does not return the paths. So have to build them
  downloaded_file_paths <- gsub(
    file.path(s3_stuff$bucket, input$projectId),
    "./input", s3_stuff$files
  )

  # read the downloaded files as raw files
  expected_files <- lapply(s3_stuff$files, readBin, what = "raw")
  downloaded_files <- lapply(downloaded_file_paths, readBin, what = "raw")

  expect_identical(expected_files, downloaded_files)
})

test_that("metadata is passed over correctly in v2", {
  samples <- mock_sample_ids(n_samples = 3)
  input <- mock_input_v2(samples)
  s3_stuff <- local_create_samples(input$projectId, samples, api_version = "v2")
  pipeline_config <- list(originals_bucket = s3_stuff$bucket)

  # metadata shape does not matter here
  expected_metadata <- data.frame(
    a = seq_len(30),
    b = paste0("meta_", seq_len(30))
  )

  input$metadata <- expected_metadata

  res <- stubbed_download_user_files(input, pipeline_config)

  downloaded_metadata <- res$output$config$metadata

  expect_identical(expected_metadata, downloaded_metadata)
})


test_that("download_user_files correctly downloads compressed files in v2", {
  samples <- mock_sample_ids(n_samples = 2)
  input <- mock_input_v2(samples)
  s3_stuff <- local_create_samples(input$projectId, samples, compressed = TRUE, api_version = "v2")
  pipeline_config <- list(originals_bucket = s3_stuff$bucket)

  res <- stubbed_download_user_files(input, pipeline_config)

  # download_user_files does not return the paths. So have to build them
  downloaded_file_paths <- gsub(
    file.path(s3_stuff$bucket, input$projectId),
    "./input", s3_stuff$files
  )
  # read the downloaded files as raw files
  expected_files <- lapply(s3_stuff$files, readBin, what = "raw")
  downloaded_files <- lapply(downloaded_file_paths, readBin, what = "raw")

  expect_identical(expected_files, downloaded_files)
})
