create_samples <- function(bucket, samples) {
  fs::dir_create(bucket)
  for (name in names(samples)) {
    f <- fs::path(bucket, "test_samples", name, samples[[name]])
    fs::dir_create(dirname(f))
    fs::file_create(f)
  }
}

local_create_samples <- function(samples, env = parent.frame()) {
  bucket <- tempdir()
  create_samples(bucket, samples)
  withr::defer(fs::dir_delete(bucket), envir = env)

  bucket
}


test_that("download_user_files downloads user's files.", {

  files_10x <- c("features.tsv", "barcodes.tsv", "matrix.mtx")
  samples <- list(sample1 = files_10x, sample2 = files_10x)

  bucket <- local_create_samples(samples)

  print(list.files(bucket, recursive = T, full.names = T))
  expect_true(is.character(bucket))

})


