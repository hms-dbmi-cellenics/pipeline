mock_files <- function(counts, features, sample_dir, type = "10x") {
  features_fnames <- list("10x" = "features.tsv", "parse" = "all_genes.csv")
  barcodes_fnames <- list("10x" = "barcodes.tsv", "parse" = "cell_metadata.csv")
  matrix_fnames <- list("10x" = "matrix.mtx", "parse" = "DGE.mtx")
  file_sep <- list("10x" = "\t", "parse" = ",")
  # Parse features and barcodes files include the header, so we need to emulate it.
  col_names <- list("10x" = FALSE, "parse" = TRUE)

  # save features
  features_path <- file.path(sample_dir, features_fnames[[type]])
  write.table(features, features_path, col.names = col_names[[type]], quote = FALSE, sep = file_sep[[type]], row.names = FALSE)
  R.utils::gzip(features_path)

  # save barcodes
  barcodes <- data.frame(barcodes = colnames(counts))
  barcodes_path <- file.path(sample_dir, barcodes_fnames[[type]])
  # writeLines(barcodes, barcodes_path)
  write.table(barcodes, barcodes_path, col.names = col_names[[type]], quote = FALSE, sep = file_sep[[type]], row.names = FALSE)
  R.utils::gzip(barcodes_path)

  # save Matrix
  if (is(counts, "data.frame")) counts <- as.matrix(counts)
  sparse.mat <- Matrix::Matrix(counts, sparse = TRUE)
  matrix_path <- file.path(sample_dir, matrix_fnames[[type]])
  Matrix::writeMM(sparse.mat, matrix_path)
  R.utils::gzip(matrix_path)
}

local_experiment <- function(counts, features, experiment_dir, sample_dir, env = parent.frame(), type = "10x") {
  sample_path <- file.path(experiment_dir, sample_dir)
  dir.create(sample_path, recursive = T)
  mock_files(counts, features, sample_path, type)
  withr::defer(unlink(experiment_dir, recursive = T), envir = env)
}
