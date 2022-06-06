mock_cellranger_files <- function(counts, features, sample_dir) {

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
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )
}

local_cellranger_experiment <- function(counts, features, experiment_dir, sample_dir, env = parent.frame()) {

  sample_path <- file.path(experiment_dir, sample_dir)
  dir.create(sample_path, recursive = T)
  mock_cellranger_files(counts, features, sample_path)
  withr::defer(unlink(experiment_dir, recursive = T), envir = env)

}

mock_rhapsody_matrix <- function(counts, sample_dir) {
  counts$Gene <- rownames(counts)
  counts <- tidyr::pivot_longer(counts, -Gene,
    names_to = "barcode",
    values_to = "DBEC_Adjusted_Molecules"
  )

  counts$Cell_Index <- as.integer(factor(counts$barcode))
  counts$RSEC_Adjusted_Molecules <- counts$DBEC_Adjusted_Molecules + 5

  matrix_path <- file.path(sample_dir, "expression_matrix.st")

  # prepend some of that nice header
  header <- c(
    "####################",
    "## BD Targeted Multiplex Rhapsody Analysis Pipeline Version 1.9.1",
    "## Analysis Date: 2020-09-29 23:41:40",
    "## Sample: SampleMultiplexDemo",
    "## Reference: BD_Rhapsody_Immune_Response_Panel_Hs.fasta",
    "## Sample Tags Version: Single-Cell Multiplex Kit - Human",
    "####################"
  )
  writeLines(header, matrix_path)
  write.table(counts,
    file = matrix_path,
    append = TRUE,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )

  matrix_path
}


local_rhapsody_experiment <- function(samples, env = parent.frame()) {
  # calls creates_samples but makes them "local" (in withr speech), deleting
  # created stuff after the test finishes.
  bucket <- "./input"
  dir.create(bucket)
  files <- c()

  for (sample in samples) {
    sample_path <- file.path(bucket, sample$name)
    dir.create(sample_path)
    files <- c(files, mock_rhapsody_matrix(sample$counts, sample_path))
  }

  withr::defer(unlink(bucket, recursive = TRUE), envir = env)
  files
}


test_that("format_annot keeps unique rows", {
  annot_list <- list(
    sample1 = data.frame(ENSID = 1:5, SYMBOL = paste0("gene", 1:5)),
    sample2 = data.frame(ENSID = 1:5, SYMBOL = paste0("gene", 1:5))
  )

  annot <- pipeline:::format_annot(annot_list)

  expect_s3_class(annot, "data.frame")
  expect_true(nrow(annot) == nrow(annot_list$sample1))
})


test_that("format_annot deduplicates name column", {
  annot_list <- list(
    sample1 = data.frame(ENSID = 1:6, SYMBOL = paste0("gene", c(1, 1:5)))
  )

  annot <- pipeline:::format_annot(annot_list)

  expect_true(length(annot$name) == length(unique(annot$name)))
})


test_that("format_annot removes duplicated input (Ensembl IDs) column", {
  ensids <- c(1, 1:4)
  annot_list <- list(
    sample1 = data.frame(ENSID = ensids, SYMBOL = paste0("gene", 1:5))
  )

  annot <- pipeline:::format_annot(annot_list)

  expect_equal(
    length(unique(ensids)),
    length(annot$input)
  )
})


test_that("load_user_files loads a 10x count matrix", {
  counts <- mock_counts()
  features <- data.frame(
    ensid = paste0("ENSFAKE", seq_len(nrow(counts))),
    symbol = row.names(counts)
  )

  outdir <- tempdir()
  sample <- "sample_a"
  sample_dir <- file.path(outdir, sample)
  dir.create(sample_dir)

  mock_cellranger_files(counts, features, sample_dir)

  prev_out <- list(config = list(samples = sample, input = list(type = "10x")))
  out <- load_user_files(NULL, NULL, prev_out, outdir)$output

  expect_true("counts_list" %in% names(out))
  expect_true(sample %in% names(out$counts_list))

  expect_s4_class(out$counts_list[[1]], "dgCMatrix")
  unlink(sample_dir, recursive = TRUE)
})


test_that("load_user_files generates feature annotation for 10x data", {
  counts <- mock_counts()
  features <- data.frame(
    ensid = paste0("ENSFAKE", seq_len(nrow(counts))),
    symbol = row.names(counts)
  )

  outdir <- tempdir()
  sample <- "sample_a"
  sample_dir <- file.path(outdir, sample)
  dir.create(sample_dir)

  mock_cellranger_files(counts, features, sample_dir)

  prev_out <- list(config = list(samples = sample, input = list(type = "10x")))
  out <- load_user_files(NULL, NULL, prev_out, outdir)$output

  expect_true("annot" %in% names(out))
  expect_true(
    all(c("input", "name", "original_name") %in% colnames(out$annot))
  )

  unlink(sample_dir, recursive = TRUE)
})


test_that("load_user_files deduplicates gene symbols for 10x data", {
  counts <- mock_counts()

  symbols <- row.names(counts)
  symbols[1:5] <- "DUPLICATED"

  features <- data.frame(
    ensid = paste0("ENSFAKE", seq_len(nrow(counts))),
    symbol = symbols
  )

  outdir <- tempdir()
  sample <- "sample_a"
  sample_dir <- file.path(outdir, sample)
  dir.create(sample_dir)

  mock_cellranger_files(counts, features, sample_dir)

  prev_out <- list(config = list(samples = sample, input = list(type = "10x")))
  annot <- load_user_files(NULL, NULL, prev_out, outdir)$output$annot

  # unique gene names is same as number of gene names
  expect_length(unique(annot$name), length(symbols))

  # unique original names is same as unique gene names
  expect_length(unique(annot$original_name), length(unique(symbols)))
  unlink(sample_dir, recursive = TRUE)
})


test_that("load_user_files uses appropriate feature columns for 10x data", {
  counts <- mock_counts()

  symbols <- row.names(counts)
  features <- data.frame(
    ensid = paste0("ENSFAKE", seq_len(nrow(counts))),
    symbol = symbols
  )

  outdir <- tempdir()
  sample <- "sample_a"
  sample_dir <- file.path(outdir, sample)
  dir.create(sample_dir)

  mock_cellranger_files(counts, features, sample_dir)

  prev_out <- list(config = list(samples = sample, input = list(type = "10x")))
  out <- load_user_files(NULL, NULL, prev_out, outdir)$output

  # ensembl ids are counts row names
  expect_equal(
    features$ensid,
    row.names(out$counts_list[[1]])
  )

  # ensembl ids are in column 'input' of annot
  expect_equal(out$annot$input, features$ensid)

  # symbols are in column 'name' of annot
  expect_equal(out$annot$name, symbols)

  unlink(sample_dir, recursive = TRUE)
})


test_that("load_user_files loads 10x multisample experiments", {
  counts <- mock_counts()

  symbols <- row.names(counts)
  ensids <- paste0("ENSFAKE", seq_len(nrow(counts)))

  features <- data.frame(ensid = ensids, symbol = symbols)

  # most overlapping, some unique
  ensids2 <- ensids
  symbols2 <- symbols
  ensids2[1:20] <- paste0(ensids2[1:20], "_2")
  symbols2[1:20] <- paste0(symbols2[1:20], "_2")

  features2 <- data.frame(ensid = ensids2, symbols = symbols2)

  outdir <- tempdir()
  samples <- c("sample_a", "samble_b")
  sample_dirs <- file.path(outdir, samples)
  dir.create(sample_dirs[1])
  dir.create(sample_dirs[2])

  mock_cellranger_files(counts, features, sample_dirs[1])
  mock_cellranger_files(counts, features2, sample_dirs[2])

  prev_out <- list(config = list(samples = samples, input = list(type = "10x")))
  out <- load_user_files(NULL, NULL, prev_out, outdir)$output

  # loaded both
  expect_equal(names(out$counts_list), samples)

  # kept only unique ensembl ids
  unique_ensids <- unique(c(ensids, ensids2))
  expect_length(out$annot$input, length(unique_ensids))

  # removed duplicated ensembl ids
  expect_lt(
    length(unique_ensids),
    length(c(ensids, ensids2))
  )

  unlink(sample_dirs, recursive = TRUE)
})


test_that("load_user_files reads rhapsody files", {
  samples <- list(sample_1 = list(name = "sample_1", counts = mock_counts()))

  files <- local_rhapsody_experiment(samples)

  prev_out <- list(config = list(samples = names(samples), input = list(type = "rhapsody")))
  input_dir <- "./input"

  res <- load_user_files(NULL, NULL, prev_out, input_dir)

  expect_true("counts_list" %in% names(res$output))
  expect_true(names(samples) %in% names(res$output$counts_list))
  expect_s4_class(res$output$counts_list[[1]], "dgCMatrix")
})

test_that("read_rhapsody_files reads a rhapsody matrix", {
  samples <- list(sample_1 = list(name = "sample_1", counts = mock_counts()))

  files <- local_rhapsody_experiment(samples)

  config <- list(samples = names(samples))
  input_dir <- "./input"

  res <- read_rhapsody_files(config, input_dir)

  expect_true("counts_list" %in% names(res))
  expect_true(names(samples) %in% names(res$counts_list))
  expect_s4_class(res$counts_list[[1]], "dgCMatrix")
})


test_that("parse_rhapsody_matrix reads a rhapsody matrix", {
  samples <- list(sample_1 = list(name = "sample_1", counts = mock_counts()))

  files <- local_rhapsody_experiment(samples)

  config <- list(samples = names(samples))
  input_dir <- "./input"

  res <- parse_rhapsody_matrix(config, input_dir)

  expect_true("counts_list" %in% names(res))
  expect_true(names(samples) %in% names(res$counts_list))
  expect_s4_class(res$counts_list[[1]], "dgCMatrix")
})


test_that("parse_rhapsody_matrix keeps the counts where it counts (correct gene-cell-value)", {
  samples <- list(sample_1 = list(name = "sample_1", counts = mock_counts()))
  files <- local_rhapsody_experiment(samples)
  config <- list(samples = names(samples))
  input_dir <- "./input"

  # read original table and get vector of expected values (originals)
  original <- data.table::fread(file.path(input_dir, names(samples), "expression_matrix.st"))
  expected_values <- original$DBEC_Adjusted_Molecules

  res <- parse_rhapsody_matrix(config, input_dir)

  # create row and column indices, to cbind and use matrix indexing to get values
  # as a vector. Given that Simple Triplet sparse matrices basically contain vectors
  # of row and column indices, we just transform them to ints, keeping order.
  row_idx <- as.integer(factor(original$Gene))
  col_idx <- match(original$Cell_Index, unique(original$Cell_Index))

  values <- res$counts_list$sample_1[cbind(row_idx, col_idx)]

  expect_equal(values, expected_values)
})


test_that("read_10x_files returns error if files missing", {
  counts <- mock_counts()
  features <- data.frame(
    ensid = paste0("ENSFAKE", seq_len(nrow(counts))),
    symbol = row.names(counts)
  )

  outdir <- tempdir()
  sample <- "sample_a"
  sample_dir <- file.path(outdir, sample)
  dir.create(sample_dir)

  mock_cellranger_files(counts, features, sample_dir)

  prev_out <- list(config = list(samples = sample, input = list(type = "10x")))

  files <- c("features.tsv.gz", "barcodes.tsv.gz", "matrix.mtx.gz")

  # remove files one by one renaming
  for (file in files) {
    file.rename(file.path(sample_dir, file), file.path(sample_dir, "blah"))
    expect_error(load_user_files(NULL, NULL, prev_out, outdir), "file missing")
    file.rename(file.path(sample_dir, "blah"), file.path(sample_dir, file))
  }

  unlink(sample_dir, recursive = TRUE)
})


test_that("parse_rhapsody_matrix returns error if files are missing", {
  samples <- list(sample_1 = list(name = "sample_1", counts = mock_counts()))

  files <- local_rhapsody_experiment(samples)

  config <- list(samples = names(samples))
  input_dir <- "./input"

  # remove file
  unlink(files[1])

  expect_error(
    parse_rhapsody_matrix(config, input_dir),
    "File .* does not exist or is non-readable."
  )
})


test_that("parse_rhapsody_matrix returns error if a column is invalid", {
  samples <- list(sample_1 = list(name = "sample_1", counts = mock_counts()))

  files <- local_rhapsody_experiment(samples)

  config <- list(samples = names(samples))
  input_dir <- "./input"

  # remove file
  tab <- data.table::fread(files[1])
  tab[, DBEC_Adjusted_Molecules := paste0("corrupt_", DBEC_Adjusted_Molecules)]
  data.table::fwrite(tab, file = files[1])

  expect_error(
    parse_rhapsody_matrix(config, input_dir)
  )
})


test_that("parse_rhapsody_matrix uses RSEC if DBEC corrected counts are missing", {
  samples <- list(sample_1 = list(name = "sample_1", counts = mock_counts()))

  files <- local_rhapsody_experiment(samples)

  config <- list(samples = names(samples))
  input_dir <- "./input"

  # remove DBEC column
  original <- data.table::fread(files[1])
  original[, DBEC_Adjusted_Molecules := NULL]
  data.table::fwrite(original, file = files[1])

  # keep RSEC values
  expected_values <- original$RSEC_Adjusted_Molecules

  res <- pipeline:::parse_rhapsody_matrix(config, input_dir)

  row_idx <- as.integer(factor(original$Gene))
  col_idx <- match(original$Cell_Index, unique(original$Cell_Index))

  values <- res$counts_list$sample_1[cbind(row_idx, col_idx)]

  expect_equal(values, expected_values)
})


test_that("read_10x_files removes rows with empty feature names both in count matrix and annotation if present and < 0.1%", {
  # mock count matrix replicating it 10 times to mock a matrix with < 0.1% of empty features
  counts <- mock_counts()[rep(seq_len(nrow(mock_counts())), each = 10), ]
  rownames(counts)[2] <- ""
  rownames(counts)[3] <- ".1"

  features <- data.frame(
    ensid = paste0("ENSFAKE", seq_len(nrow(counts))),
    symbol = row.names(counts)
  )
  features[2:3, 1:2] <- ""

  experiment_dir <- "./experiment_1"
  sample <- "sample_a"

  local_cellranger_experiment(counts, features, experiment_dir, sample)

  prev_out <- list(config = list(samples = sample, input = list(type = "10x")))

  out <- load_user_files(NULL, NULL, prev_out, experiment_dir)$output

  counts_list <- out$counts_list
  annot <- out$annot

  expect_equal(length(which(rownames(counts_list[[1]]) == "")), 0)
  expect_equal(length(which(annot[, 1] == "")), 0)
})


test_that("read_10x_files removes single row with empty feature names both in count matrix and annotation if present and < 0.1%", {
  # mock count matrix replicating it 10 times to mock a matrix with < 0.1% of empty features
  counts <- mock_counts()[rep(seq_len(nrow(mock_counts())), each = 10), ]
  rownames(counts)[2] <- ""

  features <- data.frame(
    ensid = paste0("ENSFAKE", seq_len(nrow(counts))),
    symbol = row.names(counts)
  )
  features[2, 1:2] <- ""

  experiment_dir <- "./experiment_1"
  sample <- "sample_a"

  local_cellranger_experiment(counts, features, experiment_dir, sample)

  prev_out <- list(config = list(samples = sample, input = list(type = "10x")))

  out <- load_user_files(NULL, NULL, prev_out, experiment_dir)$output
  counts_list <- out$counts_list
  annot <- out$annot

  expect_equal(length(which(rownames(counts_list[[1]]) == "")), 0)
  expect_equal(length(which(annot[, 1] == "")), 0)
})


test_that("read_10x_files doesn't remove any rows with empty feature names both in count matrix and annotation if present and >= 0.1%", {
  counts <- mock_counts()
  rownames(counts)[2] <- ""
  rownames(counts)[3] <- ".1"

  features <- data.frame(
    ensid = paste0("ENSFAKE", seq_len(nrow(counts))),
    symbol = row.names(counts)
  )
  features[2:3, 1:2] <- ""

  experiment_dir <- "./experiment_1"
  sample <- "sample_a"

  local_cellranger_experiment(counts, features, experiment_dir, sample)

  prev_out <- list(config = list(samples = sample, input = list(type = "10x")))

  out <- load_user_files(NULL, NULL, prev_out, experiment_dir)$output

  counts_list <- out$counts_list
  annot <- out$annot

  expect_equal(nrow(counts), nrow(counts_list[[1]]))
  # expect_equal(nrow(features), nrow(annot))  # decomment this line and delete the following line when make unique will be added in format_annot [BIOMAGE-1817]
  expect_equal(length(which(annot[, 1] == "")), 1)
})


test_that("read_10x_files doesn't remove any rows if no rows with empty rownames are present", {
  counts <- mock_counts()

  features <- data.frame(
    ensid = paste0("ENSFAKE", seq_len(nrow(counts))),
    symbol = row.names(counts)
  )

  experiment_dir <- "./experiment_1"
  sample <- "sample_a"

  local_cellranger_experiment(counts, features, experiment_dir, sample)

  prev_out <- list(config = list(samples = sample, input = list(type = "10x")))

  out <- load_user_files(NULL, NULL, prev_out, experiment_dir)$output

  counts_list <- out$counts_list
  annot <- out$annot

  expect_equal(length(which(rownames(counts_list[[1]]) == "")), 0)
  expect_equal(length(which(annot[, 1] == "")), 0)
  expect_equal(nrow(counts), nrow(counts_list[[1]]))
  expect_equal(nrow(features), nrow(annot))
})
