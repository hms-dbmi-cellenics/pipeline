mock_rhapsody_matrix <- function(counts, sample_dir, with_abseq) {
  counts$Gene <- rownames(counts)
  counts <- tidyr::pivot_longer(counts, -Gene,
    names_to = "barcode",
    values_to = "DBEC_Adjusted_Molecules"
  )

  counts$Cell_Index <- as.integer(factor(counts$barcode))
  counts$RSEC_Adjusted_Molecules <- counts$DBEC_Adjusted_Molecules + 5

  if (with_abseq) {
    counts <- counts |>
      dplyr::rename(Bioproduct = Gene) |>
      dplyr::mutate(Bioproduct = ifelse(
        startsWith(Bioproduct, "M"),
        paste0(Bioproduct, "_pAbO"),
        Bioproduct
      ))
  }

  matrix_path <- file.path(sample_dir, file_names[["rhapsody"]])

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

mock_counts <- function() {
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )
}

local_rhapsody_experiment <- function(samples, with_abseq = FALSE, env = parent.frame()) {
  # calls creates_samples but makes them "local" (in withr speech), deleting
  # created stuff after the test finishes.
  bucket <- "./input"
  dir.create(bucket)
  files <- c()

  for (sample in samples) {
    sample_path <- file.path(bucket, sample$name)
    dir.create(sample_path)
    files <- c(files, mock_rhapsody_matrix(sample$counts, sample_path, with_abseq))
  }

  withr::defer(unlink(bucket, recursive = TRUE), envir = env)
  files
}

mock_prev_out <- function(samples, include_abseq) {

  sample_options <- list()
  for (sample in names(samples)) {
    sample_options[[sample]] = list(includeAbSeq = include_abseq)
  }

  prev_out <- list(
    config = list(
      samples = names(samples),
      input = list(type = "rhapsody"),
      sampleOptions = sample_options
    )
  )

  return(prev_out)
}


test_that("load_user_files reads rhapsody files", {
  samples <- list(sample_1 = list(name = "sample_1", counts = mock_counts()))

  files <- local_rhapsody_experiment(samples)


  prev_out <- mock_prev_out(samples, include_abseq = TRUE)
  input_dir <- "./input"

  res <- load_user_files(NULL, NULL, prev_out, input_dir)

  expect_true("counts_list" %in% names(res$output))
  expect_true(names(samples) %in% names(res$output$counts_list))
  expect_s4_class(res$output$counts_list[[1]], "dgCMatrix")
})

test_that("read_rhapsody_files reads a rhapsody matrix", {
  samples <- list(sample_1 = list(name = "sample_1", counts = mock_counts()))

  files <- local_rhapsody_experiment(samples)

  prev_out <- mock_prev_out(samples, include_abseq = TRUE)
  config <- prev_out$config
  input_dir <- "./input"

  res <- read_rhapsody_files(config, input_dir)

  expect_true("counts_list" %in% names(res))
  expect_true(names(samples) %in% names(res$counts_list))
  expect_s4_class(res$counts_list[[1]], "dgCMatrix")
})


test_that("parse_rhapsody_matrix reads a rhapsody matrix", {
  samples <- list(sample_1 = list(name = "sample_1", counts = mock_counts()))

  files <- local_rhapsody_experiment(samples)

  prev_out <- mock_prev_out(samples, include_abseq = TRUE)
  config <- prev_out$config
  input_dir <- "./input"

  res <- parse_rhapsody_matrix(config, input_dir)

  expect_true("counts_list" %in% names(res))
  expect_true(names(samples) %in% names(res$counts_list))
  expect_s4_class(res$counts_list[[1]], "dgCMatrix")
})


test_that("parse_rhapsody_matrix keeps the counts where it counts (correct gene-cell-value)", {
  samples <- list(sample_1 = list(name = "sample_1", counts = mock_counts()))
  files <- local_rhapsody_experiment(samples)
  prev_out <- mock_prev_out(samples, include_abseq = TRUE)
  config <- prev_out$config
  input_dir <- "./input"

  # read original table and get vector of expected values (originals)
  original_path <- file.path(input_dir, names(samples), file_names[["rhapsody"]])
  original <- data.table::fread(original_path)
  expected_values <- original$DBEC_Adjusted_Molecules

  res <- parse_rhapsody_matrix(config, input_dir)

  # create row and column indices, to cbind and use matrix indexing to get
  # values as a vector. Given that Simple Triplet sparse matrices basically
  # contain vectors of row and column indices, we just transform them to ints,
  # keeping order.
  row_idx <- as.integer(factor(original$Gene))
  col_idx <- match(original$Cell_Index, unique(original$Cell_Index))

  values <- res$counts_list$sample_1[cbind(row_idx, col_idx)]

  expect_equal(values, expected_values)
})

test_that("parse_rhapsody_matrix returns error if files are missing", {
  samples <- list(sample_1 = list(name = "sample_1", counts = mock_counts()))

  files <- local_rhapsody_experiment(samples)

  prev_out <- mock_prev_out(samples, include_abseq = TRUE)
  config <- prev_out$config
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

  prev_out <- mock_prev_out(samples, include_abseq = TRUE)
  config <- prev_out$config
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

  prev_out <- mock_prev_out(samples, include_abseq = TRUE)
  config <- prev_out$config
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


test_that("parse_rhapsody_matrix filters out AbSeq genes if asked to", {
  samples <- list(sample_1 = list(name = "sample_1", counts = mock_counts()))

  files <- local_rhapsody_experiment(samples, with_abseq = TRUE)
  input_dir <- "./input"

  prev_out <- mock_prev_out(samples, include_abseq = TRUE)
  config <- prev_out$config

  with_abseq <- pipeline:::parse_rhapsody_matrix(config, input_dir)

  prev_out <- mock_prev_out(samples, include_abseq = FALSE)
  config <- prev_out$config

  without_abseq <- pipeline:::parse_rhapsody_matrix(config, input_dir)
  abseq_pattern <- "(p_?ab_?o)$"

  for (sample in names(samples)) {
    expect_true(any(grepl(abseq_pattern, rownames(with_abseq$counts_list[[sample]]), ignore.case = TRUE)))
    expect_false(any(grepl(abseq_pattern, rownames(without_abseq$counts_list[[sample]]), ignore.case = TRUE)))

    expect_gte(nrow(with_abseq$counts_list[[sample]]), nrow(without_abseq$counts_list[[sample]]))
  }

  # there are abseq genes in annot
  expect_true(any(grepl(abseq_pattern, with_abseq$annot$input, ignore.case = TRUE)))
  expect_true(any(grepl(abseq_pattern, with_abseq$annot$original_name, ignore.case = TRUE)))
  expect_true(any(grepl(abseq_pattern, with_abseq$annot$name, ignore.case = TRUE)))

  # no abseq genes in annot if filtered
  expect_false(any(grepl(abseq_pattern, without_abseq$annot$input, ignore.case = TRUE)))
  expect_false(any(grepl(abseq_pattern, without_abseq$annot$original_name, ignore.case = TRUE)))
  expect_false(any(grepl(abseq_pattern, without_abseq$annot$name, ignore.case = TRUE)))



})
