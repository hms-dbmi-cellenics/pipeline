source("mock-gem2s-input-files.R")

# Function to generate generic gene IDs and names
mock_features <- function(num_rows) {
  # Generate generic gene IDs and names
  gene_ids <- paste0("ENSMUSG", seq(1:num_rows))
  gene_names <- paste0("Gene", seq(1:num_rows))

  # Create the mock data frame
  annotations <- data.frame(
    gene_id = gene_ids,
    gene_name = gene_names,
    genome = rep("mm38", num_rows)
  )

  return(annotations)
}

mock_counts <- function(num_genes, num_cells) {
  counts <-
    Matrix::Matrix(
      rnbinom(
        num_genes * num_cells,
        mu = 0.3,
        size = 5
      ),
      ncol = num_cells,
      byrow = TRUE
    )

  colnames(counts) <- make.unique(paste0(sample(1:96, num_cells, replace = TRUE), sample(1:96, num_cells, replace = TRUE), "_", sample(1:96, num_cells, replace = TRUE)))

  return(counts)
}

test_that("load_user_files loads a parse count matrix", {
  counts <- mock_counts(100, 100)
  features <- mock_features(100)

  experiment_dir <- "./experiment_1"
  sample <- "sample_a"

  local_experiment(counts, features, experiment_dir, sample, type = "parse")

  prev_out <- list(config = list(samples = sample, input = list(type = "parse")))
  out <- load_user_files(NULL, NULL, prev_out, experiment_dir)$output

  expect_true("counts_list" %in% names(out))
  expect_true(sample %in% names(out$counts_list))

  expect_s4_class(out$counts_list[[1]], "dgCMatrix")
})

test_that("load_user_files generates feature annotation for parse data", {
  counts <- mock_counts(100, 100)
  features <- mock_features(100)

  experiment_dir <- "./experiment_1"
  sample <- "sample_a"

  local_experiment(counts, features, experiment_dir, sample, type = "parse")

  prev_out <- list(config = list(samples = sample, input = list(type = "parse")))
  out <- load_user_files(NULL, NULL, prev_out, experiment_dir)$output

  expect_true("annot" %in% names(out))
  expect_true(
    all(c("input", "name", "original_name") %in% colnames(out$annot))
  )
})

test_that("load_user_files uses appropriate feature columns for parse data", {
  counts <- mock_counts(100, 100)
  features <- mock_features(100)

  experiment_dir <- "./experiment_1"
  sample <- "sample_a"

  local_experiment(counts, features, experiment_dir, sample, type = "parse")


  prev_out <- list(config = list(samples = sample, input = list(type = "parse")))
  out <- load_user_files(NULL, NULL, prev_out, experiment_dir)$output

  # ensembl ids are counts row names
  expect_equal(
    features$gene_id,
    row.names(out$counts_list[[1]])
  )

  # ensembl ids are in column 'input' of annot
  expect_equal(out$annot$input, features$gene_id)

  # symbols are in column 'name' of annot
  expect_equal(out$annot$name, features$gene_name)
})

test_that("load_user_files loads 10x multisample experiments", {
  features <- mock_features(100)

  experiment_dir <- "./experiment_1"
  samples <- c("sample_a", "sample_b")

  local_experiment(mock_counts(100, 100), features, experiment_dir, samples[1], type = "parse")
  local_experiment(mock_counts(100, 100), features, experiment_dir, samples[2], type = "parse")


  prev_out <- list(config = list(samples = samples, input = list(type = "parse")))
  out <- load_user_files(NULL, NULL, prev_out, experiment_dir)$output

  # loaded both
  expect_equal(names(out$counts_list), samples)

  expect_equal(out$annot$name, features$gene_name)
})

test_that("read_parse_files returns error if files missing", {
  counts <- mock_counts(100, 100)
  features <- mock_features(100)

  experiment_dir <- "./experiment_1"
  sample <- "sample_a"

  local_experiment(counts, features, experiment_dir, sample, type = "parse")

  prev_out <- list(config = list(samples = sample, input = list(type = "parse")))

  files <- c("cell_metadata.csv.gz", "DGE.mtx.gz", "all_genes.csv.gz")

  sample_dir <- file.path(experiment_dir, sample)

  # remove files one by one renaming
  for (file in files) {
    file.rename(file.path(sample_dir, file), file.path(sample_dir, "blah"))
    expect_error(supressWarnings(load_user_files(NULL, NULL, prev_out, experiment_dir), "Cannot find"))
    file.rename(file.path(sample_dir, "blah"), file.path(sample_dir, file))
  }
})

test_that("read_parse_annotations makes features unique", {
  counts <- mock_counts(100, 100)
  features <- mock_features(100)

  features$gene_id[[2]] <- features$gene_id[[1]]

  experiment_dir <- "./experiment_1"
  sample <- "sample_a"

  local_experiment(counts, features, experiment_dir, sample, type = "parse")

  prev_out <- list(config = list(samples = sample, input = list(type = "parse")))

  out <- load_user_files(NULL, NULL, prev_out, experiment_dir)$output

  expect_true(out$annot$input[[2]] != out$annot$input[[1]] & rownames(out$counts_list$sample_a)[[2]] != rownames(out$counts_list$sample_a)[[1]])
})
