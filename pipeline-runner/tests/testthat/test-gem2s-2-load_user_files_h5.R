local_h5_experiment <-
  function(experiment_dir, sample_dir, env = parent.frame()) {
    sample_path <- file.path(experiment_dir, sample_dir)
    dir.create(sample_path, recursive = T)
    mock_10x_h5_matrix(5, 10, sample_path)
    R.utils::gzip(paste0(sample_path, "/matrix.h5"),paste0(sample_path, "/matrix.h5.gz"))
    withr::defer(unlink(experiment_dir, recursive = T), envir = env)

  }

mock_10x_h5_matrix <- function(num_cells, num_genes, sample_path) {
  # Generate example data in sparse matrix representation
  data_matrix <-
    matrix(rpois(num_genes * num_cells, lambda = 5),
           nrow = num_genes,
           ncol = num_cells)
  sparse_matrix <- Matrix::Matrix(data_matrix, sparse = TRUE)
  indices_matrix <- sparse_matrix@i
  indptr_matrix <- sparse_matrix@p
  shape <- sparse_matrix@Dim


  # Generate example gene names and barcodes
  gene_names <- paste0("Gene", 1:num_genes)
  gene_ids <- paste0("ENS", 1:num_genes)
  barcodes <- paste0("Cell", 1:num_cells)
  feature_types <- rep("Gene Expression", num_genes)

  # Define file path and dataset names
  file_path <- paste0(sample_path, "/matrix.h5")
  features_dataset_name <- "features/"
  main_slot <- "matrix/"

  rhdf5::h5createFile(file_path)

  # Create main groups
  rhdf5::h5createGroup(file_path, main_slot)
  rhdf5::h5createGroup(file_path, paste0(main_slot, features_dataset_name))

  # Create data, indices, and indptr datasets
  rhdf5::h5write(sparse_matrix@x, file_path, paste0(main_slot, "data"))
  rhdf5::h5write(indices_matrix, file_path, paste0(main_slot, "indices"))
  rhdf5::h5write(indptr_matrix, file_path, paste0(main_slot, "indptr"))

  # Create shape dataset
  rhdf5::h5write(shape, file_path, paste0(main_slot, "shape"))

  # Create gene_names dataset and write the values
  rhdf5::h5write(gene_names, file_path, paste0(main_slot, "gene_names"))

  # Create features dataset and write the values
  rhdf5::h5write(gene_names,
                 file_path,
                 paste0(main_slot, features_dataset_name, "name"))
  rhdf5::h5write(gene_ids,
                 file_path,
                 paste0(main_slot, features_dataset_name, "id"))
  rhdf5::h5write(
    feature_types,
    file_path,
    paste0(main_slot, features_dataset_name, "feature_type")
  )

  # Create barcodes dataset and write the values
  rhdf5::h5write(barcodes, file_path, paste0(main_slot, "barcodes"))
}


test_that("load_user_files loads an h5 matrix", {
  experiment_dir <- "./experiment_1"
  sample <- "sample_a"

  local_h5_experiment(experiment_dir, sample)

  prev_out <- list(config = list(samples = sample, input = list
                                 (type = "10x_h5")))
  out <- load_user_files(NULL, NULL, prev_out, experiment_dir)$output

  expect_true("counts_list" %in% names(out))
  expect_true(sample %in% names(out$counts_list))

  expect_s4_class(out$counts_list[[1]], "dgCMatrix")
})
