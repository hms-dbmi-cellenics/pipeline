mock_scdata <- function(data_dir) {
  data("pbmc_small", package = 'SeuratObject')
  pbmc_small$seurat_clusters <- pbmc_small$RNA_snn_res.1
  saveRDS(pbmc_small, file.path(data_dir, 'r.rds'))
  return(pbmc_small)
}


test_that("load_seurat works", {
  input_dir <- tempdir()
  data_dir <- file.path(input_dir, 'pbmc_small')
  dir.create(data_dir)
  mock_scdata(data_dir)
  prev_out <- list(config = list(samples = 'pbmc_small'))
  res <- load_seurat(input = NULL, pipeline_config = NULL, prev_out = prev_out, input_dir = input_dir)
  expect_type(res, 'list')
  expect_s4_class(res$output$scdata, 'Seurat')

  # clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("load_seurat has the correct structure", {
  # setup
  input_dir <- tempdir()
  data_dir <- file.path(input_dir, 'pbmc_small')
  dir.create(data_dir)
  orig_scdata <- mock_scdata(data_dir)
  prev_out <- list(config = list(samples = 'pbmc_small'))

  res <- load_seurat(input = NULL, pipeline_config = NULL, prev_out = prev_out, input_dir = input_dir)
  scdata <- res$output$scdata

  # misc has dispersions and annotations
  expect_setequal(names(scdata@misc), c('gene_dispersion', 'gene_annotations'))

  # dispersions have expected columns
  dispersions <- scdata@misc$gene_dispersion
  expect_setequal(colnames(dispersions), c('mean', 'variance', 'variance.standardized', 'SYMBOL', 'ENSEMBL'))

  # annotations have expected columns
  annot <- scdata@misc$gene_annotations
  expect_setequal(colnames(annot), c('input', 'name', 'original_name'))

  # only the default reduction is present
  expect_equal(SeuratObject::DefaultDimReduc(orig_scdata), names(scdata@reductions))

  # all metadata is preserved
  expect_equal(scdata@meta.data, orig_scdata@meta.data)

  # clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("load_seurat fails if the file isn't a .rds file", {
  # setup
  input_dir <- tempdir()
  data_dir <- file.path(input_dir, 'pbmc_small')
  dir.create(data_dir)
  orig_scdata <- mock_scdata(data_dir)

  # overwrite with mtcars csv
  write.csv(mtcars, file.path(data_dir, 'r.rds'))

  prev_out <- list(config = list(samples = 'pbmc_small'))
  expect_error(
    load_seurat(input = NULL, pipeline_config = NULL, prev_out = prev_out, input_dir = input_dir),
    regexp = 'ERROR_SEURAT_RDS'
  )

  # clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("load_seurat fails if there is no RNA assay", {
  # setup
  input_dir <- tempdir()
  data_dir <- file.path(input_dir, 'pbmc_small')
  dir.create(data_dir)
  orig_scdata <- mock_scdata(data_dir)

  # rename RNA assay
  names(orig_scdata@assays) <- 'blah'
  saveRDS(orig_scdata, file.path(data_dir, 'r.rds'))

  prev_out <- list(config = list(samples = 'pbmc_small'))
  expect_error(
    load_seurat(input = NULL, pipeline_config = NULL, prev_out = prev_out, input_dir = input_dir),
    regexp = 'ERROR_SEURAT_COUNTS'
    )

  # clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("load_seurat fails if there is no seurat_clusters column in metadata", {
  # setup
  input_dir <- tempdir()
  data_dir <- file.path(input_dir, 'pbmc_small')
  dir.create(data_dir)
  orig_scdata <- mock_scdata(data_dir)

  # remove seurat_clusters
  orig_scdata$seurat_clusters <- NULL
  saveRDS(orig_scdata, file.path(data_dir, 'r.rds'))

  prev_out <- list(config = list(samples = 'pbmc_small'))
  expect_error(
    load_seurat(input = NULL, pipeline_config = NULL, prev_out = prev_out, input_dir = input_dir),
    regexp = 'ERROR_SEURAT_CLUSTERS'
    )

  # clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("load_seurat fails if there is no reduction", {
  # setup
  input_dir <- tempdir()
  data_dir <- file.path(input_dir, 'pbmc_small')
  dir.create(data_dir)
  orig_scdata <- mock_scdata(data_dir)

  # remove reductions
  orig_scdata@reductions$tsne <- NULL
  orig_scdata@reductions$pca <- NULL
  saveRDS(orig_scdata, file.path(data_dir, 'r.rds'))

  prev_out <- list(config = list(samples = 'pbmc_small'))
  expect_error(
    load_seurat(input = NULL, pipeline_config = NULL, prev_out = prev_out, input_dir = input_dir),
    regexp = 'ERROR_SEURAT_REDUCTION'
    )

  # clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("load_seurat fails if there is no pca reduction", {
  # setup
  input_dir <- tempdir()
  data_dir <- file.path(input_dir, 'pbmc_small')
  dir.create(data_dir)
  orig_scdata <- mock_scdata(data_dir)

  # remove reductions
  orig_scdata@reductions$pca <- NULL
  saveRDS(orig_scdata, file.path(data_dir, 'r.rds'))

  prev_out <- list(config = list(samples = 'pbmc_small'))
  expect_error(
    load_seurat(input = NULL, pipeline_config = NULL, prev_out = prev_out, input_dir = input_dir),
    regexp = 'ERROR_SEURAT_REDUCTION'
    )

  # clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("load_seurat fails if there is inappropriate logcounts data", {
  # setup
  input_dir <- tempdir()
  data_dir <- file.path(input_dir, 'pbmc_small')
  dir.create(data_dir)
  orig_scdata <- mock_scdata(data_dir)

  # sparse matrix with wrong dimensions
  bad_data <- Matrix::sparseMatrix(1, 1, x = 1)
  expect_error(test_user_sparse_mat(bad_data), NA)

  orig_scdata@assays$RNA@data <- bad_data
  saveRDS(orig_scdata, file.path(data_dir, 'r.rds'))

  prev_out <- list(config = list(samples = 'pbmc_small'))
  expect_error(
    load_seurat(input = NULL, pipeline_config = NULL, prev_out = prev_out, input_dir = input_dir),
    regexp = 'ERROR_SEURAT_LOGCOUNTS'
    )

  # clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("load_seurat fails if there if HVFInfo can't be called", {
  # setup
  input_dir <- tempdir()
  data_dir <- file.path(input_dir, 'pbmc_small')
  dir.create(data_dir)
  orig_scdata <- mock_scdata(data_dir)


  # mess up meta.features (used by HVFInfo)
  colnames(orig_scdata@assays$RNA@meta.features) <- paste0('blah.', colnames(orig_scdata@assays$RNA@meta.features))
  saveRDS(orig_scdata, file.path(data_dir, 'r.rds'))

  prev_out <- list(config = list(samples = 'pbmc_small'))
  expect_error(
    load_seurat(input = NULL, pipeline_config = NULL, prev_out = prev_out, input_dir = input_dir),
    regexp = 'ERROR_SEURAT_HVFINFO'
    )

  # clean up
  unlink(data_dir, recursive = TRUE)
})

