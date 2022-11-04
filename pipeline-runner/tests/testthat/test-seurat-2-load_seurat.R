mock_scdata <- function(data_dir) {
  data("pbmc_small", package = 'SeuratObject')
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
    load_seurat(input = NULL, pipeline_config = NULL, prev_out = prev_out, input_dir = input_dir)
    )

  # clean up
  unlink(data_dir, recursive = TRUE)
})

