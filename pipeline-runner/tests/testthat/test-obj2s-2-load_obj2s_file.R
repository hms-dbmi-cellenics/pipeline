mock_scdata <- function(data_dir, obj2s_type, remove_reductions = NULL) {
  data("pbmc_small", package = 'SeuratObject')
  pbmc_small$seurat_clusters <- pbmc_small$RNA_snn_res.1

  for (red in remove_reductions) {
    pbmc_small@reductions[[red]] <- NULL
  }

  if (obj2s_type == 'sce_object') {
    pbmc_small <- Seurat::as.SingleCellExperiment(pbmc_small)

  } else if (obj2s_type == 'anndata_object') {
    pbmc_small <- Seurat::as.SingleCellExperiment(pbmc_small)

    # setup X with counts under '/X'
    SingleCellExperiment::logcounts(pbmc_small) <- NULL
    SummarizedExperiment::assayNames(pbmc_small, 'counts') <- 'X'

    # change reduction names to match AnnData convention
    red_names <- SingleCellExperiment::reducedDimNames(pbmc_small)
    if (length(red_names))
      SingleCellExperiment::reducedDimNames(pbmc_small) <- paste0('X_', red_names)

    # write to h5ad
    anndata <- reticulate::import('anndata')
    adata <- zellkonverter::SCE2AnnData(pbmc_small)

    anndata$AnnData$write_h5ad(adata, file.path(data_dir, 'adata.h5ad'))
    return(pbmc_small)
  }

  saveRDS(pbmc_small, file.path(data_dir, 'r.rds'))
  return(pbmc_small)
}

obj2s_types <- c('seurat_object', 'sce_object', 'anndata_object')


test_that("load_obj2s_file works", {
  input_dir <- tempdir()
  data_dir <- file.path(input_dir, 'pbmc_small')
  dir.create(data_dir)

  for (obj2s_type in obj2s_types) {
    mock_scdata(data_dir, obj2s_type)
    prev_out <- list(config = list(samples = 'pbmc_small', input = list(type = obj2s_type)))
    res <- load_obj2s_file(input = NULL, pipeline_config = NULL, prev_out = prev_out, input_dir = input_dir)
    expect_type(res, 'list')
    expect_s4_class(res$output$scdata, 'Seurat')
  }

  # clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("load_obj2s_file has the correct structure", {
  # setup
  input_dir <- tempdir()
  data_dir <- file.path(input_dir, 'pbmc_small')
  dir.create(data_dir)

  for (obj2s_type in obj2s_types) {

    orig_scdata <- mock_scdata(data_dir, obj2s_type)
    mock_scdata(data_dir, obj2s_type)
    prev_out <- list(config = list(samples = 'pbmc_small', input = list(type = obj2s_type)))
    res <- load_obj2s_file(input = NULL, pipeline_config = NULL, prev_out = prev_out, input_dir = input_dir)
    scdata <- res$output$scdata

    # misc has dispersions and annotations
    expect_setequal(names(scdata@misc), c('gene_dispersion', 'gene_annotations'))

    # dispersions have expected columns
    dispersions <- scdata@misc$gene_dispersion
    expect_setequal(colnames(dispersions), c('mean', 'variance', 'variance.standardized', 'SYMBOL', 'ENSEMBL'))

    # annotations have expected columns
    annot <- scdata@misc$gene_annotations
    expect_setequal(colnames(annot), c('input', 'name', 'original_name'))

    # only the default and 'pca' reduction is present
    expect_equal(c('tsne', 'pca'), names(scdata@reductions))

    if (obj2s_type == 'seurat_object') {
      # all metadata is preserved and active.ident is added
      expect_equal(colnames(scdata@meta.data), c(colnames(orig_scdata@meta.data), 'active.ident'))

    } else if (obj2s_type %in% c('sce_object', 'anndata_object')) {
      # all metadata is preserved
      expect_equal(colnames(scdata@meta.data), colnames(orig_scdata@colData))
    }
  }

  # clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("load_obj2s_file fails if the file isn't a .rds file", {
  # setup
  input_dir <- tempdir()
  data_dir <- file.path(input_dir, 'pbmc_small')
  dir.create(data_dir)


  for (obj2s_type in obj2s_types) {

    dataset_fname <- switch(
      obj2s_type,
      'seurat_object' = 'r.rds',
      'sce_object' = 'r.rds',
      'anndata_object' = 'adata.h5ad'
    )

    mock_scdata(data_dir, obj2s_type)
    # overwrite with mtcars csv
    write.csv(mtcars, file.path(data_dir, dataset_fname))

    prev_out <- list(config = list(samples = 'pbmc_small', input = list(type = obj2s_type)))
    expect_error(
      load_obj2s_file(input = NULL, pipeline_config = NULL, prev_out = prev_out, input_dir = input_dir),
      regexp = 'ERROR_OBJ2S_READ'
    )
  }

  # clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("load_obj2s_file fails if there is no RNA assay for seurat_object", {
  # setup
  input_dir <- tempdir()
  data_dir <- file.path(input_dir, 'pbmc_small')
  dir.create(data_dir)
  orig_scdata <- mock_scdata(data_dir, obj2s_type = 'seurat_object')

  # rename RNA assay
  names(orig_scdata@assays) <- 'blah'
  saveRDS(orig_scdata, file.path(data_dir, 'r.rds'))

  prev_out <- list(config = list(samples = 'pbmc_small', input = list(type = 'seurat_object')))
  expect_error(
    load_obj2s_file(input = NULL, pipeline_config = NULL, prev_out = prev_out, input_dir = input_dir),
    regexp = 'ERROR_OBJ2S_COUNTS'
    )

  # clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("load_obj2s_file does not fail if there is no seurat_clusters column in metadata", {
  # setup
  input_dir <- tempdir()
  data_dir <- file.path(input_dir, 'pbmc_small')
  dir.create(data_dir)
  orig_scdata <- mock_scdata(data_dir, 'seurat_object')

  # remove seurat_clusters
  orig_scdata$seurat_clusters <- NULL
  saveRDS(orig_scdata, file.path(data_dir, 'r.rds'))

  prev_out <- list(config = list(samples = 'pbmc_small', input = list(type = 'seurat_object')))
  expect_error(
    load_obj2s_file(input = NULL, pipeline_config = NULL, prev_out = prev_out, input_dir = input_dir),
    NA
    )

  # clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("load_obj2s_file fails if there is no reduction", {
  # setup
  input_dir <- tempdir()
  data_dir <- file.path(input_dir, 'pbmc_small')
  dir.create(data_dir)

  for (obj2s_type in obj2s_types) {
    scdata <- mock_scdata(data_dir, obj2s_type, remove_reductions = c('pca', 'tsne'))
    prev_out <- list(config = list(samples = 'pbmc_small', input = list(type = obj2s_type)))

    expect_error(
      load_obj2s_file(input = NULL, pipeline_config = NULL, prev_out = prev_out, input_dir = input_dir),
      regexp = 'ERROR_OBJ2S_REDUCTION'
      )
  }

  # clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("load_obj2s_file does not fail if there is no pca reduction", {
  # setup
  input_dir <- tempdir()
  data_dir <- file.path(input_dir, 'pbmc_small')
  dir.create(data_dir)

  for (obj2s_type in obj2s_types) {
    scdata <- mock_scdata(data_dir, obj2s_type, remove_reductions = 'pca')
    prev_out <- list(config = list(samples = 'pbmc_small', input = list(type = obj2s_type)))

    expect_error(
      load_obj2s_file(input = NULL, pipeline_config = NULL, prev_out = prev_out, input_dir = input_dir),
      NA
    )
  }

  # clean up
  unlink(data_dir, recursive = TRUE)
})


test_that("load_obj2s_file generates HVFInfo if it is not present", {
  # setup
  input_dir <- tempdir()
  data_dir <- file.path(input_dir, 'pbmc_small')
  dir.create(data_dir)
  orig_scdata <- mock_scdata(data_dir, 'seurat_object')


  # mess up meta.features (used by HVFInfo)
  colnames(orig_scdata@assays$RNA@meta.features) <- paste0('blah.', colnames(orig_scdata@assays$RNA@meta.features))
  saveRDS(orig_scdata, file.path(data_dir, 'r.rds'))

  prev_out <- list(config = list(samples = 'pbmc_small', input = list(type = 'seurat_object')))
  expect_error(
    res <- load_obj2s_file(input = NULL, pipeline_config = NULL, prev_out = prev_out, input_dir = input_dir),
    NA)

  scdata <- res$output$scdata
  disp_colnames <- c('mean', 'variance', 'variance.standardized', 'ENSEMBL', 'SYMBOL')
  expect_true(methods::is(scdata@misc$gene_dispersion, 'data.frame'))
  expect_equal(colnames(scdata@misc$gene_dispersion), disp_colnames)

  # clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("load_obj2s_file identifies and log-transforms counts stored in data assay of seurat_object", {
  # setup
  input_dir <- tempdir()
  data_dir <- file.path(input_dir, 'pbmc_small')
  dir.create(data_dir)
  orig_scdata <- mock_scdata(data_dir, 'seurat_object')

  # replace logcounts with counts
  orig_logcounts <- orig_scdata[['RNA']]$data
  orig_counts <- orig_scdata[['RNA']]$counts
  orig_scdata[['RNA']]$data <- orig_counts

  # checking setup
  expect_true(max(orig_counts) > 100)
  expect_false(max(orig_logcounts) > 100)
  expect_identical(orig_scdata[['RNA']]$data, orig_counts)

  saveRDS(orig_scdata, file.path(data_dir, 'r.rds'))
  prev_out <- list(config = list(samples = 'pbmc_small', input = list(type = 'seurat_object')))
  res <- load_obj2s_file(input = NULL, pipeline_config = NULL, prev_out = prev_out, input_dir = input_dir)
  scdata <- res$output$scdata

  # logcounts were log-transformed
  counts <- scdata[['RNA']]$counts
  expect_identical(orig_counts, counts)
  expect_identical(scdata[['RNA']]$data, Seurat::NormalizeData(orig_counts))

  # clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("load_obj2s_file calculates logcounts for sce_object if not present", {
  # setup
  input_dir <- tempdir()
  data_dir <- file.path(input_dir, 'pbmc_small')
  dir.create(data_dir)
  orig_scdata <- mock_scdata(data_dir, 'sce_object')

  # remove logcounts
  SingleCellExperiment::logcounts(orig_scdata) <- NULL

  saveRDS(orig_scdata, file.path(data_dir, 'r.rds'))
  prev_out <- list(config = list(samples = 'pbmc_small', input = list(type = 'sce_object')))
  res <- load_obj2s_file(input = NULL, pipeline_config = NULL, prev_out = prev_out, input_dir = input_dir)
  scdata <- res$output$scdata

  # logcounts are log-transformed counts
  orig_counts <- SingleCellExperiment::counts(orig_scdata)
  counts <- scdata[['RNA']]$counts
  expect_identical(orig_counts, counts)
  expect_identical(scdata[['RNA']]$data, Seurat::NormalizeData(orig_counts))

  # clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("load_obj2s_file works when default reducion for seurat_object is different than umap or tsne", {
  # setup
  input_dir <- tempdir()
  data_dir <- file.path(input_dir, 'pbmc_small')
  dir.create(data_dir)
  orig_scdata <- mock_scdata(data_dir, 'seurat_object')
  withr::defer(unlink(data_dir, recursive = TRUE))

  # simulate the condition where renaming is needed
  names(orig_scdata@reductions)[names(orig_scdata@reductions) == "tsne"] <- "tsne.ref"
  orig_scdata@reductions$tsne <- orig_scdata@reductions$tsne.ref
  saveRDS(orig_scdata, file.path(data_dir, 'r.rds'))

  prev_out <- list(config = list(samples = 'pbmc_small', input = list(type = 'seurat_object')))
  res <- load_obj2s_file(input = NULL, pipeline_config = NULL, prev_out = prev_out, input_dir = input_dir)
  scdata <- res$output$scdata

  expect_true('tsne' %in% names(scdata@reductions))
  expect_equal('tsne', SeuratObject::DefaultDimReduc(scdata))
})

test_that("load_obj2s_file works when reduction for sce_object is a string match for umap or tsne", {
  # setup
  input_dir <- tempdir()
  data_dir <- file.path(input_dir, 'pbmc_small')
  dir.create(data_dir)
  orig_scdata <- mock_scdata(data_dir, 'sce_object')
  withr::defer(unlink(data_dir, recursive = TRUE))

  # simulate the condition where renaming is needed
  SingleCellExperiment::reducedDimNames(orig_scdata) <- c('PCA', 'TSNE.REF')
  saveRDS(orig_scdata, file.path(data_dir, 'r.rds'))

  prev_out <- list(config = list(samples = 'pbmc_small', input = list(type = 'sce_object')))
  res <- load_obj2s_file(input = NULL, pipeline_config = NULL, prev_out = prev_out, input_dir = input_dir)
  scdata <- res$output$scdata

  expect_true('tsne' %in% names(scdata@reductions))
  expect_equal('tsne', SeuratObject::DefaultDimReduc(scdata))
})

test_that("update_reduction_name correctly renames dimensionality reductions", {
  # setup
  input_dir <- tempdir()
  data_dir <- file.path(input_dir, 'pbmc_small')
  dir.create(data_dir)
  scdata <- mock_scdata(data_dir, 'seurat_object')
  withr::defer(unlink(data_dir, recursive = TRUE))

  # simulate the condition where renaming is needed
  names(scdata@reductions)[names(scdata@reductions) == 'tsne'] <- 'tsne.ref'
  scdata@reductions$tsne <- scdata@reductions$tsne.ref

  red_name <- SeuratObject::DefaultDimReduc(scdata)
  red_match <- grep('umap|tsne', red_name, value = TRUE)
  is_umap <- grepl('umap', red_match)
  is_tsne <- grepl('tsne', red_match)
  new_red_name <- ifelse(is_umap, 'umap', ifelse(is_tsne, 'tsne', NA))

  updated_scdata <- update_reduction_name(scdata, red_name, new_red_name)
  expect_equal(SeuratObject::DefaultDimReduc(updated_scdata), new_red_name)
  expect_true('tsne.ori' %in% names(updated_scdata@reductions))
})

test_that("load_obj2s_file works with an Assay5 RNA for seurat_object", {
  input_dir <- tempdir()
  data_dir <- file.path(input_dir, 'pbmc_small')
  dir.create(data_dir)
  scdata <- mock_scdata(data_dir, 'seurat_object')

  # convert to Assay5
  scdata[['RNA']] <- as(scdata[['RNA']], 'Assay5')
  saveRDS(scdata, file.path(data_dir, 'r.rds'))

  prev_out <- list(config = list(samples = 'pbmc_small', input = list(type = 'seurat_object')))
  res <- load_obj2s_file(input = NULL, pipeline_config = NULL, prev_out = prev_out, input_dir = input_dir)
  expect_type(res, 'list')
  expect_s4_class(res$output$scdata, 'Seurat')
  expect_s4_class(res$output$scdata[['RNA']], 'Assay5')

  # clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("load_obj2s_file works with split Assay5 RNA for seurat_object", {
  input_dir <- tempdir()
  data_dir <- file.path(input_dir, 'pbmc_small')
  dir.create(data_dir)
  scdata <- mock_scdata(data_dir, 'seurat_object')

  # convert to Assay5 and split
  scdata[['RNA']] <- as(scdata[['RNA']], 'Assay5')
  scdata$batch <- sample(c("batchA", "batchB", "batchC"), ncol(scdata), replace = TRUE)
  scdata[["RNA"]] <- split(scdata[["RNA"]], f = scdata$batch)

  # was split
  expect_setequal(names(scdata[["RNA"]]@layers),
                  c("data.batchB", "data.batchC", "data.batchA", "scale.data", "counts.batchB", "counts.batchC", "counts.batchA"))

  saveRDS(scdata, file.path(data_dir, 'r.rds'))

  prev_out <- list(config = list(samples = 'pbmc_small', input = list(type = 'seurat_object')))
  res <- load_obj2s_file(input = NULL, pipeline_config = NULL, prev_out = prev_out, input_dir = input_dir)
  expect_type(res, 'list')
  expect_s4_class(res$output$scdata, 'Seurat')
  expect_s4_class(res$output$scdata[['RNA']], 'Assay5')

  # is joined back
  expect_setequal(names(res$output$scdata[['RNA']]@layers),
                  c('counts', 'data'))

  # clean up
  unlink(data_dir, recursive = TRUE)
})
