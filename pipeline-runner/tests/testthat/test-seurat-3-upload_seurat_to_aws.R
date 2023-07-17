mock_scdata <- function() {
  data("pbmc_small", package = 'SeuratObject')
  return(pbmc_small)
}

mock_put_object_in_s3 <- function(pipeline_config, bucket, object, key) {
  return(NULL)
}

mock_upload_matrix_to_s3 <- function(pipeline_config, experiment_id, data) {
  return('object_key')
}

test_that("upload_seurat_to_aws completes successfully", {

  mockery::stub(upload_seurat_to_aws, 'put_object_in_s3', mock_put_object_in_s3)
  mockery::stub(upload_seurat_to_aws, 'upload_matrix_to_s3', mock_upload_matrix_to_s3)

  scdata <- mock_scdata()
  samples <- rep(c('A', 'B', 'C', 'D'), each = ncol(scdata)/4)
  scdata$samples <- samples
  scdata$seurat_clusters <- rep(letters[1:8], length.out = ncol(scdata))

  input <- list(experimentId = '1234')
  prev_out <- list(scdata = scdata, config = list())

  expect_error(upload_seurat_to_aws(input, NULL, prev_out), NA)
})


test_that("mock_metadata can be used to add columns that are missing", {

  scdata <- mock_scdata()
  metadata_cols <- list('percent.mt' = 0, 'doublet_scores' = 0, 'doublet_class' = 'singlet')

  expect_true(all(!names(metadata_cols) %in% colnames(scdata@meta.data)))

  scdata <- mock_metadata(scdata, metadata_cols)
  expect_true(all(names(metadata_cols) %in% colnames(scdata@meta.data)))
  expect_true(all(scdata$percent.mt == 0))
  expect_true(all(scdata$doublet_scores == 0))
  expect_true(all(scdata$doublet_class == 'singlet'))
})

test_that("mock_metadata skips over columns that are present", {

  scdata <- mock_scdata()
  metadata_cols <- list('percent.mt' = 0, 'doublet_scores' = 0, 'doublet_class' = 'singlet')

  expect_true(all(!names(metadata_cols) %in% colnames(scdata@meta.data)))
  scdata$percent.mt <- 20

  scdata <- mock_metadata(scdata, metadata_cols)
  expect_true(all(names(metadata_cols) %in% colnames(scdata@meta.data)))
  expect_true(all(scdata$percent.mt == 20))
  expect_true(all(scdata$doublet_scores == 0))
  expect_true(all(scdata$doublet_class == 'singlet'))
})


test_that("add_samples_col labels all cells as one sample if no sample or samples column", {

  scdata <- mock_scdata()

  expect_null(scdata@meta.data$sample)
  expect_null(scdata@meta.data$samples)

  scdata <- add_samples_col(scdata)
  expect_true(all(scdata$samples == 'NA'))
})

test_that("add_samples_col uses existing 'samples' or 'sample' metadata column", {

  scdata <- mock_scdata()
  samples1 <- rep(c('A', 'B'), length.out = ncol(scdata))

  # will add 'sample' if present
  scdata$sample <- samples1
  scdata <- add_samples_col(scdata)
  expect_equal(scdata@meta.data$samples, samples1)

  # will leave 'samples' if present
  scdata <- mock_scdata()
  scdata$samples <- samples1
  scdata <- add_samples_col(scdata)
  expect_equal(scdata@meta.data$samples, samples1)

  # will prefer 'samples' or 'sample'
  scdata <- mock_scdata()
  scdata$sample <- samples1
  samples2 <- rep(c('C', 'D'), length.out = ncol(scdata))
  scdata$samples <- samples2
  scdata <- add_samples_col(scdata)
  expect_equal(scdata@meta.data$samples, samples2)
})


test_that("format_seurat adds requires metadata to a SeuratObject", {

  scdata <- mock_scdata()
  scdata <- format_seurat(scdata, '1234')

  # added samples
  expect_true(all(scdata$samples == 'NA'))

  # added 0-indexed cells_id
  expect_true(all(scdata$cells_id == seq(0, ncol(scdata)-1)))

  # added misc
  expect_equal(scdata@misc$experimentId, '1234')
  expect_setequal(names(scdata@misc), c('experimentId', 'color_pool', 'ingestionDate'))

  # added required metadata columns
  expect_true(all(c('percent.mt', 'doublet_scores', 'doublet_class') %in% colnames(scdata@meta.data)))
})


test_that("find_group_columns finds all group columns that are superset of sample identities", {

  scdata <- mock_scdata()
  samples <- rep(c('A', 'B', 'C', 'D'), each = ncol(scdata)/4)
  scdata$samples <- samples

  # ignores column that is just duplicate of samples
  scdata$samples_copy <- as.numeric(as.factor(samples))
  expect_length(find_group_columns(scdata), 0)

  # finds column that has two samples per group
  scdata$group1 <- ifelse(samples %in% c('A', 'B'), 'AB', 'CD')
  expect_equal(find_group_columns(scdata), 'group1')

  # finds column that has three samples in one group
  scdata$group2 <- ifelse(samples %in% c('A', 'B', 'C'), 'ABC', 'D')
  expect_equal(find_group_columns(scdata), c('group1', 'group2'))

  # ignores column that has all samples together (trivial)
  scdata$group3 <- 'ABCD'
  expect_equal(find_group_columns(scdata), c('group1', 'group2'))

  # ignores columns that divide a sample between two+ groups
  scdata$group4 <- rep(LETTERS[1:8], each = ncol(scdata)/8)
  expect_equal(find_group_columns(scdata), c('group1', 'group2'))
})


test_that("add_metadata_to_input adds group metadata to input list", {

  scdata <- mock_scdata()
  samples <- rep(c('A', 'B', 'C', 'D'), each = ncol(scdata)/4)
  scdata$samples <- samples

  # add column that has two samples per group
  scdata$group1 <- ifelse(samples %in% c('A', 'B'), 'AB', 'CD')
  # add column that has three samples in one group
  scdata$group2 <- ifelse(samples %in% c('A', 'B', 'C'), 'ABC', 'D')

  input <- add_metadata_to_input(scdata, input = list())

  expect_equal(input$metadata$group1, c('AB', 'AB', 'CD', 'CD'))
  expect_equal(input$metadata$group2, c('ABC', 'ABC', 'ABC', 'D'))
})


test_that("add_samples_to_input adds sampleNames and sampleIds to input", {

  scdata <- mock_scdata()
  sample_names <- c('A', 'B', 'C', 'D')
  scdata$samples <- rep(sample_names, each = ncol(scdata)/4)

  input <- add_samples_to_input(scdata, input = list())

  expect_equal(input$sampleNames, sample_names)
  expect_length(input$sampleIds, 4)
})


test_that("change_sample_names_to_ids substitutes samples column for values in sampleIds input using order based on sampleNames input", {

  scdata <- mock_scdata()
  sample_names <- c('A', 'B', 'C', 'D')
  samples <- rep(sample_names, each = ncol(scdata)/4)
  scdata$samples <- samples

  input <- list(sampleNames = sample_names, sampleIds = letters[1:4])
  scdata <- change_sample_names_to_ids(scdata, input)
  expect_equal(scdata@meta.data$samples, tolower(samples))

  # uses reversed order from input$sampleNames
  scdata$samples <- samples
  input$sampleNames <- rev(input$sampleNames)
  scdata <- change_sample_names_to_ids(scdata, input)
  expect_equal(scdata@meta.data$samples, rev(tolower(samples)))
})
