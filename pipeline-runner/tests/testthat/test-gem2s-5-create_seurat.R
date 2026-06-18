mock_doublet_scores <- function(counts) {
  doublet_scores <- runif(ncol(counts))
  doublet_class <- ifelse(doublet_scores < 0.8, "singlet", "doublet")

  data.frame(
    row.names = colnames(counts),
    barcodes = colnames(counts),
    doublet_class = doublet_class,
    doublet_scores = doublet_scores
  )
}

mock_prev_out <- function(samples = "sample_a", counts = NULL, use_bpcells = FALSE) {
  if (is.null(counts)) {
    counts <- DropletUtils:::simCounts()
    colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  }

  counts <- as(counts, "dgCMatrix")
  eout <- DropletUtils::emptyDrops(counts)

  counts_list <- list()
  edrops <- list()
  doublet_scores <- list()

  for (sample in samples) {
    counts_list[[sample]] <- maybe_bpcells(
      counts,
      withr::local_tempfile(.local_envir = parent.frame()),
      use_bpcells
    )
    edrops[[sample]] <- eout
    doublet_scores[[sample]] <- mock_doublet_scores(counts)
  }

  list(
    counts_list = counts_list,
    edrops = edrops,
    doublet_scores = doublet_scores,
    annot = data.frame(name = row.names(counts), input = row.names(counts)),
    config = list(name = "project name")
  )
}



test_that("construct_metadata works with syntactically invalid column names", {
  sample <- "hello"
  counts <- mock_counts()

  config <- list(
    samples = "hello",
    metadata = list(
      "TRUE" = list("a"),
      "1first" = list("b"),
      "with space" = list("c"),
      "with-dash" = list("d")
    )
  )

  metadata <- construct_metadata(counts, sample, config)
  valid_names <- make.names(colnames(metadata), unique = TRUE)

  # changes column names of metadata
  expect_equal(colnames(metadata), valid_names)

  # stores correct values in metadata
  expect_true(all(metadata$TRUE. == "a"))
  expect_true(all(metadata$X1first == "b"))
  expect_true(all(metadata$with.space == "c"))
  expect_true(all(metadata$with.dash == "d"))
})


test_that("create_seurat works with bpcells", {
  prev_out <- mock_prev_out(use_bpcells = TRUE)
  out <- expect_no_error(
    create_seurat(NULL, NULL, prev_out)$output
  )
  scdata <- out$scdata_list[[1]]

  expect_s4_class(scdata, "Seurat")
})

test_that("create_seurat persists the declared technology onto the object", {
  # lets the worker (which has no pipeline config) dispatch on technology
  prev_out <- mock_prev_out()
  prev_out$config$input <- list(type = "xenium")

  out <- create_seurat(NULL, NULL, prev_out)$output
  scdata <- out$scdata_list[[1]]

  expect_equal(scdata@misc$technology, "xenium")
})


test_that("create_seurat works without emptyDrops result", {
  prev_out <- mock_prev_out()
  prev_out$edrops$sample_a <- NULL

  expect_message(create_seurat(NULL, NULL, prev_out), "^emptyDrops .+ skipping")
})

test_that("create_seurat fails with missing items in prev_out", {
  prev_out <- mock_prev_out()

  for (item in names(prev_out)) {
    obj <- prev_out[[item]]

    prev_out[[item]] <- NULL
    expect_error(create_seurat(NULL, NULL, prev_out))

    prev_out[[item]] <- obj
  }
})

test_that("create_seurat adds mitochondrial percentage", {

  # without mitcondrial genes - gets set to 0
  prev_out <- mock_prev_out()
  out <- create_seurat(NULL, NULL, prev_out)$output
  scdata <- out$scdata_list[[1]]

  expect_true(all(scdata$percent.mt == 0))

  # with mitochondrial genes - some not zero
  counts <- DropletUtils:::simCounts()
  row.names(counts)[1:10] <- paste0("mt-", 1:10)
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  prev_out <- mock_prev_out(counts = counts)

  out <- create_seurat(NULL, NULL, prev_out)$output
  scdata <- out$scdata_list[[1]]

  expect_false(all(scdata$percent.mt == 0))

  # case insensitive
  row.names(counts)[1:10] <- paste0("MT-", 1:10)
  prev_out <- mock_prev_out(counts = counts)

  out <- create_seurat(NULL, NULL, prev_out)$output
  scdata <- out$scdata_list[[1]]

  expect_false(all(scdata$percent.mt == 0))
})

test_that("create_seurat adds mitochondrial percentage with colons", {

    # without mitcondrial genes - gets set to 0
    prev_out <- mock_prev_out()
    out <- create_seurat(NULL, NULL, prev_out)$output
    scdata <- out$scdata_list[[1]]

    expect_true(all(scdata$percent.mt == 0))

    # with mitochondrial genes - some not zero
    counts <- DropletUtils:::simCounts()
    row.names(counts)[1:10] <- paste0("mt:", 1:10)
    colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
    prev_out <- mock_prev_out(counts = counts)

    out <- create_seurat(NULL, NULL, prev_out)$output
    scdata <- out$scdata_list[[1]]

    expect_false(all(scdata$percent.mt == 0))

    # case insensitive
    row.names(counts)[1:10] <- paste0("MT:", 1:10)
    prev_out <- mock_prev_out(counts = counts)

    out <- create_seurat(NULL, NULL, prev_out)$output
    scdata <- out$scdata_list[[1]]

    expect_false(all(scdata$percent.mt == 0))
})

test_that("create_seurat adds emptyDrops_FDR to SeuratObject", {
  prev_out <- mock_prev_out()
  out <- create_seurat(NULL, NULL, prev_out)$output
  scdata <- out$scdata_list[[1]]

  expect_true("emptyDrops_FDR" %in% colnames(scdata@meta.data))
})

test_that("create_seurat adds doublet_scores and doublet_class to SeuratObject", {
  prev_out <- mock_prev_out()
  out <- create_seurat(NULL, NULL, prev_out)$output
  scdata <- out$scdata_list[[1]]

  expect_true(all(c("doublet_scores", "doublet_class") %in% colnames(scdata@meta.data)))
})

test_that("create_seurat works with multiple samples", {
  prev_out <- mock_prev_out(samples = c("a", "b"))
  scdata_list <- create_seurat(NULL, NULL, prev_out)$output$scdata_list

  # scdata_list names are right
  expect_true(all(c("a", "b") %in% names(scdata_list)))

  # samples column added to SeuratObjects
  expect_true(all(scdata_list[["a"]]$samples == "a"))
  expect_true(all(scdata_list[["b"]]$samples == "b"))
})


test_that("add_segmentations returns scdata unchanged when segmentations is NULL", {
  scdata <- mock_spatial_scdata(ncells = 36, sample_id = "s1")
  out <- add_segmentations(scdata, NULL, "s1")
  expect_identical(out, scdata)
})

test_that("add_spatial_local_outliers writes <metric>_z and <metric>_log columns", {
  scdata <- mock_spatial_scdata(ncells = 64, sample_id = "s1")
  out <- add_spatial_local_outliers(scdata)

  md_cols <- colnames(out@meta.data)
  # z-score columns for all three spatial metrics
  expect_true(all(c("nCount_RNA_z", "nFeature_RNA_z", "percent.mt_z") %in% md_cols))
  # log columns only for the log-scaled metrics (UMI, num genes), not percent.mt
  expect_true(all(c("nCount_RNA_log", "nFeature_RNA_log") %in% md_cols))
  expect_false("percent.mt_log" %in% md_cols)

  # one z-score per cell, all finite, row-aligned
  expect_equal(nrow(out@meta.data), ncol(out))
  expect_true(all(is.finite(out@meta.data$nCount_RNA_z)))
  # log column matches log1p of the metric
  expect_equal(out@meta.data$nCount_RNA_log, log1p(out@meta.data$nCount_RNA))
})

test_that("add_spatial_local_outliers returns non-spatial scdata unchanged", {
  prev_out <- mock_prev_out()
  scdata <- create_seurat(NULL, NULL, prev_out)$output$scdata_list[[1]]
  before <- colnames(scdata@meta.data)

  out <- add_spatial_local_outliers(scdata)
  # no tissue coordinates -> no z columns added
  expect_identical(out, scdata)
  expect_false(any(grepl("_z$", colnames(out@meta.data))))
  expect_equal(colnames(out@meta.data), before)
})

test_that("create_seurat does not exclude genes without counts", {
  counts <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )
  counts["NOT-EXPRESSED", ] = 0
  counts <- as(as.matrix(counts), "dgCMatrix")

  prev_out <- mock_prev_out(counts = counts)

  # check that have genes with 0 counts
  counts <- prev_out$counts_list[[1]]
  counts_per_gene <- Matrix::rowSums(counts)
  expect_equal(sum(counts_per_gene == 0), 1)

  out <- create_seurat(NULL, NULL, prev_out)$output
  scdata <- out$scdata_list[[1]]

  # gene is still there
  expect_true("NOT-EXPRESSED" %in% row.names(scdata))
})
