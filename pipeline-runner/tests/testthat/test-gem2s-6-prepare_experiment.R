mock_counts <- function() {
  read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )
}

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

mock_prev_out <- function(samples = "sample_a", counts = NULL) {
  if (is.null(counts)) {
    counts <- DropletUtils:::simCounts()
    colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  }

  eout <- DropletUtils::emptyDrops(counts)

  counts_list <- list()
  edrops <- list()
  doublet_scores <- list()

  for (sample in samples) {
    counts_list[[sample]] <- counts
    edrops[[sample]] <- eout
    doublet_scores[[sample]] <- mock_doublet_scores(counts)
  }

  # as passed to create_seurat
  prev_out <- list(
    counts_list = counts_list,
    edrops = edrops,
    doublet_scores = doublet_scores,
    annot = data.frame(name = row.names(counts), input = row.names(counts)),
    config = list(name = "project name")
  )

  # call create_seurat to get prev_out to pass to prepare_experiment
  create_seurat(NULL, NULL, prev_out)$output
}



test_that("prepare_experiment merges multiple SeuratObjects", {
  prev_out <- mock_prev_out(samples = c("a", "b", "c"))
  scdata_list <- prev_out$scdata_list

  task_out <- expect_warning(prepare_experiment(NULL, NULL, prev_out)$output)

  scdata <- task_out$scdata

  expect_equal(ncol(scdata), sum(sapply(scdata_list, ncol)))
})

test_that("prepare_experiment shuffles cells after merge", {
  prev_out <- mock_prev_out(samples = c("a", "b", "c"))
  scdata_list <- prev_out$scdata_list

  task_out <- expect_warning(prepare_experiment(NULL, NULL, prev_out)$output)

  scdata <- task_out$scdata
  merged_scdatas <- expect_warning(merge_scdatas(scdata_list))

  set.seed(pipeline::gem2s$random.seed)
  shuffle_mask <- sample(colnames(merged_scdatas))

  expect_equal(ncol(scdata), ncol(merged_scdatas))
  expect_equal(merged_scdatas$samples[shuffle_mask],scdata$samples)
  expect_true(all(shuffle_mask == colnames(scdata)))
  expect_false(all(shuffle_mask == colnames(merged_scdatas)))
  expect_true(all(merged_scdatas$samples[1:ncol(scdata_list[[1]])]=="a"))
  expect_false(all(scdata$samples[1:ncol(scdata_list[[1]])]=="a"))
})

test_that("prepare_experiment ensures gene_annotations are indexed the same as scdata", {
  prev_out <- mock_prev_out()

  # shuffle gene order of annot
  annot <- prev_out$annot
  prev_out$annot <- annot[sample(nrow(annot)), ]

  scdata <- prepare_experiment(NULL, NULL, prev_out)$output$scdata

  expect_equal(row.names(scdata), scdata@misc$gene_annotations$input)
})

test_that("prepare_experiment adds 0 indexed cell_ids and other metadata to scdata", {
  prev_out <- mock_prev_out()
  input <- list(experimentId = "1234")
  scdata <- prepare_experiment(input, NULL, prev_out)$output$scdata

  added_to_misc <- c("gene_annotations", "color_pool", "experimentId", "ingestionDate")
  expect_true(all(added_to_misc %in% names(scdata@misc)))

  added_ids <- unname(scdata$cells_id)
  expected_ids <- seq(0, ncol(scdata) - 1)
  expect_equal(added_ids, expected_ids)
})

test_that("prepare_experiment generates qc_config that matches snapshot", {
  prev_out <- mock_prev_out()
  input <- list(experimentId = "1234")
  task_out <- prepare_experiment(input, NULL, prev_out)$output

  expect_snapshot(str(task_out$qc_config))
})

test_object <- function() {
  prev_out <- mock_prev_out(samples = c("a", "b", "c"))
  scdata_list <- prev_out$scdata_list

  task_out <- expect_warning(prepare_experiment(NULL, NULL, prev_out)$output)

  scdata <- task_out$scdata

  test_that("prepare_experiment creates a valid Seurat object", {
    expect_type(scdata, "S4")
    expect_true(nrow(scdata) == 100, "Seurat")
    expect_true(scdata@active.assay == "RNA")
  })

  test_that("prepare_experiment properly populates the misc slot", {
    expect_type(scdata@misc, "list")
    misc <- scdata@misc
    expect_true("gene_annotations" %in% names(misc))
    expect_true("color_pool" %in% names(misc))
    expect_type(misc$color_pool, "character")
    expect_true(all(misc$gene_annotations$input == rownames(scdata)))
    expect_equal(sum(duplicated(misc$gene_annotations$name)), 0)
    # check that in duplicated positions (including the first) we have the gene id instead of the name.
  })

  test_that("prepare_experiment properly populates the metadata slot", {
    md <- scdata@meta.data
    annotations <- scdata@misc$gene_annotations

    expect_type(md, "list")
    expect_true(nrow(md) == ncol(scdata))
    expect_true("barcode" %in% colnames(md))
    expect_true(all(colnames(scdata) %in% rownames(md)))
    expect_true("samples" %in% colnames(md))

    if (any(grepl("^mt-", annotations$name, ignore.case = T))) {
      expect_true("percent.mt" %in% colnames(md))
    }

    expect_true("doublet_scores" %in% colnames(md))
    expect_true("cells_id" %in% colnames(md))
    expect_true("samples" %in% colnames(md))

    test_that("Cell ids are assigned correctly", {
      cellNumber <- ncol(scdata@assays$RNA@data)
      expect_equal(md$cells_id, 0:(cellNumber - 1))
    })

    test_that("Mitochondrial percentage is correct", {
      expect_true(max(md$percent.mt) <= 100)
      expect_true(min(md$percent.mt) >= 0)
      #Verify that we have percent mt and not fraction
      expect_true(max(md$percent.mt) > 1 || all(md$percent.mt == 0))
    })
  })
}

test_that("prepare_experiment creates a valid Seurat object", {
  test_object()
})
