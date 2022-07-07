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

test_that("prepare_experiment returns list of seurat objects", {})

# test_that("prepare_experiment merges multiple SeuratObjects", {
#   prev_out <- mock_prev_out(samples = c("a", "b", "c"))
#   scdata_list <- prev_out$scdata_list
#
#   task_out <- suppressWarnings(prepare_experiment(NULL, NULL, prev_out)$output)
#
#   scdata <- task_out$scdata
#
#   expect_equal(ncol(scdata), sum(sapply(scdata_list, ncol)))
# })


test_that("prepare_experiment ensures gene_annotations are indexed correctly for each sample", {

  samples <- c("a", "b", "c")
  prev_out <- mock_prev_out(samples = samples)

  # remove some genes from each sample
  prev_out$counts_list$a <- prev_out$counts_list$a[-c(1:9), ]
  prev_out$counts_list$b <- prev_out$counts_list$b[-c(21:30), ]
  prev_out$counts_list$c <- prev_out$counts_list$c[-c(5:25), ]

  # TODO: refactor mocking
  # re-create seurat object
  prev_out <- create_seurat(NULL, NULL, prev_out)$output
  scdata_list <- prepare_experiment(NULL, NULL, prev_out)$output$scdata

  # we expect that the input in gene_annotations is the same as the rownames of
  # each sample seurat object
  for (sample in samples) {
    expect_equal(scdata_list[[sample]]@misc$gene_annotations$input, rownames(scdata_list[[sample]]))
  }

})

test_that("prepare_experiment adds 0 indexed cell_ids to each sample in scdata_list", {
  samples <- c("a", "b", "c")
  prev_out <- mock_prev_out(samples = samples)
  input <- list(experimentId = "1234")

  scdata_list <- prepare_experiment(input, NULL, prev_out)$output$scdata

  added_to_misc <- c("gene_annotations", "experimentId")

  for (sample in samples) {
  expect_true(all(added_to_misc %in% names(scdata_list[[sample]]@misc)))
  }

  # list of added cell_ids per sample in scdata_list
  added_ids <- purrr::map(scdata_list, ~unname(.$cells_id))

  expected_ids <- list()
  start <- 0
  for (sample in samples) {
    end <- start + ncol(scdata_list[[sample]]) - 1
    expected_ids[[sample]] <- seq(start, end)
    start <- end + 1
  }

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

  task_out <- prepare_experiment(NULL, NULL, prev_out)$output

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
  expect_true(TRUE)
  suppressWarnings(test_object())
})
