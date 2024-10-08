mock_counts <- function() {
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )
  pbmc_raw <- as(as.matrix(pbmc_raw), 'dgCMatrix')
  return(pbmc_raw)
}


mock_doublet_scores <- function(counts) {
  set.seed(1)
  doublet_scores <- runif(ncol(counts))
  doublet_class <- ifelse(doublet_scores < 0.8, "singlet", "doublet")

  data.frame(
    row.names = colnames(counts),
    barcodes = colnames(counts),
    doublet_class = doublet_class,
    doublet_scores = doublet_scores
  )
}

mock_prev_out <- function(samples = "sample_a", counts = NULL, prev_out_config = NULL) {
  if (is.null(counts)) {
    set.seed(1)
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

  if (!is.null(prev_out_config)) {
    prev_out$qc_config <- prev_out_config
  }

  # call create_seurat to get prev_out to pass to prepare_experiment
  create_seurat(NULL, NULL, prev_out)$output
}

mock_gem2s_input <- function(){
  input <- list(input = list(type="10X"), experimentId = "1234")
}

mock_subset_input <- function(){
  input <- list(
    experimentId = "1234",
    taskName = "prepareExperiment",
    processName = "subset",
    ignoreSslCert = FALSE,
    experimentName = "Subset of parentExperiment",
    authJWT = "Bearar somehash",
    input = list(type="10X")
  )
}

test_that("prepare_experiment ensures gene_annotations are indexed correctly for each sample", {

  samples <- c("a", "b", "c")
  prev_out <- mock_prev_out(samples = samples)

  input <- mock_gem2s_input()

  # remove some genes from each sample
  prev_out$counts_list$a <- prev_out$counts_list$a[-c(1:9), ]
  prev_out$counts_list$b <- prev_out$counts_list$b[-c(21:30), ]
  prev_out$counts_list$c <- prev_out$counts_list$c[-c(5:25), ]

  # re-create seurat object
  prev_out <- create_seurat(NULL, NULL, prev_out)$output
  scdata_list <- prepare_experiment(input, NULL, prev_out)$output$scdata

  # we expect that the input in gene_annotations is the same as the rownames of
  # each sample seurat object
  for (sample in samples) {
    expect_equal(scdata_list[[sample]]@misc$gene_annotations$input, rownames(scdata_list[[sample]]))
  }

})


test_that("add_metadata_to_samples adds 0 indexed cell_ids to each sample in scdata_list", {
  samples <- c("a", "b", "c")
  prev_out <- mock_prev_out(samples = samples)
  experiment_id <- "1234"

  scdata_list <- prev_out$scdata_list
  annot <- prev_out$annot
  scdata_list <- add_metadata_to_samples(scdata_list, annot, experiment_id)


  added_to_misc <- c("gene_annotations", "experimentId")

  for (sample in samples) {
    expect_true(all(added_to_misc %in% names(scdata_list[[sample]]@misc)))
  }

  # list of added cell_ids per sample in scdata_list
  added_ids <- purrr::map(scdata_list, ~unname(.$cells_id))

  set.seed(RANDOM_SEED)
  total_cells <- sum(sapply(scdata_list, ncol))
  cell_ids <- 0:(total_cells-1)
  start <- 0
  expected_ids <- list()
  for (sample in samples) {
    sample_size <- ncol(scdata_list[[sample]])
    idxs <- sample(seq_along(cell_ids), sample_size)
    expected_ids[[sample]] <- cell_ids[idxs]
    # remove the selected cell ids for next samples
    cell_ids <- cell_ids[-idxs]
  }
  expect_equal(added_ids, expected_ids)
})




test_that("prepare_experiment generates qc_config that matches snapshot", {
  prev_out <- mock_prev_out()
  input <- mock_gem2s_input()
  task_out <- prepare_experiment(input, NULL, prev_out)$output

  expect_snapshot(str(task_out$qc_config))
})


test_that("prepare_experiment creates a list of valid Seurat objects", {

  samples <- c("a", "b", "c")
  prev_out <- mock_prev_out(samples )
  scdata_list <- prev_out$scdata_list
  input <- mock_gem2s_input()

  task_out <- prepare_experiment(input, NULL, prev_out)$output
  scdata_list <- task_out$scdata_list


  expect_type(scdata_list, "list")
  for (sample in samples) {
    scdata <- scdata_list[[sample]]
    expect_type(scdata, "S4")
    expect_true(nrow(scdata) == 100, "Seurat")
    expect_true(scdata@active.assay == "RNA")
  }
})


test_that("prepare_experiment properly populates the misc slot", {

  samples <- c("a", "b", "c")
  prev_out <- mock_prev_out(samples )
  scdata_list <- prev_out$scdata_list

  input <- mock_gem2s_input()

  task_out <- prepare_experiment(input, NULL, prev_out)$output
  scdata_list <- task_out$scdata_list

  for (sample in samples) {

    scdata <- scdata_list[[sample]]
    expect_type(scdata@misc, "list")
    misc <- scdata@misc
    expect_true("gene_annotations" %in% names(misc))
    expect_true(all(misc$gene_annotations$input == rownames(scdata)))
    expect_equal(sum(duplicated(misc$gene_annotations$name)), 0)
    # check that in duplicated positions (including the first) we have the gene id instead of the name.
  }

})


test_that("prepare_experiment properly populates the metadata slot", {
  samples <- c("a", "b", "c")
  prev_out <- mock_prev_out(samples)
  scdata_list <- prev_out$scdata_list

  input <- mock_gem2s_input()

  task_out <- prepare_experiment(input, NULL, prev_out)$output
  scdata_list <- task_out$scdata_list

  for (sample in samples) {
    scdata <- scdata_list[[sample]]

    metadata <- scdata@meta.data
    annotations <- scdata@misc$gene_annotations

    expect_type(metadata, "list")
    expect_true(nrow(metadata) == ncol(scdata))
    expect_true("barcode" %in% colnames(metadata))
    expect_true(all(colnames(scdata) %in% rownames(metadata)))
    expect_true("samples" %in% colnames(metadata))

    if (any(grepl("^mt-", annotations$name, ignore.case = T))) {
      expect_true("percent.mt" %in% colnames(metadata))
    }

    expect_true("doublet_scores" %in% colnames(metadata))
    expect_true("cells_id" %in% colnames(metadata))
    expect_true("samples" %in% colnames(metadata))

  }

})


test_that("Mitochondrial percentage is correct", {
  samples <- c("a", "b", "c")
  prev_out <- mock_prev_out(samples)
  scdata_list <- prev_out$scdata_list

  input <- mock_gem2s_input()

  task_out <- prepare_experiment(input, NULL, prev_out)$output
  scdata_list <- task_out$scdata_list

  for (sample in samples) {
    scdata <- scdata_list[[sample]]

    metadata <- scdata@meta.data


    expect_true(max(metadata$percent.mt) <= 100)
    expect_true(min(metadata$percent.mt) >= 0)
    #Verify that we have percent mt and not fraction
    expect_true(max(metadata$percent.mt) > 1 || all(metadata$percent.mt == 0))

  }
})

test_that("Skips qc config creation if running a subset experiment (it is already created in prev_out)", {
  # If the config is already created in a previous step (we are running subset)
  # we will pass that config in prev_out

  samples <- c("a", "b", "c")
  prev_out <- mock_prev_out(samples = samples, prev_out_config = c('mocked'))

  input <- mock_subset_input()

  # re-create seurat object
  prev_out <- create_seurat(NULL, NULL, prev_out)$output
  scdata_list <- prepare_experiment(input, NULL, prev_out)

  expect_true(prev_out$qc_config == c('mocked'))
})
