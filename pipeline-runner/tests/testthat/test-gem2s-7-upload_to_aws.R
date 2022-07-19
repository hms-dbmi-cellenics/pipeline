set.seed(1)

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

mock_scdata_list = function (config) {
  prev_out <- mock_prev_out(config)
  scdata_list <- prev_out$scdata_list

  task_out <- prepare_experiment(NULL, NULL, prev_out)$output
  scdata_list <- task_out$scdata_list
}

mock_input <- function(metadata = NULL) {
  input <- list(
    name = "project name",
    sampleNames = list("a", "b", "c"),
    sampleIds = list("123abc", "123def", "123ghi"),
    metadata = metadata
  )

  return(input)
}

mock_config <- function(input) {
  config <- list(
    name = input$name,
    samples = input$sampleIds,
    metadata = input$metadata
  )

  return (config)
}

mock_prev_out <- function(config, counts = NULL) {
  samples <- config$samples

  if (is.null(counts)) {
    counts <- DropletUtils:::simCounts()
    colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  }

  eout <- DropletUtils::emptyDrops(counts)

  counts_list <- list()
  edrops <- list()
  doublet_scores <- list()

  for (sampleId in samples) {
    counts_list[[sampleId]] <- counts
    edrops[[sampleId]] <- eout
    doublet_scores[[sampleId]] <- mock_doublet_scores(counts)
  }

  # as passed to create_seurat
  prev_out <- list(
    counts_list = counts_list,
    edrops = edrops,
    doublet_scores = doublet_scores,
    annot = data.frame(name = row.names(counts), input = row.names(counts)),
    config = config
  )
  # call create_seurat to get prev_out to pass to prepare_experiment
  return(create_seurat(NULL, NULL, prev_out)$output)
}

check_metadata_cell_ids <- function(metadata_key, metadata_value, sample_keys, cell_sets) {
  keys <- sapply(cell_sets$cellSets, `[[`, "key")

  metadata_set <- cell_sets$cellSets[[which(keys == metadata_key)]]
  metadata_name <- sapply(metadata_set$children, `[[`, "name")

  sample_set <- cell_sets$cellSets[[which(keys == "sample")]]
  sample_names <- sapply(sample_set$children, `[[`, "name")

  metadata_cell_ids <- unlist(metadata_set$children[[which(metadata_name == metadata_value)]]$cellId)

  sample_cell_sets = purrr::keep(sample_set$children, \(x) x[["key"]] %in% sample_keys)
  sample_cell_ids = unlist(lapply(sample_cell_sets, `[[`, "cellIds"))

  expect_equal(metadata_cell_ids, sample_cell_ids)
}

test_that("get_cell_sets creates scratchpad and sample sets if no metadata", {
  input <- mock_input()
  config <- mock_config(input)

  scdata_list <- mock_scdata_list(config)

  cell_sets <- get_cell_sets(scdata_list, input)
  keys <- sapply(cell_sets$cellSets, `[[`, "key")

  expect_setequal(keys, c("scratchpad", "sample"))
})


test_that("get_cell_sets adds correct cell ids for each sample", {
  input <- mock_input()
  config <- mock_config(input)
  scdata_list <- mock_scdata_list(config)

  cell_sets <- get_cell_sets(scdata_list, input)
  sets_key <- sapply(cell_sets$cellSets, `[[`, "key")

  sample_sets <- cell_sets$cellSets[[which(sets_key == "sample")]]
  samples_key <- sapply(sample_sets$children, `[[`, "key")

  for (sample_id in config$samples) {
    sample_cells <- sample_sets$children[[which(samples_key == sample_id)]]$cellIds
    expected_cells <- unname(scdata_list[[sample_id]]$cells_id)

    expect_equal(sample_cells, expected_cells)
  }
})

test_that("get_cell_sets adds a single metadata column", {
  metadata <- list(Group = list("Hello", "WT2", "WT2"))
  input <- mock_input(metadata)
  config <- mock_config(input)
  scdata_list <- mock_scdata_list(config)

  cell_sets <- get_cell_sets(scdata_list, input)

  keys <- sapply(cell_sets$cellSets, `[[`, "key")
  # has sample key as one of the keys
  expect_setequal(keys, c("scratchpad", "sample", "Group"))

  # Check that each sample/metadata intersection contains the correct cell ids
  check_metadata_cell_ids("Group", "WT2", c("123def", "123ghi"), cell_sets)
  check_metadata_cell_ids("Group", "Hello", c("123abc"), cell_sets)
})

test_that("get_cell_sets uses user-supplied syntactically invalid metadata column names", {
  metadata <- list("TRUE" = list("Hello", "WT2", "WT2"))
  input <- mock_input(metadata)
  config <- mock_config(input)
  scdata_list <- mock_scdata_list(config)

  cell_sets <- get_cell_sets(scdata_list, input)

  # has sample key as one of the keys
  keys <- sapply(cell_sets$cellSets, `[[`, "key")
  expect_setequal(keys, c("scratchpad", "sample", "TRUE"))

  check_metadata_cell_ids("TRUE", "WT2", c("123def", "123ghi"), cell_sets)
  check_metadata_cell_ids("TRUE", "Hello", c("123abc"), cell_sets)
})


test_that("get_cell_sets adds two metadata columns", {
  metadata <- list(Group1 = list("Hello", "WT2", "WT2"), Group2 = list("WT", "WT", "WTA"))
  input <- mock_input(metadata)
  config <- mock_config(input)
  scdata_list <- mock_scdata_list(config)

  cell_sets <- get_cell_sets(scdata_list, input)

  # have as keys
  keys <- sapply(cell_sets$cellSets, `[[`, "key")
  expect_setequal(keys, c("scratchpad", "sample", "Group1", "Group2"))

  check_metadata_cell_ids("Group1", "Hello", c("123abc"), cell_sets)
  check_metadata_cell_ids("Group1", "WT2", c("123def", "123ghi"), cell_sets)

  check_metadata_cell_ids("Group2", "WT", c("123abc", "123def"), cell_sets)
  check_metadata_cell_ids("Group2", "WTA", c( "123ghi"), cell_sets)
})


test_that("get_cell_sets uses unique colors for each cell set", {
  metadata <- list(Group1 = list("Hello", "WT2", "WT2"), Group2 = list("WT", "WT", "WTA"))
  input <- mock_input(metadata)
  config <- mock_config(input)
  scdata_list <- mock_scdata_list(config)

  cell_sets <- get_cell_sets(scdata_list, input)

  flat_cell_sets <- unlist(cell_sets)
  colors <- flat_cell_sets[grepl("[.]color", names(flat_cell_sets))]
  colors <- unname(colors)

  expect_equal(unique(colors), colors)
})


test_that("get_cell_sets without metadata matches snapshot", {
  input <- mock_input()
  config <- mock_config(input)
  scdata_list <- mock_scdata_list(config)

  cell_sets <- get_cell_sets(scdata_list, input)
  expect_snapshot(str(cell_sets))
})


test_that("get_cell_sets with two metadata groups matches snapshot", {
  metadata <- list(Group1 = list("Hello", "WT2", "WT2"), Group2 = list("WT", "WT", "WT124"))
  input <- mock_input(metadata)
  config <- mock_config(input)
  scdata_list <- mock_scdata_list(config)

  cell_sets <- get_cell_sets(scdata_list, input)

  expect_snapshot(str(cell_sets))
})
