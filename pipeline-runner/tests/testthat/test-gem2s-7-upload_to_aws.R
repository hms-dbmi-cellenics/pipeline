mock_config <- function(metadata = NULL) {
  config <- list(
    sampleNames = list("WT1", "WT2"),
    sampleIds = list("123abc", "123def"),
    metadata = metadata
  )

  return(config)
}


mock_scdata <- function(config) {
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )

  # construct metadata
  # add samples
  samples <- unlist(config$sampleIds)

  rest <- config$metadata
  keys <- c("samples", names(rest))
  metadata <- data.frame(row.names = colnames(pbmc_raw), samples = rep(samples, each = 40))

  # Add "metadata" if exists in config
  if (!is.null(rest)) {
    rest <- lapply(rest, unlist)
    rest <- data.frame(rest, row.names = samples, check.names = FALSE)
    metadata[names(rest)] <- rest[metadata$samples, ]
  }

  # make syntactically valid
  colnames(metadata) <- make.names(colnames(metadata), unique = TRUE)

  scdata <- Seurat::CreateSeuratObject(counts = pbmc_raw, meta.data = metadata)

  scdata$cells_id <- seq(0, ncol(scdata) - 1)
  return(scdata)
}


test_that("get_cell_sets creates scratchpad and sample sets if no metadata", {
  config <- mock_config()
  scdata <- mock_scdata(config)

  cell_sets <- get_cell_sets(scdata, config)
  keys <- sapply(cell_sets$cellSets, `[[`, "key")

  expect_setequal(keys, c("scratchpad", "sample"))
})


test_that("get_cell_sets adds correct cell ids for each sample", {
  config <- mock_config()
  scdata <- mock_scdata(config)

  cell_sets <- get_cell_sets(scdata, config)
  sets_key <- sapply(cell_sets$cellSets, `[[`, "key")

  sample_sets <- cell_sets$cellSets[[which(sets_key == "sample")]]
  samples_key <- sapply(sample_sets$children, `[[`, "key")

  for (sample_id in config$sampleIds) {
    sample_cells <- sample_sets$children[[which(samples_key == sample_id)]]$cellIds
    expected_cells <- unname(scdata$cells_id)[scdata$samples == sample_id]

    expect_equal(sample_cells, expected_cells)
  }
})


test_that("get_cell_sets adds correct cell ids for each sample", {
  config <- mock_config()
  scdata <- mock_scdata(config)

  cell_sets <- get_cell_sets(scdata, config)
  sets_key <- sapply(cell_sets$cellSets, `[[`, "key")

  sample_sets <- cell_sets$cellSets[[which(sets_key == "sample")]]
  samples_key <- sapply(sample_sets$children, `[[`, "key")

  # ids are correct for each child
  for (sample_id in config$sampleIds) {
    sample_cells <- sample_sets$children[[which(samples_key == sample_id)]]$cellIds
    expected_cells <- unname(scdata$cells_id)[scdata$samples == sample_id]

    expect_equal(sample_cells, expected_cells)
  }
})


test_that("get_cell_sets adds a single metadata column", {
  metadata <- list(Group = list("Hello", "WT2"))
  config <- mock_config(metadata)
  scdata <- mock_scdata(config)

  cell_sets <- get_cell_sets(scdata, config)

  # have it as a key
  keys <- sapply(cell_sets$cellSets, `[[`, "key")
  expect_setequal(keys, c("scratchpad", "sample", "Group"))

  group_set <- cell_sets$cellSets[[which(keys == "Group")]]
  group_names <- sapply(group_set$children, `[[`, "name")


  # cell ids are correct for each child
  for (group_name in group_names) {
    group_cells <- group_set$children[[which(group_names == group_name)]]$cellIds
    expected_cells <- unname(scdata$cells_id)[scdata$Group == group_name]

    expect_equal(group_cells, expected_cells)
  }
})


test_that("get_cell_sets uses user-supplied syntactically invalid metadata column names", {
  metadata <- list("TRUE" = list("Hello", "WT2"))
  config <- mock_config(metadata)
  scdata <- mock_scdata(config)

  cell_sets <- get_cell_sets(scdata, config)

  # have TRUE as a key
  keys <- sapply(cell_sets$cellSets, `[[`, "key")
  expect_true("TRUE" %in% keys)

  group_set <- cell_sets$cellSets[[which(keys == "TRUE")]]
  group_names <- sapply(group_set$children, `[[`, "name")

  # cell ids are correct for each child
  for (group_name in group_names) {
    group_cells <- group_set$children[[which(group_names == group_name)]]$cellIds
    expected_cells <- unname(scdata$cells_id)[scdata$TRUE. == group_name]

    expect_equal(group_cells, expected_cells)
  }
})


test_that("get_cell_sets adds two metadata columns", {
  metadata <- list(Group1 = list("Hello", "WT2"), Group2 = list("WT", "WT"))
  config <- mock_config(metadata)
  scdata <- mock_scdata(config)

  cell_sets <- get_cell_sets(scdata, config)

  # have as keys
  keys <- sapply(cell_sets$cellSets, `[[`, "key")
  expect_setequal(keys, c("scratchpad", "sample", "Group1", "Group2"))

  # check that Group2 has all cells
  group2_set <- cell_sets$cellSets[[which(keys == "Group2")]]
  group2_cells <- group2_set$children[[1]]$cellIds
  expect_equal(group2_cells, unname(scdata$cells_id))
})


test_that("get_cell_sets uses unique colors for each cell set", {
  metadata <- list(Group1 = list("Hello", "WT2"), Group2 = list("WT", "WT"))
  config <- mock_config(metadata)
  scdata <- mock_scdata(config)

  cell_sets <- get_cell_sets(scdata, config)

  flat_cell_sets <- unlist(cell_sets)
  colors <- flat_cell_sets[grepl("[.]color", names(flat_cell_sets))]
  colors <- unname(colors)

  expect_equal(unique(colors), colors)
})


test_that("get_cell_sets without metadata matches snapshot", {
  config <- mock_config()
  scdata <- mock_scdata(config)

  cell_sets <- get_cell_sets(scdata, config)
  expect_snapshot(str(cell_sets))
})


test_that("get_cell_sets with two metadata groups matches snapshot", {
  metadata <- list(Group1 = list("Hello", "WT2"), Group2 = list("WT", "WT"))

  config <- mock_config(metadata)
  scdata <- mock_scdata(config)

  scdata@misc <- list(metadata_lookups = c(Group1 = "Group1", Group2 = "Group2"))
  cell_sets <- get_cell_sets(scdata, config)

  expect_snapshot(str(cell_sets))
})
