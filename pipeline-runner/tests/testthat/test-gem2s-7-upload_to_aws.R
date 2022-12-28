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

mock_scdata_list <- function(config) {
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
    metadata = metadata,
    experimentId = "mock_experiment_id",
    projectId = "mock_experiment_id"
  )

  return(input)
}

mock_config <- function(input) {
  config <- list(
    name = input$name,
    samples = input$sampleIds,
    metadata = input$metadata
  )

  return(config)
}


mock_parsed_cellsets <- function(scdata_list) {
  samples <- unlist(lapply(scdata_list, function(x) {
    x$samples
  }))
  cells_ids_by_sample <- unlist(lapply(scdata_list, function(x) {
    x$cells_id
  }))
  sample_id <- list()
  for (i in 1:length(table(samples))) {
    sample_id[[i]] <- rep(paste(paste(rep(i, 8), collapse = ""),
      paste(rep(i, 4), collapse = ""),
      paste(rep(i, 4), collapse = ""),
      paste(rep(i, 4), collapse = ""),
      paste(rep(i, 12), collapse = ""),
      sep = "-"
    ), table(samples)[i])
  }
  key <- c(
    rep("louvain-0", length(samples) / 2), rep("louvain-1", length(samples) / 2),
    rep("074f2946-be51-445c-8222-be953f7da981", length(samples) / 3), rep("73e4d06e-b6b6-42a9-90ae-16586f41b7e9", 2 * length(samples) / 3),
    unlist(sample_id),
    rep("Metadata-1", length(samples))
  )
  name <- c(
    rep("Cluster 0", length(samples) / 2), rep("Cluster 1", length(samples) / 2),
    rep("Custom-A", length(samples) / 3), rep("Custom-B", 2 * length(samples) / 3),
    samples,
    rep("Metadata-A", 2 * length(samples) / 3), rep("Metadata-B", length(samples) / 3)
  )
  ordered_cell_ids <- sort(cells_ids_by_sample)
  cell_id <- c(
    ordered_cell_ids,
    ordered_cell_ids[1:(length(ordered_cell_ids) / 3)], ordered_cell_ids[(length(ordered_cell_ids) / 3 + 1):(length(ordered_cell_ids))],
    cells_ids_by_sample,
    ordered_cell_ids[(length(ordered_cell_ids) / 3 + 1):(length(ordered_cell_ids))], ordered_cell_ids[1:(length(ordered_cell_ids) / 3)]
  )
  type <- c(
    rep("cluster", length(cells_ids_by_sample)), rep("scratchpad", length(cells_ids_by_sample)),
    rep("sample", length(cells_ids_by_sample)), rep("metadata", length(cells_ids_by_sample))
  )

  parsed_cellsets <- data.table(key, "name" = name, "cell_id" = cell_id, "type" = type)

  return(parsed_cellsets)
}


mock_prev_out <- function(config, counts = NULL) {
  samples <- config$samples

  if (is.null(counts)) {
    set.seed(1)
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

mock_subset_data <- function(scdata_list, cell_ids) {
  subset_scdata <- list()
  sample_ids <- c("new-123abc", "new-123def", "new-123ghi")
  for (i in 1:length(scdata_list)) {

    sample_scdata <- subset_ids(scdata_list[[i]], cell_ids)
    sample_scdata@meta.data$parent_samples <- sample_scdata@meta.data$samples
    sample_scdata@meta.data$samples <- sample_ids[i]

    subset_scdata[[i]] <- sample_scdata
  }

  names(subset_scdata) <- sample_ids
  return(subset_scdata)
}

mock_sample_id_map <- function() {
  sample_id_map <- list("11111111-1111-1111-1111-111111111111" = "new-123abc",
                        "22222222-2222-2222-2222-222222222222" = "new-123def",
                        "33333333-3333-3333-3333-333333333333" = "new-123ghi")
  return(sample_id_map)
}

check_metadata_cell_ids <- function(metadata_key, metadata_value, sample_keys, cell_sets) {
  keys <- sapply(cell_sets$cellSets, `[[`, "key")

  metadata_set <- cell_sets$cellSets[[which(keys == metadata_key)]]
  metadata_name <- sapply(metadata_set$children, `[[`, "name")

  sample_set <- cell_sets$cellSets[[which(keys == "sample")]]
  sample_names <- sapply(sample_set$children, `[[`, "name")

  metadata_cell_ids <- unlist(metadata_set$children[[which(metadata_name == metadata_value)]]$cellId)

  sample_cell_sets <- purrr::keep(sample_set$children, \(x) x[["key"]] %in% sample_keys)
  sample_cell_ids <- unlist(lapply(sample_cell_sets, `[[`, "cellIds"))

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
  check_metadata_cell_ids("Group2", "WTA", c("123ghi"), cell_sets)
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


test_that("upload_to_aws tries to upload the correct files to aws", {
  metadata <- list(Group1 = list("Hello", "WT2", "WT2"), Group2 = list("WT", "WT", "WT124"))
  input <- mock_input(metadata)
  config <- mock_config(input)
  scdata_list <- mock_scdata_list(config)

  paths <- setup_test_paths()

  pipeline_config <- mock_pipeline_config()

  prev_out <- list(
    config = config,
    counts_list = list(),
    annot = list(),
    doublet_scores = list(),
    scdata_list = scdata_list,
    qc_config = list("mock_qc_config"),
    disable_qc_filters = FALSE
  )

  res <- stubbed_upload_to_aws(input, pipeline_config, prev_out)

  # cellsets file
  expect_snapshot_file(
    file.path(pipeline_config$cell_sets_bucket, input$experimentId),
    name = "cellsets.json"
  )

  # raw sample seurat objects, test that they exist where upload_to_aws puts them
  for (sample_id in prev_out$config$samples) {
    expect_true(
      file.exists(file.path(
        pipeline_config$source_bucket,
        input$experimentId,
        sample_id,
        "r.rds"
      ))
    )
  }

  # cleanup
  withr::defer(unlink(pipeline_config$cell_sets_bucket, recursive = TRUE))
  withr::defer(unlink(pipeline_config$source_bucket, recursive = TRUE))
  withr::defer(unlink(file.path(paths$mock_data, "temp"), recursive = TRUE))
})


test_that("get_subset_cell_sets filters out louvain clusters from parent cellset", {
  input <- mock_input()
  config <- mock_config(input)
  scdata_list <- mock_scdata_list(config)
  disable_qc_filters <- TRUE

  cell_ids_to_keep <- c(4, 8, 15, 16, 23, 42)
  subset_scdata <- mock_subset_data(scdata_list, cell_ids_to_keep)
  mock_parsed_cellset <- mock_parsed_cellsets(scdata_list)

  prev_out <- mock_prev_out(config)
  prev_out$parent_cellsets <- mock_parsed_cellset
  sample_id_map <- mock_sample_id_map()
  prev_out$sample_id_map <- sample_id_map
  cell_sets <- get_subset_cell_sets(subset_scdata, input, prev_out, disable_qc_filters)

  sets_key <- sapply(cell_sets$cellSets, `[[`, "key")
  expect_equal(integer(0), which(sets_key == "louvain"))
})


test_that("get_subset_cell_sets produces a cellset with correct cell_ids", {
  input <- mock_input()
  config <- mock_config(input)
  scdata_list <- mock_scdata_list(config)
  disable_qc_filters <- TRUE

  cell_ids_to_keep <- c(4, 8, 15, 16, 23, 42)
  subset_scdata <- mock_subset_data(scdata_list, cell_ids_to_keep)
  mock_parsed_cellset <- mock_parsed_cellsets(scdata_list)

  prev_out <- mock_prev_out(config)
  prev_out$parent_cellsets <- mock_parsed_cellset
  sample_id_map <- mock_sample_id_map()
  prev_out$sample_id_map <- sample_id_map
  cell_sets <- get_subset_cell_sets(subset_scdata, input, prev_out, disable_qc_filters)

  sets_key <- sapply(cell_sets$cellSets, `[[`, "key")
  sample_sets <- cell_sets$cellSets[[which(sets_key == "sample")]]

  subset_cell_ids <- as.integer(unlist(lapply(subset_scdata, function(x) {
    x$cells_id
  })))
  cellset_cell_ids <- unlist(lapply(sample_sets$children, function(x) {
    x$cellIds
  }))

  expect_setequal(subset_cell_ids, cellset_cell_ids)
})


test_that("get_subset_cell_sets produces a cellset with correct new sample ids", {
  input <- mock_input()
  config <- mock_config(input)
  scdata_list <- mock_scdata_list(config)
  disable_qc_filters <- TRUE

  cell_ids_to_keep <- c(4, 8, 15, 16, 23, 42)
  subset_scdata <- mock_subset_data(scdata_list, cell_ids_to_keep)
  mock_parsed_cellset <- mock_parsed_cellsets(scdata_list)

  prev_out <- mock_prev_out(config)
  prev_out$parent_cellsets <- mock_parsed_cellset
  sample_id_map <- mock_sample_id_map()
  prev_out$sample_id_map <- sample_id_map
  cell_sets <- get_subset_cell_sets(subset_scdata, input, prev_out, disable_qc_filters)

  sets_key <- sapply(cell_sets$cellSets, `[[`, "key")
  sample_sets <- cell_sets$cellSets[[which(sets_key == "sample")]]

  subset_sample_ids <- unique(as.character(unlist(lapply(subset_scdata, function(x) {
    x$samples
  }))))
  cellset_sample_ids <- unlist(lapply(sample_sets$children, function(x) {
    x$key
  }))
  expect_equal(subset_sample_ids, cellset_sample_ids)
})


test_that("get_subset_cell_sets produces a cellset with correct cell_ids for each metadata group present in the parent experiment", {
  input <- mock_input()
  config <- mock_config(input)
  scdata_list <- mock_scdata_list(config)
  disable_qc_filters <- TRUE

  cell_ids_to_keep <- c(4, 8, 15, 16, 23, 42)
  subset_scdata <- mock_subset_data(scdata_list, cell_ids_to_keep)
  mock_parsed_cellset <- mock_parsed_cellsets(scdata_list)

  prev_out <- mock_prev_out(config)
  prev_out$parent_cellsets <- mock_parsed_cellset
  sample_id_map <- mock_sample_id_map()
  prev_out$sample_id_map <- sample_id_map

  cell_sets <- get_subset_cell_sets(subset_scdata, input, prev_out, disable_qc_filters)

  metadata_sets <- cell_sets$cellSets[[3]]
  metadata_key <- sapply(metadata_sets, `[[`, "key")

  for (i in 1:length(metadata_key)) {
    metadata_names <- sapply(metadata_sets[[i]]$children, `[[`, "name")
    for (j in 1:length(metadata_names)) {
      cellset_cell_ids <- metadata_sets[[i]]$children[[j]]$cellIds
      parent_cellset_metadata <- mock_parsed_cellset[cell_id %in% cellset_cell_ids][type == "metadata"]
      expect_equal(unique(parent_cellset_metadata$name), metadata_names[[j]])
    }
  }
})


test_that("get_subset_cell_sets produces a cellset with correct cell_ids for each scratchpad present in the parent experiment", {
  input <- mock_input()
  config <- mock_config(input)
  scdata_list <- mock_scdata_list(config)
  disable_qc_filters <- TRUE

  cell_ids_to_keep <- c(4, 8, 15, 16, 23, 42)
  subset_scdata <- mock_subset_data(scdata_list, cell_ids_to_keep)
  mock_parsed_cellset <- mock_parsed_cellsets(scdata_list)

  prev_out <- mock_prev_out(config)
  prev_out$parent_cellsets <- mock_parsed_cellset
  sample_id_map <- mock_sample_id_map()
  prev_out$sample_id_map <- sample_id_map

  cell_sets <- get_subset_cell_sets(subset_scdata, input, prev_out, disable_qc_filters)

  sets_key <- sapply(cell_sets$cellSets, `[[`, "key")
  scratchpad_sets <- cell_sets$cellSets[[which(sets_key == "scratchpad")]]
  scratchpad_key <- sapply(scratchpad_sets$children, `[[`, "key")

  for (i in 1:length(scratchpad_key)) {
    cellset_cell_ids <- scratchpad_sets$children[[i]]$cellIds
    scratchpad_names <- sapply(scratchpad_sets$children, `[[`, "name")
    parent_cellset_scratchpad <- mock_parsed_cellset[cell_id %in% cellset_cell_ids][type == "scratchpad"]
    expect_equal(unique(parent_cellset_scratchpad$key), scratchpad_key[[i]])
  }
})


test_that("filter_parent_cellsets only keeps cell_ids_to_keep", {

  input <- mock_input()
  config <- mock_config(input)
  scdata_list <- mock_scdata_list(config)
  cell_ids_to_keep <- c(4, 8, 15, 16, 23, 42)

  subset_scdata <- mock_subset_data(scdata_list, cell_ids_to_keep)
  parent_cellsets <- mock_parsed_cellsets(scdata_list)

  res <- filter_parent_cellsets(parent_cellsets, cell_ids_to_keep)

  expect_setequal(unique(res[,cell_id]), cell_ids_to_keep)
})


test_that("filter_parent_cellsets keeps all other variables equal to the parents' value", {

  input <- mock_input()
  config <- mock_config(input)
  scdata_list <- mock_scdata_list(config)
  cell_ids_to_keep <- c(4, 8, 15, 16, 23, 42)

  subset_scdata <- mock_subset_data(scdata_list, cell_ids_to_keep)
  parent_cellsets <- mock_parsed_cellsets(scdata_list)

  expected_parent_cellsets <- parent_cellsets[cell_id %in% cell_ids_to_keep & type != "cluster"]

  res <- filter_parent_cellsets(parent_cellsets, cell_ids_to_keep)

  expect_setequal(names(res), names(expected_parent_cellsets))
  expect_equal(res, expected_parent_cellsets)
})


test_that("build_scratchpad_cellsets builds single cellset when subsetting removes all but one", {

  input <- mock_input()
  config <- mock_config(input)
  color_pool <- get_color_pool()
  scdata_list <- mock_scdata_list(config)
  cell_ids_to_keep <- c(4, 8, 15, 16, 23, 42)

  subset_scdata <- mock_subset_data(scdata_list, cell_ids_to_keep)
  parent_cellsets <- mock_parsed_cellsets(scdata_list)
  subset_cellsets <- filter_parent_cellsets(parent_cellsets,cell_ids_to_keep)

  # remove one cell only from scratchpad cellsets to mimic real situation
  subset_cellsets <- subset_cellsets[!(type == "scratchpad" & cell_id == 16)]


  res <- build_scratchpad_cellsets(color_pool, subset_cellsets)

  expect_named(res, c("key", "name", "rootNode", "children", "type"))
  expect_named(res$children[[1]], c("key", "name", "color", "cellIds"))

  expect_equal(res$key, "scratchpad")
  expect_equal(res$children[[1]]$key, unique(subset_cellsets[type == "scratchpad", key]))
  expect_equal(res$children[[1]]$name, unique(subset_cellsets[type == "scratchpad", name]))
  expect_setequal(res$children[[1]]$cellIds, subset_cellsets[type == "scratchpad", cell_id])
})


test_that("build_scratchpad_cellsets builds multiple cellsets", {

  input <- mock_input()
  config <- mock_config(input)
  color_pool <- get_color_pool()
  scdata_list <- mock_scdata_list(config)
  cell_ids_to_keep <- c(4, 8, 15, 16, 23, 42, 23019:23027)

  subset_scdata <- mock_subset_data(scdata_list, cell_ids_to_keep)
  parent_cellsets <- mock_parsed_cellsets(scdata_list)
  subset_cellsets <- filter_parent_cellsets(parent_cellsets,cell_ids_to_keep)

  # remove a couple cells from scratchpad cellsets to mimic real situation
  subset_cellsets <- subset_cellsets[!(type == "scratchpad" & cell_id %in% c(16,23024))]


  res <- build_scratchpad_cellsets(color_pool, subset_cellsets)

  expect_equal(length(res$children), nrow(unique(subset_cellsets[type == "scratchpad"], by = "key")))

  for(cs in res$children) {
    expect_setequal(cs$cellIds, subset_cellsets[(type == "scratchpad" & key == cs$key & name == cs$name), cell_id])
  }

  })
