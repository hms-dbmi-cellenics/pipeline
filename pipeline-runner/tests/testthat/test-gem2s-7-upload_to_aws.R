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

mock_prepare_experiment <- function(config, use_bpcells = FALSE) {
  prev_out <- mock_create_seurat(config, use_bpcells = use_bpcells)
  input <- mock_input()

  prepare_experiment(input, NULL, prev_out)$output
}

mock_input <- function(metadata = NULL) {
  input <- list(
    name = "project name",
    sampleNames = list("a", "b", "c"),
    sampleIds = list("123abc", "123def", "123ghi"),
    metadata = metadata,
    experimentId = "mock_experiment_id",
    projectId = "mock_experiment_id",
    input = list(type = "10x")
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
    rep("metadata_var-1-value_A", 2 * length(samples) / 3), rep("metadata_var-1-value_B", length(samples) / 3)
  )
  name <- c(
    rep("Cluster 0", length(samples) / 2), rep("Cluster 1", length(samples) / 2),
    rep("Custom-A", length(samples) / 3), rep("Custom-B", 2 * length(samples) / 3),
    samples,
    rep("value_A", 2 * length(samples) / 3), rep("value_B", length(samples) / 3)
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


mock_create_seurat <- function(config, counts = NULL, use_bpcells = FALSE) {
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
  matrix_dir_list <- list()

  for (sampleId in samples) {
    sample_counts <- counts

    if (use_bpcells) {
      matrix_dir <- file.path(tempdir(), paste0(sampleId, "_matrix_dir"))
      unlink(matrix_dir, recursive = TRUE)
      sample_counts <- counts_to_bpcells(sample_counts, matrix_dir)
      matrix_dir_list[[sampleId]] <- matrix_dir
    }
    counts_list[[sampleId]] <- sample_counts
    edrops[[sampleId]] <- eout
    doublet_scores[[sampleId]] <- mock_doublet_scores(counts)
  }

  # as passed to create_seurat
  prev_out <- list(
    counts_list = counts_list,
    matrix_dir_list = matrix_dir_list,
    edrops = edrops,
    doublet_scores = doublet_scores,
    annot = data.frame(name = row.names(counts), input = row.names(counts)),
    config = config
  )
  # call create_seurat to get prev_out to pass to prepare_experiment
  return(create_seurat(NULL, NULL, prev_out)$output)
}

mock_subset_data <- function(scdata_list, cell_ids, mock_parsed_cellset) {
  subset_scdata <- list()
  sample_ids <- c("new-123abc", "new-123def", "new-123ghi")
  for (i in 1:length(scdata_list)) {
    parent_cell_ids <- scdata_list[[i]]$cells_id
    parent_metadata <- mock_parsed_cellset[cell_id %in% parent_cell_ids, ][type == "metadata", ]
    parent_metadata <- parent_metadata[order(match(cell_id, parent_cell_ids))]
    user_metadata <- extract_subset_user_metadata(mock_parsed_cellset)
    # add metadata from parent cellset to the seurat object
    for (j in length(names(user_metadata))) {
      valid_metadata_name <- make.names(c("samples", names(user_metadata)[1]), unique = TRUE)[-1]
      scdata_list[[i]]@meta.data[[valid_metadata_name]] <- parent_metadata[grepl(names(user_metadata)[1], key), name]
    }

    sample_scdata <- subset_ids(scdata_list[[i]], cell_ids)
    sample_scdata@meta.data$parent_samples <- sample_scdata@meta.data$samples
    sample_scdata@meta.data$samples <- sample_ids[i]

    subset_scdata[[i]] <- sample_scdata
  }

  names(subset_scdata) <- sample_ids
  return(subset_scdata)
}

mock_sample_id_map <- function() {
  sample_id_map <- list(
    "11111111-1111-1111-1111-111111111111" = "new-123abc",
    "22222222-2222-2222-2222-222222222222" = "new-123def",
    "33333333-3333-3333-3333-333333333333" = "new-123ghi"
  )
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

  scdata_list <- mock_prepare_experiment(config)$scdata_list

  cell_sets <- get_cell_sets(scdata_list, input)
  keys <- sapply(cell_sets$cellSets, `[[`, "key")

  expect_setequal(keys, c("scratchpad", "sample"))
})


test_that("get_cell_sets adds correct cell ids for each sample", {
  input <- mock_input()
  config <- mock_config(input)
  scdata_list <- mock_prepare_experiment(config)$scdata_list

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
  scdata_list <- mock_prepare_experiment(config)$scdata_list

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
  scdata_list <- mock_prepare_experiment(config)$scdata_list

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
  scdata_list <- mock_prepare_experiment(config)$scdata_list

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
  metadata <- list(
    Group1 = list("Hello", "WT2", "WT2"),
    Group2 = list("WT", "WT", "WTA")
  )
  input <- mock_input(metadata)
  config <- mock_config(input)
  scdata_list <- mock_prepare_experiment(config)$scdata_list

  cell_sets <- get_cell_sets(scdata_list, input)

  flat_cell_sets <- unlist(cell_sets)
  colors <- flat_cell_sets[grepl("[.]color", names(flat_cell_sets))]
  colors <- unname(colors)

  expect_equal(unique(colors), colors)
})


test_that("get_cell_sets without metadata matches snapshot", {
  input <- mock_input()
  config <- mock_config(input)
  scdata_list <- mock_prepare_experiment(config)$scdata_list

  cell_sets <- get_cell_sets(scdata_list, input)
  expect_snapshot(str(cell_sets))
})


test_that("get_cell_sets with two metadata groups matches snapshot", {
  metadata <- list(
    Group1 = list("Hello", "WT2", "WT2"),
    Group2 = list("WT", "WT", "WT124")
  )
  input <- mock_input(metadata)
  config <- mock_config(input)
  scdata_list <- mock_prepare_experiment(config)$scdata_list

  cell_sets <- get_cell_sets(scdata_list, input)

  expect_snapshot(str(cell_sets))
})


test_that("get_cell_sets converts numeric metadata values to strings", {
  # Test case for issue: obj2s Seurat uploads with numeric sample-level metadata
  # should have string names in cellsets, not numeric values
  metadata <- list(Seq_Batch = list(1, 2, 2))
  input <- mock_input(metadata)
  config <- mock_config(input)
  scdata_list <- mock_prepare_experiment(config)$scdata_list

  cell_sets <- get_cell_sets(scdata_list, input)

  # Get the Seq_Batch cellset
  seq_batch_cellset <- NULL
  for (cs in cell_sets$cellSets) {
    if (cs$key == "Seq_Batch") {
      seq_batch_cellset <- cs
      break
    }
  }

  expect_true(!is.null(seq_batch_cellset))
  expect_equal(seq_batch_cellset$type, "metadataCategorical")

  # Test that all children have string names, not numeric
  for (child in seq_batch_cellset$children) {
    expect_type(child$name, "character")
    expect_true(child$name %in% c("1", "2"))
  }
})


test_that("upload_to_aws tries to upload the correct files to aws", {
  metadata <- list(
    Group1 = list("Hello", "WT2", "WT2"),
    Group2 = list("WT", "WT", "WT124")
  )
  input <- mock_input(metadata)
  config <- mock_config(input)
  scdata_list <- mock_prepare_experiment(config)$scdata_list

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

test_that("upload_to_aws tries to upload the correct files to aws with bpcells", {
  metadata <- list(
    Group1 = list("Hello", "WT2", "WT2"),
    Group2 = list("WT", "WT", "WT124")
  )
  input <- mock_input(metadata)
  config <- mock_config(input)
  prev_out <- mock_prepare_experiment(config, use_bpcells = TRUE)
  paths <- setup_test_paths()
  pipeline_config <- mock_pipeline_config()

  # cleanup
  withr::defer(unlink(pipeline_config$cell_sets_bucket, recursive = TRUE))
  withr::defer(unlink(pipeline_config$source_bucket, recursive = TRUE))
  withr::defer(unlink(file.path(paths$mock_data, "temp"), recursive = TRUE))

  prev_out <- list(
    config = config,
    counts_list = list(),
    annot = list(),
    doublet_scores = list(),
    scdata_list = prev_out$scdata_list,
    matrix_dir_list = prev_out$matrix_dir_list,
    qc_config = list("mock_qc_config"),
    disable_qc_filters = FALSE
  )

  res <- stubbed_upload_to_aws(input, pipeline_config, prev_out)

  # raw sample seurat objects, test that they exist where upload_to_aws puts them
  for (sample_id in prev_out$config$samples) {
    print(list.files(file.path(
      pipeline_config$source_bucket,
      input$experimentId,
      sample_id
    )))
    expect_true(
      file.exists(file.path(
        pipeline_config$source_bucket,
        input$experimentId,
        sample_id,
        "r.rds"
      ))
    )
    expect_true(
      file.exists(file.path(
        pipeline_config$source_bucket,
        input$experimentId,
        sample_id,
        "matrix_dir.tar.zst"
      ))
    )
  }

  # cellsets file
  expect_snapshot_file(
    file.path(pipeline_config$cell_sets_bucket, input$experimentId),
    name = "cellsets.json"
  )

})


test_that("get_subset_cell_sets filters out louvain clusters from parent cellset", {
  input <- mock_input()
  config <- mock_config(input)
  scdata_list <- mock_prepare_experiment(config)$scdata_list
  disable_qc_filters <- TRUE

  cell_ids_to_keep <- c(4, 8, 15, 16, 23, 42)
  mock_parsed_cellset <- mock_parsed_cellsets(scdata_list)
  subset_scdata <- mock_subset_data(
    scdata_list, cell_ids_to_keep,
    mock_parsed_cellset
  )

  prev_out <- mock_create_seurat(config)
  prev_out$parent_cellsets <- mock_parsed_cellset
  sample_id_map <- mock_sample_id_map()
  prev_out$sample_id_map <- sample_id_map
  cell_sets <- get_subset_cell_sets(
    subset_scdata,
    input, 
    prev_out,
    disable_qc_filters
  )

  sets_key <- sapply(cell_sets$cellSets, `[[`, "key")
  expect_equal(integer(0), which(sets_key == "louvain"))
})


test_that("get_subset_cell_sets produces a cellset with correct cell_ids", {
  input <- mock_input()
  config <- mock_config(input)
  scdata_list <- mock_prepare_experiment(config)$scdata_list
  disable_qc_filters <- TRUE

  cell_ids_to_keep <- c(4, 8, 15, 16, 23, 42)
  mock_parsed_cellset <- mock_parsed_cellsets(scdata_list)
  subset_scdata <- mock_subset_data(
    scdata_list,
    cell_ids_to_keep,
    mock_parsed_cellset
  )

  prev_out <- mock_create_seurat(config)
  prev_out$parent_cellsets <- mock_parsed_cellset
  sample_id_map <- mock_sample_id_map()
  prev_out$sample_id_map <- sample_id_map
  cell_sets <- get_subset_cell_sets(
    subset_scdata,
    input,
    prev_out,
    disable_qc_filters
  )

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
  scdata_list <- mock_prepare_experiment(config)$scdata_list
  disable_qc_filters <- TRUE

  cell_ids_to_keep <- c(4, 8, 15, 16, 23, 42)
  mock_parsed_cellset <- mock_parsed_cellsets(scdata_list)
  subset_scdata <- mock_subset_data(
    scdata_list,
    cell_ids_to_keep,
    mock_parsed_cellset
  )

  prev_out <- mock_create_seurat(config)
  prev_out$parent_cellsets <- mock_parsed_cellset
  sample_id_map <- mock_sample_id_map()
  prev_out$sample_id_map <- sample_id_map
  cell_sets <- get_subset_cell_sets(
    subset_scdata,
    input,
    prev_out,
    disable_qc_filters
  )

  sets_key <- sapply(cell_sets$cellSets, `[[`, "key")
  sample_sets <- cell_sets$cellSets[[which(sets_key == "sample")]]

  subset_sample_ids <- unique(
    as.character(unlist(lapply(subset_scdata, function(x) {
      x$samples
    })))
  )
  cellset_sample_ids <- unlist(lapply(sample_sets$children, function(x) {
    x$key
  }))
  expect_equal(subset_sample_ids, cellset_sample_ids)
})


test_that("get_subset_cell_sets produces a cellset with correct cell_ids for each metadata group present in the parent experiment", {
  input <- mock_input()
  config <- mock_config(input)
  scdata_list <- mock_prepare_experiment(config)$scdata_list
  disable_qc_filters <- TRUE

  cell_ids_to_keep <- c(4, 8, 15, 16, 23, 42)
  mock_parsed_cellset <- mock_parsed_cellsets(scdata_list)
  subset_scdata <- mock_subset_data(
    scdata_list,
    cell_ids_to_keep,
    mock_parsed_cellset
  )

  prev_out <- mock_create_seurat(config)
  prev_out$parent_cellsets <- mock_parsed_cellset
  sample_id_map <- mock_sample_id_map()
  prev_out$sample_id_map <- sample_id_map

  cell_sets <- get_subset_cell_sets(
    subset_scdata,
    input,
    prev_out,
    disable_qc_filters
  )

  metadata_sets <- cell_sets$cellSets[[2]]
  metadata_key <- sapply(metadata_sets$children, `[[`, "key")

  for (i in seq_along(metadata_key)) {
    metadata_names <- metadata_sets$children[[i]]$name
    for (j in seq_along(metadata_names)) {
      cellset_cell_ids <- unlist(metadata_sets$children[[j]]$cellIds)
      parent_cellset_metadata <-
        mock_parsed_cellset[cell_id %in% cellset_cell_ids][type == "metadata"]
      expect_equal(unique(parent_cellset_metadata$name), metadata_names[[j]])
    }
  }
})


test_that("get_subset_cell_sets produces a cellset with correct cell_ids for each scratchpad present in the parent experiment", {
  input <- mock_input()
  config <- mock_config(input)
  scdata_list <- mock_prepare_experiment(config)$scdata_list
  disable_qc_filters <- TRUE

  cell_ids_to_keep <- c(4, 8, 15, 16, 23, 42)
  mock_parsed_cellset <- mock_parsed_cellsets(scdata_list)
  subset_scdata <- mock_subset_data(
    scdata_list,
    cell_ids_to_keep,
    mock_parsed_cellset
  )

  prev_out <- mock_create_seurat(config)
  prev_out$parent_cellsets <- mock_parsed_cellset
  sample_id_map <- mock_sample_id_map()
  prev_out$sample_id_map <- sample_id_map

  cell_sets <- get_subset_cell_sets(
    subset_scdata,
    input,
    prev_out,
    disable_qc_filters
  )

  sets_key <- sapply(cell_sets$cellSets, `[[`, "key")
  scratchpad_sets <- cell_sets$cellSets[[which(sets_key == "scratchpad")]]
  scratchpad_key <- sapply(scratchpad_sets$children, `[[`, "key")

  for (i in seq_along(scratchpad_key)) {
    cellset_cell_ids <- scratchpad_sets$children[[i]]$cellIds
    scratchpad_names <- sapply(scratchpad_sets$children, `[[`, "name")
    parent_cellset_scratchpad <-
      mock_parsed_cellset[cell_id %in% cellset_cell_ids][type == "scratchpad"]
    expect_equal(unique(parent_cellset_scratchpad$key), scratchpad_key[[i]])
  }
})


test_that("filter_parent_cellsets only keeps cell_ids_to_keep", {
  input <- mock_input()
  config <- mock_config(input)
  scdata_list <- mock_prepare_experiment(config)$scdata_list
  cell_ids_to_keep <- c(4, 8, 15, 16, 23, 42)

  parent_cellsets <- mock_parsed_cellsets(scdata_list)
  subset_scdata <- mock_subset_data(
    scdata_list,
    cell_ids_to_keep,
    parent_cellsets
  )

  res <- filter_parent_cellsets(parent_cellsets, cell_ids_to_keep)

  expect_setequal(unique(res[, cell_id]), cell_ids_to_keep)
})


test_that("filter_parent_cellsets keeps all other variables equal to the parents' value", {
  input <- mock_input()
  config <- mock_config(input)
  scdata_list <- mock_prepare_experiment(config)$scdata_list
  cell_ids_to_keep <- c(4, 8, 15, 16, 23, 42)

  parent_cellsets <- mock_parsed_cellsets(scdata_list)
  subset_scdata <- mock_subset_data(
    scdata_list,
    cell_ids_to_keep,
    parent_cellsets
  )

  expected_parent_cellsets <-
    parent_cellsets[cell_id %in% cell_ids_to_keep & type != "cluster"]

  res <- filter_parent_cellsets(parent_cellsets, cell_ids_to_keep)

  expect_setequal(names(res), names(expected_parent_cellsets))
  expect_equal(res, expected_parent_cellsets)
})


test_that("build_scratchpad_cellsets builds single cellset when subsetting removes all but one", {
  input <- mock_input()
  config <- mock_config(input)
  color_pool <- get_color_pool()
  scdata_list <- mock_prepare_experiment(config)$scdata_list
  cell_ids_to_keep <- c(4, 8, 15, 16, 23, 42)

  parent_cellsets <- mock_parsed_cellsets(scdata_list)
  subset_scdata <- mock_subset_data(
    scdata_list,
    cell_ids_to_keep,
    parent_cellsets
  )
  subset_cellsets <- filter_parent_cellsets(parent_cellsets, cell_ids_to_keep)

  # remove one cell only from scratchpad cellsets to mimic real situation
  subset_cellsets <- subset_cellsets[!(type == "scratchpad" & cell_id == 16)]


  res <- build_scratchpad_cellsets(color_pool, subset_cellsets)

  expect_named(res, c("key", "name", "rootNode", "children", "type"))
  expect_named(res$children[[1]], c("key", "name", "color", "cellIds"))

  expect_equal(res$key, "scratchpad")
  expect_equal(
    res$children[[1]]$key,
    unique(subset_cellsets[type == "scratchpad", key])
  )
  expect_equal(
    res$children[[1]]$name,
    unique(subset_cellsets[type == "scratchpad", name])
  )
  expect_setequal(
    res$children[[1]]$cellIds,
    subset_cellsets[type == "scratchpad", cell_id]
  )
})


test_that("build_scratchpad_cellsets builds multiple cellsets", {
  input <- mock_input()
  config <- mock_config(input)
  color_pool <- get_color_pool()
  scdata_list <- mock_prepare_experiment(config)$scdata_list
  cell_ids_to_keep <- c(4, 8, 15, 16, 23, 42, 23019:23027)

  parent_cellsets <- mock_parsed_cellsets(scdata_list)
  subset_scdata <- mock_subset_data(
    scdata_list,
    cell_ids_to_keep,
    parent_cellsets
  )
  subset_cellsets <- filter_parent_cellsets(parent_cellsets, cell_ids_to_keep)

  # remove a couple cells from scratchpad cellsets to mimic real situation
  subset_cellsets <-
    subset_cellsets[!(type == "scratchpad" & cell_id %in% c(16, 23024))]


  res <- build_scratchpad_cellsets(color_pool, subset_cellsets)

  expect_equal(
    length(res$children),
    nrow(unique(subset_cellsets[type == "scratchpad"], by = "key"))
  )

  for (cs in res$children) {
    expect_setequal(
      cs$cellIds,
      subset_cellsets[
        (type == "scratchpad" & key == cs$key & name == cs$name),
        cell_id
      ]
    )
  }
})


test_that("extract_subset_user_metadata extracts subset cellsets correctly when regex characters present", {
  input <- mock_input()
  config <- mock_config(input)
  color_pool <- get_color_pool()
  scdata_list <- mock_prepare_experiment(config)$scdata_list

  parent_cellsets <- mock_parsed_cellsets(scdata_list)

  # create new metadata cellsets with some annoying symbols
  new_parent_cellsets <- data.table::copy(parent_cellsets[type == "metadata"])
  new_parent_cellsets[
    1:.N / 2,
    `:=`(key = "metadata_var-2-value_A+", name = "value_A+")
  ]
  new_parent_cellsets[
    (.N / 2):.N,
    `:=`(key = "metadata_var-2-value+_B-", name = "value+_B-")
  ]
  parent_cellsets <- data.table::rbindlist(
    list(parent_cellsets, new_parent_cellsets)
  )

  res <- extract_subset_user_metadata(parent_cellsets)

  expect_equal(names(res), c("metadata_var-1", "metadata_var-2"))
  expect_equal(res[[1]], c("value_A", "value_B"))
  expect_equal(res[[2]], c("value_A+", "value+_B-"))
})

# Tests for tar_matrix_dir function
test_that("tar_matrix_dir creates tarfile with correct name", {
  # Create a temporary directory with a test file
  base_temp <- withr::local_tempdir()
  sample_id <- paste0("sample_", sample(100000, 1))
  matrix_dir <- file.path(base_temp, sample_id, "matrix_dir")
  dir.create(matrix_dir, showWarnings = FALSE, recursive = TRUE)

  # Create a test file in the directory
  test_file <- file.path(matrix_dir, "test.txt")
  writeLines("test content", test_file)

  # Call tar_matrix_dir
  tarfile <- tar_matrix_dir(sample_id, matrix_dir)

  expect_true(file.exists(tarfile))
  expect_match(tarfile, paste0(sample_id, "_matrix_dir\\.tar\\.zst$"))
})

test_that("tar_zstd creates compressed tarfile", {
  # Create a temporary directory with test files
  base_temp <- withr::local_tempdir()
  test_id <- paste0("tar_test_", sample(100000, 1))
  test_dir <- file.path(base_temp, test_id)
  dir.create(test_dir, showWarnings = FALSE)

  test_file <- file.path(test_dir, "test.txt")
  writeLines("test content", test_file)

  tarfile <- file.path(base_temp, paste0(test_id, ".tar.zst"))

  # Change to base directory and create tar
  current_dir <- getwd()
  tryCatch({
    setwd(base_temp)
    tar_zstd(tarfile, files = test_id)
  }, finally = {
    setwd(current_dir)
  })

  expect_true(file.exists(tarfile))
})

test_that("untar_zstd extracts compressed tarfile", {
  # Create a temporary directory with test files
  base_temp <- withr::local_tempdir()
  test_id <- paste0("untar_test_", sample(100000, 1))
  test_dir <- file.path(base_temp, test_id)
  dir.create(test_dir, showWarnings = FALSE)

  test_file <- file.path(test_dir, "test.txt")
  writeLines("test content", test_file)

  # Create tarfile
  tarfile <- file.path(base_temp, paste0(test_id, "_extract.tar.zst"))
  current_dir <- getwd()
  tryCatch({
    setwd(base_temp)
    tar_zstd(tarfile, files = test_id)
  }, finally = {
    setwd(current_dir)
  })

  # Create extraction directory
  extract_dir <- file.path(base_temp, paste0(test_id, "_extracted"))
  dir.create(extract_dir, showWarnings = FALSE)

  # Extract tarfile
  untar_zstd(tarfile, exdir = extract_dir)

  # Verify extraction
  extracted_file <- file.path(extract_dir, test_id, "test.txt")
  expect_true(file.exists(extracted_file))
})


# ─── Xenium spatial upload (gem2s-7) ────────────────────────────────────────

test_that("get_image_scale bypasses the scale lookup for Xenium (returns NULL)", {
  # a Xenium FOV has no @image; the scaleless branch must return before reading it
  scdata <- mock_spatial_scdata(ncells = 16, technology = "xenium")
  expect_null(get_image_scale(Seurat::Images(scdata), scdata))
})


test_that("get_polygon_coords reads segmentation coords with no ScaleFactors multiply", {
  scdata <- mock_spatial_scdata(ncells = 9, sample_id = "s1")
  cells <- colnames(scdata)[1:3]

  fake_coords <- data.frame(
    x = c(1, 2, 3), y = c(4, 5, 6), cell = cells, stringsAsFactors = FALSE
  )
  gtc <- mockery::mock(fake_coords)
  scalefactors_called <- FALSE

  mockery::stub(get_polygon_coords, "Seurat::GetTissueCoordinates", gtc)
  mockery::stub(
    get_polygon_coords, "Seurat::ScaleFactors",
    function(...) {
      scalefactors_called <<- TRUE
      stop("ScaleFactors must not be used for FOV polygon coords")
    }
  )

  res <- get_polygon_coords(Seurat::Images(scdata), scdata, scale = NULL)

  # the boundary is selected via which = "segmentations", scale passed through
  args <- mockery::mock_args(gtc)[[1]]
  expect_equal(args$which, "segmentations")
  expect_null(args$scale)

  # no pixel->coord scaling (coords are already microns)
  expect_false(scalefactors_called)

  # cells_id attached + arranged for the polygon upload
  expect_true("cells_id" %in% colnames(res))
})


test_that("upload_to_aws (xenium) skips the image upload, uploads polygons, cleans the segmentations bucket", {
  input <- mock_input()
  input$input <- list(type = "xenium")
  config <- mock_config(input)
  config$input <- list(type = "xenium")

  scdata_list <- mock_prepare_experiment(config)$scdata_list
  scdata_list <- add_tissue_coords(scdata_list) # attach an FOV per sample

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

  removed_buckets <- c()
  image_uploads <- 0
  polygon_uploads <- 0

  mockery::stub(upload_to_aws, "put_object_in_s3", stub_put_object_in_s3)
  mockery::stub(upload_to_aws, "put_object_in_s3_multipart", stub_put_object_in_s3_multipart)
  mockery::stub(upload_to_aws, "tempdir", stub_tempdir)
  mockery::stub(upload_to_aws, "remove_bucket_folder", function(pipeline_config, bucket, ...) {
    removed_buckets <<- c(removed_buckets, bucket)
  })
  mockery::stub(upload_to_aws, "upload_image_to_s3", function(...) image_uploads <<- image_uploads + 1)
  mockery::stub(upload_to_aws, "upload_polygons_to_s3", function(...) polygon_uploads <<- polygon_uploads + 1)
  mockery::stub(upload_to_aws, "upload_molecules_to_s3", function(...) NULL)
  mockery::stub(upload_to_aws, "get_polygon_coords", function(...) data.frame(x = c(1, 2, 3), y = c(1, 2, 3)))

  upload_to_aws(input, pipeline_config, prev_out)

  # no tissue image for Xenium; one polygon upload per sample
  expect_equal(image_uploads, 0)
  expect_equal(polygon_uploads, length(scdata_list))

  # the segmentation OME-Zarr bucket is included in the re-run cleanup (leak fix)
  expect_true(pipeline_config$spatial_segmentations_bucket %in% removed_buckets)

  # the molecule artifact bucket is also cleaned on re-run (same leak class)
  expect_true(pipeline_config$spatial_molecules_bucket %in% removed_buckets)

  withr::defer(unlink(pipeline_config$cell_sets_bucket, recursive = TRUE))
  withr::defer(unlink(pipeline_config$source_bucket, recursive = TRUE))
})


test_that("upload_to_aws (xenium) builds the molecule artifact for every sample", {
  input <- mock_input()
  input$input <- list(type = "xenium")
  config <- mock_config(input)
  config$input <- list(type = "xenium")

  scdata_list <- mock_prepare_experiment(config)$scdata_list
  scdata_list <- add_tissue_coords(scdata_list)

  # transcripts.parquet is a required Xenium input, so every sample carries a
  # molecules frame onto the object
  for (sample_name in names(scdata_list)) {
    scdata_list[[sample_name]]@misc$transcripts <- data.frame(
      x = runif(20), y = runif(20), gene = sample(c("Gad1", "Sst"), 20, TRUE)
    )
  }

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

  molecule_samples <- c()

  mockery::stub(upload_to_aws, "put_object_in_s3", stub_put_object_in_s3)
  mockery::stub(upload_to_aws, "put_object_in_s3_multipart", stub_put_object_in_s3_multipart)
  mockery::stub(upload_to_aws, "tempdir", stub_tempdir)
  mockery::stub(upload_to_aws, "remove_bucket_folder", function(...) NULL)
  mockery::stub(upload_to_aws, "upload_polygons_to_s3", function(...) NULL)
  mockery::stub(upload_to_aws, "get_polygon_coords", function(...) data.frame(x = c(1, 2, 3), y = c(1, 2, 3)))
  mockery::stub(
    upload_to_aws, "upload_molecules_to_s3",
    function(pipeline_config, input, experiment_id, transcripts, sample_id, ...) {
      molecule_samples <<- c(molecule_samples, sample_id)
    }
  )

  upload_to_aws(input, pipeline_config, prev_out)

  # built for every xenium sample
  expect_setequal(molecule_samples, names(scdata_list))

  withr::defer(unlink(pipeline_config$cell_sets_bucket, recursive = TRUE))
  withr::defer(unlink(pipeline_config$source_bucket, recursive = TRUE))
})


test_that("upload_molecules_to_s3 filters QV, writes one Feather per gene + meta.json, uploads + registers", {
  set.seed(7)
  # 6 distinct genes; qv spans the Q20 threshold so the filter drops rows
  transcripts <- data.frame(
    x = runif(300, 0, 100),
    y = runif(300, 0, 100),
    gene = sample(paste0("Gene", 1:6), 300, replace = TRUE),
    qv = runif(300, 5, 40)
  )

  pipeline_config <- mock_pipeline_config()
  input <- list(authJWT = "mock_jwt")

  uploaded <- list()
  registered <- list()

  mockery::stub(upload_molecules_to_s3, "tempdir", stub_tempdir)
  mockery::stub(
    upload_molecules_to_s3, "put_object_in_s3_multipart",
    function(pipeline_config, bucket, object, key) {
      uploaded <<- list(bucket = bucket, object = object, key = key)
    }
  )
  mockery::stub(
    upload_molecules_to_s3, "create_sample_file",
    function(api_url, experiment_id, sample_id, file_type, file_size, sample_file_id, overwrite_existing, auth_jwt) {
      registered <<- list(
        file_type = file_type, sample_file_id = sample_file_id,
        sample_id = sample_id, overwrite_existing = overwrite_existing
      )
    }
  )

  kept <- transcripts[transcripts$qv >= 20, ]
  expected_kept <- nrow(kept)
  ngenes <- length(unique(kept$gene))

  upload_molecules_to_s3(
    pipeline_config, input, "mock_experiment_id", transcripts, "s1",
    overwrite_existing = TRUE
  )

  molecules_dir <- file.path(stub_tempdir(), "s1.molecules.bygene")

  # meta.json: gene-partitioned shape + the dense dictionary covering every kept gene
  meta <- jsonlite::read_json(file.path(molecules_dir, "meta.json"), simplifyVector = FALSE)
  expect_equal(meta$version, 2)
  expect_equal(meta$qvThreshold, 20)
  # no quadtree fields any more
  expect_null(meta$maxDepth)
  expect_null(meta$tileSize)
  expect_null(meta$pointColumns)
  expect_false(is.null(meta$rootExtent))
  expect_length(meta$genes, ngenes)

  # dense 0-based codes (colour is assigned in the UI from the code, not baked here)
  meta_codes <- vapply(meta$genes, function(g) g$code, numeric(1))
  expect_equal(sort(meta_codes), 0:(ngenes - 1))
  expect_true(all(vapply(meta$genes, function(g) is.null(g$color), logical(1))))

  # QV filter applied before partitioning: per-gene nPoints sum to the kept total
  expect_equal(sum(vapply(meta$genes, function(g) g$nPoints, numeric(1))), expected_kept)

  # each gene's entry is a Feather with x/y only (float32), nrow == nPoints, in the
  # kept frame. float32 rounding => compare with a tolerance, not exact equality.
  for (g in meta$genes) {
    entry_path <- file.path(molecules_dir, g$entry)
    expect_equal(g$entry, paste0(g$code, ".feather"))
    expect_true(file.exists(entry_path))
    frame <- as.data.frame(arrow::read_feather(entry_path))
    expect_equal(colnames(frame), c("x", "y"))
    expect_equal(nrow(frame), g$nPoints)
    expect_equal(sort(frame$x), sort(kept$x[kept$gene == g$gene]), tolerance = 1e-4)
  }

  # uploaded to the molecules bucket + registered as molecules_by_gene
  expect_equal(uploaded$bucket, pipeline_config$spatial_molecules_bucket)
  expect_equal(uploaded$key, registered$sample_file_id)
  expect_equal(registered$file_type, "molecules_by_gene")
  expect_equal(registered$sample_id, "s1")
  expect_true(registered$overwrite_existing)
})


test_that("upload_molecules_to_s3 records qvThreshold null when no qv column", {
  set.seed(8)
  transcripts <- data.frame(
    x = runif(50, 0, 10),
    y = runif(50, 0, 10),
    gene = sample(c("A", "B"), 50, replace = TRUE)
  )

  pipeline_config <- mock_pipeline_config()
  input <- list(authJWT = "mock_jwt")

  mockery::stub(upload_molecules_to_s3, "tempdir", stub_tempdir)
  mockery::stub(upload_molecules_to_s3, "put_object_in_s3_multipart", function(...) NULL)
  mockery::stub(upload_molecules_to_s3, "create_sample_file", function(...) NULL)

  upload_molecules_to_s3(
    pipeline_config, input, "mock_experiment_id", transcripts, "s2",
    overwrite_existing = TRUE
  )

  meta_path <- file.path(stub_tempdir(), "s2.molecules.bygene", "meta.json")
  meta <- jsonlite::read_json(meta_path, simplifyVector = FALSE)
  # no qv column -> threshold recorded as null (not silently applied)
  expect_null(meta$qvThreshold)
  # both genes present, all 50 molecules partitioned across them
  expect_length(meta$genes, 2)
  expect_equal(sum(vapply(meta$genes, function(g) g$nPoints, numeric(1))), 50)
})
