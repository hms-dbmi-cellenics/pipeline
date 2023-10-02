mock_scdata <- function() {
  data("pbmc_small", package = "SeuratObject", envir = environment())
  # shuffle cell ids, as we do in the platform
  set.seed(1)
  pbmc_small$cells_id <- sample(0:(ncol(pbmc_small) - 1))
  pbmc_small$samples <- rep_len(paste0("sample_", 1:4), ncol(pbmc_small))
  pbmc_small@misc$gene_annotations <- data.frame(
    input = paste0("ENSG", seq_len(nrow(pbmc_small))),
    name = row.names(pbmc_small),
    row.names = paste0("ENSG", seq_len(nrow(pbmc_small)))
  )
  # we have ~500 colors in the color pool
  pbmc_small@misc$color_pool <- mock_color_pool(500)
  return(pbmc_small)
}

mock_cellset_object <- function(n_cells, n_clusters) {

  if (n_clusters == 0) {
    return(data.frame(cluster = integer(), cell_ids = integer()))
  }

  # cell_ids with no replacement to avoid repeated cell_ids
  data.frame(cluster = sample(1:n_clusters, size = n_cells, replace = T),
             cell_ids = sample(1:(2*n_cells), size = n_cells))
}

mock_color_pool <- function(n) {
  paste0("color_", 1:n)
}


mock_cl_metadata <- function(scdata) {
  barcode <-  rownames(scdata@meta.data)
  samples <- rep_len(paste0("sample_", 1:4), length(barcode))
  cell_type <- rep_len(paste0("cell_type_", 1:10), length(barcode))
  group_var <- rep_len(paste0("group_", 1:2), length(barcode))
  redundant_group_var <- paste0("red_group_", group_var)
  continuous_var <- rnorm(length(barcode))

  data.table::data.table(barcode, samples, cell_type, group_var, redundant_group_var, continuous_var)
}


local_mock_cl_metadata_table <- function(mock_cl_metadata, experiment_id, env = parent.frame()) {
  # calls mock_cl_metadata but makes it "local" (in withr speech), deleting
  # created stuff after the test finishes.
  bucket <- file.path(".", bucket_list$cl_metadata_bucket)

  dir.create(bucket)
  withr::defer(unlink(bucket, recursive = TRUE), envir = env)

  data.table::fwrite(mock_cl_metadata, file.path(bucket, experiment_id), sep = "\t")
}


mock_config <- function() {
  config <-
    list(clusteringSettings = list(
      method = "louvain",
      methodSettings = list(louvain = list(resolution = "0.8"))
    ),
    aws_config = list(aws = "mock_aws"),
    metadataS3Path = file.path(".", bucket_list$cl_metadata_bucket, "mock_experiment_id")
    )

  return(config)
}

test_that("runClusters returns correct keys", {
  algos <- c("louvain", "leiden")
  data <- mock_scdata()
  expected_keys <- c("cluster", "cell_ids")
  resolution <- 0.8

  for (algo in algos) {
    res <- runClusters(algo, resolution, data)
    expect_equal(names(res), expected_keys)
  }
})

test_that("runClusters returns one value per cell", {
  algos <- c("louvain", "leiden")
  data <- mock_scdata()
  expected_n_cells <- ncol(data)
  resolution <- 0.8


  for (algo in algos) {
    res <- runClusters(algo, resolution, data)
    n_cells <- nrow(res)

    expect_equal(n_cells, expected_n_cells)
  }
})

test_that("runClusters returns same order of barcodes as seurat object", {
  algos <- c("louvain", "leiden")
  data <- mock_scdata()
  expected_barcodes <- colnames(data)
  resolution <- 0.8


  for (algo in algos) {
    res <- runClusters(algo, resolution, data)
    barcodes <- rownames(res)
    expect_equal(barcodes, expected_barcodes)
  }
})

test_that("runClusters returns at least one cluster", {
  algos <- c("louvain", "leiden")
  data <- mock_scdata()
  resolution <- 0.8

  for (algo in algos) {
    res <- runClusters(algo, resolution, data)
    n_clusters <- length(unique(res$cluster))
    expect_gte(n_clusters, 1)
  }
})

test_that("runClusters uses active.reduction in misc slot", {
  algos <- c("louvain", "leiden")
  data <- mock_scdata()
  resolution <- 0.8

  # get error if no PCA/SNN graph
  blah_reduction <- data@reductions$pca
  data@reductions$pca <- NULL
  data@graphs$RNA_snn <- NULL

  for (algo in algos) {
    expect_error(runClusters(algo, resolution, data), "Cannot find 'pca'")
  }

  # will use active.reduction to get SNN graph
  data@reductions$blah_reduction <- blah_reduction
  data@misc$active.reduction <- "blah_reduction"

  for (algo in algos) {
    expect_error(runClusters(algo,resolution, data), NA)
  }
})

with_fake_http(
  test_that("replace_cell_class_through_api sends patch request", {
    expect_PATCH(
      replace_cell_class_through_api(
        list(),
        "api_url",
        "experiment_id",
        "cell_set_key",
        "auth",
        FALSE
    )
    )
  })
)

with_fake_http(
 test_that("replace_cell_class_through_api diables SSL certificate checking if disabled", {

    expect_PATCH(
      replace_cell_class_through_api(
        list(),
        "api_url",
        "experiment_id",
        "cell_set_key",
        "auth",
        TRUE
    )
    )

    # set random config to get the previously used config
    old_config <- httr::set_config(httr::config(ssl_verifyhost = 0L))
    expect_equal(old_config$options$ssl_verifypeer, 0)

  })
)

test_that("format_cluster_cellsets returns correct items", {
  n_clusters <- 5
  cell_sets <- mock_cellset_object(1000, n_clusters)
  color_pool <- mock_color_pool(n_clusters)
  algos <- c("louvain", "leiden")

  types <- c(
    "key" = "character",
    "name" = "character",
    "rootNode" = "logical",
    "type" = "character",
    "children" = "list"
  )

  for (algo in algos) {
    res <- format_cluster_cellsets(cell_sets, algo, color_pool)

    # expect all keys present
    expect_setequal(names(res), names(types))

    for (item in names(types)) {
      expect_type(res[[item]], types[[item]])
    }
  }
})


test_that("format_cluster_cellsets correctly formats a cellset object", {

  n_clusters <- 5
  cell_sets <- mock_cellset_object(1000, n_clusters)
  color_pool <- mock_color_pool(n_clusters)
  algos <- c("louvain", "leiden")

  expected_items <- c("key", "name", "rootNode", "type", "color", "cellIds")

  for (algo in algos) {
    res <- format_cluster_cellsets(cell_sets, algo, color_pool)


    # each children has the expected items
    for (cluster in res$children) {
      expect_setequal(names(cluster), expected_items)
    }

    # each cluster has correct content
    for (cluster in res$children) {
      expect_true(startsWith(cluster$key, algo))
      expect_true(startsWith(cluster$name, "Cluster"))
      expect_type(cluster$cellIds, "integer")
      expect_gte(length(cluster$cellIds), 1)
    }
  }
})


test_that("format_cluster_cellsets result has correct number of clusters",{
  n_clusters <- c(0, 1, 4, 6)
  algos <- c("louvain", "leiden")

  for (algo in algos) {
    for (n in n_clusters) {
      cell_sets <- mock_cellset_object(100, n)
      color_pool <- mock_color_pool(n)

      res <- format_cluster_cellsets(cell_sets, algo, color_pool)

      # number of clusters is the number of children elements
      expect_equal(length(res$children), n)

    }}
})

test_that("format_cluster_cellsets returns empty children on empty cellset", {
  n_clusters <- 5
  cell_sets <- mock_cellset_object(0, n_clusters)

  color_pool <- mock_color_pool(n_clusters)
  algos <- c("louvain", "leiden")

  for (algo in algos) {
    res <- format_cluster_cellsets(cell_sets, algo, color_pool)
    expect_equal(length(res$children), 0)
  }
})



test_that("embed_and_cluster works", {
  scdata <- mock_scdata()

  config <- list(
    clusteringSettings = list(
      method = "louvain",
      methodSettings = list(louvain = list(resolution = 0.8))
    ),
    embeddingSettings = list(
      method = "umap",
      methodSettings = list(
        tsne = list(
          learningRate = 738.75,
          perplexity = 30
        ),
        umap = list(
          distanceMetric = "cosine",
          minimumDistance = 0.3
        )
      )
    )
  )

  sample_id <- "mock_sample_id"
  cells_id <- "83abbeea-b664-499a-8e6e-d2dbae4c60a9"
  task_name <- "configureEmbedding"
  ignore_ssl_cert <- FALSE

  cell_sets_bucket <- "./mock_data/cell_sets_bucket"
  stubbed_embed_and_cluster(
    scdata,
    config,
    sample_id,
    cells_id,
    task_name,
    ignore_ssl_cert
)

  expect_snapshot_file(file.path(cell_sets_bucket, "cluster_cellsets.json"),
                       name = "cluster_cell_sets.json")
  withr::defer(unlink(cell_sets_bucket, recursive = TRUE))

})


test_that("runClusters does not crash with less than 10 dimensions available", {
  algos <- c("louvain", "leiden")
  scdata <- mock_scdata()
  expected_keys <- c("cluster", "cell_ids")
  resolution <- 0.8

  # remove all pre-existing reductions and calculate low-PC PCA
  scdata <- Seurat::DietSeurat(scdata, scale.data = T)
  scdata <- suppressWarnings(Seurat::RunPCA(scdata, assay = "RNA", npcs = 2, verbose = F))

  for (algo in algos) {
    res <- runClusters(algo, resolution, scdata)
    expect_equal(names(res), expected_keys)
  }
})


test_that("getClusters uses the default value of 10 if there are enough PCs available",{
  algos <- c("louvain", "leiden")
  scdata <- mock_scdata()
  resolution <- 0.8

  # remove all pre-existing reductions and calculate low-PC PCA
  scdata <- Seurat::DietSeurat(scdata, scale.data = T)
  scdata@commands <- list()
  scdata <- suppressWarnings(Seurat::RunPCA(scdata, assay = "RNA", npcs = 20, verbose = F))

  for (algo in algos) {
    clustered_scdata <- getClusters(algo, resolution, scdata)
    if (algo == "louvain") expect_equal(clustered_scdata@commands$FindNeighbors.RNA.pca$dims, 1:10)
    # difficult to test in leiden, so test internal state as proxy
    if (algo == "leiden") expect_true("seurat_clusters" %in% names(clustered_scdata@meta.data))
  }
})



test_that("get_cell_id_barcode_map returns correct table with barcode and cell ids", {
  scdata <- mock_scdata()
  res <- get_cell_id_barcode_map(scdata)

  expect_named(res, c("barcode", "cells_id", "samples"))
  expect_type(res$barcode, "character")
  expect_type(res$cells_id, "integer")
  expect_type(res$samples, "character")
  expect_equal(nrow(res), ncol(scdata))
})



test_that("make_cl_metadata_table correctly joins cell ids to metadata table", {
  scdata <- mock_scdata()
  cl_meta <- mock_cl_metadata(scdata)

  cell_id_barcode_map <- get_cell_id_barcode_map(scdata)
  res <- make_cl_metadata_table(cl_meta, cell_id_barcode_map)

  expect_named(res, c(names(cl_meta), "cells_id"))
  expect_equal(rownames(scdata@meta.data), res$barcode)
  expect_equal(scdata@meta.data$cells_id, res$cells_id)
})


test_that("make_cl_metadata_table uses composite primary key for join if samples in user-supplied cl metadata", {
  scdata <- mock_scdata()
  cl_meta <- mock_cl_metadata(scdata)

  cell_id_barcode_map <- get_cell_id_barcode_map(scdata)

  # tables should be identical save for the duplicated barcode
  expected <- make_cl_metadata_table(cl_meta, cell_id_barcode_map)
  expected$barcode[1] <- expected$barcode[2]

  # artificially duplicate a barcode
  cl_meta$barcode[1] <- cl_meta$barcode[2]
  cell_id_barcode_map$barcode[1] <- cell_id_barcode_map$barcode[2]

  res <- make_cl_metadata_table(cl_meta, cell_id_barcode_map)

  expect_equal(res, expected)
})


test_that("make_cl_metadata_table joins on barcode only if no samples column in user-supplied cl metadata", {
  scdata <- mock_scdata()
  cl_meta <- mock_cl_metadata(scdata)

  # remove samples column from cell-level metadata
  cl_meta[, samples:= NULL]

  cell_id_barcode_map <- get_cell_id_barcode_map(scdata)
  res <- make_cl_metadata_table(cl_meta, cell_id_barcode_map)

  # artificially duplicate a barcode
  cl_meta$barcode[1] <- cl_meta$barcode[2]
  cell_id_barcode_map$barcode[1] <- cell_id_barcode_map$barcode[2]

  # the join will return 2 extra rows (when there are 2 dups) if the primary key
  # is the barcode only
  res2 <- make_cl_metadata_table(cl_meta, cell_id_barcode_map)

  expect_true(nrow(res2) == nrow(res) + 2)
})


test_that("make_cl_metadata_cellclass makes correctly formatted cellsets", {
  scdata <- mock_scdata()
  cl_meta <- mock_cl_metadata(scdata)

  # remove samples column from cell-level metadata
  cl_meta[, samples:= NULL]

  cell_id_barcode_map <- get_cell_id_barcode_map(scdata)
  cl_metadata_table <- make_cl_metadata_table(cl_meta, cell_id_barcode_map)

  var_to_cellset <- "cell_type"
  cellset_type <- "CLM"

  res <- make_cl_metadata_cellclass(var_to_cellset,
                                  cellset_type,
                                  cl_metadata_table,
                                  mock_color_pool(20))

  cell_class_names <- c("key", "name", "rootNode", "type", "children")
  expect_named(res, cell_class_names)

  # cellsets have the same keys as cell classes except children, color and cellIds
  purrr::walk(res$children, expect_named, c(cell_class_names[-5], "color", "cellIds"))

})


test_that("detect_variable_types correctly detects variable types", {
  scdata <- mock_scdata()
  cl_meta <- mock_cl_metadata(scdata)

  cell_id_barcode_map <- get_cell_id_barcode_map(scdata)
  cl_metadata <-
    make_cl_metadata_table(cl_meta, cell_id_barcode_map)

  expected_var_types <-
    list(
      CLM = c("cell_type"),
      CLMPerSample = c("group_var", "redundant_group_var")
    )

  res <- detect_variable_types(cl_metadata)

  expect_equal(res, expected_var_types)
})


test_that("detect_variable_types removes samples variable", {
  scdata <- mock_scdata()
  cl_meta <- mock_cl_metadata(scdata)

  cell_id_barcode_map <- get_cell_id_barcode_map(scdata)
  cl_metadata <-
    make_cl_metadata_table(cl_meta, cell_id_barcode_map)

  res <- detect_variable_types(cl_metadata)

  expect_false("samples" %in% names(res))
})


test_that("detect_variable_types removes continuous variables", {
  scdata <- mock_scdata()
  cl_meta <- mock_cl_metadata(scdata)

  cell_id_barcode_map <- get_cell_id_barcode_map(scdata)
  cl_metadata <-
    make_cl_metadata_table(cl_meta, cell_id_barcode_map)

  res <- detect_variable_types(cl_metadata)

  expect_false("continuous_var" %in% names(res))
})


test_that("detect_variable_types removes high-cardinality variables", {
  scdata <- mock_scdata()
  cl_meta <- mock_cl_metadata(scdata)

  cell_id_barcode_map <- get_cell_id_barcode_map(scdata)
  cl_metadata <-
    make_cl_metadata_table(cl_meta, cell_id_barcode_map)

  res <- detect_variable_types(cl_metadata)

  # cells_id is explicitly excluded from cellset creation downstream, but it is
  # a high cardinality variable, so useful testing
  expect_false("cells_id" %in% names(res))
})


test_that("download_cl_metadata_file loads cl_metadata tables correctly", {
  config <- mock_config()
  scdata <- mock_scdata()
  cl_metadata <- mock_cl_metadata(scdata)

  local_mock_cl_metadata_table(cl_metadata, "mock_experiment_id")
  res <- stubbed_download_cl_metadata_file(config)
  withr::defer(unlink(file.path(".", basename(config$metadataS3Path))))

  expect_s3_class(res, "data.table")
  expect_named(res, names(cl_metadata))
  expect_equal(nrow(res), nrow(cl_metadata))
  expect_equal(ncol(res), ncol(cl_metadata))
})


test_that("make_cl_metadata_cellsets makes cell-level metadata cellsets.", {
  config <- mock_config()
  scdata <- mock_scdata()
  cl_metadata <- mock_cl_metadata(scdata)

  local_mock_cl_metadata_table(cl_metadata, "mock_experiment_id")

  res <- stubbed_make_cl_metadata_cellsets(scdata, config)
  withr::defer(unlink(file.path(".", basename(config$metadataS3Path))))

  expect_equal(length(res), 3)
  expect_equal(length(res[[1]]$children), length(unique(cl_metadata$cell_type)))
  expect_equal(length(res[[2]]$children), length(unique(cl_metadata$group_var)))
  expect_equal(length(res[[3]]$children), length(unique(cl_metadata$redundant_group_var)))

  cell_class_names <- c("key", "name", "rootNode", "type", "children")
  purrr::walk(res, expect_named, cell_class_names)

  # cellsets have the same keys as cell classes except children, color and cellIds
  for (i in seq_along(res)) {
    purrr::walk(res[[i]]$children, expect_named, c(cell_class_names[-5], "color", "cellIds"))
  }
})


with_fake_http(
  test_that("replace_cl_metadata_through_api sends patch request", {
    expect_PATCH(
      replace_cl_metadata_through_api(
        list(),
        "api_url",
        "experiment_id",
        "auth",
        FALSE
      )
    )
  })
)
