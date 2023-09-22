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
  redundant_group_var <- paste0("red_group_", samples)
  continuous_var <- rnorm(length(barcode))

  data.table::data.table(barcode, samples, cell_type, group_var, redundant_group_var, continuous_var)
}


mock_config <- function() {
  config <-
    list(clusteringSettings = list(
      method = "louvain",
      methodSettings = list(louvain = list(resolution = "0.8"))
    ),
    aws_config = list(aws = "mock_aws"),
    metadataS3Path = "mock_s3_path"
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
  test_that("update_clusters_through_api sends patch request", {
    expect_PATCH(
      update_clusters_through_api(
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
 test_that("update_clusters_through_api diables SSL certificate checking if disabled", {

    expect_PATCH(
      update_clusters_through_api(
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


test_that("make_cl_metadata_table joins on samples and barcode if samples column in user-supplied cl metadata", {
  scdata <- mock_scdata()
  cl_meta <- mock_cl_metadata(scdata)

  cell_id_barcode_map <- get_cell_id_barcode_map(scdata)
  res <- make_cl_metadata_table(cl_meta, cell_id_barcode_map)

  expect_true("samples" %in% names(res))
})


test_that("make_cl_metadata_table joins barcode only if no samples column in user-supplied cl metadata", {
  scdata <- mock_scdata()
  cl_meta <- mock_cl_metadata(scdata)

  # remove samples column from cell-level metadata
  cl_meta[, samples:= NULL]

  cell_id_barcode_map <- get_cell_id_barcode_map(scdata)
  res <- make_cl_metadata_table(cl_meta, cell_id_barcode_map)

  expect_true(!"samples" %in% names(res))
})


test_that("make_cl_metadata_cellset makes correctly formatted cellsets", {

})
