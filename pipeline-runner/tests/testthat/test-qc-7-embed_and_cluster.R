mock_req <- function(type = "louvain") {
  req <- list(body =
                list(
                  config = list(
                    resolution = 2,
                    apiUrl = "mock_api_url",
                    authJwt = "mock_auth"

                  ),
                  type = type
                ))
}

mock_scdata <- function() {
  data("pbmc_small", package = "SeuratObject", envir = environment())
  pbmc_small$cells_id <- 0:(ncol(pbmc_small) - 1)
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
  test_that("update_sets_through_api sends patch request", {
    expect_PATCH(
      update_sets_through_api(
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
 test_that("update_sets_through_api diables SSL certificate checking if disabled", {

    expect_PATCH(
      update_sets_through_api(
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

test_that("format_cell_sets_object returns correct items", {
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
    res <- format_cell_sets_object(cell_sets, algo, color_pool)

    # expect all keys present
    expect_setequal(names(res), names(types))

    for (item in names(types)) {
      expect_type(res[[item]], types[[item]])
    }
  }
})


test_that("format_cell_sets_object correctly formats a cellset object", {

  n_clusters <- 5
  cell_sets <- mock_cellset_object(1000, n_clusters)
  color_pool <- mock_color_pool(n_clusters)
  algos <- c("louvain", "leiden")

  expected_items <- c("key", "name", "rootNode", "type", "color", "cellIds")

  for (algo in algos) {
    res <- format_cell_sets_object(cell_sets, algo, color_pool)


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


test_that("format_cell_sets_object result has correct number of clusters",{
  n_clusters <- c(0, 1, 4, 6)
  algos <- c("louvain", "leiden")

  for (algo in algos) {
    for (n in n_clusters) {
      cell_sets <- mock_cellset_object(100, n)
      color_pool <- mock_color_pool(n)

      res <- format_cell_sets_object(cell_sets, algo, color_pool)

      # number of clusters is the number of children elements
      expect_equal(length(res$children), n)

    }}
})

test_that("format_cell_sets_object returns empty children on empty cellset", {
  n_clusters <- 5
  cell_sets <- mock_cellset_object(0, n_clusters)

  color_pool <- mock_color_pool(n_clusters)
  algos <- c("louvain", "leiden")

  for (algo in algos) {
    res <- format_cell_sets_object(cell_sets, algo, color_pool)
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



