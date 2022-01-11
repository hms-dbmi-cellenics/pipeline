# IMPORTANT: these tests are also present in the worker. If update, change both.

mock_req <- function(type = "louvain") {
    req <- list(
        body =
            list(
                config = list(
                    resolution = 2
                ),
                type = type
            )
    )
}

mock_scdata <- function() {
    data("pbmc_small", package = "SeuratObject", envir = environment())
    pbmc_small$cells_id <- 0:(ncol(pbmc_small) - 1)
    pbmc_small@misc$gene_annotations <- data.frame(
        input = paste0("ENSG", seq_len(nrow(pbmc_small))),
        name = row.names(pbmc_small),
        row.names = paste0("ENSG", seq_len(nrow(pbmc_small)))
    )
    return(pbmc_small)
}

test_that("runClusters returns correct keys", {
    algos <- c("louvain", "leiden")
    data <- mock_scdata()
    expected_keys <- c("cluster", "cell_ids")

    for (algo in algos) {
        req <- mock_req(type = algo)
        res <- runClusters(req, data)
        expect_equal(names(res), expected_keys)
    }
})

test_that("runClusters returns one value per cell", {
    algos <- c("louvain", "leiden")
    data <- mock_scdata()
    expected_n_cells <- ncol(data)

    for (algo in algos) {
        req <- mock_req(type = algo)
        res <- runClusters(req, data)
        n_cells <- nrow(res)

        expect_equal(n_cells, expected_n_cells)
    }
})

test_that("runClusters orders barcodes correctly", {
    algos <- c("louvain", "leiden")
    data <- mock_scdata()
    expected_barcodes <- colnames(data)

    for (algo in algos) {
        req <- mock_req(type = algo)
        res <- runClusters(req, data)
        barcodes <- rownames(res)
        expect_equal(barcodes, expected_barcodes)
    }
})

test_that("runClusters returns at least one cluster", {
    algos <- c("louvain", "leiden")
    data <- mock_scdata()

    for (algo in algos) {
        req <- mock_req(type = algo)
        res <- runClusters(req, data)
        n_clusters <- length(unique(res$cluster))
        expect_gte(n_clusters, 1)
    }
})


test_that("runClusters uses active.reduction in misc slot", {

    algos <- c("louvain", "leiden")
    data <- mock_scdata()

    # get error if no PCA/SNN graph
    blah_reduction <- data@reductions$pca
    data@reductions$pca <- NULL
    data@graphs$RNA_snn <- NULL

    for (algo in algos) {
        req <- mock_req(type = algo)
        expect_error(runClusters(req, data), "Cannot find 'pca'")
    }

    # will use active.reduction to get SNN graph
    data@reductions$blah_reduction <- blah_reduction
    data@misc$active.reduction <- 'blah_reduction'

    for (algo in algos) {
        req <- mock_req(type = algo)
        expect_error(runClusters(req, data), NA)
    }
})
