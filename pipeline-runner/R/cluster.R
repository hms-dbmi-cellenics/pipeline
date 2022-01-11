# IMPORTANT: this function/utilities mirrored in worker. If update, change both.

#' Get Clusters
#'
#' @param clustering_method
#' @param resolution
#' @param data
#'
#' @return
#' @export
#'
#' @examples
runClusters <- function(req, data) {
    resol <- req$body$config$resolution
    type <- req$body$type

    data <- getClusters(type, resol, data)
    res_col <- paste0(data@active.assay, "_snn_res.", toString(resol))
    # In the meta data slot the clustering is stored with the resolution
    # used to calculate it
    # RNA_snn_res.#resolution
    df <- data.frame(
        cluster = data@meta.data[, res_col],
        cell_ids = data@meta.data$cells_id
    )
    # get the cell barcodes as rownames
    rownames(df) <- rownames(data@meta.data)
    return(df)
}

#' Compute clusters and return object with clusters
#'
#' @param algorithm
#' @param resolution
#' @param data
#'
#' @return
#'
#' @examples
getClusters <- function(type, resolution, data) {
    res_col <- paste0(data@active.assay, "_snn_res.", toString(resolution))
    algorithm <- list("louvain" = 1, "leiden" = 4)[[type]]

    # use the reduction from data integration for nearest neighbors graph
    if ("active.reduction" %in% names(data@misc)) {
        active.reduction <- data@misc[["active.reduction"]]
    } else {
        active.reduction <- "pca"
    }

    if (type == "leiden") {

        # emulate FindClusters, which overwrites seurat_clusters slot and meta.data column
        g <- getSNNiGraph(data, active.reduction)
        clus_res <- igraph::cluster_leiden(g, "modularity", resolution_parameter = resolution)
        clusters <- clus_res$membership
        names(clusters) <- clus_res$names
        clusters <- clusters[colnames(data)]
        data$seurat_clusters <- data@meta.data[, res_col] <- factor(clusters - 1)

    } else {

        graph.name <- paste0(Seurat::DefaultAssay(data), "_snn")
        if (!graph.name %in% names(data)) {
            data <- Seurat::FindNeighbors(data, annoy.metric = "cosine", verbose = FALSE, reduction = active.reduction)
        }
        data <- Seurat::FindClusters(data, resolution = resolution, verbose = FALSE, algorithm = algorithm)
    }

    return(data)
}


#' Get and Convert SNN Graph object into igraph object
#'
#' This is used to facilitate leiden clustering.
#'
#' @param data \code{Seurat} object
#'
#' @return boolean indicating if SNN Graph object exists
#'
getSNNiGraph <- function(data, active.reduction) {

    # check to see if we already have Seurat SNN Graph object
    snn_name <- paste0(data@active.assay, "_snn")

    # if doesn't exist, run SNN
    if (!snn_name %in% names(data)) data <- Seurat::FindNeighbors(data, reduction = active.reduction)

    # convert Seurat Graph object to igraph
    # similar to https://github.com/joshpeters/westerlund/blob/46609a68855d64ed06f436a6e2628578248d3237/R/functions.R#L85
    adj_matrix <- Matrix::Matrix(as.matrix(data@graphs[[snn_name]]), sparse = TRUE)
    g <- igraph::graph_from_adjacency_matrix(adj_matrix,
                                             mode = "undirected",
                                             weighted = TRUE
    )
    return(g)
}
