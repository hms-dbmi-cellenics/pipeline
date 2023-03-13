#' Perform geometric sketching
#'
#' See https://github.com/brianhie/geosketch
#'
#' @param object Seurat object
#' @param dims Number of dimensions of PC space to use
#' @param num_cells Number of desired cells
#'
#' @return Seurat object downsampled to desired number of cells
#' @export
#'

run_geosketch <- function(scdata, dims, perc_num_cells, reduction = "pca") {

  num_cells <- round(ncol(scdata) * perc_num_cells / 100)

  # use the min of what the user wants and what is available
  dims <- min(ncol(scdata@reductions[[reduction]]) - 1, dims)

  message("Geosketching to ", num_cells, " cells")

  if (!exists("geosketch")) {
    geosketch <- reticulate::import("geosketch")
  }
  stopifnot(
    "The requested reduction is not present in the Seurat object." = reduction %in% names(scdata@reductions),
    "The number of cells is lower that the number of dimensions." = ncol(scdata@reductions[[reduction]]) >= dims
  )

  embeddings <- scdata@reductions[[reduction]]@cell.embeddings[, 1:dims]
  index <- unlist(geosketch$gs(embeddings, as.integer(num_cells), one_indexed = TRUE))
  sketch <- scdata[, index]
  Seurat::DefaultAssay(sketch) <- "RNA"
  sketch@misc[["active.reduction"]] <- reduction

  return(list(scdata = scdata, sketch = sketch))
}


#' Learn integration transformation from sketches
#'
#' Uses the integrated sketches to learn the integration transformation and
#' apply it to the whole dataset
#'
#' @param scdata Seurat object
#' @param scdata_sketch Sketched Seurat object
#' @param scdata_sketch_integrated Sketched integrated Seurat object
#' @param method Reduction method
#' @param dims Number of dimensions of PC space to use
#'
#' @return Integrated Seurat object with original number of cells
#' @export
#'
learn_from_sketches <- function(scdata, scdata_sketch, scdata_sketch_integrated, dims, reduction = "pca") {
  # get embeddings from splitted Seurat object
  active_reduction <- scdata_sketch_integrated@misc[["active.reduction"]]
  embeddings_orig <- list(scdata@reductions[[reduction]]@cell.embeddings[, 1:dims])
  embeddings_sketch <- list(scdata_sketch@reductions[[reduction]]@cell.embeddings[, 1:dims])
  embeddings_sketch_int <- list(scdata_sketch_integrated@reductions[[active_reduction]]@cell.embeddings[, 1:dims])

  # use python script to learn integration from sketches and apply to whole dataset
  geosketch_script_path <- "/src/pipeline-runner/inst/python/learn-apply-transformation.py"

  reticulate::source_python(geosketch_script_path)
  learned_int <- apply_transf(embeddings_orig, embeddings_sketch, embeddings_sketch_int)
  rownames(learned_int[[1]]) <- colnames(scdata)

  scdata[[active_reduction]] <- Seurat::CreateDimReducObject(
    embeddings = learned_int[[1]],
    key = paste0(active_reduction, "_"),
    assay = Seurat::DefaultAssay(scdata)
  )

  scdata@misc[["active.reduction"]] <- active_reduction
  scdata@misc[["geosketch"]] <- TRUE

  return(scdata)
}
