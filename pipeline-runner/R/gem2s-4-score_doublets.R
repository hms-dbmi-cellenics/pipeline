#' Compute doublet score for count matrices
#'
#' @inheritParams download_cellranger
#' @param prev_out 'output' slot from call to \code{run_emptydrops}
#'
#' @return \code{prev_out} with added slot 'doublet_scores' containing result of call to
#'   \code{\link{compute_doublet_scores}} for each sample.
#'
#' @export
#'
score_doublets <- function(input, pipeline_config, prev_out) {
  message("Calculating probability of droplets being doublets...")

  counts_list <- prev_out$counts_list
  samples <- names(counts_list)

  scores <- list()
  for (sample in samples) {
    message("\tSample --> ", sample, "...")
    scores[[sample]] <- compute_doublet_scores(counts_list[[sample]])
  }

  prev_out$doublet_scores <- scores
  res <- list(
    data = list(),
    output = prev_out)

  message("\tStep 4 completed.")
  return(res)
}


#' Compute doublets scores per sample.
#'
#' @param scdata Raw sparse matrix with the counts for one sample.
#'
#' @return data.frame with doublet scores and assigned classes
#'
compute_doublet_scores <- function(counts) {

  set.seed(0)
  sce <- scDblFinder::scDblFinder(counts)
  dbl.df <- data.frame(
    row.names = colnames(sce),
    barcodes = colnames(sce),
    doublet_class = sce$scDblFinder.class,
    doublet_scores = sce$scDblFinder.score
  )

  return(dbl.df)
}
