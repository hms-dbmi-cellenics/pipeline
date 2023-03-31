#' Compute doublet score for count matrices
#'
#' @inheritParams download_user_files
#' @param prev_out 'output' slot from call to \code{run_emptydrops}
#'
#' @return \code{prev_out} with added slot 'doublet_scores' containing result of call to
#'   \code{\link{compute_doublet_scores}} for each sample.
#'
#' @export
#'
score_doublets <- function(input, pipeline_config, prev_out) {
  message("Calculating probability of droplets being doublets...")

  # NOTE: edrops is not required
  check_prev_out(prev_out, c("config", "counts_list", "annot"))

  edrops_list <- prev_out$edrops
  counts_list <- prev_out$counts_list
  samples <- names(counts_list)

  scores <- list()
  for (sample in samples) {
    message("\nSample --> ", sample)
    sample_counts <- counts_list[[sample]]
    sample_edrops <- edrops_list[[sample]]

    # scDblFinder expects empty droplets to be removed
    if (!is.null(sample_edrops)) {
      keep <- which(sample_edrops$FDR <= gem2s$max.edrops.fdr)
      sample_counts <- sample_counts[, keep]
    }

    scores[[sample]] <- get_doublet_scores(sample_counts)

  }

  prev_out$doublet_scores <- scores
  res <- list(
    data = list(),
    output = prev_out
  )

  message("\nScoring doublets step complete.")
  return(res)
}


#' Compute doublets scores per sample.
#'
#' @param sample_counts Sparse matrix with the counts for one sample.
#'
#' @return data.frame with doublet scores and assigned classes
#'
compute_sample_doublet_scores <- function(sample_counts) {
  set.seed(RANDOM_SEED)
  sce <- scDblFinder::scDblFinder(sample_counts)
  doublet_res <- data.frame(
    row.names = colnames(sce),
    barcodes = colnames(sce),
    doublet_class = sce$scDblFinder.class,
    doublet_scores = sce$scDblFinder.score
  )

  return(doublet_res)
}


get_doublet_scores <- function(sample_counts, max_attempts = 5) {
  # also filter low UMI as per scDblFinder:::.checkSCE()
  ntot <- Matrix::colSums(sample_counts)

  # retry increasing the minimum counts in case of low sparsity in the sample
  retry <- NULL
  attempt <- 1
  while (is.null(retry) && attempt <= max_attempts) {
    message("\nTrying to score doublets, attempt: ", attempt)
    # make the threshold stricter in every attempt
    empty_cells_mask <- ntot > (200 * attempt)
    try({
      scores <- compute_sample_doublet_scores(sample_counts[, empty_cells_mask])
      retry <- "not null"
    })
    attempt <- attempt + 1
  }

  return(scores)
}
