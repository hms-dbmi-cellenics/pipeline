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

  technology <- input$input$type
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

    # TODO: Pass also parse_kit when available from the UI
    scores[[sample]] <- get_doublet_scores(sample_counts, technology = technology)

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
#' @param technology Technology used to generate the data (e.g. 10x, Parse).
#' @param parse_kit Kit used to generate the data (specific to Parse data).
#'
#' @return data.frame with doublet scores and assigned classes
#'
compute_sample_doublet_scores <- function(sample_counts, technology, parse_kit = "WT") {
  set.seed(RANDOM_SEED)

  dbr <- NULL
  if (technology == "parse") {
    dbr <- switch(parse_kit,
      "mini" = 0.046,
      "WT" = 0.034,
      "mega" = 0.064,
      stop("Invalid parse kit value: ", parse_kit)
    )
  }
  sce <- scDblFinder::scDblFinder(sample_counts, dbr = dbr)

  doublet_res <- data.frame(
    row.names = colnames(sce),
    barcodes = colnames(sce),
    doublet_class = sce$scDblFinder.class,
    doublet_scores = sce$scDblFinder.score
  )

  return(doublet_res)
}


get_doublet_scores <- function(sample_counts, max_attempts = 5, technology = "10x") {
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
      # TODO: Pass also parse_kit when available from the UI
      scores <- compute_sample_doublet_scores(sample_counts[, empty_cells_mask], technology)
      retry <- "not null"
    })
    attempt <- attempt + 1
  }

  return(scores)
}
