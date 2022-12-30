#' Run emptyDrops from DropletUtils
#'
#' @inheritParams download_user_files
#' @param prev_out 'output' slot from call to \code{load_user_files}
#'
#' @return \code{prev_out} with added slot 'edrops' containing result of call to
#'   \code{\link[DropletUtils]{emptyDrops}} for each sample.
#'
#' @export
#'
run_emptydrops <- function(input, pipeline_config, prev_out) {
  message("Testing if droplets are empty...")
  check_prev_out(prev_out, c("config", "counts_list", "annot"))

  # destructure previous output
  counts_list <- prev_out$counts_list
  samples <- names(counts_list)

  edrops <- list()
  for (sample in samples) {
    message("\nSample --> ", sample)
    edrops[[sample]] <- compute_sample_edrops(counts_list[[sample]])
  }

  prev_out$edrops <- edrops
  res <- list(
    data = list(),
    output = prev_out
  )

  message("\nRunning of emptydrops step complete.")
  return(res)
}


#' Calculate empty drops scores for sample
#'
#' @param sample_counts dgCMatrix with counts for one sample.
#'
#' @return data.frame with edrops scores
#' @export
#'
compute_sample_edrops <- function(sample_counts) {
  # check if filtered
  num_empty_drops <- sum(Matrix::colSums(sample_counts) < gem2s$max.empty.counts)

  if (num_empty_drops < gem2s$max.empty.drops) {
    message("Detected sample as filtered --> Skipping emptyDrops.")
    sample_edrops <- NULL
  } else {
    set.seed(RANDOM_SEED)
    sample_edrops <- DropletUtils::emptyDrops(sample_counts)
  }

  return(sample_edrops)
}
