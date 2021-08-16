#' Run emptyDrops from DropletUtils
#'
#' @inheritParams download_cellranger
#' @param prev_out 'output' slot from call to \code{load_cellranger}
#'
#' @return \code{prev_out} with added slot 'edrops' containing result of call to
#'   \code{\link[DropletUtils]{emptyDrops}} for each sample.
#'
#' @export
#'
run_emptydrops <- function(input, pipeline_config, prev_out) {
  message("Testing if droplets are empty...")

  # destructure previous output
  check_prev_out(prev_out, 'counts_list')
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
    output = prev_out)

  message("\nRunning of emptydrops step complete.")
  return(res)
}


#' @param sample_counts dgCMatrix with counts for one sample.
compute_sample_edrops <- function(sample_counts) {
  # check if filtered
  nempty <- sum(Matrix::colSums(sample_counts) < gem2s$max.empty.counts)

  if (nempty < gem2s$max.empty.drops) {
    message("Detected sample as filtered --> Skipping emptyDrops.")
    sample_edrops <- NULL
  } else {
    sample_edrops <- DropletUtils::emptyDrops(sample_counts)
  }

  return(sample_edrops)
}
