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
  counts_list <- prev_out$counts_list
  samples <- names(counts_list)

  edrops <- list()
  for (sample in samples) {
    message("\nSample --> ", sample)
    edrops[[sample]] <- compute_emptydrops(counts_list[[sample]])
  }

  prev_out$edrops <- edrops
  res <- list(
    data = list(),
    output = prev_out)

  message("\nStep 3 completed.")
  return(res)
}


#' @param counts dgCMatrix with counts for one sample.
compute_emptydrops <- function(counts) {
  # check if filtered
  nempty <- sum(Matrix::colSums(counts) < gem2s$max.empty.counts)

  if (nempty < gem2s$max.empty.drops) {
    message("Detected sample as filtered --> Skipping emptyDrops.")
    out <- NULL
  } else {
    out <- DropletUtils::emptyDrops(counts)
  }

  return(out)
}
