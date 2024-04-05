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

determine_dynamic_threshold <- function(sample_counts, adjust = 2, min_distance_from_peak = 0.5) {
  umi_counts <- Matrix::colSums(sample_counts)
  log_umi_counts <- log10(umi_counts + 1)

  # Calculate the kernel density estimation on the log10-transformed data with smoothing
  density_est <- density(log_umi_counts, adjust = adjust)
  valleys <- find_valleys(density_est, min_distance_from_peak)
  thresholds <- 10^valleys - 1
  return(min(thresholds))
}

find_valleys <- function(density_est, min_distance_from_peak) {
  gradient <- diff(density_est$y)
  sign_change <- diff(sign(gradient))
  potential_valleys <- which(sign_change == 2) + 1

  peak_location <- which.max(density_est$y)
  peak_x_value <- density_est$x[peak_location]

  # Identify valid valleys based on distance from peak
  valid_valleys <- potential_valleys[abs(density_est$x[potential_valleys] - peak_x_value) > min_distance_from_peak]

  return(density_est$x[valid_valleys])
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
    return(NULL)
  }

  set.seed(RANDOM_SEED)
  tryCatch({
      dynamic_threshold <- determine_dynamic_threshold(sample_counts)
      message("EmptyDrops lower = ", dynamic_threshold)
      sample_edrops <- DropletUtils::emptyDrops(sample_counts, lower = dynamic_threshold)
      return(sample_edrops)
  }, error = function(e) {
      message("Number of cells in sample too low for emptyDrops --> Skipping emptyDrops.")
      return(NULL)
  })
}
