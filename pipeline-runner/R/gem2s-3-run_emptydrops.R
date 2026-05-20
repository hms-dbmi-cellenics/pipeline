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


  # ~3M droplets uses 15GB of RAM
  # calculate nworkers from sample with most cells
  # ensure will not exceed 85% of RAM on 60GB machine
  nworkers <- get_edrops_nworkers(counts_list)
  message("Number of workers: ", nworkers)

  edrops <- BiocParallel::bplapply(
    setNames(samples, samples),
    function(sample) {
      message("\nSample --> ", sample)
      compute_sample_edrops(counts_list[[sample]])
    },
    BPPARAM = BiocParallel::MulticoreParam(workers = nworkers)
  )

  edrops <- Filter(Negate(is.null), edrops)

  prev_out$edrops <- edrops
  res <- list(
    data = list(),
    output = prev_out
  )

  message("\nRunning of emptydrops step complete.")
  return(res)
}

get_edrops_nworkers <- function(counts_list) {
  sample_ndrops <- sapply(counts_list, ncol)
  nsamples <- length(counts_list)
  max_ndrops <- max(sample_ndrops)
  est_ram_per_worker <- (max_ndrops / 3e6) * 15
  max_workers <- floor(0.85 * BATCH_POD_MEMORY / est_ram_per_worker)

  # use at most ncpus workers
  min(nsamples, max_workers, BATCH_POD_CPUS)
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
  num_empty_drops <- sum(
    Matrix::colSums(sample_counts) < gem2s$max.empty.counts
  )

  if (num_empty_drops < gem2s$max.empty.drops) {
    message("Detected sample as filtered --> Skipping emptyDrops.")
    return(NULL)
  }

  if (methods::is(sample_counts, "IterableMatrix")) {
    # emptyDrops doesn't support IterableMatrix (bpcells)
    sample_counts <- as(sample_counts, "dgCMatrix")
  }

  set.seed(RANDOM_SEED)
  tryCatch({
    sample_edrops <- DropletUtils::emptyDrops(sample_counts)
    return(sample_edrops)
  }, error = function(e) {
    message(
      "Number of cells in sample too low for emptyDrops --> skipping."
    )
    return(NULL)
  })
}
