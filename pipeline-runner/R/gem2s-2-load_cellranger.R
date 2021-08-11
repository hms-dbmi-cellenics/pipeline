#' Read input folder of 10x data
#'
#' @inheritParams download_cellranger
#' @param prev_out list with experiment configuration settings
#'
#' @return list with 'output' slot containing \itemize{
#'   \item{"counts_list"}{named list of dgCMatrix per sample}
#'   \item{"annot"}{data.frame from reading features.tsv.gz for each sample}
#' }
#' @export
#'
load_cellranger <- function(input, pipeline_config, prev_out) {
  message("Loading cellranger output ...")

  # destructure previous output
  config <- prev_out$config

  output <- c(prev_out, call_read10x(config))

  res <- list(
    data = list(),
    output = output)

  message("\nStep 2 completed.")
  return(res)
}

#' Calls Read10X
#'
#' Cellranger outputs from V2 and V3 kits were renamed to look like V3 (features.tsv.gz).
#'
#' @param config experiment settings.
#'
call_read10x <- function(config) {
  counts_list <- list()
  annot_list <- list()

  samples <- config$samples
  message("Samples to include in the analysis:\n- ", paste(samples, collapse = "\n- "))
  message("Loading 10x data set from input folder.")


  for (sample in samples) {
    sample_dir <- file.path("/input", sample)
    sample_fpaths <- list.files(sample_dir)
    annot_fpath <- file.path(sample_dir, 'features.tsv.gz')

    message("\nSample --> ", sample)
    message("Reading files from ", sample_dir, ' --> ', paste(sample_fpaths, collapse = ' - '))

    counts <- Seurat::Read10X(sample_dir, gene.column = 1)
    annot <- read.delim(annot_fpath, header = FALSE)

    message(
      sprintf("Sample %s has %s genes and %s droplets.", sample, nrow(counts), ncol(counts))
    )

    counts_list[[sample]] <- counts
    annot_list[[sample]] <- annot
  }
  annot <- unique(do.call("rbind", annot_list))
  annot <- annot[, c(1, 2)]
  colnames(annot) <- c("input", "name")

  return(list(counts_list = counts_list, annot = annot))
}
