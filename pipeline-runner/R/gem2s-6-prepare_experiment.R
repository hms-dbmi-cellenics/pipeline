#' Prepare experiment for upload to AWS
#'
#'  1) Merges the samples for the current experiment
#'  2) Adds metadata: cellsId, color_pool, and gene annotation
#'  3) Preparing QC configuration
#'
#' @inheritParams download_user_files
#' @param prev_out  'output' slot from call to \code{create_seurat}
#'
#' @return prev_out \code{prev_out} with added slots 'scdata' containing merged
#'   \code{SeuratObject} and 'qc_config' containing default config for QC steps.
#'
#' @export
#'
prepare_experiment <- function(input, pipeline_config, prev_out) {
  message("Preparing experiment ...")

  check_names <- c("config", "counts_list", "annot", "doublet_scores", "scdata_list")
  check_prev_out(prev_out, check_names)

  scdata_list <- prev_out$scdata_list
  samples <- names(scdata_list)

  message("Total cells:", sum(sapply(scdata_list, ncol)))

  scdata_list <- add_metadata_to_samples(scdata_list, prev_out$annot, input$experimentId)
  prev_out$scdata_list <- scdata_list

  # construct default QC config and update prev out
  message("Constructing default QC configuration...")
  any_filtered <- !(length(prev_out$edrops) == length(samples))
  prev_out$qc_config <- construct_qc_config(scdata_list, any_filtered)

  res <- list(
    data = list(),
    output = prev_out
  )

  message("\nPreparation for AWS upload step complete.")
  return(res)
}

# TODO add description
add_metadata_to_samples <- function(scdata_list, annot, experiment_id) {
  # cell_ids need to be consecutive among the samples and each sample can have different sizes
  # so we need to keep track of start IDs
  start <- 0
  for (sample in names(scdata_list)) {
    end <- start + (ncol(scdata_list[[sample]]) - 1)

    # select only the annotations of the current sample
    annot <- annot[match(rownames(scdata_list[[sample]]), annot$input), ]
    scdata_list[[sample]]@misc[["gene_annotations"]] <- annot

    # add the experiment ID so it's available later
    scdata_list[[sample]]@misc[["experimentId"]] <- experiment_id

    # generate sample cell IDs
    scdata_list[[sample]]$cells_id <- start:end
    start <- end + 1
  }

  return(scdata_list)
}
