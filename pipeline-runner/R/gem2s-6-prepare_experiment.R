#' Prepare experiment for upload to AWS
#'
#'  1) Adds metadata: cellsId, color_pool, and gene annotation
#'  2) Prepares QC configuration
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

  check_names <- c("config", "edrops", "annot", "scdata_list", "disable_qc_filters")
  check_prev_out(prev_out, check_names)

  scdata_list <- prev_out$scdata_list
  disable_qc_filters <- prev_out$disable_qc_filters
  samples <- names(scdata_list)

  message("Total cells:", sum(sapply(scdata_list, ncol)))

  scdata_list <- add_metadata_to_samples(scdata_list, prev_out$annot, input$experimentId)
  prev_out$scdata_list <- scdata_list

  # construct default QC config and update prev out
  message("Constructing default QC configuration...")
  any_filtered <- !(length(prev_out$edrops) == length(samples))
  prev_out$qc_config <- construct_qc_config(scdata_list, any_filtered, disable_qc_filters)

  res <- list(
    data = list(),
    output = prev_out
  )

  message("\nPreparation for AWS upload step complete.")
  return(res)
}

add_metadata_to_samples <- function(scdata_list, annot, experiment_id) {
  # cell ids will be generated at random in order to shuffle samples. it is done
  # here because merging samples in QC means that shuffling the cells will not
  # result in a shuffled cell_ids
  set.seed(RANDOM_SEED)
  total_cells <- sum(sapply(scdata_list, ncol))
  cell_ids <- 0:(total_cells-1)

  # we need to iterate always in the same order because we are assigning the cell ids
  # at random and we want the relation cell <-> ID to be always the same to avoid
  # any changes in downstream analysis (like in the clustering)
  scdata_list <- order_by_size(scdata_list)
  for (sample in names(scdata_list)) {
    sample_size <- ncol(scdata_list[[sample]])

    # select only the annotations of the current sample
    sample_annotations_idx <- match(rownames(scdata_list[[sample]]), annot$input)
    sample_annot <- annot[sample_annotations_idx, ]
    scdata_list[[sample]]@misc[["gene_annotations"]] <- sample_annot

    # add the experiment ID so it's available later
    scdata_list[[sample]]@misc[["experimentId"]] <- experiment_id

    # sample cell ids to shuffle them
    idxs <- sample(seq_along(cell_ids), sample_size)
    scdata_list[[sample]]$cells_id <- cell_ids[idxs]
    # remove the selected cell ids for next samples
    cell_ids <- cell_ids[-idxs]
  }

  return(scdata_list)
}
