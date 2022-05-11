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

  message("Merging Seurat Objects...")
  # saveRDS(scdata_list, "/debug/scdata_list.rds")
  # saveRDS(prev_out, "/debug/prev_out.rds")
  # saveRDS(pipeline_config, "/debug/pipeline_config.rds")
  # sum(sapply(scdata_list, ncol))r
  # scdata <- merge_scdatas(scdata_list)
  # scdata_list <- add_metadata(scdata_list, prev_out$annot, input$experimentId)
  prev_out$scdata <- scdata_list

  # construct default QC config and update prev out
  message("Constructing default QC configuration...")
  any_filtered <- !(length(prev_out$edrops) == length(samples))
  # prev_out$qc_config <- construct_qc_config(scdata_list, any_filtered)

  res <- list(
    data = list(),
    output = prev_out
  )

  message("\nPreparation for AWS upload step complete.")
  return(res)
}

merge_scdatas <- function(scdata_list) {
  if (length(scdata_list) == 1) {
    scdata <- scdata_list[[1]]
  } else {
    scdata <- merge(scdata_list[[1]], y = scdata_list[-1])
  }

  return(scdata)
}

add_metadata <- function(scdata, annot, experiment_id) {

  # Ensure index by rownames in scdata
  annot <- annot[match(rownames(scdata), annot$input), ]
  scdata@misc[["gene_annotations"]] <- annot

  message("Storing cells id...")
  # Keeping old version of ids starting from 0
  scdata$cells_id <- 0:(ncol(scdata) - 1)

  message("Storing color pool...")
  # We store the color pool in a slot in order to be able to access it during configureEmbedding
  scdata@misc[["color_pool"]] <- get_color_pool()
  scdata@misc[["experimentId"]] <- experiment_id
  scdata@misc[["ingestionDate"]] <- Sys.time()

  return(scdata)
}
