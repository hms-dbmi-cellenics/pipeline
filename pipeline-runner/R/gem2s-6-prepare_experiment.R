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
  message(sum(sapply(scdata_list, ncol)))

  # scdata <- merge_scdatas(scdata_list)

  #If subsetting all cells, Seurat will not reorder the cells in the object. We need to subset to [-1] and [1] and merge to shuffle.
  set.seed(gem2s$random.seed)
  # shuffle_mask <- sample(colnames(scdata))
  # scdata <- merge(scdata[,shuffle_mask[1]],scdata[,shuffle_mask[-1]])

  # scdata <- add_metadata(scdata, prev_out$annot, input$experimentId)
  # prev_out$scdata <- scdata
  scdata_list <- add_metadata_to_each(scdata_list, prev_out$annot, input$experimentId)
  # saveRDS(scdata_list, '/debug/scdata_list.metadata.rds')
  prev_out$scdata_list <- scdata_list

  # construct default QC config and update prev out
  message("Constructing default QC configuration...")
  any_filtered <- !(length(prev_out$edrops) == length(samples))
  # NEEDED ?
  # scdata_list$samples <- samples
  prev_out$qc_config <- construct_qc_config(scdata_list, any_filtered)

  res <- list(
    data = list(),
    output = prev_out
  )

  message("\nPreparation for AWS upload step complete.")
  return(res)
}

#' Merge scdatas: merge rds files in input list
#'
#' @param scdata_list
#'
#' @return
#' @export
#'
#' @examples
merge_scdatas <- function(scdata_list) {
  if (length(scdata_list) == 1) {
    scdata <- scdata_list[[1]]
  } else {
    scdata <- merge(scdata_list[[1]], y = scdata_list[-1])
  }

  return(scdata)
}

add_metadata_to_each <- function(scdata_list, annot, experiment_id) {

  message("add_metadata_to_each")
  message("names(scdata_list): ", names(scdata_list))
  for (sample in names(scdata_list)) {
  # for (scdata in scdata_list) {
    # Ensure index by rownames in scdata
    annot <- annot[match(rownames(scdata_list[[sample]]), annot$input), ]
    scdata_list[[sample]]@misc[["gene_annotations"]] <- annot

    message("Storing cells id...")
    # Keeping old version of ids starting from 0
    scdata_list[[sample]]$cells_id <- 0:(ncol(scdata_list[[sample]]) - 1)

    message("Storing color pool...")
    # We store the color pool in a slot in order to be able to access it during configureEmbedding
    scdata_list[[sample]]@misc[["color_pool"]] <- get_color_pool()
    scdata_list[[sample]]@misc[["experimentId"]] <- experiment_id
    scdata_list[[sample]]@misc[["ingestionDate"]] <- Sys.time()
  }

  return(scdata_list)
}
