#' Prepare experiment for upload to AWS
#'
#'  1) Merges the samples for the current experiment
#'  2) Adds metadata: cellsId, color_pool, and gene annotation
#'  3) Preparing QC configuration
#'
#' @inheritParams download_cellranger
#' @param prev_out  'output' slot from call to \code{create_seurat}
#'
#' @return prev_out \code{prev_out} with added slots 'scdata' containing merged
#'   \code{SeuratObject} and 'qc_config' containing default config for QC steps.
#'
#' @export
#'
prepare_experiment <- function(input, pipeline_config, prev_out) {
  message("Preparing experiment ...")
  check_prev_out(prev_out, c('scdata_list', 'annot'))

  scdata_list <- prev_out$scdata_list
  edrops <- prev_out$edrops
  samples <- names(scdata_list)

  message("Merging Seurat Objects...")
  if (length(scdata_list) == 1) {
    scdata <- scdata_list[[1]]
  } else {
    scdata <- merge(scdata_list[[1]], y = scdata_list[-1])
  }

  message("Deduplicating gene annotations...")
  annot <- prev_out$annot

  # add ENSEMBL ID for genes that are duplicated (geneNameDuplicated-ENSEMBL)
  # original name kept in 'original_name' column
  gname <- annot$name
  annot$original_name <- gname
  is.dup <- duplicated(gname) | duplicated(gname, fromLast = TRUE)

  #We need to convert the gene inputs from _ to - bc when we create the Seurat object we do this, and the match would return NA values if any of the inputs still has _.
  annot$input <- gsub('_', '-', annot$input)
  annot$name[is.dup] <- paste(gname[is.dup], annot$input[is.dup], sep = " - ")

  # Ensure index by rownames in scdata
  annot <- annot[match(rownames(scdata), annot$input), ]
  rownames(annot) <- annot$input

  scdata@misc[["gene_annotations"]] <- annot

  message("Storing cells id...")
  # Keeping old version of ids starting from 0
  scdata$cells_id <- 0:(ncol(scdata) - 1)

  message("Storing color pool...")
  # We store the color pool in a slot in order to be able to access it during configureEmbedding
  scdata@misc[["color_pool"]] <- get_color_pool()
  scdata@misc[["experimentId"]] <- input$experimentId
  scdata@misc[["ingestionDate"]] <- Sys.time()

  # construct default QC config and update prev out
  message("Constructing default QC configuration...")
  any_filtered <- !(length(edrops) == length(samples))
  prev_out$scdata <- scdata
  prev_out$qc_config <- construct_qc_config(scdata, any_filtered)

  res <- list(
    data = list(),
    output = prev_out)

  message("\nPreperation for AWS upload step complete.")
  return(res)
}

