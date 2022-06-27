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
#'
myFun <- function(n = 5000) {
    a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
    paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

prepare_experiment <- function(input, pipeline_config, prev_out) {
  message("Preparing experiment ...")
  check_names <- c("config", "counts_list", "annot", "doublet_scores", "scdata_list")
  check_prev_out(prev_out, check_names)

  scdata_list <- prev_out$scdata_list
  samples <- names(scdata_list)

  message("Merging Seurat Objects...")
  x <- myFun(100000000) # 7.2GB
  Sys.sleep(10)
  message("Merging Seurat Objects...2")
  x <- myFun(100000000) # 7.2GB
  Sys.sleep(10)
  message("Merging Seurat Objects...3")
  x <- myFun(100000000) # 7.2GB
  Sys.sleep(10)
  message("Merging Seurat Objects...4")
  x <- myFun(100000000) # 7.2GB
  Sys.sleep(10)
  message("Merging Seurat Objects...5")
  x <- myFun(100000000) # 7.2GB
  Sys.sleep(10)
  message("Merging Seurat Objects...6")
  x <- myFun(100000000) # 7.2GB
  scdata <- merge_scdatas(scdata_list)

  #If subsetting all cells, Seurat will not reorder the cells in the object. We need to subset to [-1] and [1] and merge to shuffle.
  set.seed(gem2s$random.seed)
  shuffle_mask <- sample(colnames(scdata))
  scdata <- merge(scdata[,shuffle_mask[1]],scdata[,shuffle_mask[-1]])

  scdata <- add_metadata(scdata, prev_out$annot, input$experimentId)
  prev_out$scdata <- scdata

  # construct default QC config and update prev out
  message("Constructing default QC configuration...")
  any_filtered <- !(length(prev_out$edrops) == length(samples))
  prev_out$qc_config <- construct_qc_config(scdata, any_filtered)

  res <- list(
    data = list(),
    output = prev_out
  )

  message("\nPreperation for AWS upload step complete.")
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
