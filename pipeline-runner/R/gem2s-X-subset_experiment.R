#' Extract cellset type as data.frame
#'
#' Gets the cellsets list and converts it to tidy tibble, keeping only the
#' the required cellset type
#'
#' @param cellsets
#' @param cellset_type
#'
#' @return
#' @export
#'
parse_cellsets <- function(cellsets, cellset_type) {

  cellsets$cellSets |>
    dplyr::filter(key == cellset_type) %>%
    .$children |>
    as.data.frame() |>
    tibble::as_tibble() |>
    dplyr::select(key, name, cellIds)
}


#' Filters cellsets, getting vector of cell_ids to keep
#'
#' @param cellsets_df data.frame
#' @param cellset_keys character
#'
#' @return
#' @export
#'
get_cell_ids <- function(cellsets_df, cellset_keys) {

  cellsets_df |>
    dplyr::filter(key %in% cellset_keys) |>
    tidyr::unnest(cellIds) |>
    dplyr::pull(cellIds) |>
    unique()

}



create_subset_experiment <- function(input, pipeline_config) {

  parent_experiment_id <- input$parentExperimentId
  subset_experiment_id <- input$subsetExperimentId
  cellset_keys <- input$cellSetKeys

  # load parent processed scdata and cellsets
  s3 <- paws::s3(config = pipeline_config$aws_config)
  parent_scdata <- load_processed_scdata(s3, pipeline_config, parent_experiment_id)
  parent_cellsets <- load_cellsets(s3, pipeline_config, parent_experiment_id)

  cell_ids_to_keep <- get_cell_sets(parent_cellsets, cellset_keys)

  sample_id_mapping <- input$sampleIdMapping

  # subset seurat object
  scdata <- subset_ids(scdata, cell_ids_to_keep)

  # add subset experiment name to the subset seurat object
  scdata$project <- input$name

  # add new sample_ids, keep originals in a new variable
  scdata$parent_samples <- scdata$samples
  scdata$samples <- sample_id_mapping[match(parent_samples, sample_id_mapping)]

  # split by sample
  scdata_list <- Seurat::SplitObject(scdata, split.by = "samples")

  prev_out$scdata_list <- scdata_list
  prev_out$annot <- scdata@misc
  res <- list(
    data = list(),
    output = prev_out
  )

  message("\nSubsetting of Seurat object step complete.")
  return(res)
}
