# #' create a subset experiment
# #'
# #' This is the first step of a subset pipeline, which basically takes the parent
# #' experiment ID and cellset keys to keep as input, extracts the cell ids to keep
# #' and subsets and slims down the parent seurat object.
# #'
# #' @param input list containing:
# #'   - parentExperimentId character
# #'   - subsetExperimentId character
# #'   - cellSetKeys character vector of cellset keys to subset
# #'   - experimentName character
# #'   - parentProcessingConfig The processingConfig of the parent experiment
# #' @param pipeline_config list
# #' @param prev_out list, ignored because this is the first step in the subset pipeline
# #'
# #' @return list containing scdata_list, annotations and sample_id_map
# #' @export
# #'
# subset_seurat <- function(input, pipeline_config, prev_out = NULL) {
#   parent_data <- load_parent_experiment_data(input, pipeline_config)

#   subset_scdata <- subset_experiment(input, parent_data)
#   sample_id_map <- create_sample_id_map(unique(subset_scdata$samples))
#   subset_scdata <- add_subset_metadata(input, subset_scdata, sample_id_map)

#   subset_scdata_list <- Seurat::SplitObject(subset_scdata, split.by = "samples")

#   # TODO: remove from here and refactor all pipeline.
#   config <- list(
#     name = input$experimentName,
#     samples = sample_id_map$subset_sample_id
#   )

#   # structure step output
#   res <- list(
#     data = list(
#       sampleIdMap = sample_id_map
#       ),
#     output = list(
#       scdata_list = subset_scdata_list,
#       annot = subset_scdata@misc$gene_annotations,
#       edrops = NULL,
#       sample_id_map = sample_id_map,
#       config = config,
#       disable_qc_filters = TRUE,
#       parent_cellsets = parent_data$cellsets,
#       qc_config = generate_subset_config(input$parentProcessingConfig, sample_id_map)
#     )
#   )

#   message("\nSubsetting of Seurat object step complete.")
#   return(res)
# }

copy_s3_objects <- function(input, pipeline_config, prev_out = NULL) {
  parent_data <- load_from_experiment_data(input, pipeline_config)

  parent_scdata <- parent_data$parent_scdata
  parent_cellsets <- parent_data$parent_cellsets

  
}


#' load parent experiment data
#'
#' Loads the processed rds and cellsets file from the parent experiment from s3.
#'
#' @param input list of input parameters
#' @param pipelne_config list of pipeline parameters
#'
#' @return list with scdata and parsed cellsets
#' @export
#'
load_from_experiment_data <- function(input, pipeline_config) {

  str("inputDebug")
  str(input)

  message("inputDebug")
  message(input)

  # load parent processed scdata and cellsets
  s3 <- paws::s3(config = pipeline_config$aws_config)
  parent_scdata <-
    load_processed_scdata(s3, pipeline_config, input$fromExperimentId)
  parent_cellsets <-
    parse_cellsets(load_cellsets(s3, pipeline_config, input$fromExperimentId))

  str("inputsampleIdsMapDebug")
  str(input$sampleIdsMap)

  return(list(scdata = parent_scdata, cellsets = parent_cellsets))
}

#' Add new sample ids to the subset Seurat Object
#'
#' @param scdata Seurat Object
#' @param sample_id_map data.table of parent/subset sample id map
#'
#' @return SeuratObject with new sample ids
#' @export
#'
add_new_sample_ids <- function(scdata, sample_id_map) {
  sample_map_idx <- match(subset_scdata$samples, names(sample_id_map))
  str("sample_map_idxDebug")
  str(sample_map_idx)
  subset_scdata$samples <- unname(unlist(sample_id_map[sample_map_idx]))
  return(subset_scdata)
}
