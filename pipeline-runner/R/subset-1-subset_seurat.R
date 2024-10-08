library(uuid)

#' Translate processing config
#'
#' Auxiliary function for subset_seurat
#'
#' @param parentProcessingConfig The processingConfig of the parent experiment
#' @param sample_id_map Maps parent sample ids to the subset sample ids
#'  (only the parent samples that have cells after the subset are included)
#'
#' @return Processing config for the subset experiment
#'
generate_subset_config <- function(parent_processing_config, sample_id_map) {
  steps <- names(parent_processing_config)
  # dataIntegration and configureEmbedding don't have sample-based config
  # so we don't need to run over it
  steps_to_filter <- steps[!steps %in% c("dataIntegration", "configureEmbedding")]

  processing_config <- parent_processing_config

  for (step_key in steps_to_filter) {
    new_step_config <- processing_config[[step_key]]

    # Keep only the config of samples that survived the subset
    new_step_config <- new_step_config[names(new_step_config) %in% names(sample_id_map)]

    # Translate each of the sample ids to their subset counterpart
    for (sample_id in names(new_step_config)) {
      subset_sample_id <- sample_id_map[[sample_id]]

      # Set config on new sample id
      new_step_config[[subset_sample_id]] <- new_step_config[[sample_id]]

      # All filters on subset exps are disabled by default
      new_step_config[[subset_sample_id]]$enabled <- FALSE

      # Clean up old sample config
      new_step_config[[sample_id]] <- NULL
    }

    # Replace old step config with the new step config we created
    processing_config[[step_key]] <- new_step_config
  }

  return(processing_config)
}

#' create a subset experiment
#'
#' This is the first step of a subset pipeline, which basically takes the parent
#' experiment ID and cellset keys to keep as input, extracts the cell ids to keep
#' and subsets and slims down the parent seurat object.
#'
#' @param input list containing:
#'   - parentExperimentId character
#'   - subsetExperimentId character
#'   - cellSetKeys character vector of cellset keys to subset
#'   - experimentName character
#'   - parentProcessingConfig The processingConfig of the parent experiment
#' @param pipeline_config list
#' @param prev_out list, ignored because this is the first step in the subset pipeline
#'
#' @return list containing scdata_list, annotations and sample_id_map
#' @export
#'
subset_seurat <- function(input, pipeline_config, prev_out = NULL) {
  parent_data <- load_parent_experiment_data(input, pipeline_config)

  subset_scdata <- subset_experiment(input, parent_data)
  sample_id_map <- create_sample_id_map(unique(subset_scdata$samples))
  subset_scdata <- add_subset_metadata(input, subset_scdata, sample_id_map)

  subset_scdata_list <- Seurat::SplitObject(subset_scdata, split.by = "samples")

  # TODO: remove from here and refactor all pipeline.
  config <- list(
    name = input$experimentName,
    samples = sample_id_map$subset_sample_id
  )

  # structure step output
  res <- list(
    data = list(
      sampleIdMap = sample_id_map
      ),
    output = list(
      scdata_list = subset_scdata_list,
      annot = subset_scdata@misc$gene_annotations,
      edrops = NULL,
      sample_id_map = sample_id_map,
      config = config,
      disable_qc_filters = TRUE,
      parent_cellsets = parent_data$cellsets,
      qc_config = generate_subset_config(input$parentProcessingConfig, sample_id_map)
    )
  )

  message("\nSubsetting of Seurat object step complete.")
  return(res)
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
load_parent_experiment_data <- function(input, pipeline_config) {
  # load parent processed scdata and cellsets
  s3 <- paws::s3(config = pipeline_config$aws_config)
  parent_scdata <- load_processed_scdata(s3, pipeline_config, input$parentExperimentId)
  parent_cellsets <- parse_cellsets(load_cellsets(s3, pipeline_config, input$parentExperimentId))

  return(list(scdata = parent_scdata, cellsets = parent_cellsets))
}


#' Remove all unnecessary data from the parent seurat object
#'
#' Seurat::DietSeurat is not able to remove certain slots from a seurat object.
#' This function also removes elements from the misc slot which are not necessary
#'
#' @param scdata SeuratObject
#'
#' @return leaner SeuratObject
#' @export
#'
diet_scdata <- function(scdata) {
  lean_scdata <- Seurat::CreateSeuratObject(
    counts = scdata@assays$RNA$counts,
    meta.data = scdata@meta.data,
    min.cells = 0,
    min.features = 0
  )

  lean_scdata@misc <- list(
    gene_annotations = scdata@misc$gene_annotations,
    parent_experimentId = scdata@misc$experimentId
  )

  return(lean_scdata)
}


#' Subset seurat object by the input cellset keys
#'
#' This function takes the cellset keys sent by the API, extracts the cell_ids
#' that belong to them, subsets the seurat object and removes all unnecessary
#' data from it.
#'
#' @param input list of input parameters, containing cellSetKeys to subset
#' @param parent list containing parent scdata and parsed cellsets
#'
#' @return subset seurat object
#' @export
#'
subset_experiment <- function(input, parent_data) {
  # subset seurat object, remove unnecesary data
  cell_ids_to_keep <- unique(parent_data$cellsets[key %in% input$cellSetKeys, cell_id])
  scdata <- subset_ids(parent_data$scdata, cell_ids_to_keep)
  scdata <- filter_low_cell_samples(scdata)
  scdata <- diet_scdata(scdata)
  return(scdata)
}


#' Remove samples with low cell numbers
#'
#' When subsetting, there's an increased likelihood of having samples with low
#' cell numbers. This breaks the pipeline in different places, especially during
#' integration. This function removes samples with cell numbers below an empirically
#' determined threshold
#'
#' @param scdata Seurat Object
#' @param min_cells integer of minimum number of cells for a sample to be kept
#'
#' @return Seurat object with samples containing at least min_cells
#' @export
#'
filter_low_cell_samples <- function(scdata, min_cells = MIN_CELLS_IN_SAMPLE) {
  cells_per_sample <- table(scdata$samples)
  samples_to_keep <- names(cells_per_sample)[cells_per_sample > min_cells]
  scdata <- subset(scdata, samples %in% samples_to_keep)

  return(scdata)
}


#' generate a sample id mapping for remaining samples after subset
#'
#' New sample ids must be created, but the number of samples depends on which
#' cells have been subset by the user. Sample Ids that belong to the parent
#' experiment are also kept, which is useful for the addition of the new subclusters
#' to the parent experiment.
#'
#' @param parent_sample_id character vector of unique parent sample ids
#'
#' @return data.table with sample id map
#' @importFrom uuid UUIDgenerate
#' @export
#'
create_sample_id_map <- function(parent_sample_id) {
  subset_sample_id <- UUIDgenerate(n = length(parent_sample_id))

  sample_id_map <- as.list(subset_sample_id)
  names(sample_id_map) <- parent_sample_id

  return(sample_id_map)
}


#' Add new sample ids to the subset Seurat Object
#'
#' @param scdata Seurat Object
#' @param sample_id_map data.table of parent/subset sample id map
#'
#' @return SeuratObject with new sample ids
#' @export
#'
add_new_sample_ids <- function(subset_scdata, sample_id_map) {
  sample_map_idx <- match(subset_scdata$parent_samples, names(sample_id_map))
  subset_scdata$samples <- unname(unlist(sample_id_map[sample_map_idx]))
  return(subset_scdata)
}


#' add experiment level metadata to subset seurat object
#'
#' @param input list of input params, containing the experimentId
#' @param scdata seurat object
#' @param sample_id_map list with mapping between sample_ids from
#'  parent and subset experiments
#'
#' @return scdata with additional metadata
#' @export
#'
add_subset_metadata <- function(input, subset_scdata, sample_id_map) {
  # add new sample_ids, keep originals in a new variable
  subset_scdata$parent_samples <- subset_scdata$samples
  subset_scdata <- add_new_sample_ids(subset_scdata, sample_id_map)
  subset_scdata@misc$experimentId <- input$experimentId

  return(subset_scdata)
}
