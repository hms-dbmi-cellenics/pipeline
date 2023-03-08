#' Constructs default QC configuration
#'
#' This function returns the default parameters used during QC as a nested list.
#' It is sent to the API, which in turn saves it as a jsonb object in the
#' PostgreSQL database.
#'
#' For each pipeline step, it customizes and adds the QC parameters for all samples
#' in the running experiment.
#'
#' @param scdata_list list of Seurat objects
#' @param unfiltered_samples character vector of unfiltered sample ids
#' @param disable_qc_filters bool indicating if filters should be disabled.
#'
#' @return list of QC configuration parameters
#'
construct_qc_config <-
  function(scdata_list,
           disable_qc_filters,
           unfiltered_samples) {
    samples <- names(scdata_list)

    config_classifier <-
      add_custom_config_per_sample(
        customize_classifier_config,
        processing_config_template[["classifier"]],
        scdata_list,
        disable_qc_filters,
        unfiltered_samples
      )

    config_cell_size <-
      add_custom_config_per_sample(
        customize_cellsize_config,
        processing_config_template[["cell_size"]],
        scdata_list,
        disable_qc_filters
      )

    config_mitochondrial <-
      add_custom_config_per_sample(
        customize_mitochondrial_config,
        processing_config_template[["mitochondrial"]],
        scdata_list,
        disable_qc_filters
      )

    config_genes_vs_umis <-
      add_custom_config_per_sample(
        customize_genes_vs_umis_config,
        processing_config_template[["genes_vs_umis"]],
        scdata_list,
        disable_qc_filters
      )


    config_doublet <-
      add_custom_config_per_sample(
        customize_doublet_config,
        processing_config_template[["doublet"]],
        scdata_list,
        disable_qc_filters
      )

    config_data_integration <-
      processing_config_template[["data_integration"]]

    config_embedding_clustering <-
      get_embedding_config(scdata_list, processing_config_template[["embedding_clustering"]])

    # combine config for all steps
    config <- list(
      cellSizeDistribution = config_cell_size,
      mitochondrialContent = config_mitochondrial,
      classifier = config_classifier,
      numGenesVsNumUmis = config_genes_vs_umis,
      doubletScores = config_doublet,
      dataIntegration = config_data_integration,
      configureEmbedding = config_embedding_clustering
    )

    return(config)
  }


customize_classifier_config <-
  function(scdata,
           config,
           sample_name,
           disable_qc_filters,
           unfiltered_samples) {
    config$enabled <- sample_name %in% unfiltered_samples && !disable_qc_filters
    config$prefiltered <- !(sample_name %in% unfiltered_samples)

    return(config)
  }


customize_cellsize_config <-
  function(scdata,
           config,
           sample_name,
           disable_qc_filters,
           unfiltered_samples) {
    minCellSize <- generate_default_values_cellSizeDistribution(scdata, config)
    config$filterSettings$minCellSize <- minCellSize
    return(config)
  }


customize_mitochondrial_config <-
  function(scdata,
           config,
           sample_name,
           disable_qc_filters,
           unfiltered_samples) {
    default_max_fraction <- generate_default_values_mitochondrialContent(scdata, config)
    config$filterSettings$methodSettings$absoluteThreshold$maxFraction <-
      default_max_fraction

    return(config)
  }


customize_doublet_config <-
  function(scdata,
           config,
           sample_name,
           disable_qc_filters,
           unfiltered_samples) {
    probabilityThreshold <- generate_default_values_doubletScores(scdata)
    config$filterSettings$probabilityThreshold <- probabilityThreshold

    return(config)
  }


customize_genes_vs_umis_config <-
  function(scdata,
           config,
           sample_name,
           disable_qc_filters,
           unfiltered_samples) {
    # Sensible values are based on the function "gene.vs.molecule.cell.filter"
    # from the pagoda2 package
    p.level <- min(0.001, 1 / ncol(scdata))
    regression_type <- config$filterSettings$regressionType
    config$filterSettings$regressionTypeSettings[[regression_type]]$p.level <- p.level

    return(config)
  }


get_embedding_config <- function(scdata_list, config) {
  # tsne parameters depend on number of cells in sample
  default_perplexity <- config$embeddingSettings$methodSettings$tsne$perplexity
  default_learning_rate <- config$embeddingSettings$methodSettings$tsne$learningRate

  config$embeddingSettings$methodSettings$tsne <- list(
    perplexity = min(
      default_perplexity,
      min(vapply(scdata_list, ncol, integer(1))) / 100),
    learningRate = max(
      default_learning_rate,
      min(vapply(scdata_list, ncol, integer(1))) / 12)
  )

  return(config)
}


#' Customize configuration for each sample in an experiment
#'
#' Takes the step config template and the required data and calls the corresponding
#' step customization function, which takes care of changing the QC parameters to
#' sensible values for each sample in the running experiment.
#'
#' Values that apply to all experiments are defined in the `data-raw/processing_config_template.R`
#' file.
#'
#' @param customize_template_config function - step customization function
#' @param config_template list - template of step configuration parameters
#' @param scdata_list list - with Seurat objects
#' @param disable_qc_filters logical
#' @param unfiltered_samples character vector
#'
#' @return list of customized QC parameters for each sample
#' @export
#'
add_custom_config_per_sample <-
  function(customize_template_config,
           config_template,
           scdata_list,
           disable_qc_filters = FALSE,
           unfiltered_samples = NA) {
    config <- list()
    for (sample_name in names(scdata_list)) {
      # subset the Seurat object list to a single sample
      sample_scdata <- scdata_list[[sample_name]]

      # run the function to generate config for a sample
      sample_config <-
        customize_template_config(
          sample_scdata,
          config_template,
          sample_name,
          disable_qc_filters,
          unfiltered_samples
        )

      sample_config$enabled <-
        sample_config$enabled && !disable_qc_filters

      # update sample config thresholds
      config[[sample_name]] <- sample_config
    }

    return(config)
  }
