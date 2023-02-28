#' Constructs default QC configuration
#'
#' This function returns the default parameters used during QC as a nested list.
#' It is sent to the API, which in turn saves it as a jsonb object in the
#' PostgreSQL database.
#'
#' @param scdata_list list of seurat objects
#' @param unfiltered_samples character of unfiltered sample ids
#' @param disable_qc_filters bool indicating if the data derives from the
#' subsetting of another experiment
#'
#' @return list of QC configuration parameters
#'
construct_qc_config <- function(scdata_list, unfiltered_samples, disable_qc_filters) {
  samples <- names(scdata_list)

  config_classifier <-
    add_custom_config_per_sample(get_classifier_config,
                                 processing_config_template[["classifier"]],
                                 scdata_list,
                                 disable_qc_filters,
                                 unfiltered_samples)

  config_cell_size <-
    add_custom_config_per_sample(get_cellsize_config,
                                 processing_config_template[["cell_size"]],
                                 scdata_list,
                                 disable_qc_filters)

  config_mitochondrial <-
    add_custom_config_per_sample(get_sample_mitochondrial_config,
                                 processing_config_template[["mitochondrial"]],
                                 scdata_list,
                                 disable_qc_filters)

  config_genes_vs_umis <-
    add_custom_config_per_sample(get_gene_umi_config,
                                 processing_config_template[["genes_vs_umis"]],
                                 scdata_list,
                                 disable_qc_filters)


  config_doublet <-
    add_custom_config_per_sample(get_dblscore_config,
                                 processing_config_template[["doublet"]],
                                 scdata_list,
                                 disable_qc_filters)

  config_data_integration <- processing_config_template[["data_integration"]]

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


get_classifier_config <- function(scdata, config) {
  return(config)
}

get_cellsize_config <- function(scdata, config) {
  minCellSize <- generate_default_values_cellSizeDistribution(scdata, config)
  config$filterSettings$minCellSize <- minCellSize
  return(config)
}

get_sample_mitochondrial_config <- function(scdata, config) {
  default_max_fraction <- generate_default_values_mitochondrialContent(scdata, config)
  config$filterSettings$methodSettings$absoluteThreshold$maxFraction <- default_max_fraction

  return(config)
}


# threshold for doublet score is the max score given to a singlet (above score => doublets)
get_dblscore_config <- function(scdata, config) {
  probabilityThreshold <- generate_default_values_doubletScores(scdata)

  config$filterSettings$probabilityThreshold <- probabilityThreshold

  return(config)
}


get_gene_umi_config <- function(scdata, config) {
  # Sensible values are based on the function "gene.vs.molecule.cell.filter" from the pagoda2 package
  p.level <- min(0.001, 1 / ncol(scdata))
  config$filterSettings$regressionTypeSettings[[config$filterSettings$regressionType]]$p.level <- p.level

  return(config)
}

get_embedding_config <- function(scdata_list, config) {
  # tsne parameters depend on number of cells in sample
  default_perplexity <- config$embeddingSettings$methodSettings$tsne$perplexity
  default_learning_rate <- config$embeddingSettings$methodSettings$tsne$learningRate

  config$embeddingSettings$methodSettings$tsne <- list(
    perplexity = min(default_perplexity, min(vapply(scdata_list, ncol, integer(1))) / 100),
    learningRate = max(default_learning_rate, min(vapply(scdata_list, ncol, integer(1))) / 12)
  )

  return(config)
}


add_custom_config_per_sample <- function(generate_sample_config, default_config, scdata_list, disable_qc_filters = FALSE, unfiltered_samples = NA) {
  # We update the config file, so to be able to access the raw config we create a copy
  raw_config <- default_config
  config <- list()
  for (sample in names(scdata_list)) {
    # subset the Seurat object list to a single sample
    sample_data <- scdata_list[[sample]]

    # run the function to generate config for a sample
    sample_config <- generate_sample_config(sample_data, raw_config)

    if (rlang::has_name(sample_config, "prefiltered")) {
      # only change these values for the classifier config
      sample_config$enabled <- sample %in% unfiltered_samples && !disable_qc_filters
      sample_config$prefiltered <- !(sample %in% unfiltered_samples)
    }

    if (rlang::has_name(sample_config, "enabled")) {
      sample_config$enabled <- sample_config$enabled && !disable_qc_filters
    }

    # update sample config thresholds
    config[[sample]] <- sample_config
  }

  return(config)
}
