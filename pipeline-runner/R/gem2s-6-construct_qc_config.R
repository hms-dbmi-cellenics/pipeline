#' Constructs default QC configuration
#'
#' This function returns the default parameters used during QC as a nested list.
#' It is sent to the API, which in turn saves it as a jsonb object in the
#' PostgreSQL database.
#'
#' @param scdata_list list of seurat objects
#' @param any_filtered bool indicating if barcodes were filtered by emptyDrops
#' @param subset_experiment bool indicating if the data derives from the
#' subsetting of another experiment
#'
#' @return list of QC configuration parameters
#'
construct_qc_config <- function(scdata_list, any_filtered, subset_experiment) {
  samples <- names(scdata_list)

  # classifier
  classifier_config_to_duplicate <- list(
    enabled = !any_filtered && !subset_experiment,
    prefiltered = any_filtered,
    auto = TRUE,
    filterSettings = list(FDR = 0.01)
  )

  config.classifier <- duplicate_config_per_sample(classifier_config_to_duplicate, samples)


  # cell size
  default_cellSizeDistribution_config <- list(
    enabled = FALSE,
    auto = TRUE,
    filterSettings = list(minCellSize = 1080, binStep = 200)
  )

  config.cellSizeDistribution <- add_custom_config_per_sample(get_cellsize_config, default_cellSizeDistribution_config, scdata_list)

  # mito
  default_mitochondrialContent_config <- list(
    enabled = !subset_experiment,
    auto = TRUE,
    filterSettings = list(
      method = "absoluteThreshold",
      methodSettings = list(
        absoluteThreshold = list(
          maxFraction = 0.1,
          binStep = 0.3
        )
      )
    )
  )

  config.mitochondrialContent <- add_custom_config_per_sample(get_sample_mitochondrial_config, default_mitochondrialContent_config, scdata_list)

  # ngenes vs umis
  default_numGenesVsNumUmis_config <- list(
    enabled = !subset_experiment,
    auto = TRUE,
    filterSettings = list(
      regressionType = "linear",
      regressionTypeSettings = list(
        "linear" = list(p.level = 0.001),
        "spline" = list(p.level = 0.001)
      )
    )
  )

  config.numGenesVsNumUmis <- add_custom_config_per_sample(get_gene_umi_config, default_numGenesVsNumUmis_config, scdata_list)


  # doublet scores
  default_doubletScores_config <- list(
    enabled = !subset_experiment,
    auto = TRUE,
    filterSettings = list(
      probabilityThreshold = 0.5,
      binStep = 0.02
    )
  )

  config.doubletScores <- add_custom_config_per_sample(get_dblscore_config, default_doubletScores_config, scdata_list)

  # data integration
  config.dataIntegration <- list(
    dataIntegration = list(
      method = "harmony",
      methodSettings = list(
        seuratv4 = list(numGenes = 2000, normalisation = "logNormalize"),
        unisample = list(numGenes = 2000, normalisation = "logNormalize"),
        harmony = list(numGenes = 2000, normalisation = "logNormalize"),
        fastmnn = list(numGenes = 2000, normalisation = "logNormalize")
      )
    ),
    dimensionalityReduction = list(
      method = "rpca",
      numPCs = NULL,
      excludeGeneCategories = list()
    )
  )


  # embedding
  config.configureEmbedding <- list(
    embeddingSettings = list(
      method = "umap",
      methodSettings = list(
        umap = list(
          minimumDistance = 0.3,
          distanceMetric = "cosine"
        ),
        tsne = list(
          perplexity = min(30, ncol(scdata_list) / 100),
          learningRate = max(200, ncol(scdata_list) / 12)
        )
      )
    ),
    clusteringSettings = list(
      method = "louvain",
      methodSettings = list(louvain = list(resolution = 0.8))
    )
  )

  # combine config for all steps
  config <- list(
    cellSizeDistribution = config.cellSizeDistribution,
    mitochondrialContent = config.mitochondrialContent,
    classifier = config.classifier,
    numGenesVsNumUmis = config.numGenesVsNumUmis,
    doubletScores = config.doubletScores,
    dataIntegration = config.dataIntegration,
    configureEmbedding = config.configureEmbedding
  )

  return(config)
}


get_cellsize_config <- function(scdata_list, config) {
  minCellSize <- generate_default_values_cellSizeDistribution(scdata_list, config)
  config$filterSettings$minCellSize <- minCellSize
  return(config)
}

get_sample_mitochondrial_config <- function(scdata_list.sample, config) {
  config.sample <- list(
    enabled = config$enabled,
    auto = TRUE,
    filterSettings = list(
      method = "absoluteThreshold",
      methodSettings = list()
    )
  )

  config.sample$filterSettings$methodSettings$absoluteThreshold <- list(
    maxFraction = generate_default_values_mitochondrialContent(scdata_list.sample, config.sample),
    binStep = 0.3
  )

  return(config.sample)
}


# threshold for doublet score is the max score given to a singlet (above score => doublets)
get_dblscore_config <- function(scdata_list, config) {
  probabilityThreshold <- max(scdata_list$doublet_scores[scdata_list$doublet_class == "singlet"], na.rm = TRUE)
  config$filterSettings$probabilityThreshold <- probabilityThreshold

  return(config)
}


get_gene_umi_config <- function(scdata_list, config) {
  # Sensible values are based on the function "gene.vs.molecule.cell.filter" from the pagoda2 package
  p.level <- min(0.001, 1 / ncol(scdata_list))
  config$filterSettings$regressionTypeSettings[[config$filterSettings$regressionType]]$p.level <- p.level

  return(config)
}


duplicate_config_per_sample <- function(step_config, samples) {
  config <- list()
  for (sample in unique(samples)) {
    config[[sample]] <- step_config
    config[[sample]]$defaultFilterSettings <- step_config$filterSettings
  }

  return(config)
}


add_custom_config_per_sample <- function(generate_sample_config, default_config, scdata_list) {
  # We update the config file, so to be able to access the raw config we create a copy
  raw_config <- default_config
  config <- list()
  for (sample in names(scdata_list)) {
    # subset the Seurat object list to a single sample
    sample_data <- scdata_list[[sample]]

    # run the function to generate config for a sample
    sample_config <- generate_sample_config(sample_data, raw_config)

    # update sample config thresholds
    config[[sample]] <- sample_config

    # add auto settings
    config[[sample]]$defaultFilterSettings <- sample_config$filterSettings
  }

  return(config)
}
