# constructs default QC configuration for merged SeuratObject
construct_qc_config <- function(scdata, any_filtered) {
  samples <- scdata$samples

  # classifier
  config.classifier <- list(
    enabled = !any_filtered,
    prefiltered = any_filtered,
    auto = TRUE,
    filterSettings = list(FDR = 0.01)
  )

  classifier_config_to_duplicate <- list(
    enabled = !any_filtered,
    auto = TRUE,
    filterSettings = list(FDR = 0.01)
  )

  config.classifier <- duplicate_config_per_sample(classifier_config_to_duplicate, config.classifier, samples)


  # cell size
  config.cellSizeDistribution <- list(
    enabled = FALSE,
    auto = TRUE,
    filterSettings = list(minCellSize = 1080, binStep = 200)
  )

  config.cellSizeDistribution <- add_custom_config_per_sample(get_cellsize_config, config.cellSizeDistribution, scdata)


  # mito
  config.mitochondrialContent <- list(
    enabled = TRUE,
    auto = TRUE,
    filterSettings = list(
      method = "absolute_threshold",
      methodSettings = list(
        absolute_threshold = list(
          maxFraction = 0.1,
          binStep = 0.05
        )
      )
    )
  )

  mitochondrial_config_to_duplicate <- list(
    auto = TRUE,
    filterSettings = list(
      method = "absolute_threshold",
      methodSettings = list(absolute_threshold = list(maxFraction = 0.1, binStep = 0.05))
    )
  )

  config.mitochondrialContent <- duplicate_config_per_sample(mitochondrial_config_to_duplicate, config.mitochondrialContent, samples)


  # ngenes vs umis
  config.numGenesVsNumUmis <- list(
    enabled = TRUE,
    auto = TRUE,
    filterSettings = list(
      regressionType = "linear",
      regressionTypeSettings = list(
        "linear" = list(p.level = 0.001),
        "spline" = list(p.level = 0.001)
      )
    )
  )

  config.numGenesVsNumUmis <- add_custom_config_per_sample(get_gene_umi_config, config.numGenesVsNumUmis, scdata)


  # doublet scores
  config.doubletScores <- list(
    enabled = TRUE,
    auto = TRUE,
    filterSettings = list(
      probabilityThreshold = 0.5,
      binStep = 0.05
    )
  )

  config.doubletScores <- add_custom_config_per_sample(get_dblscore_config, config.doubletScores, scdata)

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
          perplexity = min(30, ncol(scdata) / 100),
          learningRate = max(200, ncol(scdata) / 12)
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


get_cellsize_config <- function(scdata, config) {
  minCellSize <- generate_default_values_cellSizeDistribution(scdata, config)
  config$filterSettings$minCellSize <- minCellSize
  return(config)
}

# threshold for doublet score is the max score given to a singlet (above score => doublets)
get_dblscore_config <- function(scdata, config) {
  probabilityThreshold <- max(scdata$doublet_scores[scdata$doublet_class == "singlet"], na.rm = TRUE)
  config$filterSettings$probabilityThreshold <- probabilityThreshold

  return(config)
}


get_gene_umi_config <- function(scdata, config) {
  # Sensible values are based on the function "gene.vs.molecule.cell.filter" from the pagoda2 package
  p.level <- min(0.001, 1 / ncol(scdata))
  config$filterSettings$regressionTypeSettings[[config$filterSettings$regressionType]]$p.level <- p.level

  return(config)
}



duplicate_config_per_sample <- function(step_config, config, samples) {
  for (sample in unique(samples)) {
    config[[sample]] <- step_config
    config[[sample]]$defaultFilterSettings <- step_config$filterSettings
  }

  return(config)
}

add_custom_config_per_sample <- function(step_fn, config, scdata) {

  # We update the config file, so to be able to access the raw config we create a copy
  config.raw <- config

  samples <- scdata$samples

  for (sample in unique(samples)) {
    # Downsample the seurat object to a unisample experiment
    scdata_sample <- scdata[, samples %in% sample]
    # Run the step fun with the unisample experiment and keep the config result
    result_config <- step_fn(scdata_sample, config.raw)
    # Update config with the unisample thresholds
    config[[sample]] <- result_config

    # Add auto settings
    config[[sample]]$defaultFilterSettings <- result_config$filterSettings
  }

  return(config)
}
