

check_prev_out <- function(prev_out, check_names) {
  for (check_name in check_names) {
    if (check_name %in% names(prev_out)) next()

    stop(check_name, " is missing from prev_out but is required.")
  }
}
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
          binStep = 0.3
        )
      )
    )
  )

  config.mitochondrialContent <- add_custom_config_per_sample(get_sample_mitochondrial_config, config.mitochondrialContent, scdata)

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
      binStep = 0.02
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

get_sample_mitochondrial_config <- function(scdata.sample, config) {

  config.sample <- list(
    auto = TRUE,
    filterSettings = list(
      method = "absolute_threshold",
      methodSettings = list()
    )
  )

  config.sample$filterSettings$methodSettings$absolute_threshold <- list(
    maxFraction = generate_default_values_mitochondrialContent(scdata.sample, config.sample),
    binStep = 0.3
  )

  return(config.sample)
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

add_custom_config_per_sample <- function(generate_sample_config, config, scdata) {

  # We update the config file, so to be able to access the raw config we create a copy
  config.raw <- config

  samples <- scdata$samples

  for (sample in unique(samples)) {
    # subset the Seurat object to a single sample
    scdata.sample <- scdata[, samples %in% sample]

    # run the function to generate config for a sample
    config.sample <- generate_sample_config(scdata.sample, config.raw)

    # update sample config thresholds
    config[[sample]] <- config.sample

    # add auto settings
    config[[sample]]$defaultFilterSettings <- config.sample$filterSettings
  }

  return(config)
}

generate_default_values_cellSizeDistribution <- function(scdata, config) {
  # `CalculateBarcodeInflections` including inflection point calculation
  threshold.low <- if (ncol(scdata) <= 200) NULL else 100
  scdata_tmp <- Seurat::CalculateBarcodeInflections(
    object = scdata,
    barcode.column = "nCount_RNA",
    group.column = "samples",
    threshold.low = threshold.low
  )
  # returned is both the rank(s) as well as inflection point
  # extracting only inflection point(s)
  sample_subset <- Seurat::Tool(scdata_tmp, slot = "Seurat::CalculateBarcodeInflections")$inflection_points
  return(sample_subset$nCount_RNA)
}


get_cellsize_config <- function(scdata, config) {
  minCellSize <- generate_default_values_cellSizeDistribution(scdata, config)
  config$filterSettings$minCellSize <- minCellSize
  return(config)
}


add_custom_config_per_sample <- function(generate_sample_config, config, scdata) {

  # We update the config file, so to be able to access the raw config we create a copy
  config.raw <- config

  samples <- scdata$samples

  for (sample in unique(samples)) {
    # subset the Seurat object to a single sample
    scdata.sample <- scdata[, samples %in% sample]

    # run the function to generate config for a sample
    config.sample <- generate_sample_config(scdata.sample, config.raw)

    # update sample config thresholds
    config[[sample]] <- config.sample

    # add auto settings
    config[[sample]]$defaultFilterSettings <- config.sample$filterSettings
  }

  return(config)
}


add_custom_config_per_sample <- function(generate_sample_config, config, scdata) {

  # We update the config file, so to be able to access the raw config we create a copy
  config.raw <- config

  samples <- scdata$samples

  for (sample in unique(samples)) {
    # subset the Seurat object to a single sample
    scdata.sample <- scdata[, samples %in% sample]

    # run the function to generate config for a sample
    config.sample <- generate_sample_config(scdata.sample, config.raw)

    # update sample config thresholds
    config[[sample]] <- config.sample

    # add auto settings
    config[[sample]]$defaultFilterSettings <- config.sample$filterSettings
  }

  return(config)
}

generate_default_values_mitochondrialContent <- function(scdata, config) {
  if (config$filterSettings$method == "absolute_threshold") {
    is.out <- 20

    threshold <- ifelse(!sum(is.out),
                        max(scdata$percent.mt),
                        min(scdata$percent.mt[is.out])) / 100
  }

  return(threshold)
}

get_sample_mitochondrial_config <- function(scdata.sample, config) {

  config.sample <- list(
    auto = TRUE,
    filterSettings = list(
      method = "absolute_threshold",
      methodSettings = list()
    )
  )

  config.sample$filterSettings$methodSettings$absolute_threshold <- list(
    maxFraction = generate_default_values_mitochondrialContent(scdata.sample, config.sample),
    binStep = 0.3
  )

  return(config.sample)
}

prepare_experiment <- function() {
  message("Preparing experiment ...")
  prev_out <- readRDS('prev_out.rds')
  pipeline_config <- readRDS('pipeline_config.rds')
  check_names <- c("config", "counts_list", "annot", "doublet_scores", "scdata_list")
  input <- {}
  input$experimentId <- 'randomEID'
  check_prev_out(prev_out, check_names)

  scdata_list <- prev_out$scdata_list
  prev_out$scdata_list <- NULL
  samples <- names(scdata_list)

  message("Merging Seurat Objects...")
  message(sum(sapply(scdata_list, ncol)))
#   saveRDS(scdata_list, "/debug/scdata_list.rds")
#   saveRDS(prev_out, "/debug/prev_out.rds")
#   saveRDS(pipeline_config, "/debug/pipeline_config.rds")
    message("Merging ", length(scdata_list))
  if (length(scdata_list) == 1) {
    scdata <- scdata_list[[1]]
  } else {
    # scdata <- merge(scdata_list[[1]], y = scdata_list[-1])
    scdata <- merge(scdata_list[[1]], y = scdata_list[[2]])
    scdata_list <- scdata_list[3:length(scdata_list)]
    while(length(scdata_list) > 0)
      scdata <- merge(scdata, scdata_list[[1]])
      scdata_list <- scdata_list[2:length(scdata_list)]
      qc()
  }

  message("Merging done")
  message("Adding metadata...")
  scdata <- add_metadata(scdata, prev_out$annot, input$experimentId)
  message("Added metadata...")
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

# merge_scdatas <- function(scdata_list) {
#   message("Merging ", length(scdata_list))
#   if (length(scdata_list) == 1) {
#     scdata <- scdata_list[[1]]
#   } else {
#     # scdata <- merge(scdata_list[[1]], y = scdata_list[-1])
#     scdata <- merge(scdata_list[[1]], y = scdata_list[[2]])
#     scdata_list <- scdata_list[3:length(scdata_list)]
#     while(length(scdata_list) > 0)
#       scdata <- merge(scdata, scdata_list[[1]])
#   }

#   message("Merging done")
#   return(scdata)
# }

add_metadata <- function(scdata, annot, experiment_id) {

  # Ensure index by rownames in scdata
  annot <- annot[match(rownames(scdata), annot$input), ]
  scdata@misc[["gene_annotations"]] <- annot

  message("Storing cells id...")
  # Keeping old version of ids starting from 0
  scdata$cells_id <- 0:(ncol(scdata) - 1)

  message("Storing color pool...")
  # We store the color pool in a slot in order to be able to access it during configureEmbedding
  # scdata@misc[["color_pool"]] <- get_color_pool()
  scdata@misc[["experimentId"]] <- experiment_id
  scdata@misc[["ingestionDate"]] <- Sys.time()

  return(scdata)
}

res <- prepare_experiment()