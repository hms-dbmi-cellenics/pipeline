processing_config_template <- list(
  "classifier" = list(
    enabled = NA,
    prefiltered = NA,
    auto = TRUE,
    filterSettings = list(FDR = 0.01)
  ),

  "cell_size" = list(
    enabled = FALSE,
    auto = TRUE,
    filterSettings = list(minCellSize = 1080, binStep = 200)
  ),

  "mitochondrial" = list(
    enabled = TRUE,
    auto = TRUE,
    filterSettings = list(
      method = "absoluteThreshold",
      methodSettings = list(absoluteThreshold = list(
        maxFraction = 0.1,
        binStep = 0.3
      ))
    )
  ),

  "genes_vs_umis" = list(
    enabled = TRUE,
    auto = TRUE,
    filterSettings = list(
      regressionType = "linear",
      regressionTypeSettings = list(
        "linear" = list(p.level = 0.001),
        "spline" = list(p.level = 0.001)
      )
    )
  ),

  "doublet" = list(
    enabled = TRUE,
    auto = TRUE,
    filterSettings = list(probabilityThreshold = 0.5,
                          binStep = 0.02)
  ),

  "data_integration" = list(
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
  ),

  "embedding_clustering" = list(
    embeddingSettings = list(
      method = "umap",
      methodSettings = list(
        umap = list(minimumDistance = 0.3,
                    distanceMetric = "cosine"),
        tsne = list(
          perplexity = 30,
          learningRate = 200
        )
      )
    ),
    clusteringSettings = list(method = "louvain",
                              methodSettings = list(louvain = list(resolution = 0.8)))
  )
)
