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

  # Spatial (Visium HD, Xenium) QC metric filters. method is "localOutlier"
  # (z-score cutoff vs the spatial neighbourhood) or "absoluteThreshold" (raw
  # metric bounds; NULL = no bound). Direction is fixed per filter in code.
  # The absolute defaults are filled in per sample by construct_qc_config.
  "spatial_umi_outlier" = list(
    enabled = TRUE,
    auto = TRUE,
    filterSettings = list(
      method = "localOutlier",
      methodSettings = list(
        localOutlier = list(cutoff = 3, binStep = 0.2),
        absoluteThreshold = list(min = NULL, max = NULL, binStep = 0.05)
      )
    )
  ),

  "spatial_num_genes_outlier" = list(
    enabled = TRUE,
    auto = TRUE,
    filterSettings = list(
      method = "localOutlier",
      methodSettings = list(
        localOutlier = list(cutoff = 3, binStep = 0.2),
        absoluteThreshold = list(min = NULL, max = NULL, binStep = 0.05)
      )
    )
  ),

  "spatial_mito_outlier" = list(
    enabled = TRUE,
    auto = TRUE,
    filterSettings = list(
      method = "localOutlier",
      methodSettings = list(
        localOutlier = list(cutoff = 3, binStep = 0.2),
        absoluteThreshold = list(max = NULL, binStep = 1)
      )
    )
  ),

  # Xenium count-density filter (OSTA imaging QC): 3-MAD lower outlier on
  # log(counts / cell area), or control probe + codeword counts above
  # maxControlCounts when filterControls is TRUE.
  "spatial_count_density" = list(
    enabled = TRUE,
    auto = TRUE,
    filterSettings = list(
      nmads = 3,
      filterControls = TRUE,
      maxControlCounts = 0,
      binStep = 0.05
    )
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
    ),
    # qc-6-integrate_scdata.R will treat below downsampling config
    # as default and use geosketch if warrented + update config
    downsampling = list(
      method = "default"
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
