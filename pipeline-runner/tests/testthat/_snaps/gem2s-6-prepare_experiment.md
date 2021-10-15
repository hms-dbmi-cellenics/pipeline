# prepare_experiment generates qc_config that matches snapshot

    Code
      str(task_out$qc_config)
    Output
      List of 7
       $ cellSizeDistribution:List of 4
        ..$ enabled       : logi FALSE
        ..$ auto          : logi TRUE
        ..$ filterSettings:List of 2
        .. ..$ minCellSize: num 1080
        .. ..$ binStep    : num 200
        ..$ sample_a      :List of 4
        .. ..$ enabled              : logi FALSE
        .. ..$ auto                 : logi TRUE
        .. ..$ filterSettings       :List of 2
        .. .. ..$ minCellSize: num 10
        .. .. ..$ binStep    : num 200
        .. ..$ defaultFilterSettings:List of 2
        .. .. ..$ minCellSize: num 10
        .. .. ..$ binStep    : num 200
       $ mitochondrialContent:List of 4
        ..$ enabled       : logi TRUE
        ..$ auto          : logi TRUE
        ..$ filterSettings:List of 2
        .. ..$ method        : chr "absolute_threshold"
        .. ..$ methodSettings:List of 1
        .. .. ..$ absolute_threshold:List of 2
        .. .. .. ..$ maxFraction: num 0.1
        .. .. .. ..$ binStep    : num 0.05
        ..$ sample_a      :List of 3
        .. ..$ auto                 : logi TRUE
        .. ..$ filterSettings       :List of 2
        .. .. ..$ method        : chr "absolute_threshold"
        .. .. ..$ methodSettings:List of 1
        .. .. .. ..$ absolute_threshold:List of 2
        .. .. .. .. ..$ maxFraction: num 0.1
        .. .. .. .. ..$ binStep    : num 0.05
        .. ..$ defaultFilterSettings:List of 2
        .. .. ..$ method        : chr "absolute_threshold"
        .. .. ..$ methodSettings:List of 1
        .. .. .. ..$ absolute_threshold:List of 2
        .. .. .. .. ..$ maxFraction: num 0.1
        .. .. .. .. ..$ binStep    : num 0.05
       $ classifier          :List of 5
        ..$ enabled       : logi TRUE
        ..$ prefiltered   : logi FALSE
        ..$ auto          : logi TRUE
        ..$ filterSettings:List of 1
        .. ..$ FDR: num 0.01
        ..$ sample_a      :List of 4
        .. ..$ enabled              : logi TRUE
        .. ..$ auto                 : logi TRUE
        .. ..$ filterSettings       :List of 1
        .. .. ..$ FDR: num 0.01
        .. ..$ defaultFilterSettings:List of 1
        .. .. ..$ FDR: num 0.01
       $ numGenesVsNumUmis   :List of 4
        ..$ enabled       : logi TRUE
        ..$ auto          : logi TRUE
        ..$ filterSettings:List of 2
        .. ..$ regressionType        : chr "gam"
        .. ..$ regressionTypeSettings:List of 1
        .. .. ..$ gam:List of 1
        .. .. .. ..$ p.level: num 0.001
        ..$ sample_a      :List of 4
        .. ..$ enabled              : logi TRUE
        .. ..$ auto                 : logi TRUE
        .. ..$ filterSettings       :List of 2
        .. .. ..$ regressionType        : chr "gam"
        .. .. ..$ regressionTypeSettings:List of 1
        .. .. .. ..$ gam:List of 1
        .. .. .. .. ..$ p.level: num 0.00013
        .. ..$ defaultFilterSettings:List of 2
        .. .. ..$ regressionType        : chr "gam"
        .. .. ..$ regressionTypeSettings:List of 1
        .. .. .. ..$ gam:List of 1
        .. .. .. .. ..$ p.level: num 0.00013
       $ doubletScores       :List of 4
        ..$ enabled       : logi TRUE
        ..$ auto          : logi TRUE
        ..$ filterSettings:List of 2
        .. ..$ probabilityThreshold: num 0.5
        .. ..$ binStep             : num 0.05
        ..$ sample_a      :List of 4
        .. ..$ enabled              : logi TRUE
        .. ..$ auto                 : logi TRUE
        .. ..$ filterSettings       :List of 2
        .. .. ..$ probabilityThreshold: num 0.8
        .. .. ..$ binStep             : num 0.05
        .. ..$ defaultFilterSettings:List of 2
        .. .. ..$ probabilityThreshold: num 0.8
        .. .. ..$ binStep             : num 0.05
       $ dataIntegration     :List of 2
        ..$ dataIntegration        :List of 2
        .. ..$ method        : chr "harmony"
        .. ..$ methodSettings:List of 4
        .. .. ..$ seuratv4 :List of 2
        .. .. .. ..$ numGenes     : num 2000
        .. .. .. ..$ normalisation: chr "logNormalize"
        .. .. ..$ unisample:List of 2
        .. .. .. ..$ numGenes     : num 2000
        .. .. .. ..$ normalisation: chr "logNormalize"
        .. .. ..$ harmony  :List of 2
        .. .. .. ..$ numGenes     : num 2000
        .. .. .. ..$ normalisation: chr "logNormalize"
        .. .. ..$ fastmnn  :List of 2
        .. .. .. ..$ numGenes     : num 2000
        .. .. .. ..$ normalisation: chr "logNormalize"
        ..$ dimensionalityReduction:List of 3
        .. ..$ method               : chr "rpca"
        .. ..$ numPCs               : num 30
        .. ..$ excludeGeneCategories: list()
       $ configureEmbedding  :List of 2
        ..$ embeddingSettings :List of 2
        .. ..$ method        : chr "umap"
        .. ..$ methodSettings:List of 2
        .. .. ..$ umap:List of 2
        .. .. .. ..$ minimumDistance: num 0.3
        .. .. .. ..$ distanceMetric : chr "cosine"
        .. .. ..$ tsne:List of 2
        .. .. .. ..$ perplexity  : num 30
        .. .. .. ..$ learningRate: num 640
        ..$ clusteringSettings:List of 2
        .. ..$ method        : chr "louvain"
        .. ..$ methodSettings:List of 1
        .. .. ..$ louvain:List of 1
        .. .. .. ..$ resolution: num 0.8

