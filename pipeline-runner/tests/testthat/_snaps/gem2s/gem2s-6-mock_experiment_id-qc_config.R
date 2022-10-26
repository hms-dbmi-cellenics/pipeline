qc_config <-
list(cellSizeDistribution = list(mock_sample_2_id = list(enabled = FALSE, 
    auto = TRUE, filterSettings = list(minCellSize = 17, binStep = 200), 
    defaultFilterSettings = list(minCellSize = 17, binStep = 200)), 
    mock_sample_1_id = list(enabled = FALSE, auto = TRUE, filterSettings = list(
        minCellSize = 27, binStep = 200), defaultFilterSettings = list(
        minCellSize = 27, binStep = 200))), mitochondrialContent = list(
    mock_sample_2_id = list(enabled = TRUE, auto = TRUE, filterSettings = list(
        method = "absoluteThreshold", methodSettings = list(absoluteThreshold = list(
            maxFraction = 0.52054794520547942, binStep = 0.29999999999999999))), 
        defaultFilterSettings = list(method = "absoluteThreshold", 
            methodSettings = list(absoluteThreshold = list(maxFraction = 0.52054794520547942, 
                binStep = 0.29999999999999999)))), mock_sample_1_id = list(
        enabled = TRUE, auto = TRUE, filterSettings = list(method = "absoluteThreshold", 
            methodSettings = list(absoluteThreshold = list(maxFraction = 0.60256410256410253, 
                binStep = 0.29999999999999999))), defaultFilterSettings = list(
            method = "absoluteThreshold", methodSettings = list(
                absoluteThreshold = list(maxFraction = 0.60256410256410253, 
                  binStep = 0.29999999999999999))))), classifier = list(
    mock_sample_2_id = list(enabled = TRUE, prefiltered = FALSE, 
        auto = TRUE, filterSettings = list(FDR = 0.01), defaultFilterSettings = list(
            FDR = 0.01)), mock_sample_1_id = list(enabled = TRUE, 
        prefiltered = FALSE, auto = TRUE, filterSettings = list(
            FDR = 0.01), defaultFilterSettings = list(FDR = 0.01))), 
    numGenesVsNumUmis = list(mock_sample_2_id = list(enabled = TRUE, 
        auto = TRUE, filterSettings = list(regressionType = "linear", 
            regressionTypeSettings = list(linear = list(p.level = 0.001), 
                spline = list(p.level = 0.001))), defaultFilterSettings = list(
            regressionType = "linear", regressionTypeSettings = list(
                linear = list(p.level = 0.001), spline = list(
                  p.level = 0.001)))), mock_sample_1_id = list(
        enabled = TRUE, auto = TRUE, filterSettings = list(regressionType = "linear", 
            regressionTypeSettings = list(linear = list(p.level = 0.001), 
                spline = list(p.level = 0.001))), defaultFilterSettings = list(
            regressionType = "linear", regressionTypeSettings = list(
                linear = list(p.level = 0.001), spline = list(
                  p.level = 0.001))))), doubletScores = list(
        mock_sample_2_id = list(enabled = TRUE, auto = TRUE, 
            filterSettings = list(probabilityThreshold = 0.97920405864715576, 
                binStep = 0.02), defaultFilterSettings = list(
                probabilityThreshold = 0.97920405864715576, binStep = 0.02)), 
        mock_sample_1_id = list(enabled = TRUE, auto = TRUE, 
            filterSettings = list(probabilityThreshold = 0.83960545063018799, 
                binStep = 0.02), defaultFilterSettings = list(
                probabilityThreshold = 0.83960545063018799, binStep = 0.02))), 
    dataIntegration = list(dataIntegration = list(method = "harmony", 
        methodSettings = list(seuratv4 = list(numGenes = 2000, 
            normalisation = "logNormalize"), unisample = list(
            numGenes = 2000, normalisation = "logNormalize"), 
            harmony = list(numGenes = 2000, normalisation = "logNormalize"), 
            fastmnn = list(numGenes = 2000, normalisation = "logNormalize"))), 
        dimensionalityReduction = list(method = "rpca", numPCs = NULL, 
            excludeGeneCategories = list())), configureEmbedding = list(
        embeddingSettings = list(method = "umap", methodSettings = list(
            umap = list(minimumDistance = 0.29999999999999999, 
                distanceMetric = "cosine"), tsne = list(perplexity = 30, 
                learningRate = 200))), clusteringSettings = list(
            method = "louvain", methodSettings = list(louvain = list(
                resolution = 0.80000000000000004)))))
