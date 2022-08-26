send_output_to_api_mock_data <- list(
  pipeline_config = list(
    cluster_env = "development",
    sandbox_id = "default",
    aws_account_id = "000000000000",
    aws_region = "eu-west-1",
    pod_name = "local",
    activity_arn = "arn:aws:states:eu-west-1:000000000000:activity:pipeline-development-201b0f11-310e-4a60-962f-524c870ce409",
    api_url = "http://host.docker.internal:3000",
    api_version = "v2",
    debug_config = list(
      step = "",
      path = "/mocked/path"
    ),
    aws_config = list(
      region = "eu-west-1",
      endpoint = "http://host.docker.internal:4566",
      credentials = list(
        creds = list(
          access_key_id = "mock-access-key",
          secret_access_key = "mock-secure-acces-key"
        )
      )
    ),
    originals_bucket = "biomage-originals-development-000000000000",
    source_bucket = "biomage-source-development-000000000000",
    processed_bucket = "processed-matrix-development-000000000000",
    results_bucket = "worker-results-development-000000000000",
    cells_id_bucket = "biomage-filtered-cells-development-000000000000",
    plot_data_bucket = "plots-tables-development-000000000000",
    cell_sets_bucket = "cell-sets-development-000000000000",
    debug_bucket = "biomage-pipeline-debug-development-000000000000",
    sns_topic = "arn:aws:sns:eu-west-1:000000000000:work-results-mocked-default-v2"
  ),
  input = list(
    experimentId = "eb6eb7453062ff1018d7c6f136cb713b",
    taskName = "configureEmbedding",
    processName = "qc",
    config = list(
      embeddingSettings = list(
        method = "umap",
        methodSettings = list(
          tsne = list(
            perplexity = 9.18,
            learningRate = 200
          ),
          umap = list(
            distanceMetric = "cosine",
            minimumDistance = 0.3
          )
        )
      ),
      clusteringSettings = list(
        method = "louvain",
        methodSettings = list(
          louvain = list(
            resolution = 0.8
          )
        )
      )
    ),
    sampleUuid = "",
    uploadCountMatrix = FALSE,
    authJWT = 'mockedAuthJWT'
  ),
  plot_data_keys = list(),
  output = list(
    config = list(
      embeddingSettings = list(
        method = "umap",
        methodSettings = list(
          tsne = list(
            perplexity = 9.18,
            learningRate = 200
          ),
          umap = list(
            distanceMetric = "cosine",
            minimumDistance = 0.3
          )
        )
      ),
      clusteringSettings = list(
        method = "louvain",
        methodSettings = list(
          louvain = list(
            resolution = 0.8
          )
        )
      )
    ),
    plotData = list()
  )
)