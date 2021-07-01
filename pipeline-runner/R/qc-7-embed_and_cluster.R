# STEP 7. Compute embedding

# Compute embedding step where we run dimensional reduction technniques such as t-SNE and  UMAP. Moreover, the cluster analysis is
# done also in this step.

embed_and_cluster <- function(scdata, config, sample_id, task_name = 'configureEmbedding') {

  plots <- list()

  # the result object will have to conform to this format: {data, config, plotData : {plot1, plot2}}
  result <- list(
    data = scdata,
    config = config,
    plotData = plots
  )

  return(result)
}
