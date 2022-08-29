# STEP 7. Compute embedding

# Compute embedding step where we run dimensional reduction technniques such as t-SNE and  UMAP. Moreover, the cluster analysis is
# done also in this step.

embed_and_cluster <-
  function(scdata,
           config,
           sample_id,
           cells_id,
           task_name = "configureEmbedding") {

    message("starting clusters")
    clustering_method <- config$clusteringSettings$method
    methodSettings <-
      config$clusteringSettings$methodSettings[[clustering_method]]
    message("Running clustering")
    cellSets <-
      runClusters(clustering_method, methodSettings$resolution, scdata)
    message("formatting cellsets")
    formated_cell_sets <-
      format_cell_sets_object(cellSets, clustering_method, scdata@misc$color_pool)
    message("updating through api")

    update_sets_through_api(
      formated_cell_sets,
      config$api_url,
      scdata@misc$experimentId,
      clustering_method,
      config$auth_JWT
    )

    # the result object will have to conform to this format:
    # {data, config, plotData : {plot1, plot2}}
    result <- list(
      data = scdata,
      new_ids = cells_id,
      config = config,
      plotData = list()
    )

    return(result)
  }

format_cell_sets_object <-
  function(cell_sets, clustering_method, color_pool) {
    name <- paste0(clustering_method, " clusters")

    # careful with capital l on type for the key.
    cell_sets_object <-
      list(
        key = clustering_method,
        name = name,
        rootNode = TRUE,
        type = "cellSets",
        children = list()
      )
    for (i in sort(unique(cell_sets$cluster))) {
      cells <- cell_sets[cell_sets$cluster == i, "cell_ids"]
      new_set <- list(
        key = paste0(clustering_method, "-", i),
        name = paste0("Cluster ", i),
        rootNode = FALSE,
        type = "cellSets",
        color = color_pool[1],
        cellIds = unname(cells)
      )
      color_pool <- color_pool[-1]
      cell_sets_object$children <-
        append(cell_sets_object$children, list(new_set))
    }
    return(cell_sets_object)
  }

update_sets_through_api <-
  function(cell_sets_object,
           api_url,
           experiment_id,
           cell_set_key,
           auth_JWT) {

    httr_query <- paste0("$[?(@.key == \"", cell_set_key, "\")]")

    httr::PATCH(
      paste0(api_url, "/v2/experiments/", experiment_id, "/cellSets"),
      body = list(list(
        "$match" = list(query = httr_query, value = list("$remove" = TRUE))
      ),
      list("$prepend" = cell_sets_object)),
      encode = "json",
      httr::add_headers("Content-Type" = "application/boschni-json-merger+json",
                        "Authorization" = auth_JWT)
    )
  }
