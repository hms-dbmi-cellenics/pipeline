#' run clustering
#'
#'
#' @param scdata seurat object
#' @param config list with clustering parameters
#' @param sample_id character
#' @param cells_id list of cell ids that passed all filters so far
#' @param task_name character
#'
#' @return list
#' @export
#'
embed_and_cluster <-
  function(scdata,
           config,
           sample_id,
           cells_id,
           task_name = "configureEmbedding",
           ignore_ssl_cert = FALSE) {

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
      config$auth_JWT,
      ignore_ssl_cert
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

#' Ensure is list in json
#'
#' When sending responses as json, Vectors of length 0 or 1 are converted to
#' null and scalar (respectively) Using as.list fixes this, however, long R
#' lists take a VERY long time to be converted to JSON.
#' This function deals with the problematic cases, leaving vector as a vector
#' when it isnt a problem.
#'
#' @param vector
#'
#' @export
#'
ensure_is_list_in_json <- function(vector) {
  if (length(vector) <= 1) {
    return(as.list(vector))
  } else {
    return(vector)
  }
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
        cellIds = ensure_is_list_in_json(unname(cells))
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
           auth_JWT,
           ignore_ssl_cert
    ) {

    httr_query <- paste0("$[?(@.key == \"", cell_set_key, "\")]")

    if (ignore_ssl_cert) {
        httr::set_config(httr::config(ssl_verifypeer = 0L))
    }

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
