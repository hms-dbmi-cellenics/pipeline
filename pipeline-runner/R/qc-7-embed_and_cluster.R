# STEP 7. Compute embedding

# Compute embedding step where we run dimensional reduction technniques such as t-SNE and  UMAP. Moreover, the cluster analysis is
# done also in this step.

embed_and_cluster <- function(scdata, config, sample_id, task_name = 'configureEmbedding') {
  message("starting clusters")
  type <- config$clusteringSettings$method
  methodSettings <- config$clusteringSettings$methodSettings[[type]]
  message("Running clustering")
  res <- runClusters(type,methodSettings$resolution,scdata)
  message("formatting cellsets")
  formated_cell_sets <- format_cell_sets_object(res,type,scdata@misc$color_pool)
  message("udating through api")
  update_sets_through_api(formated_cell_sets,config$api_url,scdata@misc$experimentId,type)

  # the result object will have to conform to this format: {data, config, plotData : {plot1, plot2}}
  result <- list(
    data = scdata,
    config = config,
    plotData = list()
  )

  return(result)
}

format_cell_sets_object <- function (cell_sets,type,color_pool) {
  key <- type
  name <- paste0(type," clusters")

  #careful with capital l on type for the key.
  cell_sets_object <- list(key=key,name=name,rootNode=TRUE,type="cellSets",children=list())
  for(i in sort(unique(cell_sets$cluster))){
    cells <- cell_sets[cell_sets$cluster==i,"cell_ids"]
    new_set <- list(
      key=paste0(key,"-",i),
      name=paste0("Cluster ",i),
      rootNode=FALSE,
      type="cellSets",
      color=color_pool[1],
      cellIds=unname(cells)
    )
    color_pool <- color_pool[-1]
    cell_sets_object$children <- append(cell_sets_object$children,list(new_set))
  }
  return(cell_sets_object)
}

update_sets_through_api <- function(cell_sets_object,api_url,experiment_id,cell_set_key){
  httr::PATCH(
    paste0(api_url,"/v1/experiments/",experiment_id,"/cellSets"),
    body = list(list("$match"=list(query=paste0("$[?(@.key == ",cell_set_key,")]"),"$remove"=TRUE)),list("$prepend"=cell_sets_object)),
    encode="json",
    httr::add_headers(
      "Content-Type"= "application/boschni-json-merger+json",
      "Authorization"= ""))
}


