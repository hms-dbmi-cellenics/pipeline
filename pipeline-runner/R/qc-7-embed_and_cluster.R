# STEP 7. Compute embedding

# Compute embedding step where we run dimensional reduction technniques such as t-SNE and  UMAP. Moreover, the cluster analysis is
# done also in this step.

embed_and_cluster <- function(scdata, config, sample_id, task_name = 'configureEmbedding') {
  message("starting clusters")
  req <- list()
  type <- config$clusteringSettings$method
  methodSettings <- config$clusteringSettings$methodSettings[[type]]
  req$body <- list(type=type, config=list(resolution=methodSettings$resolution))
  message("Running clustering")
  res <- runClusters(req,scdata)
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
  for(i in unique(cell_sets$cluster)){
    cells <- cell_sets[cell_sets$cluster==i,"cell_ids"]
    new_set <- list(
      key=paste0("key","-",i),
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

update_sets_through_api <- function(cell_sets_object,apiUrl,expId,cellSetKey){
  message(expId)
  httr::PATCH(
    paste0("http://localhost:3000","/v1/experiments/",expId,"/cellSets"),
    body = list("$match"=list(query=paste0("$[?(@.key == ",cellSetKey,")]"),"$remove"=TRUE),"$prepend"=cell_sets_object),
    encode="json",
    httr::add_headers(
      "Content-Type"= "application/boschni-json-merger+json",
      "Authorization"= "Bearer eyJraWQiOiJ2TytRZ1lud0lnOU5pT2Y3azJTNEFEY2xvaDBwVlNUbkNNdDJMOU8xQU1RPSIsImFsZyI6IlJTMjU2In0.eyJhdF9oYXNoIjoiN0tFQUFkeVE1dVpBcDhCNHFuZDRzUSIsInN1YiI6ImI3MjM0YmVjLTE4OGQtNDAzMS1iYWRmLWQ2NGM1NDg5ZjgyYSIsImVtYWlsX3ZlcmlmaWVkIjp0cnVlLCJpc3MiOiJodHRwczpcL1wvY29nbml0by1pZHAuZXUtd2VzdC0xLmFtYXpvbmF3cy5jb21cL2V1LXdlc3QtMV9zVlFMMTVZeHUiLCJjb2duaXRvOnVzZXJuYW1lIjoiYjcyMzRiZWMtMTg4ZC00MDMxLWJhZGYtZDY0YzU0ODlmODJhIiwiYXVkIjoiMWR1N2FoOXByaXNvN3A1b2F2aGEyaTluamgiLCJldmVudF9pZCI6IjViZDE3MzdlLTJhMWMtNGM2MC1hNTU1LTg0OTYxMzcwN2Q4NSIsInRva2VuX3VzZSI6ImlkIiwiYXV0aF90aW1lIjoxNjMxMTE5MTMzLCJuYW1lIjoiT2xpdmVyIiwiZXhwIjoxNjMxMTMwNzg1LCJpYXQiOjE2MzExMjcxODUsImVtYWlsIjoib2xpdmVyQGJpb21hZ2UubmV0In0.S5rpbxhnOcuyjzKir6Rbu3lIsWGyThxo0Y2fc__aS1WNQLPGOwiTJ46tpxIFp877ir_ltRemVnP07IE_bKaiuJ9Bn_0KqHhiU2ka2Wp8L925-uKuqJqFLs4e7_-EKdVUThg91x370fLxaO3pfYfOezN9kMPjmvjBi3OzyAKai8-xg-dv9tjJH-YR7XaeB2dxexPUyiH3bMjsqdz1BjBCqRqj7ZQ9MeOYTTcGtzFHV4-uXof0TcfuqNV3DsgwTlV9COmYSTRjRfgljDOxcnGqpCO9vQhxIqdvMUkIDUPjLI46PQjJwlIe2Fg5y_p42qNC44Fgf4EEgqC14NhneluLRg"))
}


