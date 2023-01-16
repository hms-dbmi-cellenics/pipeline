prepare_scdata <- function(scdata_list, config, cells_id) {
  # pre-process

  scdata_list <- order_by_size(scdata_list)
  message("Started create_scdata")
  scdata <- create_scdata(scdata_list, cells_id)
  message("Finished create_scdata")

  exclude_groups <- config$dimensionalityReduction$excludeGeneCategories
  # remove cell cycle genes if needed
  if (length(exclude_groups) > 0) {
    scdata <- remove_genes(scdata, exclude_groups)
  }

  return(scdata)
}


integrate_scdata <- function(scdata_list, config, sample_id, cells_id, task_name = "dataIntegration", use_geosketch = FALSE, perc_num_cells = 5) {

  scdata <- prepare_scdata(scdata_list, config, cells_id)

  message("Started data integration")

  # get method and settings
  method <- config$dataIntegration$method
  npcs <- config$dimensionalityReduction$numPCs

  if (length(unique(scdata$samples) == 1)) {
    method <- "unisample"
    message("Only one sample detected or method is non integrate.")
  }

  # we need RNA assay to compute the integrated matrix
  Seurat::DefaultAssay(scdata) <- "RNA"

  # integrate
  integration_function <- get(paste0("run_", method))
  scdata <- integration_function(scdata, config, use_geosketch)

  message("Finished data integration")

  # get npcs from the PCA
  npcs <- length(scdata_integrated@reductions$pca)
  message("\nSet config numPCs to npcs used in last UMAP call: ", npcs, "\n")
  config$dimensionalityReduction$numPCs <- npcs

  var_explained <- get_explained_variance(scdata_integrated)

  # This same numPCs will be used throughout the platform.
  scdata_integrated@misc[["numPCs"]] <- config$dimensionalityReduction$numPCs

  scdata_integrated <- colorObject(scdata_integrated)

  plots <- generate_elbow_plot_data(scdata_integrated, config, task_name, var_explained)

  # the result object will have to conform to this format: {data, config, plotData : {plot1, plot2}}
  result <- list(
    data = scdata_integrated,
    new_ids = cells_id,
    config = config,
    plotData = plots
  )

  return(result)

}


