#' STEP 3. Mitochondrial content filter
#'
#' This is a simplest filter that looks at a threshold value for the mitochondrial content.
#'
#' To separate cells with low MT-content from the ones that have a high MT-content what makes us think that are dead.
#' This can be a useful first guess. The settings for such a filter can also contain a simple "probabilityThreshold" setting.
#'
#' @param config list containing the following information
#'          - enable: true/false. Refering to apply or not the filter.
#'          - auto: true/false. 'True' indicates that the filter setting need to be changed depending on some sensible value (it requires
#'          to call generate_default_values_mitochondrialContent)
#'          - filterSettings: slot with thresholds
#'                  - method: String. Method to be used {absoluteThreshold}
#'                  - methodSettings: List with the method as key and contain all the filterSettings for this specific method.
#'                          * absoluteThreshold: based on a cut-off threshold
#'                                  - maxFraction: Float. maximun pct MT-content that we considere for a alive cell
#'                                  - binStep: Float. Bin size for the histogram
#'                          * we are supposed to add more methods ....
#' @export
#' @return a list with the filtered seurat object by mitochondrial content, the config and the plot values
#'
filter_high_mito <- function(scdata_list, config, sample_id, cells_id, task_name = "mitochondrialContent", num_cells_to_downsample = 6000) {
  sample_cell_ids <- cells_id[[sample_id]]

  if (length(sample_cell_ids) == 0) {
    return(list(data = scdata_list[[sample_id]], new_ids = cells_id, config = config, plotData = list()))
  }

  sample_data <- subset_ids(scdata_list[[sample_id]], sample_cell_ids)

  # Check if the experiment has MT-content
  if (!"percent.mt" %in% colnames(scdata_list[[sample_id]]@meta.data)) {
    message("Warning! No MT-content has been computed for this experiment!")
    guidata <- list()
    return(list(data = scdata_list[[sample_id]], config = config, plotData = guidata))
  }

  max_fraction <- config$filterSettings$methodSettings[[config$filterSettings$method]]$maxFraction

  plot_data <- generate_mito_plot_data(sample_data)

  if (exists("auto", where = config)) {
    if (as.logical(toupper(config$auto))) {
      max_fraction <- generate_default_values_mitochondrialContent(sample_data, config)
    }
  }


  config$filterSettings$methodSettings[[config$filterSettings$method]]$maxFraction <- max_fraction
  
  # Assign updated config to global env so that it can be accessed if there is an error
  config_key <- paste0("config-", task_name, "-", sample_id)
  assign(config_key, config, envir = globalenv())

  if (as.logical(toupper(config$enabled))) {
    remaining_ids <- sample_data@meta.data$cells_id[sample_data$percent.mt <= max_fraction * 100]
  } else {
    remaining_ids <- sample_cell_ids
  }


  # Downsample plotData
  num_cells_to_downsample <- downsample_plotdata(ncol(sample_data), num_cells_to_downsample)

  set.seed(RANDOM_SEED)
  cells_position_to_keep <- sample(1:ncol(sample_data), num_cells_to_downsample, replace = FALSE)
  cells_position_to_keep <- sort(cells_position_to_keep)
  plot1_data <- plot_data$plot1_data[cells_position_to_keep]
  plot2_data <- plot_data$plot2_data[cells_position_to_keep]

  guidata <- list()

  # plot 1: histgram of MT-content
  guidata[[generate_gui_uuid(sample_id, task_name, 0)]] <- plot1_data
  # plot 2: scatterplot
  guidata[[generate_gui_uuid(sample_id, task_name, 1)]] <- plot2_data

  # populate with filter statistics
  filter_stats <- list(
    before = calc_filter_stats(sample_data),
    after = calc_filter_stats(subset_ids(sample_data, remaining_ids))
  )

  guidata[[generate_gui_uuid(sample_id, task_name, 2)]] <- filter_stats

  cells_id[[sample_id]] <- remaining_ids

  result <- list(
    data = scdata_list,
    new_ids = cells_id,
    config = config,
    plotData = guidata
  )
  return(result)
}


generate_default_values_mitochondrialContent <- function(scdata, config) {
  if (config$filterSettings$method == "absoluteThreshold") {
    is.out <- scuttle::isOutlier(scdata$percent.mt, type = "higher")

    threshold <- ifelse(!sum(is.out),
                        max(scdata$percent.mt),
                        min(scdata$percent.mt[is.out])) / 100
  }

  return(threshold)
}

generate_mito_plot_data <- function(scdata) {
  plot1_data <- lapply(unname(scdata$percent.mt), function(x) {
    c("fracMito" = x)
  })
  plot2_data <- unname(purrr::map2(scdata$percent.mt, scdata$nCount_RNA, function(x, y) {
    c("cellSize" = y, "fracMito" = x)
  }))

  plot_data <- list(plot1_data = plot1_data, plot2_data = plot2_data)
  return(plot_data)
}
