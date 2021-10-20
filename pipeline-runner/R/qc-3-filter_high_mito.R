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
#'                  - method: String. Method to be used {absolute_threshold}
#'                  - methodSettings: List with the method as key and contain all the filterSettings for this specific method.
#'                          * absolute_threshold: based on a cut-off threshold
#'                                  - maxFraction: Float. maximun pct MT-content that we considere for a alive cell
#'                                  - binStep: Float. Bin size for the histogram
#'                          * we are supposed to add more methods ....
#' @export
#' @return a list with the filtered seurat object by mitochondrial content, the config and the plot values
#'
filter_high_mito <- function(scdata, config, sample_id, task_name = "mitochondrialContent", num_cells_to_downsample = 6000) {
  # Check if the experiment has MT-content
  if (!"percent.mt" %in% colnames(scdata@meta.data)) {
    # This return value breaks the step, no config data is returned!!!!
    message("Warning! No MT-content has been computed for this experiment!")
    return(scdata)
  }

  maxFraction <- config$filterSettings$methodSettings[[config$filterSettings$method]]$maxFraction

  # Computing the subset for the current sample
  meta <- scdata@meta.data
  barcode_names_this_sample <- rownames(meta)[meta$samples == sample_id]
  if (length(barcode_names_this_sample) == 0) {
    return(list(data = scdata, config = config, plotData = list()))
  }
  sample_subset <- subset(scdata, cells = barcode_names_this_sample)
  if (exists("auto", where = config)) {
    if (as.logical(toupper(config$auto))) {
      maxFraction <- generate_default_values_mitochondrialContent(sample_subset, config)
    }
  }

  if (as.logical(toupper(config$enabled))) {
    # extract cell id that do not(!) belong to current sample (to not apply filter there)
    barcode_names_non_sample <- rownames(meta)[meta$samples != sample_id]
    # all barcodes that match threshold in the subset
    barcode_names_keep_current_sample <- colnames(sample_subset)[sample_subset$percent.mt <= maxFraction * 100]
    # combine the 2:
    barcodes_to_keep <- union(barcode_names_non_sample, barcode_names_keep_current_sample)

    scdata.filtered <- subset_safe(scdata, barcodes_to_keep)
  } else {
    scdata.filtered <- scdata
  }
  plot1_data <- lapply(unname(sample_subset$percent.mt), function(x) {
    c("fracMito" = x)
  })
  plot2_data <- unname(purrr::map2(sample_subset$percent.mt, sample_subset$nCount_RNA, function(x, y) {
    c("cellSize" = y, "fracMito" = x)
  }))

  config$filterSettings$methodSettings[[config$filterSettings$method]]$maxFraction <- maxFraction

  # Downsample plotData
  num_cells_to_downsample <- downsample_plotdata(ncol(sample_subset), num_cells_to_downsample)

  set.seed(123)
  cells_position_to_keep <- sample(1:ncol(sample_subset), num_cells_to_downsample, replace = FALSE)
  cells_position_to_keep <- sort(cells_position_to_keep)
  plot1_data <- plot1_data[cells_position_to_keep]
  plot2_data <- plot2_data[cells_position_to_keep]

  guidata <- list()

  # plot 1: histgram of MT-content
  guidata[[generate_gui_uuid(sample_id, task_name, 0)]] <- plot1_data

  guidata[[generate_gui_uuid(sample_id, task_name, 1)]] <- plot2_data

  # populate with filter statistics
  filter_stats <- list(
    before = calc_filter_stats(scdata, sample_id),
    after = calc_filter_stats(scdata.filtered, sample_id)
  )

  guidata[[generate_gui_uuid(sample_id, task_name, 2)]] <- filter_stats

  result <- list(
    data = scdata.filtered,
    config = config,
    plotData = guidata
  )
  stop('FAILURE!!')
  return(result)
}


# The most uses values in MT-content are between [0.1, 0.2]. There are not too much literature about how to compute
# a threshold.
# --> Absolute threshold: In order to be not too extrictive the threshold is set to 0.1
generate_default_values_mitochondrialContent <- function(scdata, config) {
  if (config$filterSettings$method == "absolute_threshold") {
    # HARDCODE
    threshold <- 0.1
  }

  return(threshold)
}
