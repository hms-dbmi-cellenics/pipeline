#' STEP 5. Doublet score filter
#'
#' Filters seurat object based on doubletScores scores
#'
#' This is a simplest filter that looks at a threshold value for the doublet scores.
#' To separate cells with low droplet score from the ones that have a high droplet score content what makes us think that the are mistakenly considered as a single cell but they are actully two or more.
#' This can be a useful first guess. The settings for such a filter can also contain a simple "probabilityThreshold" setting.
#'
#' @param config list containing the following information
#'          - enable: true/false. Refering to apply or not the filter.
#'          - auto: true/false. 'True' indicates that the filter setting need to be changed depending on some sensible value (it requires
#'          to call generate_default_values_doubletScores)
#'          - filterSettings: slot with thresholds
#'                  - probabilityThreshold: Float. cut-off for the maximun probability scores that could have a cell. The doubletScores scores have been computed
#'                  through "scrublet" [1]
#'                  - binStep: Float. Bin size for the histogram
#' @export
#' @return a list with the filtered seurat object by doublet score, the config and the plot values
#'
filter_doublets <- function(scdata, config, sample_id, cells_id,task_name = "doubletScores", num_cells_to_downsample = 6000) {
  cells_id.sample <- cells_id[[sample_id]]

  if (length(cells_id.sample) == 0) {
    guidata <- list()
    return(list(data = scdata, config = config, plotData = guidata))
  }

  scdata.sample <- subset_ids(scdata,cells_id.sample)

  # Check if the experiment has doubletScores
  if (!"doublet_scores" %in% colnames(scdata@meta.data)) {
    message("Warning! No doubletScores scores has been computed for this experiment!")
    guidata <- list()
    return(list(data = scdata, config = config, plotData = guidata))
  }

  probabilityThreshold <- config$filterSettings[["probabilityThreshold"]]

  # Check if it is required to compute sensible values. From the function 'generate_default_values_doubletScores', it is expected
  # to get a value --> probabilityThreshold.
  if (exists("auto", where = config)) {
    if (as.logical(toupper(config$auto))) {
      probabilityThreshold <- generate_default_values_doubletScores(scdata.sample)
    }
  }

  plot1_data <- generate_doublets_plot_data(scdata.sample,num_cells_to_downsample)

  # Check whether the filter is set to true or false
  if (as.logical(toupper(config$enabled))) {
    # all barcodes that match threshold in the subset data
    # treat NA doublet scores as defacto singlets
    doublet_scores <- scdata.sample$doublet_scores
    doublet_scores[is.na(doublet_scores)] <- 0
    remaining_ids <- scdata.sample@meta.data$cells_id[doublet_scores <= probabilityThreshold]
  }
  else {
    remaining_ids <- cells_id.sample
  }

  # update config
  config$filterSettings$probabilityThreshold <- probabilityThreshold

  guidata <- list()

  # plot 1: histogram of doublet scores
  guidata[[generate_gui_uuid(sample_id, task_name, 0)]] <- plot1_data

  # populate with filter statistics
  filter_stats <- list(
    before = calc_filter_stats(scdata, sample_id),
    after = calc_filter_stats(scdata.filtered, sample_id)
  )

  guidata[[generate_gui_uuid(sample_id, task_name, 1)]] <- filter_stats

  cells_id[[sample_id]] <- remaining_ids

  result <- list(
    data = scdata,
    new_ids = cells_id,
    config = config,
    plotData = guidata
  )

  return(result)
}

generate_default_values_doubletScores <- function(scdata) {
  # default doublet score based of scDblFinder classification
  is.singlet <- scdata$doublet_class == "singlet"
  threshold <- max(scdata$doublet_scores[is.singlet], na.rm = TRUE)

  return(threshold)
}

generate_doublets_plot_data <- function(scdata,num_cells_to_downsample){
  plot1_data <- lapply(unname(scdata.sample$doublet_scores), function(x) {
    c("doubletP" = x)
  })

  plot1_data <- plot1_data[get_positions_to_keep(scdata,num_cells_to_downsample)]

  return(plot1_data)
}