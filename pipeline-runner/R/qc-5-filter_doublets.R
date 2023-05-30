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
#'                  through scDblFinder
#'                  - binStep: Float. Bin size for the histogram
#' @export
#' @return a list with the filtered seurat object by doublet score, the config and the plot values
#'
filter_doublets <- function(scdata_list, config, sample_id, cells_id, task_name = "doubletScores", num_cells_to_downsample = 6000) {
  sample_cell_ids <- cells_id[[sample_id]]

  if (length(sample_cell_ids) == 0) {
    return(list(data = scdata_list[[sample_id]], new_ids = cells_id, config = config, plotData = list()))
  }

  sample_data <- subset_ids(scdata_list[[sample_id]], sample_cell_ids)

  if ("recomputeDoubletScore" %in% names(config)) {
    if (config$recomputeDoubletScore) {
    scores <- get_doublet_scores(sample_data@assays$RNA@counts)
    sample_data <- add_dblscore(sample_data, scores)
    # update doublet scores in original scdata
    scdata_list[[sample_id]] <- add_dblscore(scdata_list[[sample_id]], scores)
    }
  }

  # Check if the experiment has doubletScores
  if (!"doublet_scores" %in% colnames(scdata_list[[sample_id]]@meta.data)) {
    message("Warning! No doubletScores scores has been computed for this experiment!")
    guidata <- list()
    return(list(data = scdata_list[[sample_id]], config = config, plotData = guidata))
  }

  probability_threshold <- config$filterSettings[["probabilityThreshold"]]

  # Check if it is required to compute sensible values. From the function 'generate_default_values_doubletScores', it is expected
  # to get a probability threshold value
  if (exists("auto", where = config)) {
    if (as.logical(toupper(config$auto))) {
      probability_threshold <- generate_default_values_doubletScores(sample_data)
    }
  }

  plot1_data <- generate_doublets_plot_data(sample_data, num_cells_to_downsample)

  # update config
  config$filterSettings$probabilityThreshold <- probability_threshold

  # Assign updated config to global env so that it can be accessed if there is an error
  config_key <- paste0("config-", task_name, "-", sample_id)
  assign(config_key, config, envir = globalenv())

  # Check whether the filter is set to true or false
  if (as.logical(toupper(config$enabled))) {
    # all barcodes that match threshold in the subset data
    # treat NA doublet scores as defacto singlets
    doublet_scores <- sample_data$doublet_scores
    doublet_scores[is.na(doublet_scores)] <- 0
    remaining_ids <- sample_data@meta.data$cells_id[doublet_scores <= probability_threshold]
  } else {
    remaining_ids <- sample_cell_ids
  }

  guidata <- list()

  # plot 1: histogram of doublet scores
  guidata[[generate_gui_uuid(sample_id, task_name, 0)]] <- plot1_data

  # populate with filter statistics
  filter_stats <- list(
    before = calc_filter_stats(sample_data),
    after = calc_filter_stats(subset_ids(sample_data, remaining_ids))
  )

  guidata[[generate_gui_uuid(sample_id, task_name, 1)]] <- filter_stats

  cells_id[[sample_id]] <- remaining_ids

  result <- list(
    data = scdata_list,
    new_ids = cells_id,
    config = config,
    plotData = guidata
  )
  return(result)
}

generate_default_values_doubletScores <- function(scdata) {
  # default doublet score based of scDblFinder classification
  # threshold for is the max score given to a singlet (above score => doublets)

  is_singlet <- scdata$doublet_class == "singlet"

  # if no singlets set threshold to 0, preventing some experiments to fail
  threshold <- max(scdata$doublet_scores[is_singlet], 0.0, na.rm = TRUE)

  return(threshold)
}

generate_doublets_plot_data <- function(scdata, num_cells_to_downsample) {
  plot1_data <- lapply(unname(scdata$doublet_scores), function(x) {
    c("doubletP" = x)
  })

  plot1_data <- plot1_data[get_positions_to_keep(scdata, num_cells_to_downsample)]

  return(plot1_data)
}
