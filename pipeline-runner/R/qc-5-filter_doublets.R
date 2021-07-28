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
filter_doublets <- function(scdata, config, sample_id, task_name = "doubletScores", num_cells_to_downsample = 6000) {

  # Check if the experiment has doubletScores
  if (!"doublet_scores" %in% colnames(scdata@meta.data)) {
    message("Warning! No doubletScores scores has been computed for this experiment!")
    return(scdata)
  }
  probabilityThreshold <- config$filterSettings[["probabilityThreshold"]]

  # Check if it is required to compute sensible values. From the function 'generate_default_values_doubletScores', it is expected
  # to get a value --> probabilityThreshold.
  if (exists("auto", where = config)) {
    if (as.logical(toupper(config$auto))) {
      probabilityThreshold <- generate_default_values_doubletScores(scdata, sample_id)
    }
  }

  # extract plotting data of original data to return to plot slot later
  meta <- scdata@meta.data
  barcode_names_this_sample <- rownames(meta)[meta$samples == sample_id]
  if (length(barcode_names_this_sample) == 0) {
    return(list(data = scdata, config = config, plotData = list()))
  }
  sample_subset <- subset(scdata, cells = barcode_names_this_sample)

  # Check whether the filter is set to true or false
  if (as.logical(toupper(config$enabled))) {
    # extract cells id that do not(!) belong to current sample (to not apply filter there)
    barcode_names_non_sample <- rownames(meta)[meta$samples != sample_id, ]
    # all barcodes that match threshold in the subset data
    # treat NA doublet scores as defacto singlets
    doublet_scores <- sample_subset$doublet_scores
    doublet_scores[is.na(doublet_scores)] <- 0
    barcode_names_keep_current_sample <- colnames(sample_subset)[doublet_scores <= probabilityThreshold]
    # combine the 2:
    barcodes_to_keep <- union(barcode_names_non_sample, barcode_names_keep_current_sample)
    # Information regarding doublet score is pre-computed during the 'data-ingest'.
    scdata.filtered <- subset_safe(scdata, barcodes_to_keep)
  }
  else {
    scdata.filtered <- scdata
  }

  # update config
  config$filterSettings$probabilityThreshold <- probabilityThreshold

  plot1_data <- lapply(unname(sample_subset$doublet_scores), function(x) {
    c("doubletP" = x)
  })

  # Downsample plotData
  num_cells_to_downsample <- downsample_plotdata(ncol(sample_subset), num_cells_to_downsample)

  set.seed(123)
  cells_position_to_keep <- sample(1:ncol(sample_subset), num_cells_to_downsample, replace = FALSE)
  cells_position_to_keep <- sort(cells_position_to_keep)
  plot1_data <- plot1_data[cells_position_to_keep]

  guidata <- list()

  # plot 1: histogram of doublet scores
  guidata[[generate_gui_uuid(sample_id, task_name, 0)]] <- plot1_data

  # populate with filter statistics
  filter_stats <- list(
    before = calc_filter_stats(scdata, sample_id),
    after = calc_filter_stats(scdata.filtered, sample_id)
  )

  guidata[[generate_gui_uuid(sample_id, task_name, 1)]] <- filter_stats

  result <- list(
    data = scdata.filtered,
    config = config,
    plotData = guidata
  )
  return(result)
}

generate_default_values_doubletScores <- function(scdata, sample) {
  # default doublet score based of scDblFinder classification
  is.sample <- scdata$samples == sample
  is.singlet <- scdata$doublet_class == "singlet"
  threshold <- max(scdata$doublet_scores[is.sample & is.singlet], na.rm = TRUE)

  return(threshold)
}
