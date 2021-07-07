# STEP 1. Classifier filter

#' @description Filters seurat object based on mitochondrialContent
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
filter_emptydrops <- function(scdata, config, sample_id, task_name = 'classifier', num_cells_to_downsample = 6000) {

  tmp_sample <- sub("sample-", "", sample_id)

  before <- calc_filter_stats(scdata, tmp_sample)

  FDR <- config$filterSettings$FDR

  if (isTRUE(config$auto)) {
    FDR <- generate_default_values_classifier(scdata, config)
  }

  if (config$enabled) {

    # check if filter data is actually available
    if (is.null(scdata@meta.data$emptyDrops_FDR)) {
      message("Classify is enabled but has no classify data available: will dissable it: no filtering!")
      config$enabled <- FALSE
      guidata <- list()
    } else {
      message("Classify is enabled: filtering with FDR=", FDR)
      obj_metadata <- scdata@meta.data
      # extract plotting data of original data to return to plot slot later
      barcode_names_this_sample <- rownames(obj_metadata[grep(tmp_sample, rownames(obj_metadata)), ])
      if (length(barcode_names_this_sample) == 0) {
        return(list(data = scdata, config = config, plotData = list()))
      }
      sample_subset <- subset(scdata, cells = barcode_names_this_sample)
      message("Info: empty-drops table of FDR threshold categories (# UMIs for a given threshold interval)")
      print(table(obj_metadata$samples, cut(obj_metadata$emptyDrops_FDR, breaks = c(-Inf, 0, 0.0001, 0.01, 0.1, 0.5, 1, Inf)), useNA = "ifany"))

      # prevents filtering of NA FDRs if FDR=1
      ed_fdr <- sample_subset$emptyDrops_FDR
      ed_fdr[is.na(ed_fdr)] <- 1

      message(
        "Number of barcodes to filter for this sample: ",
        sum(ed_fdr > FDR, na.rm = TRUE), "/", length(ed_fdr)
      )

      numis <- log10(sample_subset@meta.data$nCount_RNA)

      fdr_data <- unname(purrr::map2(ed_fdr, numis, function(x, y) {
        c("FDR" = x, "log_u" = y)
      }))
      # extract cell id that do not(!) belong to current sample (to not apply filter there)
      barcode_names_non_sample <- rownames(obj_metadata[-grep(tmp_sample, rownames(obj_metadata)), ])
      # all barcodes that match threshold in the subset data
      barcode_names_keep_current_sample <- colnames(sample_subset[, ed_fdr <= FDR])
      # combine the 2:
      barcodes_to_keep <- union(barcode_names_non_sample, barcode_names_keep_current_sample)
      scdata <- subset_safe(scdata, barcodes_to_keep)
      # update config
      config$filterSettings$FDR <- FDR

      # Downsample plotData
      # Handle when the number of remaining cells is less than the number of cells to downsample
      nkeep <- downsample_plotdata(ncol(sample_subset), num_cells_to_downsample)

      set.seed(123)
      keep <- sample(1:ncol(sample_subset), nkeep, replace = FALSE)
      keep <- sort(keep)
      fdr_data <- fdr_data[keep]

      knee_data <- get_bcranks_plot_data(sample_subset, is.cellsize = FALSE)[['knee']]

      # Populate guidata list
      guidata <- list()
      guidata[[generate_gui_uuid(sample_id, task_name, 0)]] <- fdr_data
      guidata[[generate_gui_uuid(sample_id, task_name, 1)]] <- knee_data
    }
  } else {
    message("filter disabled: data not filtered!")
    guidata <- list()
  }

  # get filter stats after filtering
  filter_stats <- list(
    before = before,
    after = calc_filter_stats(scdata, tmp_sample)
  )
  guidata[[generate_gui_uuid(sample_id, task_name, 2)]] <- filter_stats

  result <- list(
    data = scdata,
    config = config,
    plotData = guidata
  )

  return(result)
}



#' @description Filters seurat object based on classifier filter using emptyDrops
#               https://rdrr.io/github/MarioniLab/DropletUtils/man/emptyDrops.html
#' @param config list containing the following information
#'          - enable: true/false. Refering to apply or not the filter.
#'          - auto: true/false. 'True' indicates that the filter setting need to be changed depending on some sensible value (it requires
#'          to call generate_default_values_classifier)
#'          - filterSettings: slot with thresholds
#'                  - minProbabiliy:
#'                  - filterThreshold:
#' @export
#' @return a list with the filtered seurat object by probabilities classifier, the config and the plot values
#'
generate_default_values_classifier <- function(scdata, config) {

  # HARDCODE
  threshold <- 0.01

  return(threshold)
}
