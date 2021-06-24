# STEP 2. Cell size distribution filter
filter_low_cellsize <- function(scdata, config, sample_id, task_name = 'cellSizeDistribution', num_cells_to_downsample = 6000) {

  tmp_sample <- sub("sample-", "", sample_id)

  minCellSize <- as.numeric(config$filterSettings$minCellSize)

  # extract plotting data of original data to return to plot slot later
  obj_metadata <- scdata@meta.data
  barcode_names_this_sample <- rownames(obj_metadata[grep(tmp_sample, rownames(obj_metadata)), ])
  if (length(barcode_names_this_sample) == 0) {
    guidata <- list()
    return(list(data = scdata, config = config, plotData = guidata))
  }
  sample_subset <- subset(scdata, cells = barcode_names_this_sample)

  # umi histogram plot
  numis <- sort(sample_subset$nCount_RNA, decreasing = TRUE)
  plot1_data <- lapply(unname(numis), function(x) {
    c("u" = x)
  })

  # barcode ranks plot data
  # unique average ranks maintain plot shape in case of downsampling
  ranks <- rank(-numis, ties.method = "average")
  dups <- duplicated(ranks)
  ranks <- ranks[!dups]
  numis <- numis[!dups]

  plot2_data <- unname(purrr::map2(log(numis), ranks, function(x, y) {
    c("log_u" = x, "rank" = y)
  }))

  # downsample plot data
  set.seed(123)
  nkeep1 <- downsample_plotdata(ncol(sample_subset), num_cells_to_downsample)
  nkeep2 <- downsample_plotdata(length(ranks), num_cells_to_downsample)

  keep1 <- sort(sample(ncol(sample_subset), nkeep1))
  keep2 <- sort(sample(length(ranks), nkeep2))

  plot1_data <- plot1_data[keep1]
  plot2_data <- plot2_data[keep2]

  # Check if it is required to compute sensible values. From the function 'generate_default_values_cellSizeDistribution', it is expected
  # to get a list with two elements {minCellSize and binStep}
  if (exists("auto", where = config)) {
    if (as.logical(toupper(config$auto))) {
      # config not really needed for this one (maybe later for threshold.low/high):
      # HARDCODE Value. threshold.low 
      # [ Parameter for function CalculateBarcodeInflections. Description: Ignore barcodes of rank below this threshold in inflection calculation]
      threshold.low <- 1e2
      # If there are less cells than the value threshold.low, the function CalculateBarcodeInflections fails. So we need to handle by not removing any cells, that is,
      # consider the minCellSize as the minimun UMIs in the dataset.
      # This should be handled in a long-term by adding a different function for computing the default value.
      if (ncol(sample_subset) < threshold.low) {
        minCellSize <- min(sample_subset$nCount_RNA)
      } else {
        minCellSize <- generate_default_values_cellSizeDistribution(sample_subset, config, threshold.low)
      }
    }
  }
  if (as.logical(toupper(config$enabled))) {

    # extract cell id that do not(!) belong to current sample (to not apply filter there)
    barcode_names_non_sample <- rownames(obj_metadata[-grep(tmp_sample, rownames(obj_metadata)), ])
    # all barcodes that match threshold in the subset data
    barcode_names_keep_current_sample <- rownames(sample_subset@meta.data[sample_subset@meta.data$nCount_RNA >= minCellSize, ])
    # combine the 2:
    barcodes_to_keep <- union(barcode_names_non_sample, barcode_names_keep_current_sample)

    scdata.filtered <- subset_safe(scdata, barcodes_to_keep)
  } else {
    scdata.filtered <- scdata
  }

  # update config
  config$filterSettings$minCellSize <- minCellSize
  # Populate data for UI
  guidata <- list()
  guidata[[generate_gui_uuid(sample_id, task_name, 0)]] <- plot1_data
  guidata[[generate_gui_uuid(sample_id, task_name, 1)]] <- plot2_data

  # Populate with filter statistics
  filter_stats <- list(
    before = calc_filter_stats(scdata, tmp_sample),
    after = calc_filter_stats(scdata.filtered, tmp_sample)
  )
  guidata[[generate_gui_uuid(sample_id, task_name, 3)]] <- filter_stats

  result <- list(
    data = scdata.filtered,
    config = config,
    plotData = guidata
  )

  return(result)
}


# cell_size_distribution_filter function
# This is a simplest filter that looks at the shape of the cell size (# of UMIs per cell) distribution and looks for some local minima, minimum of second derivative, etc.
# To separate the large cell barcodes that correspond to real cells from the tail containing 'empty droplets'.
# This can be a useful first guess. The settings for such a filter can also contain a simple "min cell size" setting.
#
#' @description Filters seurat object based on cell size distribution
#' @param config list containing the following information
#'          - enable: true/false. Refering to apply or not the filter.
#'          - auto: true/false. 'True' indicates that the filter setting need to be changed depending on some sensible value (it requires
#'          to call generate_default_values_cellSizeDistribution)
#'          - filterSettings: slot with thresholds
#'                  - minCellSize: Integer. Threshold that contain the minimun number of UMIs per cell
#'                  - binStep: Integer. Bin size for the histogram
#' @export
#' @return a list with the filtered seurat object by cell size ditribution, the config and the plot values

# CalculateBarcodeInflections calculates an adaptive inflection point ("knee")
# of the barcode distribution for each sample group. This is
# useful for determining a threshold for removing low-quality
# samples.
generate_default_values_cellSizeDistribution <- function(scdata, config, threshold.low) {
  # `CalculateBarcodeInflections` including inflection point calculation
  scdata_tmp <- CalculateBarcodeInflections(
    object = scdata,
    barcode.column = "nCount_RNA",
    group.column = "samples",
    # [HARDCODED]
    threshold.low = threshold.low,
    threshold.high = NULL
  )
  # returned is both the rank(s) as well as inflection point
  sample_subset <- Tool(scdata_tmp, slot = "CalculateBarcodeInflections")$inflection_points
  # extracting only inflection point(s)
  return(sample_subset$nCount_RNA)
}
