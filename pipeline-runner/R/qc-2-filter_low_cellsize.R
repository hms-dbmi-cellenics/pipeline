#' STEP 2. Cell size distribution filter
#'
#' cell_size_distribution_filter function
#'
#' This is a simplest filter that looks at the shape of the cell size (# of UMIs per cell) distribution and looks for some local minima, minimum of second derivative, etc.
#' To separate the large cell barcodes that correspond to real cells from the tail containing 'empty droplets'.
#' This can be a useful first guess. The settings for such a filter can also contain a simple "min cell size" setting.
#'
#' @param config list containing the following information
#'          - enable: true/false. Refering to apply or not the filter.
#'          - auto: true/false. 'True' indicates that the filter setting need to be changed depending on some sensible value (it requires
#'          to call generate_default_values_cellSizeDistribution)
#'          - filterSettings: slot with thresholds
#'                  - minCellSize: Integer. Threshold that contain the minimun number of UMIs per cell
#'                  - binStep: Integer. Bin size for the histogram
#' @export
#' @return a list with the filtered seurat object by cell size ditribution, the config and the plot values
#'
filter_low_cellsize <- function(scdata, config, sample_id, task_name = "cellSizeDistribution", num_cells_to_downsample = 6000) {
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
  plot_data <- get_bcranks_plot_data(sample_subset, num_cells_to_downsample)

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
  guidata[[generate_gui_uuid(sample_id, task_name, 0)]] <- plot_data[['knee']]
  guidata[[generate_gui_uuid(sample_id, task_name, 1)]] <- plot_data[['hist']]

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

get_bcranks_plot_data <- function(sample_subset, nmax = 6000, is.cellsize = TRUE) {

  set.seed(123)
  plot_data <- list()
  numis <- unname(sample_subset$nCount_RNA)
  ord <- order(numis, decreasing = TRUE)
  numis <- numis[ord]

  # umi histogram plot
  if (is.cellsize) {
    plot1_data <- lapply(numis, function(x) { c(u = x) })
    nkeep1 <- downsample_plotdata(ncol(sample_subset), nmax)
    keep1 <- sort(sample(ncol(sample_subset), nkeep1))
    plot1_data <- list(plot1_data[keep1])
    plot_data <- append(plot_data, plot1_data)
  }

  # barcode ranks plot data
  # unique rank/fdr
  ranks <- rank(-numis, ties.method = "average")
  fdrs <- sample_subset$emptyDrops_FDR[ord]
  fdrs[is.na(fdrs)] <- 1

  dt <- data.table::data.table(rank = ranks, fdr = fdrs, u = numis)
  dt <- dt[, .(ndrops = .N, u = u[1]), by = c('rank', 'fdr')]

  # example of what knee regions should look like
  # plot_knee_regions(dt)

  # remove fdr specific columns for cell size filter
  if (is.cellsize) dt[, c("fdr", "ndrops") := NULL]
  knee_data <- unname(as.list(data.table::transpose(dt)))
  knee_data <- lapply(knee_data, function(x) {names(x) <- names(dt); x})
  plot_data[['knee']] <- knee_data

  # umi histogram plot
  if (is.cellsize) {
    hist_data <- lapply(numis, function(x) { c(u = x) })
    nhist <- downsample_plotdata(ncol(sample_subset), nmax)
    keep <- sort(sample(ncol(sample_subset), nhist))
    hist_data <- hist_data[keep]
    plot_data[['hist']] <- hist_data
  }

  plot_data <- append(plot_data, plot2_data)
  return(plot_data)
}

plot_knee_regions <- function(dt, thresh = 0.01) {

  # fill regions
  dt$quality <- 'unknown'
  lt.thresh <- dt$fdr < thresh

  amb.first <- which(!lt.thresh)[1]
  amb.last <- tail(which(lt.thresh), 1)

  dt$quality[1:amb.first] <- 'good'
  dt$quality[amb.last:nrow(dt)] <- 'low'

  cols <- c('green', 'red', 'gray')
  res$quality <- factor(res$quality)

  # plot
  ggplot2::ggplot(dt, aes(x=rank, y=u, fill = quality)) +
    ggplot2::geom_ribbon(mapping = ggplot2::aes(x=rank, ymax=u, ymin=0, fill = quality),alpha=0.3, inherit.aes = FALSE) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::geom_line(size = 1) +
    ggplot2::scale_x_continuous(trans='log10') +
    ggplot2::scale_y_continuous(trans='log10') +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = 'none') +
    ggplot2::theme(panel.border = ggplot2::element_rect(color = 'black', size = 0.1, fill = 'transparent')) +
    ggplot2::xlab('cell rank') +
    ggplot2::ylab('log UMIs in cell')
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
