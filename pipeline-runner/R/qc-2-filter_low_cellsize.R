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
filter_low_cellsize <- function(scdata, config,sample_id,cells_id ,task_name = "cellSizeDistribution", num_cells_to_downsample = 6000) {
  cells_id.sample <- cells_id[[sample_id]]

  if (length(cells_id.sample) == 0) {
    guidata <- list()
    return(list(data = scdata, config = config, plotData = guidata))
  }

  scdata.sample <- subset_ids(scdata,cells_id.sample)

  minCellSize <- as.numeric(config$filterSettings$minCellSize)

  # extract plotting data of original data to return to plot slot later
  plot_data <- get_bcranks_plot_data(scdata.sample, num_cells_to_downsample)

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
      if (ncol(scdata.sample) < threshold.low) {
        minCellSize <- min(scdata.sample$nCount_RNA)
      } else {
        minCellSize <- generate_default_values_cellSizeDistribution(scdata.sample, config)
      }
    }
  }

  if (as.logical(toupper(config$enabled))) {
    remaining_ids <- scdata.sample@meta.data$cells_id[scdata.sample$nCount_RNA >= minCellSize]
  } else {
    remaining_ids <- cells_id.sample
  }
  # update config
  config$filterSettings$minCellSize <- minCellSize
  # Populate data for UI
  guidata <- list()
  guidata[[generate_gui_uuid(sample_id, task_name, 0)]] <- plot_data[['knee']]
  guidata[[generate_gui_uuid(sample_id, task_name, 1)]] <- plot_data[['hist']]
  # Populate with filter statistics
  filter_stats <- list(
    before = calc_filter_stats(scdata.sample),
    after = calc_filter_stats(subset_ids(scdata.sample,remaining_ids))
  )
  guidata[[generate_gui_uuid(sample_id, task_name, 3)]] <- filter_stats

  cells_id[[sample_id]] <- remaining_ids

  result <- list(
    data = scdata,
    new_ids = cells_id,
    config = config,
    plotData = guidata
  )

  return(result)
}

#' Get Barcode Ranks Plot Data
#'
#' @param sample_subset \code{SeuratObject}
#' @param nmax maximum number of value for histogram plot. Only applies if \code{is.cellsize}
#'   is \code{TRUE}.
#' @param is.cellsize is this for cellsize filter? Default (\code{TRUE}),
#'   returns histogram plot data and knee plot data without fdr and ndrops values.
#'
#' @keywords internal
#' @return Named list of plot data. Possible items are \code{'hist'} and \code{'knee'},
#'   dependending on \code{is.cellsize}.
#' @importFrom magrittr %>%
#'
get_bcranks_plot_data <- function(sample_subset, nmax = 6000, is.cellsize = TRUE) {

  set.seed(123)
  plot_data <- list()
  numis <- unname(sample_subset$nCount_RNA)
  ord <- order(numis, decreasing = TRUE)
  numis <- numis[ord]

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
  knee_data <- dt %>% purrr::transpose() %>% purrr::simplify_all()
  plot_data[['knee']] <- knee_data

  # umi histogram plot
  if (is.cellsize) {
    hist_data <- lapply(numis, function(x) { c(u = x) })
    nhist <- downsample_plotdata(ncol(sample_subset), nmax)
    keep <- sort(sample(ncol(sample_subset), nhist))
    hist_data <- hist_data[keep]
    plot_data[['hist']] <- hist_data
  }

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
generate_default_values_cellSizeDistribution <- function(scdata, config) {
  # `CalculateBarcodeInflections` including inflection point calculation
  threshold.low <- if(ncol(scdata)<=200) NULL else 100
  scdata_tmp <- Seurat::CalculateBarcodeInflections(
    object = scdata,
    barcode.column = "nCount_RNA",
    group.column = "samples",
    threshold.low = threshold.low
  )
  # returned is both the rank(s) as well as inflection point
  # extracting only inflection point(s)
  sample_subset <- Seurat::Tool(scdata_tmp, slot = "Seurat::CalculateBarcodeInflections")$inflection_points
  return(sample_subset$nCount_RNA)
}
