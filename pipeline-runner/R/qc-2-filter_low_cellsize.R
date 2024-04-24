#' STEP 2. Cell size distribution filter
#'
#' cell_size_distribution_filter function
#'
#' This is a simplest filter that looks at the shape of the cell size (# of UMIs
#' per cell) distribution and looks for some local minima, minimum of second derivative, etc.
#' To separate the large cell barcodes that correspond to real cells from the
#' tail containing 'empty droplets'. This can be a useful first guess. The settings
#'for such a filter can also contain a simple "min cell size" setting.
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
filter_low_cellsize <- function(scdata_list, config, sample_id, cells_id, task_name = "cellSizeDistribution", num_cells_to_downsample = 6000) {
  sample_cell_ids <- cells_id[[sample_id]]

  if (length(sample_cell_ids) == 0) {
    return(list(data = scdata_list[[sample_id]], new_ids = cells_id, config = config, plotData = list()))
  }

  sample_data <- subset_ids(scdata_list[[sample_id]], sample_cell_ids)

  minCellSize <- as.numeric(config$filterSettings$minCellSize)

  # extract plotting data of original data to return to plot slot later
  plot_data <- get_bcranks_plot_data(sample_data, num_cells_to_downsample)

  # Check if it is required to compute sensible values. From the function 'generate_default_values_cellSizeDistribution', it is expected
  # to get a list with two elements {minCellSize and binStep}
  if (exists("auto", where = config)) {
    if (as.logical(toupper(config$auto))) {
      # config not really needed for this one (maybe later for threshold_low/high):
      # HARDCODE Value. threshold_low
      # [ Parameter for function CalculateBarcodeInflections. Description: Ignore barcodes of rank below this threshold in inflection calculation]
      threshold_low <- 1e2
      # If there are less cells than the value threshold_low, the function CalculateBarcodeInflections fails. So we need to handle by not removing any cells, that is,
      # consider the minCellSize as the minimun UMIs in the dataset.
      # This should be handled in a long-term by adding a different function for computing the default value.
      if (ncol(sample_data) < threshold_low) {
        minCellSize <- min(sample_data$nCount_RNA)
      } else {
        minCellSize <- generate_default_values_cellSizeDistribution(sample_data, config)
      }
    }
  }

  # update config
  config$filterSettings$minCellSize <- minCellSize
  message("MIN CELL SIZE: ", minCellSize)

  # Assign updated config to global env so that it can be accessed if there is an error
  config_key <- paste0("config-", task_name, "-", sample_id)
  assign(config_key, config, envir = globalenv())

  if (as.logical(toupper(config$enabled))) {
    remaining_ids <- sample_data@meta.data$cells_id[sample_data$nCount_RNA >= minCellSize]
  } else {
    remaining_ids <- sample_cell_ids
  }

  # Populate data for UI
  guidata <- list()
  guidata[[generate_gui_uuid(sample_id, task_name, 0)]] <- plot_data[["knee"]]
  guidata[[generate_gui_uuid(sample_id, task_name, 1)]] <- plot_data[["hist"]]
  # Populate with filter statistics
  filter_stats <- list(
    before = calc_filter_stats(sample_data),
    after = calc_filter_stats(subset_ids(sample_data, remaining_ids))
  )
  guidata[[generate_gui_uuid(sample_id, task_name, 3)]] <- filter_stats

  cells_id[[sample_id]] <- remaining_ids

  result <- list(
    data = scdata_list,
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
  set.seed(RANDOM_SEED)
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
  dt <- dt[, .(ndrops = .N, u = u[1]), by = c("rank", "fdr")]

  # remove fdr specific columns for cell size filter
  if (is.cellsize) dt[, c("fdr", "ndrops") := NULL]
  knee_data <- dt %>%
    purrr::transpose() %>%
    purrr::simplify_all()
  plot_data[["knee"]] <- knee_data

  # umi histogram plot
  if (is.cellsize) {
    hist_data <- lapply(numis, function(x) {
      c(u = x)
    })
    nhist <- downsample_plotdata(ncol(sample_subset), nmax)
    keep <- sort(sample(ncol(sample_subset), nhist))
    hist_data <- hist_data[keep]
    plot_data[["hist"]] <- hist_data
  }

  return(plot_data)
}


#' cell_size_distribution_filter function
#' This is a simplest filter that looks at the shape of the cell size (# of
#' UMIs per cell) distribution and looks for some local minima, minimum of
#' second derivative, etc.  To separate the large cell barcodes that correspond
#' to real cells from the tail containing 'empty droplets'. This can be a useful
#' first guess. The settings for such a filter can also contain a simple
#' "min cell size" setting.
#'
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

  result <- calc_count_cutoff(scdata)
  inflection_point_rank <- result$cutoff_rank
  inflection_point_UMI <- result$cutoff_count

  message("INFLECTION POINT: ", inflection_point_UMI)
  return(inflection_point_UMI)
}

mean_win_smooth <- function(vals, win_frac) {
  win_side <- round(length(vals) * win_frac / 2, 0)
  zoo::rollmean(vals, max(3, win_side * 2 + 1), fill = NA, align = "center")
}

first_derivative <- function(vals, win_frac) {
  win_side <- round(length(vals) * win_frac / 2, 0)
  # Ensure there's at least one element on each side for difference calculation
  if (win_side >= 1) {
    slopes <- (vals[-(1:win_side)] - vals[-((length(vals) - win_side + 1):length(vals))]) / (2 * win_side)
  } else {
    slopes <- c(NA, diff(vals), NA)
  }
  return(slopes)
}

umb_parab_tran_func <- function(n_steps, xfrac=0.5, yfrac=0.25, center_focus = TRUE) {
  x <- seq(0, 1, length.out = n_steps)
  y <- (x * 2 - 1)^2  # Parabolic transformation centered
  y_shift <- (1 - xfrac)^2

  # Apply additional centering focus if required
  if (center_focus) {
    y_shift <- y_shift + (x - 0.5)^2 * 4 * (1 - xfrac)  # Enhance the center focus
  }

  y <- pmax(y - y_shift, 0)
  y <- ifelse(y_shift < 1, y / (1 - y_shift), y)
  1 - y * yfrac
}

# Main function to calculate barcode rank plot cutoff
calc_count_cutoff <- function(scdata,
                              min_cnt = 30, max_cells = 30000,
                              min_cells = 10,
                              cv_smooth_win = 0.02, cv_slope_win = 0.02, center_focus = TRUE,
                              min_cell_auc = 0.4, adjust_max_cells = FALSE, parse_kit = "WT") {

  barcode_counts <- scdata@meta.data$nCount_RNA
  barcodes <- rownames(scdata@meta.data)

  df <- data.table::data.table(barcode = barcodes, count = barcode_counts)
  df <- df[order(-count)]
  df[, rank := frank(-count, ties.method = "first")]

  df$log_count <- log10(df$count + 1)
  df$log_rank <- log10(df$rank)

  # Adjust maximum cells based on sample size if needed
  if (adjust_max_cells) {
    # Define the total number of wells based on the kit type
    total_wells <- unlist(PARSE_KIT_WELLS[parse_kit])
    num_samples <- length(unique(scdata@meta.data$samples))
    # Calculate the estimated number of wells per sample
    wells_per_sample <- total_wells / num_samples
    # Scale max_cells based on the estimated number of wells per sample
    max_cells <- max_cells * (wells_per_sample / total_wells)
  }

  # Apply dynamic minimum transcript counts based on data characteristics
  if (min_cell_auc > 0) {
    total_transcripts <- sum(df$count)
    cumsum_perc <- cumsum(df$count) / total_transcripts
    min_cells <- max(min_cells, which(cumsum_perc >= min_cell_auc)[1])
  }

  # Interpolate the log-transformed rank and count data
  fit <- approxfun(df$log_rank, df$log_count, method = "linear")

  # Create a sequence for interpolation
  xout <- seq(min(df$log_rank), max(df$log_rank), length.out = nrow(df))
  yout <- fit(xout)

  # Smooth the interpolated values
  y_smooth <- mean_win_smooth(yout, cv_smooth_win)

  # Calculate the first derivative of the smoothed curve
  y_slope <- first_derivative(y_smooth, cv_slope_win)

  # Apply the umbrella-like transformation
  n_steps <- length(y_slope)
  weights <- umb_parab_tran_func(n_steps, xfrac = 0.5, yfrac = 0.25, center_focus = center_focus)
  weighted_slope <- y_slope * weights

  # Determine valid analysis window
  valid_indices <- which(df$count >= min_cnt & df$rank <= max_cells & df$rank >= min_cells)
  max_slope_index <- valid_indices[which.max(weighted_slope[valid_indices])]

  # Get the cutoff point
  cutoff_rank <- df$rank[max_slope_index]
  cutoff_count <- df$count[max_slope_index]

  list(cutoff_rank = cutoff_rank, cutoff_count = cutoff_count)
}
