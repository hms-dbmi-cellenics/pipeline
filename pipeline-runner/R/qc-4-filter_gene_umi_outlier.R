#' STEP 4. Number of genes vs UMIs filter
#'
#' Eliminates cells based on a p value and a linear regression generated from numGenes vs numUmis
#'
#' This filter focuses on filter cells that are far from the behaviour of the
#' relationship between the number of genes (it measures the number of
#' genes in a cell that has at least one count) and the number of UMIs/molecules
#' (the total number of counts in a cell).
#' The cutoff bands for filtering are derived from the prediction interval.
#' Prediction interval represents the likelihood that the predicted value will be
#' between the upper and lower limits of the prediction interval.
#' Prediction intervals are similar to confidence intervals, but on top of the sampling
#' uncertainty, they also express uncertainty around a single value.
#' They must account for the uncertainty in estimating the population mean,
#' plus the random variation of the individual values.
#' Higher prediction interval means higher probability of the value to be inside
#' the range. Consequently, the size of the interval will be wider.
#' The higher the prediction level, the less stringent we are when filtering the cells.
#' Conversely, the lower the prediction level, the more stringent we are,
#' and we exclude more cells that are far from the behaviour of the relationship
#' between the number of genes and the number of UMIs/molecules.
#'
#' @param config list containing the following information
#'          - enable: true/false. Referring to apply or not the filter.
#'          - auto: true/false. 'True' indicates that the filter setting need to be changed depending on some sensible value (it requires
#'          to call generate_default_values_numGenesVsNumUmis)
#'          - filterSettings: slot with thresholds
#'              - regressionType: String. Regression to be used: {linear or spline}
#'              - regressionTypeSettings: list with the config settings for all the regression type options
#'                          - linear and spline: for each there is only one element:
#'                             - p_level: which refers to  confidence level for deviation from the main trend
#'
#' @param scdata_list list of \code{SeuratObject}
#' @param sample_id value in \code{scdata$samples} to apply filter for
#' @param task_name name of task: \code{'numGenesVsNumUmis'}
#' @param num_cells_to_downsample maximum number of cells for returned plots
#' @export
#'
#' @return a list with the filtered seurat object by numGenesVsNumUmis, the config and the plot values
#'
#'
filter_gene_umi_outlier <- function(scdata_list, config, sample_id, cells_id, task_name = "numGenesVsNumUmis", num_cells_to_downsample = 6000) {
  sample_cell_ids <- cells_id[[sample_id]]

  if (length(sample_cell_ids) == 0) {
    return(list(data = scdata_list[[sample_id]], new_ids = cells_id, config = config, plotData = list()))
  }

  sample_data <- subset_ids(scdata_list[[sample_id]], sample_cell_ids)

  type <- config$filterSettings$regressionType

  # get p_level and update in config
  # defaults from "gene.vs.molecule.cell.filter" in pagoda2
  if (safeTRUE(config$auto))
    p_level <- min(0.001, 1 / ncol(sample_data))
  else
    p_level <- config$filterSettings$regressionTypeSettings[[type]]$p.level

  p_level <- suppressWarnings(as.numeric(p_level))
  if(is.na(p_level)) stop("p_level couldn't be interpreted as a number.")

  pred_int_auto <- 1 - p_level

  config$filterSettings$regressionTypeSettings[[type]]$p.level <- p_level

  # Assign updated config to global env so that it can be accessed if there is an error
  config_key <- paste0("config-", task_name, "-", sample_id)
  assign(config_key, config, envir = globalenv())

  # regress log10 molecules vs genes
  fit.data <- data.frame(
    log_molecules = log10(sample_data$nCount_RNA),
    log_genes = log10(sample_data$nFeature_RNA),
    row.names = sample_data$cells_id
  )

  fit.data <- fit.data[order(fit.data$log_molecules), ]

  if (type == "spline") {
    fit <- lm(log_genes ~ splines::bs(log_molecules), data = fit.data)
  } else {
    fit <- MASS::rlm(log_genes ~ log_molecules, data = fit.data)
  }

  if (safeTRUE(config$enabled)) {
    # get the interval based on p_level parameter
    preds <- suppressWarnings(predict(fit, interval = "prediction", level = 1 - p_level))

    # filter outliers above/below cutoff bands
    is.outlier <- fit.data$log_genes > preds[, "upr"] | fit.data$log_genes < preds[, "lwr"]
    remaining_ids <- as.numeric(rownames(fit.data)[!is.outlier])
    remaining_ids <- remaining_ids[order(remaining_ids)]
  } else {
    remaining_ids <- sample_cell_ids
  }

  # downsample for plot data
  nkeep <- downsample_plotdata(ncol(sample_data), num_cells_to_downsample)

  set.seed(RANDOM_SEED)
  keep_rows <- sample(nrow(fit.data), nkeep)
  keep_rows <- sort(keep_rows)
  downsampled_data <- fit.data[keep_rows, ]

  # get evenly spaced predictions on downsampled data for plotting lines
  xrange <- range(downsampled_data$log_molecules)
  newdata <- data.frame(log_molecules = seq(xrange[1], xrange[2], length.out = 10))
  line_preds <- suppressWarnings(predict(fit, newdata, interval = "prediction", level = 1 - p_level))

  pred_int_values <- sort(c(seq(0, 0.99, 0.01), 0.999, 0.9999, 0.99999, 0.999999))
  pred_int_values <- c(pred_int_values, pred_int_auto)
  line_preds_list <- list()
  i <- 1
  for (pred in pred_int_values) {
    line_preds <- suppressWarnings(predict(fit, newdata, interval = "prediction", level = pred))

    line_preds <- cbind(newdata, line_preds) %>%
      dplyr::select(-fit) %>%
      dplyr::rename(lower_cutoff = lwr, upper_cutoff = upr)

    line_preds_list[[i]] <- purrr::transpose(line_preds)
    i <- i + 1
  }

  plot_data <- list(
    pointsData = purrr::transpose(downsampled_data),
    linesData = line_preds_list
  )


  # Populate with filter statistics and plot data
  filter_stats <- list(
    before = calc_filter_stats(sample_data),
    after = calc_filter_stats(subset_ids(sample_data, remaining_ids))
  )

  guidata <- list()
  guidata[[generate_gui_uuid(sample_id, task_name, 0)]] <- plot_data
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
