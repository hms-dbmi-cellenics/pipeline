#' STEP 4. Number of genes vs UMIs filter
#'
#' Eliminates cells based on a p value and a linear regression generated from numGenes vs numUmis
#'
#' This filter focuses on filter cells that are far from the behaviour of the relationship between the number of genes (it measures the number of
#' genes in a cell that has at least one count) and the number of UMIs/molecules (the total number of counts in a cell).
#'
#' @param config list containing the following information
#'          - enable: true/false. Refering to apply or not the filter.
#'          - auto: true/false. 'True' indicates that the filter setting need to be changed depending on some sensible value (it requires
#'          to call generate_default_values_numGenesVsNumUmis)
#'          - filterSettings: slot with thresholds
#'              - regressionType: String. Regression to be used: {gam}
#'              - regressionTypeSettings: list with the config settings for all the regression type options
#'                          - gam: for the gam option there is only one element:
#'                             - p.level: which refers to  confidence level for deviation from the main trend
#'
#' @param scdata \code{SeuratObject}
#' @param sample_id value in \code{scdata$samples} to apply filter for
#' @param task_name name of task: \code{'numGenesVsNumUmis'}
#' @param num_cells_to_downsample maximum number of cells for returned plots
#' @export
#'
#' @return a list with the filtered seurat object by numGenesVsNumUmis, the config and the plot values
#'
filter_gene_umi_outlier <- function(scdata, config, sample_id, task_name = "numGenesVsNumUmis", num_cells_to_downsample = 6000) {

  type <- config$filterSettings$regressionType

  # get p.level and update in config
  # defaults from "gene.vs.molecule.cell.filter" in pagoda2
  if (!is.null(config$auto) && config$auto)
    p.level <- min(0.001, 1 / ncol(scdata))
  else
    p.level <- config$filterSettings$regressionTypeSettings[[type]]$p.level

  config$filterSettings$regressionTypeSettings[[type]]$p.level <- p.level

  # defaults if not enabled
  scdata.filtered <- scdata
  plot_data <- list()

  # barcodes for this sample
  sample_barcodes <- colnames(scdata)[scdata$samples == sample_id]

  if (config$enabled && length(sample_barcodes)) {

    # subset this sample
    sample_subset <- scdata[, sample_barcodes]

    # regress log10 molecules vs genes
    data <- data.frame(
      log_molecules = log10(sample_subset$nCount_RNA),
      log_genes = log10(sample_subset$nFeature_RNA),
      row.names = colnames(sample_subset)
    )

    data <- data[order(data$log_molecules), ]
    fit <- lm(log_genes ~ splines::bs(log_molecules), data = data)

    # get the interval based on p.level parameter
    preds <- suppressWarnings(predict(fit, interval = "prediction", level = 1 - p.level))

    # filter outliers above/below cutoff bands
    is.outlier <- data$log_genes > preds[, 'upr'] | data$log_genes < preds[, 'lwr']
    outliers <- rownames(data)[is.outlier]

    keep <- setdiff(row.names(data), outliers)
    scdata.filtered <- subset_safe(scdata, keep)

    # get evenly spaced predictions for plotting lines
    xrange <- range(data$log_molecules)
    newdata <- data.frame(log_molecules = seq(xrange[1], xrange[2], length.out = 10))
    line_preds <- predict(fit, newdata, interval = "prediction", level = 1 - p.level)

    line_preds <- cbind(newdata, line_preds) %>%
      dplyr::select(-fit) %>%
      dplyr::rename(lower_cutoff = lwr, upper_cutoff = upr)

    # downsample plot data
    nkeep <- downsample_plotdata(ncol(sample_subset), num_cells_to_downsample)

    set.seed(123)
    keep_rows <- sample(nrow(data), nkeep)
    keep_rows <- sort(keep_rows)
    plot_data <- list(
      pointsData = purrr::transpose(data[keep_rows, ]),
      linesData = purrr::transpose(line_preds)
    )
  }

  # Populate with filter statistics and plot data
  filter_stats <- list(
    before = calc_filter_stats(scdata, sample_id),
    after = calc_filter_stats(scdata.filtered, sample_id)
  )

  guidata <- list()
  guidata[[generate_gui_uuid(sample_id, task_name, 0)]] <- plot_data
  guidata[[generate_gui_uuid(sample_id, task_name, 1)]] <- filter_stats

  result <- list(
    data = scdata.filtered,
    config = config,
    plotData = guidata
  )

  return(result)
}
