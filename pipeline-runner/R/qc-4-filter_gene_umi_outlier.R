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

  p.level <- as.numeric(config$filterSettings$regressionTypeSettings[[config$filterSettings$regressionType]]$p.level)
if(is.na(p.level)) stop("P-level couldnt be interpreted as a number.")
  if (as.logical(toupper(config$enabled))) {
    if (config$filterSettings$regressionType == "gam") {
      # Subsetting this sample
      meta <- scdata@meta.data
      barcode_names_this_sample <- rownames(meta)[meta$samples == sample_id]
      if (length(barcode_names_this_sample) == 0) {
        return(list(data = scdata, config = config, plotData = list()))
      }
      sample_subset <- subset(scdata, cells = barcode_names_this_sample)

      # Check if it is required to compute sensible values. Sensible values are based on the funciton "gene.vs.molecule.cell.filter" from the pagoda2 package
      if (exists("auto", where = config)) {
        if (as.logical(toupper(config$auto))) {
          p.level <- min(0.001, 1 / ncol(scdata))
        }
      }

      # We regress the molecules vs the genes. This information are stored in nCount_RNA and nFeature_RNA respectively
      df <- data.frame(molecules = sample_subset$nCount_RNA, genes = sample_subset$nFeature_RNA)
      # We take log10 following the plot from the mock-up
      df <- log10(df)
      # Rename the rows to be able to identify the valid cells
      rownames(df) <- colnames(sample_subset)
      df <- df[order(df$molecules, decreasing = FALSE), ]
      m <- MASS::rlm(genes ~ molecules, data = df)
      # Get the interval based on p.level paramter
      suppressWarnings(pb <- data.frame(predict(m,
        interval = "prediction",
        level = 1 - p.level, type = "response"
      )))
      # Define the outliers those that are below the lower confidence band and above the upper one.
      outliers <- rownames(df)[df$genes > pb$upr | df$genes < pb$lwr]

      # Keep the ones that are not outlier
      # Because we have subseted the object before creating the dataframes, and just selected the outlier values from that df
      # In the end we wil only remove values that are in the desired sample.
      scdata.filtered <- subset_safe(scdata, colnames(scdata)[!colnames(scdata) %in% outliers])
      # Similarly, when creating the plot from the df, we will just be using the values for the required sample.
      plot1_data <- unname(purrr::map2(df$molecules, df$genes, function(x, y) {
        c("log_molecules" = x, "log_genes" = y)
      }))
      plot1_data <- purrr::map2(plot1_data, unname(pb$lwr), function(x, y) {
        append(x, c("lower_cutoff" = y))
      })
      plot1_data <- purrr::map2(plot1_data, unname(pb$upr), function(x, y) {
        append(x, c("upper_cutoff" = y))
      })
      # have downsampling done

      # Downsample plotData
      num_cells_to_downsample <- downsample_plotdata(ncol(sample_subset), num_cells_to_downsample)

      set.seed(123)
      cells_position_to_keep <- sample(1:ncol(sample_subset), num_cells_to_downsample, replace = FALSE)
      cells_position_to_keep <- sort(cells_position_to_keep)
      plot1_data <- plot1_data[cells_position_to_keep]
    }
  } else {
    scdata.filtered <- scdata
    plot1_data <- list()
  }

  config$filterSettings$regressionTypeSettings[[config$filterSettings$regressionType]][["p.level"]] <- p.level

  # Scatter plot which is composed of:
  # x-axis: log_10_UMIs
  # y-axis: log_10_genes
  # bands that are conformed with the upper_cutoff and the lower_cutoff.
  guidata <- list()
  guidata[[generate_gui_uuid(sample_id, task_name, 0)]] <- plot1_data

  # Populate with filter statistics
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
