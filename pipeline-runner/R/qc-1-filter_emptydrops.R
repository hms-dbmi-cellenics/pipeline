# STEP 1. Classifier filter

#' Filter empty droplets
#'
#' Filters seurat object based on scores estimated by the classifier implemented
#'  in `DropletUtils::emptyDrops`
#'
#' @param config list containing the following information
#'          - enable: true/false. Refering to apply or not the filter.
#'          - auto: true/false. 'True' indicates that the filter setting need to be changed depending on some sensible value (it requires
#'          to call generate_default_values_mitochondrialContent)
#'          - filterSettings: slot with thresholds
#'                  - methodSettings: List with the method as key and contain all the filterSettings for this specific method.
#'                      - FDR: False Discovery Rate. Cells with FDR above this threshold will be filtered out.
#' @import data.table
#' @export
#' @return a list with the filtered seurat object by mitochondrial content, the config and the plot values
#'
filter_emptydrops <- function(scdata_list, config, sample_id, cells_id, task_name = "classifier", num_cells_to_downsample = 6000) {
  sample_cell_ids <- cells_id[[sample_id]]
  message("Number of cells IDs: ", length(sample_cell_ids))
  message("Number of cells: ", ncol(scdata_list[[sample_id]]))

  if (length(sample_cell_ids) == 0) {
    return(list(data = scdata_list[[sample_id]], new_ids = cells_id, config = config, plotData = list()))
  }

  # TODO this is probably not needed because it's the first filter and it has not been yet applied
  # so the cells_id are the ones generate in gem2s
  sample_data <- subset_ids(scdata_list[[sample_id]], sample_cell_ids)

  FDR <- config$filterSettings$FDR
  if (isTRUE(config$auto)) {
    FDR <- processing_config_template[["classifier"]]$filterSettings$FDR
  }

  # for plots and filter stats to be populated
  guidata <- list()

  # check if filter data is actually available
  if (is.null(scdata_list[[sample_id]]@meta.data$emptyDrops_FDR)) {
    message("Classify is enabled but has no classify data available: will dissable it: no filtering!")
    config$enabled <- FALSE
  }

  if (config$enabled) {
    message("Classify is enabled: filtering with FDR=", FDR)
    meta <- sample_data@meta.data

    message("Info: empty-drops table of FDR threshold categories (# UMIs for a given threshold interval)")
    print(table(meta$samples, cut(meta$emptyDrops_FDR, breaks = c(-Inf, 0, 0.0001, 0.01, 0.1, 0.5, 1, Inf)), useNA = "ifany"))

    # prevents filtering of NA FDRs if FDR=1
    ed_fdr <- sample_data$emptyDrops_FDR
    ed_fdr[is.na(ed_fdr)] <- 1

    message(
      "Number of barcodes to filter for this sample: ",
      sum(ed_fdr > FDR, na.rm = TRUE), "/", length(ed_fdr)
    )

    numis <- log10(sample_data@meta.data$nCount_RNA)

    fdr_data <- unname(purrr::map2(ed_fdr, numis, function(x, y) {
      c("FDR" = x, "log_u" = y)
    }))
    fdr_data <- fdr_data[get_positions_to_keep(sample_data, num_cells_to_downsample)]

    remaining_ids <- sample_data@meta.data$cells_id[ed_fdr <= FDR]

    # update config
    config$filterSettings$FDR <- FDR

    # Downsample plotData
    knee_data <- get_bcranks_plot_data(sample_data, is.cellsize = FALSE)[["knee"]]

    # Populate guidata list
    guidata[[generate_gui_uuid(sample_id, task_name, 0)]] <- fdr_data
    guidata[[generate_gui_uuid(sample_id, task_name, 1)]] <- knee_data
  } else {
    message("filter disabled: data not filtered!")
    # guidata is an empty list
    guidata[[generate_gui_uuid(sample_id, task_name, 0)]] <- list()
    guidata[[generate_gui_uuid(sample_id, task_name, 1)]] <- list()
    guidata[[generate_gui_uuid(sample_id, task_name, 2)]] <- list()
    remaining_ids <- sample_cell_ids
  }

  # get filter stats after filtering
  filter_stats <- list(
    before = calc_filter_stats(sample_data),
    after = calc_filter_stats(subset_ids(sample_data, remaining_ids))
  )

  guidata[[generate_gui_uuid(sample_id, task_name, 2)]] <- filter_stats

  cells_id[[sample_id]] <- remaining_ids

  result <- list(
    data = scdata_list,
    new_ids = cells_id,
    config = config,
    plotData = guidata
  )

  return(result)
}

