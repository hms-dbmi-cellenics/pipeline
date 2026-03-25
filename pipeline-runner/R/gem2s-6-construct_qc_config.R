#' Constructs default QC configuration
#'
#' This function returns the default parameters used during QC as a nested list.
#' It is sent to the API, which in turn saves it as a jsonb object in the
#' PostgreSQL database.
#'
#' For each pipeline step, it customizes and adds the QC parameters for all samples
#' in the running experiment.
#'
#' @param scdata_list list of Seurat objects
#' @param unfiltered_samples character vector of unfiltered sample ids
#'
#' @return list of QC configuration parameters
#'
construct_qc_config <- function(scdata_list, unfiltered_samples, technology) {
  samples <- names(scdata_list)
  config_classifier <-
    add_custom_config_per_sample(
      customize_classifier_config,
      processing_config_template[["classifier"]],
      scdata_list,
      unfiltered_samples
    )

  config_cell_size <-
    add_custom_config_per_sample(
      customize_cellsize_config,
      processing_config_template[["cell_size"]],
      scdata_list,
    )

  config_mitochondrial <-
    add_custom_config_per_sample(
      customize_mitochondrial_config,
      processing_config_template[["mitochondrial"]],
      scdata_list,
    )

  config_genes_vs_umis <-
    add_custom_config_per_sample(
      customize_genes_vs_umis_config,
      processing_config_template[["genes_vs_umis"]],
      scdata_list,
      technology = technology
    )


  config_doublet <-
    add_custom_config_per_sample(
      customize_doublet_config,
      processing_config_template[["doublet"]],
      scdata_list,
    )

  tryCatch({
    ncells_unfiltered <- estimate_unfiltered_cells(
      scdata_list,
      config_classifier,
      config_cell_size,
      config_mitochondrial,
      config_genes_vs_umis,
      config_doublet
    )

    message("Estimated total cells after QC filtering: ")
    print(ncells_unfiltered)

  }, error = function(e) {
    message("Error estimating total cells after QC filtering: ", e$message)
    num_cells <- NA
  })

  config_data_integration <-
    get_data_integration_config(
      scdata_list,
      processing_config_template[["data_integration"]]
    )

  config_embedding_clustering <-
    get_embedding_config(
      scdata_list,
      processing_config_template[["embedding_clustering"]])

  # combine config for all steps
  config <- list(
    cellSizeDistribution = config_cell_size,
    mitochondrialContent = config_mitochondrial,
    classifier = config_classifier,
    numGenesVsNumUmis = config_genes_vs_umis,
    doubletScores = config_doublet,
    dataIntegration = config_data_integration,
    configureEmbedding = config_embedding_clustering
  )

  return(config)
}


customize_classifier_config <-
  function(scdata,
           config,
           sample_name,
           unfiltered_samples,
           technology) {
    config$enabled <- sample_name %in% unfiltered_samples
    config$prefiltered <- !(sample_name %in% unfiltered_samples)

    return(config)
  }


customize_cellsize_config <-
  function(scdata,
           config,
           sample_name,
           unfiltered_samples,
           technology) {
    minCellSize <- generate_default_values_cellSizeDistribution(scdata, config)
    config$filterSettings$minCellSize <- minCellSize
    return(config)
  }


customize_mitochondrial_config <-
  function(scdata,
           config,
           sample_name,
           unfiltered_samples,
           technology) {
    default_max_fraction <- generate_default_values_mitochondrialContent(scdata, config)
    config$filterSettings$methodSettings$absoluteThreshold$maxFraction <-
      default_max_fraction

    return(config)
  }


customize_doublet_config <-
  function(scdata,
           config,
           sample_name,
           unfiltered_samples,
           technology) {
    probabilityThreshold <- generate_default_values_doubletScores(scdata)
    config$filterSettings$probabilityThreshold <- probabilityThreshold

    return(config)
  }

customize_genes_vs_umis_config <-
  function(scdata,
           config,
           sample_name,
           unfiltered_samples,
           technology) {
    # Sensible values are based on the function "gene.vs.molecule.cell.filter"
    # from the pagoda2 package
    p.level <- min(0.001, 1 / ncol(scdata))

    regression_type <- ifelse( technology == "parse" , "spline" , config$filterSettings$regressionType)

    config$filterSettings$regressionType <- regression_type
    config$filterSettings$regressionTypeSettings[[regression_type]]$p.level <- p.level

    return(config)
  }


get_embedding_config <- function(scdata_list, config) {
  # tsne parameters depend on number of cells in sample
  default_perplexity <- config$embeddingSettings$methodSettings$tsne$perplexity
  default_learning_rate <- config$embeddingSettings$methodSettings$tsne$learningRate

  config$embeddingSettings$methodSettings$tsne <- list(
    perplexity = min(
      default_perplexity,
      min(vapply(scdata_list, ncol, numeric(1))) / 100),
    learningRate = max(
      default_learning_rate,
      min(vapply(scdata_list, ncol, numeric(1))) / 12)
  )

  return(config)
}

get_data_integration_config <- function(scdata_list, config) {
  FDR <- processing_config_template[["classifier"]]$filterSettings$FDR
  GEOSKETCH_CELLS_THRESHOLD <- 100000

  # TODO: calculate actual filtered cells from full config (mitochondrial, doublet, genes vs umis)
  message("Calculating total cells across samples to determine if Geosketch is needed for data integration...")

  total_cells <- 0
  # Calculate total cells after estimated empty drops filtering
  for (sample_name in names(scdata_list)) {
    scdata_sample <- scdata_list[[sample_name]]

    if (!scdata_sample@tools$flag_filtered) {
      # uses the same filtering logic as qc-1-filter_emptydrops.R
      ed_fdr <- scdata_sample$emptyDrops_FDR
      ed_fdr[is.na(ed_fdr)] <- 1  # prevents filtering of NA FDRs if FDR=1
      cells_in_sample <- sum(ed_fdr <= FDR, na.rm = TRUE)

    } else {
      cells_in_sample <- ncol(scdata_sample)
    }

    total_cells <- total_cells + cells_in_sample
  }

  message("Total cells after estimated emptyDrops filtering: ", total_cells)

  # Enable geosketch if total cells exceed threshold
  if (total_cells > GEOSKETCH_CELLS_THRESHOLD) {
    message("Total cells exceed threshold for Geosketch, enabling downsampling...")
    # Calculate percentageToKeep to downsample to approximately 500000 cells
    percentageToKeep <- 500000 / total_cells * 100

    config$downsampling <- list(
      method = "geosketch",
      methodSettings = list(
        geosketch = list(percentageToKeep = percentageToKeep)
      )
    )
  }

  return(config)
}


#' Estimate number of cells remaining after QC filtering
#'
#' Predicts how many cells will remain after applying all QC filters
#' sequentially (classifier, cell size, mitochondrial, genes vs UMIs, doublet)
#'
#' @param scdata_list list of Seurat objects
#' @param classifier_config config for classifier step
#' @param cell_size_config config for cell size step
#' @param mitochondrial_config config for mitochondrial step
#' @param genes_vs_umis_config config for genes vs UMIs step
#' @param doublet_config config for doublet step
#'
#' @return list with:
#'   - per_sample: data.frame with sample_name and cells_remaining columns
#'   - total: total cells across all samples after filtering
#'
#' @keywords internal
#'
estimate_unfiltered_cells <- function(
  scdata_list,
  classifier_config,
  cell_size_config,
  mitochondrial_config,
  genes_vs_umis_config,
  doublet_config
) {
  total_remaining <- 0
  per_sample_list <- list()

  for (sample_name in names(scdata_list)) {
    scdata <- scdata_list[[sample_name]]
    current_mask <- rep(TRUE, ncol(scdata))

    # Classifier (emptyDrops FDR filter)
    if (classifier_config[[sample_name]]$enabled &&
          !is.null(scdata@meta.data$emptyDrops_FDR)) {
      ed_fdr <- scdata$emptyDrops_FDR
      ed_fdr[is.na(ed_fdr)] <- 1
      current_mask <- current_mask &
        (ed_fdr <= classifier_config[[sample_name]]$filterSettings$FDR)
    }

    # Cell size filter
    if (cell_size_config[[sample_name]]$enabled) {
      min_cell_size <- as.numeric(
        cell_size_config[[sample_name]]$filterSettings$minCellSize
      )
      current_mask <- current_mask & (scdata$nCount_RNA >= min_cell_size)
    }

    # Mitochondrial content filter
    if (mitochondrial_config[[sample_name]]$enabled &&
          "percent.mt" %in% colnames(scdata@meta.data)) {
      filter_settings <- mitochondrial_config[[sample_name]]$filterSettings
      method <- filter_settings$method
      method_settings <- filter_settings$methodSettings[[method]]
      max_fraction <- method_settings$maxFraction
      current_mask <- current_mask & (scdata$percent.mt <= max_fraction * 100)
    }

    # Genes vs UMIs filter
    if (genes_vs_umis_config[[sample_name]]$enabled) {
      type <- genes_vs_umis_config[[sample_name]]$filterSettings$regressionType
      prediction_interval <- genes_vs_umis_config[[sample_name]]$filterSettings$predictionInterval

      if (safeTRUE(genes_vs_umis_config[[sample_name]]$auto)) {
        p_level <- min(0.001, 1 / sum(current_mask))
      } else if (!is.null(prediction_interval)) {
        p_level <- 1 - as.numeric(prediction_interval)
      } else {
        p_level <- min(0.001, 1 / sum(current_mask))
      }

      fit_data <- data.frame(
        log_molecules = log10(scdata$nCount_RNA[current_mask]),
        log_genes = log10(scdata$nFeature_RNA[current_mask]),
        cell_idx = which(current_mask)
      )

      if (nrow(fit_data) > 1) {
        fit_data <- fit_data[order(fit_data$log_molecules), ]
        tryCatch({
          fit <- if (type == "spline") {
            lm(log_genes ~ splines::bs(log_molecules), data = fit_data)
          } else {
            MASS::rlm(log_genes ~ log_molecules, data = fit_data)
          }

          preds <- suppressWarnings(
            predict(fit, interval = "prediction", level = 1 - p_level)
          )
          is_outlier <-
            fit_data$log_genes > preds[, "upr"] |
            fit_data$log_genes < preds[, "lwr"]

          current_mask[fit_data$cell_idx[is_outlier]] <- FALSE
        }, error = function(e) {})
      }
    }

    # Doublet filter
    sample_doublet_config <- doublet_config[[sample_name]]
    if (sample_doublet_config$enabled &&
          !is.null(scdata@meta.data$doublet_scores)) {
      doublet_scores <- scdata$doublet_scores
      doublet_scores[is.na(doublet_scores)] <- 0
      threshold <- sample_doublet_config$filterSettings$probabilityThreshold
      current_mask <- current_mask & (doublet_scores <= threshold)
    }

    cells_remaining <- sum(current_mask)
    per_sample_list[[sample_name]] <- cells_remaining
    total_remaining <- total_remaining + cells_remaining
  }

  per_sample_df <- data.frame(
    sample_name = names(per_sample_list),
    cells_remaining = unlist(per_sample_list),
    row.names = NULL
  )

  return(list(per_sample = per_sample_df, total = total_remaining))
}


#' Customize configuration for each sample in an experiment
#'
#' Takes the step config template and the required data and calls the corresponding
#' step customization function, which takes care of changing the QC parameters to
#' sensible values for each sample in the running experiment.
#'
#' Values that apply to all experiments are defined in the `data-raw/processing_config_template.R`
#' file.
#'
#' @param customize_template_config function - step customization function
#' @param config_template list - template of step configuration parameters
#' @param scdata_list list - with Seurat objects
#' @param unfiltered_samples character vector
#'
#' @return list of customized QC parameters for each sample
#' @export
#'
add_custom_config_per_sample <-
  function(customize_template_config,
           config_template,
           scdata_list,
           unfiltered_samples = NA,
           technology = NA) {
    config <- list()
    for (sample_name in names(scdata_list)) {
      # subset the Seurat object list to a single sample
      sample_scdata <- scdata_list[[sample_name]]

      # run the function to generate config for a sample
      sample_config <-
        customize_template_config(
          sample_scdata,
          config_template,
          sample_name,
          unfiltered_samples,
          technology
        )

      # update sample config thresholds
      config[[sample_name]] <- sample_config
    }

    return(config)
  }
