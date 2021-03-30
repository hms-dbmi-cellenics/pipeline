################################################
## STEP 3. Classifier filter 
#################################################
#
#' @description Filters seurat object based on classifier filter using emptyDrops
#               https://rdrr.io/github/MarioniLab/DropletUtils/man/emptyDrops.html
#' @param config list containing the following information
#'          - enable: true/false. Refering to apply or not the filter.
#'          - auto: true/false. 'True' indicates that the filter setting need to be changed depending on some sensible value (it requires
#'          to call generate_default_values_classifier)
#'          - filterSettings: slot with thresholds
#'                  - minProbabiliy: 
#'                  - filterThreshold: 
#' @export return a list with the filtered seurat object by probabilities classifier, the config and the plot values

source('utils.r')

generate_default_values_classifier <- function(seurat_obj, config) {
   
        # HARDCODE
        threshold <- 0.1
   
    return(threshold)
}

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
#' @export return a list with the filtered seurat object by mitochondrial content, the config and the plot values


task <- function(seurat_obj, config, task_name, sample_id){
  # config$filterSettings = list(FDR=0.82, bandwidth=-1, filterThreshold=-1)
  # As of 26/3/2021
  # PlotUUID: classifierEmptyDropsPlot
  # Labels & Inputs : 
  # X-axis : cellSize, transformed with log10
  # Y-axis : FDR of emptyDrops classifier [0,1]   
  # Plot Data Schema
  # Data point :
  # { 
  #    classifierP: classifier probability,
  #    log_u: log cell size
  # }
  # Example :
  # [
  #    {"FDR": 0.994553522823595,
  #     "log_u": 4.25168687816849},
  # ...]
  print(paste("Running",task_name,"config: ",sep=" "))
  print(config)
  print(paste0("Cells per sample before filter for sample ", sample_id))
  print(table(seurat_obj$orig.ident))
  #The format of the sample_id is
  # sample-WT1
  # we need to get only the last part, in order to grep the object.
  tmp_sample <- sub("sample-","",sample_id)

  import::here(map2, .from = purrr)

  # Check wheter the filter is set to true or false
  # For some reason the last children of named lists are computed as vectors, so we can't access them as recursive objects. 
  FDR <- as.numeric(config$filterSettings[["FDR"]])
  # THIS SLOT DOESNT EXIST ON LOCAL OBJECT.
  if (is.null(seurat_obj@meta.data$emptyDrops_FDR)) {
    print(paste("Running",task_name,"config: ",sep=" "))
    print("Classify is enabled but has no classify data available: will dissable it: no filtering!")
    # should this be json, i.e. "false" instead?
    config$enabled <- FALSE
    plot1_data <- list()
  } else {

    # extract plotting data of original data to return to plot slot later
    obj_metadata <- seurat_obj@meta.data
    barcode_names_this_sample <- rownames(obj_metadata[grep(tmp_sample, rownames(obj_metadata)),]) 
    sample_subset <- subset(seurat_obj, cells = barcode_names_this_sample)

    plot1_data_1  <- seurat_obj@meta.data$emptyDrops_FDR
    # TODO: should this be log, or is the UI taking this?
    plot1_data_2  <- log10(seurat_obj@meta.data$nCount_RNA)
    # plot(= plot1_data_1,plot1_data_2)
    # plot(seurat_obj@meta.data$emptyDrops_FDR, seurat_obj@meta.data$nCount_RNA)

    plot1_data <- unname(purrr::map2(plot1_data_1,plot1_data_2,function(x,y){c("FDR"=x,"log_u"=y)}))
  }

  # Check if it is required to compute sensible values. From the function 'generate_default_values_classifier', it is expected
  # to get a list with two elements {minProbabiliy and filterThreshold}.
  if (exists('auto', where=config)){
    if(as.logical(toupper(config$auto)))
      FDR <- generate_default_values_classifier(seurat_obj, config)
  }
  # TODO: get flag from here: seurat_obj@tools$flag_filtered <- FALSE
  if(!as.logical(toupper(config$enabled))) {
    # is.cell <- meta.data$emptyDrops_FDR <= 0.01
    # sum(is.cell, na.rm=TRUE) 
    # table(Limited=meta.data$emptyDrops_Limited, Significant=is.cell)
    # is.cell2<-is.cell
    # is.cell2[is.na(is.cell2)]<-FALSE
    # sce.filt<-sce[,is.cell2]
    # Information regarding number of UMIs per cells is pre-computed during the 'CreateSeuratObject' function. 
    # this used to be the filter below which applies for all samples
    # we had to change it for the current version where we also have to 
    # keep all barcodes for all samples (filter doesn't apply there)
    # once we ensure the input is only one sample, we can revert to the line below:
    # seurat_obj <- subset(seurat_obj, subset = nCount_RNA >= minCellSize)
    # subset(x, subset, cells = NULL, features = NULL, idents = NULL, ...)
    # cells: Cell names or indices
    # extract cell id that do not(!) belong to current sample (to not apply filter there)
    barcode_names_non_sample <- rownames(obj_metadata[-grep(tmp_sample, rownames(obj_metadata)),])
    # all barcodes that match threshold in the subset data
    barcode_names_keep_current_sample <-rownames(sample_subset@meta.data[sample_subset@meta.data$emptyDrops_FDR <= FDR,])
    # combine the 2:
    barcodes_to_keep <- union(barcode_names_non_sample, barcode_names_keep_current_sample)

    seurat_obj <- subset_safe(seurat_obj,barcodes_to_keep)
    # seurat_obj.filtered <- subset(seurat_obj, subset = emptyDrops_FDR <= FDR)
  } else {
    seurat_obj <- seurat_obj
  }

  print(paste0("Cells per sample after filter for sample ", sample_id))
  print(table(seurat_obj$orig.ident))

  # update config
  config$filterSettings$FDR <- FDR

  #Populate plots list
  plots <-list()
  plots[generate_plotuuid(sample_id, task_name, 0)] <- list(plot1_data)
  # plots[generate_plotuuid(sample_id, task_name, 1)] <- list(plot2_data)
  # > head(plots[[1]])
  # [[1]]
  # FDR    log_u
  # 0.000000 3.554852

  # [[2]]
  # FDR        log_u
  # 0.0001120495 2.7176705030

  # the result object will have to conform to this format: {data, config, plotData : {plot1, plot2}}
  result <- list(
                 data = seurat_obj,
                 config = config,
                 plotData = plots
  )

  return(result)
}
