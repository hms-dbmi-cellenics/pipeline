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


task <- function(seurat_obj, config){
    # config$filterSettings = list(minProbability=0.82, bandwidth=-1, filterThreshold=-1)
    # Check wheter the filter is set to true or false
    # For some reason the last children of named lists are computed as vectors, so we can't access them as recursive objects. 
    minProbability <- config$filterSettings[["minProbability"]]
    # THIS SLOT DOESNT EXIST ON LOCAL OBJECT.
    plot1_data_1  <- seurat_obj@meta.data$emptyDrops_FDR
    plot1_data_2  <- log10(seurat_obj$nCount_RNA)

    plot1_data <- unname(purrr::map2(plot1_data_1,plot1_data_2,function(x,y){c("FDR"=x,"log_u"=y)}))
    # Check if it is required to compute sensible values. From the function 'generate_default_values_classifier', it is expected
    # to get a list with two elements {minProbabiliy and filterThreshold}.
    if (exists('auto', where=config)){
        if(as.logical(toupper(config$auto)))
            filterSettings <- generate_default_values_classifier(seurat_obj, config)
    }
    # TODO: get flag from here: seurat_obj@tools$flag_filtered <- FALSE
    if(!as.logical(toupper(config$enabled))) {
        # is.cell <- meta.data$emptyDrops_FDR <= 0.01
        # sum(is.cell, na.rm=TRUE) 
        # table(Limited=meta.data$emptyDrops_Limited, Significant=is.cell)
        # is.cell2<-is.cell
        # is.cell2[is.na(is.cell2)]<-FALSE
        # sce.filt<-sce[,is.cell2]
        seurat_obj.filtered <- subset(seurat_obj, subset = emptyDrops_FDR <= minProbability)
    } else {
        seurat_obj.filtered <- seurat_obj
    }
    # update config
    config$filterSettings$minProbability <- minProbability

    # the result object will have to conform to this format: {data, config, plotData : {plot1, plot2}}
    result <- list(
        data = seurat_obj.filtered,
        config = config,
        plotData = list(
            classifierEmptyDropsPlot = plot1_data
        )
    )

    return(result)
}
