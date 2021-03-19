################################################
## STEP 2. Mitochondrial content filter 
#################################################
# This is a simplest filter that looks at a threshold value for the mitochondrial content.
# To separate cells with low MT-content from the ones that have a high MT-content what makes us think that are dead.  
# This can be a useful first guess. The settings for such a filter can also contain a simple "probabilityThreshold" setting. 

# The most uses values in MT-content are between [0.1, 0.2]. There are not too much literature about how to compute
# a threshold. For now, we will offer two methods:
# --> Absolute threshold: In order to be not too extrictive the threshold is set to 0.1
generate_default_values_mitochondrialContent <- function(scdata, config) {
   
   if (config$filterSettings$method == "absolute_threshold")
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

task <- function(scdata, config){
    # Check if the experiment has MT-content
    if (!"percent.mt"%in%colnames(scdata@meta.data)){
        #This return value breaks the step, no config data is returned!!!!
        message("Warning! No MT-content has been computed for this experiment!")
        return(scdata)
    }
    #The absolute threshold object is an atomic vector, not a list!
    # TODO: The [[1]] is a temp fix because the maxFraction in config is for some unknown reason a list. We need a ticket and ask for a change in that.
    maxFraction <- config$filterSettings$methodSettings[[config$filterSettings$method]][["maxFraction"]][[1]]
    
    # Check if it is required to compute sensible values. From the function 'generate_default_values_doubletScores', it is expected
    # to get a list with only one element --> maxFraction.

    #WARNING
    #TODO: auto doesnt exist, we need to uncomment this when the parameter is set!!
    #if (as.logical(toupper(config$auto)))
    #    maxFraction <- generate_default_values_mitochondrialContent(scdata, config)
    
    
    # Check whether the filter is set to true or false
    if (as.logical(toupper(config$enabled)))
        # Information regarding MT-content is pre-computed during the 'data-ingest'. 
        scdata.filtered <- subset(scdata, subset = percent.mt <= maxFraction*100)
    else
        scdata.filtered <- scdata
    plot1_data <- lapply(unname(scdata$percent.mt),function(x) {c("fracMito"=x)})
    plot2_data <- unname(purrr::map2(scdata$percent.mt,scdata$nCount_RNA,function(x,y){c("cellSize"=y,"fracMito"=x)}))

    # update config
    config$filterSettings$methodSettings[[config$filterSettings$method]][["maxFraction"]] <- maxFraction
    
    # the result object will have to conform to this format: {data, config, plotData : {plot1, plot2}}
    result <- list(
        data = scdata.filtered, # scdata filter
        config = config,
        plotData = list(
            # plot 1: histgram of MT-content
            # AAACCCAAGCGCCCAT-1 AAACCCAAGGTTCCGC-1 AAACCCACAGAGTTGG-1
            #              0.161              0.198              0.284  ...
            mitochondrialFractionHistogram = plot1_data,
            # plot 2: There are two alternavitive:
            #           - Scatter plot with UMIs in the x-axis and MT-content in the y-axis
            #           --> code: plot2 = list(u=scdata$nCount_RNA.mt, "MT-content" = scdata$percent.mt)
            #           - Barplot representing in the x-axis the log10(UMIs) and in the y-axis the MT-content. This option is the one 
            #           that is shown in the mockup.
            #           --> code: plot2 = list(log_10_UMIs=log10(scdata$nCount_RNA), MT_content =mscdata$percent.mt)
            # We have decided to use the scatter plot, but I temporaly leave the other option in the comments. 
            # Q: Should we return from the R side the cells that are going to be removed? For this plot it is interesting to color the
            # cells that are going to be excluded. 
            mitochondrialFractionLogHistogram = plot2_data
        )
    )
    return(result)
}
