################################################
## STEP 7. Compute embedding
#################################################
# Compute embedding step where we run dimensional reduction technniques such as t-SNE and  UMAP. Moreover, the cluster analysis is
# done also in this step. 

############ NEED TO CHECK CONFIG SCHEMA

# nPCs ?

#"configureEmbedding": {
#    "embeddingSettings": {
#        "method": "umap",
#        "methodSettings": {
#            "umap": {
#                "minimumDistance": 0.2,
#                "distanceMetric": "euclidean"
#            },
#            "tsne": {
#                "perplexity": 30,
#                "learningRate": 200
#            }
#        }
#    },
#    "clusteringSettings": {
#        "method": "louvain",
#        "methodSettings": {
#            "louvain": {
#                "resolution": 0.5
#            }
#        }
#    }
#},

task <- function(scdata, config){

    # Check wheter the filter is set to true or false
    # Q: Can we disable the computeEmbedding?
    #if (as.logical(toupper(config$enabled)))
    scdata.embedding <- run_computeEmbedding(scdata, config)
    #else
    #    scdata.embedding <- scdata.embedding

    cells_order <- rownames(scdata.embedding@meta.data)

    # Create a generic pattern with the UMAP coordinates
    embeddingPreview <- unname(purrr::map2(
        scdata.embedding@reductions$umap@cell.embeddings[cells_order, 1],
        scdata.embedding@reductions$umap@cell.embeddings[cells_order, 2],
        function(x,y) {list("x"=x,"y"=y)}
        )
    )

    # Create plot for cellsets
    embeddingPreviewByCellSets <- purrr::map2(embeddingPreview,
        unname(scdata.embedding@active.ident[cells_order]),
        function(x,y){append(x,list("cluster"=paste("Cluster", y, sep = " ")))}
    )
    # Adding color for cellsets
    embeddingPreviewByCellSets <- purrr::map2(embeddingPreviewByCellSets,
        unname(scdata.embedding@meta.data[cells_order, "color_active_ident"]),
        function(x,y){append(x,list("col"=y))}
    )

    # Create plot for samples
    embeddingPreviewBySamples <- purrr::map2(embeddingPreview,
        unname(scdata.embedding@meta.data[cells_order, "type"]),
        function(x,y){append(x,list("sample"=y))}
    )
    
    # Adding color for samples
    embeddingPreviewBySamples <- purrr::map2(embeddingPreviewBySamples,
        unname(scdata.embedding@meta.data[cells_order, "color_samples"]),
        function(x,y){append(x,list("col"=y))}
    )

    # Create plot for MT-content
    embeddingPreviewMitochondrialContent <- purrr::map2(embeddingPreview,
        unname(scdata.embedding@meta.data[cells_order, "percent.mt"]),
        function(x,y){append(x,list("mt-content"=y))}
    )

    # Create plot for Doublet-score
    embeddingPreviewDoubletScore <- purrr::map2(embeddingPreview,
        unname(scdata.embedding@meta.data[cells_order, "doublet_scores"]),
        function(x,y){append(x,list("doublet-score"=y))}
    )

    # the result object will have to conform to this format: {data, config, plotData : {plot1, plot2}}
    result <- list(
        data = scdata.embedding,
        config = config,
        plotData = list(
                embeddingPreviewByCellSets=embeddingPreviewByCellSets,
                embeddingPreviewBySamples=embeddingPreviewBySamples,
                embeddingPreviewMitochondrialContent=embeddingPreviewMitochondrialContent,
                embeddingPreviewDoubletScore=embeddingPreviewDoubletScore
        )
    )

    return(result)

}

# This function covers
#   - Compute embedding: t-SNE or UMAP
#   - Clustering
run_computeEmbedding <- function(scdata, config){
    
    #################
    # Embedding part
    #################

    # The threshold was selected in the dataIntegration step
    # Until we update the rds in S3 where the numPCs is stored in misc, we need to handle it by HARDCODE to 30
    if ("numPCs" %in% names(scdata@misc))
        pca_nPCs <- scdata@misc[["numPCs"]]
    else
        pca_nPCs <- 30


    if (config$embeddingSettings$method=="umap"){
        
        minimumDistance <- config$embeddingSettings$methodSettings$umap$minimumDistance
        distanceMetric <- config$embeddingSettings$methodSettings$umap$distanceMetric

        message("Running embedding --> umap")
        scdata <- Seurat::RunUMAP(scdata, reduction='pca', dims = 1:pca_nPCs, verbose = F, umap.method = "uwot-learn", min.dist = minimumDistance, metric = distanceMetric)


    }

    if (config$embeddingSettings$method=="tsne"){
        
        perplexity <- config$embeddingSettings$methodSettings$tsne$perplexity
        learningRate <- config$embeddingSettings$methodSettings$tsne$learningRate
        
        if (as.logical(toupper(config$auto))){
            perplexity <- min(30, ncol(scdata)/100)
            learningRate <- max(200, ncol(scdata)/12)
        }


        message("Running embedding --> tsne")
        
        scdata <- Seurat::RunTSNE(scdata, reduction = 'pca', dims = 1:pca_nPCs, perplexity = perplexity, learning.rate = learningRate)

    }


    #################
    # Clustering part
    #################

    message("Running clustering")

    if(config$clusteringSettings$method=="louvain"){
        clustering_method <- 1 #"Louvain"
        clustering_resolution <- config$clusteringSettings$methodSettings[["louvain"]]["resolution"]

        # HARDCODE
        annoy.metric = "cosine"
        scdata <- Seurat::FindNeighbors(scdata, k.param = 20, annoy.metric = annoy.metric, verbose=FALSE)
        scdata <- Seurat::FindClusters(scdata, resolution=clustering_resolution, verbose = FALSE, algorithm = clustering_method)

    }


    ##########################
    # Coloring active ident
    ###########################
    if ("color_pool" %in% names(scdata@misc)){
        color_pool <- scdata@misc[['color_pool']]
    }else{ # THIS SHOULD BE REMOVE ONCE THE EXPERIMENT HAS BEEN UPDATED WITH THE NEW VERSION OF THE DATA-INGEST 
        color_pool <- c("#e377c2","#8c564b","#d62728","#2ca02c","#ff7f0e","#1f77b4","#f8e71c","#3957ff","#d3fe14",
        "#c9080a","#fec7f8","#0b7b3e","#0bf0e9","#c203c8","#fd9b39","#888593","#906407","#98ba7f","#fe6794","#10b0ff",
        "#ac7bff","#fee7c0","#964c63","#1da49c","#0ad811","#bbd9fd","#fe6cfe","#297192","#d1a09c","#78579e","#81ffad",
        "#739400","#ca6949","#d9bf01","#646a58","#d5097e","#bb73a9","#ccf6e9","#9cb4b6","#b6a7d4","#9e8c62","#6e83c8",
        "#01af64","#a71afd","#cfe589","#d4ccd1","#fd4109","#bf8f0e","#2f786e","#4ed1a5","#d8bb7d","#a54509","#6a9276",
        "#a4777a","#fc12c9","#606f15","#3cc4d9","#f31c4e","#73616f","#f097c6","#fc8772","#92a6fe","#875b44","#699ab3",
        "#94bc19","#7d5bf0","#d24dfe","#c85b74","#68ff57","#b62347","#994b91","#646b8c","#977ab4","#d694fd","#c4d5b5",
        "#fdc4bd","#1cae05","#7bd972","#e9700a","#d08f5d","#8bb9e1","#fde945","#a29d98","#1682fb","#9ad9e0","#d6cafe",
        "#8d8328","#b091a7","#647579","#1f8d11","#e7eafd","#b9660b","#a4a644","#fec24c","#b1168c","#188cc1","#7ab297",
        "#4468ae","#c949a6")

    }
    scdata$color_active_ident <- color_pool[as.numeric(scdata@active.ident)]

    ##########################
    # Coloring samples
    ###########################
    remaining.colors <- color_pool[-c(1:length(unique(scdata@meta.data$color_active_ident)))]
    if ("type"%in%colnames(scdata@meta.data)) # In that case we are in multisample experiment
        scdata@meta.data[, "color_samples"] <- remaining.colors[as.numeric(as.factor(scdata$type))]
    else
        scdata@meta.data[, "color_samples"] <- remaining.colors[1]

    return(scdata)
}
