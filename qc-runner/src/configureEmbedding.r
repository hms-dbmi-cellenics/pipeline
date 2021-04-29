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


source('utils.r')

task <- function(scdata, config, task_name, sample_id){

    # Check wheter the filter is set to true or false
    # Q: Can we disable the configureEmbedding?
    if (is.null(config$enabled) || as.logical(toupper(config$enabled))){
        scdata.embedding <- run_Embedding(scdata, config)
        scdata.embedding <- run_Clustering(scdata.embedding, config)

    }else{
        scdata.embedding <- scdata
    }

    plots <- list()

    # the result object will have to conform to this format: {data, config, plotData : {plot1, plot2}}
    result <- list(
        data = scdata.embedding,
        config = config,
        plotData = plots
    )

    return(result)

}

# This function covers
#   - Compute embedding: t-SNE or UMAP
run_Embedding <- function(scdata, config){
    
    #################
    # Embedding part
    #################

    # The threshold was selected in the dataIntegration step
    # Until we update the rds in S3 where the numPCs is stored in misc, we need to handle it by HARDCODE to 30
    if ("numPCs" %in% names(scdata@misc))
        pca_nPCs <- scdata@misc[["numPCs"]]
    else
        pca_nPCs <- 30


    embed_method <- config$embeddingSettings$method
    embed_settings <- config$embeddingSettings$methodSettings[[embed_method]]

    if (embed_method == "umap"){
        message("Running embedding --> umap")
        # I get same embedding moving to Seurat default of 'uwot' with return.model = TRUE
        # I think it makes sense to do this as uwot-learn is being deprecated
        # seq_len avoids e.g. ndims=c(1,0) if user selects 0pcs (a bit more defensive)
        # I'm gonna add a ticket to make sure UI can't ask for 0 or 1 PCs as this will cause error anyways
        # please delete these comment
        scdata <- Seurat::RunUMAP(scdata, 
                                reduction = 'pca', 
                                dims = seq_len(pca_nPCs), 
                                verbose = FALSE, 
                                min.dist = embed_settings$minimumDistance, 
                                metric = embed_settings$distanceMetric, 
                                return.model = TRUE)
      
    } else if (embed_method == "tsne"){
      
        if (as.logical(toupper(config$auto))){
            perplexity <- min(30, ncol(scdata)/100)
            learningRate <- max(200, ncol(scdata)/12)
          
        }  else {
            perplexity <- embed_settings$perplexity
            learningRate <- embed_settings$learningRate
        }
      
        message("Running embedding --> tsne")
      
        scdata <- Seurat::RunTSNE(scdata,
                                reduction = 'pca',
                                dims = seq_len(pca_nPCs),
                                perplexity = perplexity,
                                learning.rate = learningRate)
      
    }

    return(scdata)

}

run_Clustering <- function(scdata, config){

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

    return(scdata)
}

