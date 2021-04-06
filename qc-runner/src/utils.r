generate_plotuuid <- function(sample_uuid, task_name, plot_idx) {

  if(sample_uuid != "") {
    return(paste(sample_uuid, task_name, plot_idx, sep="-"))
  }

  return(paste(task_name, plot_idx, sep="-"))

}

#
# Subset safe allows us to attempt to subset a seurat object even with an empty list
# this function exists for the case when the user over-filters the object, and we need to return something
# that'd allow the user to realize that they are filtering all the cells, while maintaining certain seurat functionality.
# it's a questionable function and it should be questioned.
#
# IN seurat_obj: object to filter
# IN cells: cell barcodes to subset the object with
# 
subset_safe <- function(seurat_obj,cells){
  if(length(cells)>0){
    return(subset(seurat_obj, cells = cells))
  }else{
    return(subset(seurat_obj,cells = colnames(seurat_obj)[1]))
  }
}


get_sample_barcodes <- function(obj_metadata, sample_id) {
  tmp_sample <- sub("sample-","",sample_id)
  
  barcode_names_this_sample <- rownames(obj_metadata[grep(tmp_sample, rownames(obj_metadata)),]) 
  
  # sample names not currently appended to cell ids for at least unisample
  if (!length(barcode_names_this_sample)) 
    barcode_names_this_sample <- row.names(obj_metadata)[obj_metadata$samples == sample_id]
  
  return(barcode_names_this_sample)
}