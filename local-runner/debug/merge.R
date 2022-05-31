myargs = commandArgs(trailingOnly=TRUE)


number_to_merge = myargs[1]
scdata_list <- readRDS('scdata_list.rds')
y <- scdata_list[2:number_to_merge]
length(y)
scdata <- merge(scdata_list[[1]], y = scdata_list)

